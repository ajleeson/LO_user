"""
Plot surface and deep budget for all 21 inlets
"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import matplotlib.patches as mpatches
from scipy.linalg import lstsq
import math
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
import csv
import cmocean
from scipy.stats import pearsonr
from scipy.stats import ttest_ind
from scipy.stats import f_oneway
from scipy.stats import kruskal
from scipy.stats import alexandergovern
import pingouin as pg
import matplotlib.pylab as plt
import gsw
import pickle
import get_two_layer

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
jobname = 'twentyoneinlets'
year = '2017'

# Pick ONE cluster below to analyze

group = 'All inlets'
stations = 'all'

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# set up dates
startdate = year + '.01.01'
enddate = year + '.12.31'
enddate_hrly = str(int(year)+1)+'.01.01 00:00:00'

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get stations:
if stations == 'all':
    sta_dict = job_lists.get_sta_dict(jobname)
    # remove lynchcove2
    del sta_dict['lynchcove2']
    # remove shallow inlets (< 10 m deep)
    del sta_dict['hammersley']
    del sta_dict['henderson']
    del sta_dict['oak']
    del sta_dict['totten']
    del sta_dict['similk']
    del sta_dict['budd']
    del sta_dict['eld']
    del sta_dict['killsut']
else:
    sta_dict = stations

    
# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
Lfun.make_dir(out_dir)

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')[2::]
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]
# crop time vector (because we only have jan 2 - dec 30)
dates_no_crop = dates_local_daily
dates_local_daily = dates_local_daily

print('\n')

##########################################################
##            Get all variables for analysis            ##
##########################################################

print('Getting all data for analysis\n')

# create dictionaries with interface depths
interface_dict = dict()

# create dictionaries to save dataframe of budget terms for each station
bottomlay_dict = {}
surfacelay_dict = {}
DOconcen_dict = {}
dimensions_dict = {}
for i,station in enumerate(sta_dict):
    bottomlay_dict[station] = pd.DataFrame()
    surfacelay_dict[station] = pd.DataFrame()
    DOconcen_dict[station] = pd.DataFrame()
    dimensions_dict[station] = pd.DataFrame()

# get lat and lon of grid
Ldir['ds0'] = startdate
in_dir = Ldir['roms_out'] / Ldir['gtagex']
# G, S, T = zrfun.get_basic_info(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
fn0 = xr.open_dataset(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
lonr = fn0.lon_rho.values
latr = fn0.lat_rho.values
# open box extraction
box_fn = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
ds_box = xr.open_dataset(box_fn)
DX = (ds_box.pm.values)**-1
DY = (ds_box.pn.values)**-1
DA = DX*DY # get area of each grid cell in m^2

# COLLAPSE
for i,station in enumerate(sta_dict):

        # get interface depth from csv file
        with open('interface_depths.csv', 'r') as f:
            for line in f:
                inlet, interface_depth = line.strip().split(',')
                interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
        z_interface = float(interface_dict[station])

# --------------------------- get EU exchange flow terms ----------------------------------------
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_EU_exchange' / (station + '.p')
        df_exchange = pd.read_pickle(fn)
        exchange_surf_unfiltered = df_exchange['surface [kmol/s]'].values 
        exchange_deep_unfiltered = df_exchange['deep [kmol/s]'].values
        # Godin filter
        EU_surf = zfun.lowpass(exchange_surf_unfiltered, f='godin')[36:-34:24] 
        EU_deep = zfun.lowpass(exchange_deep_unfiltered, f='godin')[36:-34:24]
        exchange_color = 'cornflowerblue'

# --------------------------- get TEF exchange flow terms ----------------------------------------
        in_dir = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
        bulk = xr.open_dataset(in_dir)
        tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
        Q_p = tef_df['q_p'] # Qin [m3/s]
        Q_m = tef_df['q_m'] # Qout [m3/s]
        DO_p = tef_df['oxygen_p'] # DOin [mmol/m3]
        DO_m = tef_df['oxygen_m'] # DOout [mmol/m3]
        # convert from mmol/s to kmol/s
        TEF_surf = (Q_m.values * DO_m.values) * (1/1000) * (1/1000) # Qout * DOout
        TEF_deep = (Q_p.values * DO_p.values) * (1/1000) * (1/1000) # Qout * DOout

# ---------------------------------- get BGC terms --------------------------------------------
        bgc_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_bgc' / station
        # get months
        months = [year+'.01.01_'+year+'.01.31',
                    year+'.02.01_'+year+'.02.28',
                    year+'.03.01_'+year+'.03.31',
                    year+'.04.01_'+year+'.04.30',
                    year+'.05.01_'+year+'.05.31',
                    year+'.06.01_'+year+'.06.30',
                    year+'.07.01_'+year+'.07.31',
                    year+'.08.01_'+year+'.08.31',
                    year+'.09.01_'+year+'.09.30',
                    year+'.10.01_'+year+'.10.31',
                    year+'.11.01_'+year+'.11.30',
                    year+'.12.01_'+year+'.12.31',]
        
        # initialize arrays to save values
        photo_surf_unfiltered = []
        photo_deep_unfiltered = []
        photo_color = 'turquoise'
        cons_surf_unfiltered = [] # nitrification, respiration
        cons_deep_unfiltered = [] # nitrification, respiration, sediment oxygen demand
        cons_color = 'deeppink'
        airsea_surf_unfiltered = []
        airsea_color = 'plum'
        o2vol_surf_unfiltered = []
        o2vol_deep_unfiltered = []
        ddtDOV_color = 'gray'

        # combine all months
        for month in months:
            fn = station + '_' + month + '.p'
            df_bgc = pd.read_pickle(bgc_dir/fn)
            # conversion factor to go from mmol O2/hr to kmol O2/s
            conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
            # get photosynthesis
            photo_surf_unfiltered = np.concatenate((photo_surf_unfiltered, df_bgc['surf photo [mmol/hr]'].values * conv)) # kmol/s
            photo_deep_unfiltered = np.concatenate((photo_deep_unfiltered, df_bgc['deep photo [mmol/hr]'].values * conv)) # kmol/s
            # get consumption
            surf_cons_terms = df_bgc['surf nitri [mmol/hr]'].values + df_bgc['surf respi [mmol/hr]'].values
            deep_cons_terms = df_bgc['deep nitri [mmol/hr]'].values + df_bgc['deep respi [mmol/hr]'].values + df_bgc['SOD [mmol/hr]'].values
            cons_surf_unfiltered = np.concatenate((cons_surf_unfiltered, surf_cons_terms * conv * -1)) # kmol/s; multiply by -1 b/c loss term
            cons_deep_unfiltered = np.concatenate((cons_deep_unfiltered, deep_cons_terms * conv * -1)) # kmol/s; multiply by -1 b/c loss term
            # get air-sea gas exchange
            airsea_surf_unfiltered = np.concatenate((airsea_surf_unfiltered, df_bgc['airsea [mmol/hr]'].values * conv)) # kmol/s
            # get (DO*V)
            o2vol_surf_unfiltered = np.concatenate((o2vol_surf_unfiltered, df_bgc['surf DO*V [mmol]'].values)) # mmol
            o2vol_deep_unfiltered = np.concatenate((o2vol_deep_unfiltered, df_bgc['deep DO*V [mmol]'].values)) # mmol

        # take time derivative of (DO*V) to get d/dt (DO*V)
        ddtDOV_surf_unfiltered = np.diff(o2vol_surf_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
        ddtDOV_deep_unfiltered = np.diff(o2vol_deep_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s

        # apply Godin filter
        photo_surf = zfun.lowpass(photo_surf_unfiltered, f='godin')[36:-34:24]
        photo_deep = zfun.lowpass(photo_deep_unfiltered, f='godin')[36:-34:24]
        cons_surf   = zfun.lowpass(cons_surf_unfiltered, f='godin')[36:-34:24]
        cons_deep   = zfun.lowpass(cons_deep_unfiltered, f='godin')[36:-34:24]
        airsea_surf = zfun.lowpass(airsea_surf_unfiltered, f='godin')[36:-34:24]
        ddtDOV_surf = zfun.lowpass(ddtDOV_surf_unfiltered, f='godin')[36:-34:24]
        ddtDOV_deep = zfun.lowpass(ddtDOV_deep_unfiltered, f='godin')[36:-34:24]

# ------------------------------- get rivers and WWTPs ----------------------------------------
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_traps' / (station + '.p')
        df_traps = pd.read_pickle(fn)
        rivers_surf_unfiltered = df_traps['surface [kmol/s]'].values
        wwtps_deep_unfiltered = df_traps['deep [kmol/s]'].values

        # Pad with zeros if no rivers or WWTPs
        if rivers_surf_unfiltered.size == 0 or np.isnan(rivers_surf_unfiltered[0]):
                rivers_surf_unfiltered = np.zeros(8761)
        if wwtps_deep_unfiltered.size == 0:
                wwtps_deep_unfiltered = np.zeros(8761)
        # Godin filter
        traps_surf = zfun.lowpass(rivers_surf_unfiltered, f='godin')[36:-34:24]
        traps_deep =  zfun.lowpass(wwtps_deep_unfiltered, f='godin')[36:-34:24]
        traps_color = 'black'

# ------------------------------- get vertical exchange ----------------------------------------

        # pick color
        vertX_color = 'mediumorchid'

        # calculate error term, and ascribe to the vertical exchange
        vertX_surf_TEF = ddtDOV_surf - (TEF_surf 
                                    + photo_surf
                                    + cons_surf
                                    + airsea_surf
                                    + traps_surf)
        
        vertX_deep_TEF = ddtDOV_deep - (TEF_deep 
                                    + photo_deep
                                    + cons_deep
                                    + traps_deep)
        
        # calculate error term, and ascribe to the vertical exchange
        vertX_surf_EU = ddtDOV_surf - (EU_surf 
                                    + photo_surf
                                    + cons_surf
                                    + airsea_surf
                                    + traps_surf)
        
        vertX_deep_EU = ddtDOV_deep - (EU_deep 
                                    + photo_deep
                                    + cons_deep
                                    + traps_deep)
        
# ---------------------------------- DO concentrations --------------------------------------------

        # get average V of each layer
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
        df_V = pd.read_pickle(fn)
        # Godin filter already applied earlier in workflow
        surf_V = df_V['surface [m3]'].values[1:-1]
        deep_V = df_V['deep [m3]'].values[1:-1]

        # apply Godin filter to DO*V (mmol)
        o2vol_surf = zfun.lowpass(o2vol_surf_unfiltered, f='godin')[36:-34:24]
        o2vol_deep = zfun.lowpass(o2vol_deep_unfiltered, f='godin')[36:-34:24]

        # get average DO of each layer by calculating (DO*V)/V
        # and convert from mmol/m3 to mg/L
        o2_surf = o2vol_surf / surf_V * 32/1000
        o2_deep = o2vol_deep / deep_V * 32/1000

        # get average inflowing DO concentration (Qin concentration)
        DO_in = DO_p.values * 32/1000

        # get bottom DO
        # get section information
        # ds = xr.open_dataset('../../../../LO_output/extract/'+gtagex+
        #                     '/tef2/extractions_'+startdate+'_'+enddate+
        #                     '/'+station+'.nc')
        
        #% Get indices of inlet
        seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
        seg_df = pd.read_pickle(seg_name)
        ji_list = seg_df[station+'_p']['ji_list']
        jj_LO = [x[0] for x in ji_list] # y; lat; jj
        ii_LO = [x[1] for x in ji_list] # x; lon; ii
        # get lat and lon corresponding to ii and jj indices
        lat_LO = latr[jj_LO,0]
        lon_LO = lonr[0,ii_LO]
        # get corresponding ii and jj indices in box extraction
        lat_box_all = ds_box['lat_rho'].values[:,0]
        lon_box_all = ds_box['lon_rho'].values[0,:]
        jj = np.zeros(len(jj_LO))
        ii = np.zeros(len(ii_LO))
        for j_ind,lat in enumerate(lat_LO):
            jj[j_ind] = np.where(lat_box_all==lat)[0][0]
        for i_ind,lon in enumerate(lon_LO):
            ii[i_ind] = np.where(lon_box_all==lon)[0][0]
        # convert to array of ints
        jj = jj.astype(int)
        ii = ii.astype(int)

        # get bottom DO from my bottom DO get_DO script extraction file
        ds_oxy = xr.open_dataset(Ldir['LOo'] / 'pugetsound_DO' / 'data' / (year + '_DO_info_' + 'withStraits' + '.nc'))

        # calculate average bottom oxygen throughout section
        bottom_oxygen_alldays = np.nanmean(ds_oxy['DO_bot'].values[:,jj,ii],axis=1) # mg/L
        bottom_oxygen = bottom_oxygen_alldays[1:-1]
        # calculate minimum bottom DO throughout section
        oxygen_min_alldays = np.nanmin(ds_oxy['DO_bot'].values[:,jj,ii],axis=1) # mg/L
        oxygen_min = oxygen_min_alldays[1:-1]

        # calculate percent hypoxic volume
        hyp_thick = ds_oxy['hyp_thick'].values # m

        # get hypoxic thickness and cell area corresponding to the inlet
        hyp_thick = hyp_thick[:,jj,ii] # m
        DA_inlet = DA[jj,ii]
        # calculate hypoxic volume
        hyp_vol = np.sum(hyp_thick * DA_inlet, axis=(1)) # m^3
        # crop to 363 days
        hyp_vol = hyp_vol[1:-1]


# -------------------------- get volume and mean depth of layers -------------------------

        # get volume of layers
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
        df_V = pd.read_pickle(fn)
        # Godin filter already applied earlier in workflow
        surf_V = df_V['surface [m3]'].values[1:-1]
        deep_V = df_V['deep [m3]'].values[1:-1]

# --------------------------- get other physical dimensions ------------------------------

        # get full volume of inlet, and mean depth
        fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'vol_df_cas7_c21.p'
        vol_df = pd.read_pickle(fn)
        inlet_vol = vol_df['volume m3'].loc[station+'_p']
        inlet_area = vol_df['area m2'].loc[station+'_p']
        mean_depth = inlet_vol / inlet_area
        # estimate aspect ratio from inlet area and width at mouth
        # aspect ratio = length/width (so longer inlet has larger ratio)
        # thus, aspect ratio = area / width^2
        # since this is an approximation, we use the nominal cell width of 500 m
        fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'sect_df_cas7_c21.p'
        mouth_df = pd.read_pickle(fn)
        mouth_df = mouth_df.loc[mouth_df['sn'] == station]
        # get number of horizontal and vertical entrance points
        uv_list = mouth_df['uv'].values.tolist()
        u_count = uv_list.count('u')
        v_count = uv_list.count('v')
        mouth_width = 500 * np.sqrt(u_count**2 + v_count**2)
        est_aspect_ratio = inlet_area / mouth_width**2
        # get mouth area, for calculating fluxes 

        
# ------------------------- save data in dataframe dict -----------------------------------
        # Note: everything is in units of kmol O2 /s

        # surface layer
        # surfacelay_dict[station]['EU Exchange Flow'] = EU_surf
        surfacelay_dict[station]['TEF Exchange Flow'] = TEF_surf
        # surfacelay_dict[station]['EU Recirculation'] = EU_surf + vertX_surf_EU
        surfacelay_dict[station]['Exchange Flow & Vertical'] = TEF_surf + vertX_surf_TEF
        surfacelay_dict[station]['Rivers'] = traps_surf
        surfacelay_dict[station]['Photosynthesis'] = photo_surf
        surfacelay_dict[station]['Bio Consumption'] = cons_surf
        surfacelay_dict[station]['Air-Sea Transfer'] = airsea_surf
        # surfacelay_dict[station]['EU Vertical'] = vertX_surf_EU
        surfacelay_dict[station]['Vertical Transport'] = vertX_surf_TEF
        surfacelay_dict[station]['d/dt(DO)'] = ddtDOV_surf
        surfacelay_dict[station]['Volume'] = surf_V

        # bottom layer
        bottomlay_dict[station]['d/dt(DO)'] = ddtDOV_deep
        # bottomlay_dict[station]['EU Exchange Flow'] = EU_deep
        bottomlay_dict[station]['TEF Exchange Flow'] = TEF_deep
        # bottomlay_dict[station]['EU Recirculation'] = EU_deep + vertX_deep_EU
        bottomlay_dict[station]['Exchange Flow & Vertical'] = TEF_deep + vertX_deep_TEF
        # bottomlay_dict[station]['EU Vertical'] = vertX_deep_EU
        bottomlay_dict[station]['Vertical Transport'] = vertX_deep_TEF
        bottomlay_dict[station]['WWTPs'] = traps_deep
        bottomlay_dict[station]['Photosynthesis'] = photo_deep
        bottomlay_dict[station]['Bio Consumption'] = cons_deep
        bottomlay_dict[station]['Photosynthesis & Consumption'] = photo_deep + cons_deep
        bottomlay_dict[station]['Volume'] = deep_V
        bottomlay_dict[station]['Qin m3/s'] = Q_p.values 

# ------------------------- save DO concentrations in dataframe dict -----------------------------------

        # mg/L units

        # DO concentrations
        DOconcen_dict[station]['Surface Layer DO'] = o2_surf
        DOconcen_dict[station]['Deep Layer DO'] = o2_deep
        DOconcen_dict[station]['Bottom Sigma DO'] = bottom_oxygen
        # DOconcen_dict[station]['Minimum Bottom Layer DO'] = oxygen_min
        DOconcen_dict[station]['Qin DO'] = DO_in
        DOconcen_dict[station]['percent hypoxic volume'] = hyp_vol/[inlet_vol] * 100

# ------------------------ save inlet dimensions in dataframe dict ----------------------------------------

        dimensions_dict[station]['Inlet volume'] = [inlet_vol] # m3
        dimensions_dict[station]['Mean depth'] = [mean_depth] # m
        # dimensions_dict[station]['Mouth width'] = mouth_width
        # dimensions_dict[station]['L/W aspect ratio'] = [est_aspect_ratio] # dimensionless

        # print keys
        if i == 0:
            print(list(surfacelay_dict[station].keys()))
            print(list(DOconcen_dict[station].keys()))
            print(list(dimensions_dict[station].keys()))


# write to pickle files
with open('deeplay_dict.pickle', 'wb') as handle:
    pickle.dump(bottomlay_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('shallowlay_dict.pickle', 'wb') as handle:
    pickle.dump(surfacelay_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('dimensions_dict.pickle', 'wb') as handle:
    pickle.dump(dimensions_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

# ##########################################################
# ##              Plot DO budget of every inlet           ##
# ##########################################################

# DO_budget = False

# residual = True # recirculation (combined horizontal and vertical exchange)
# show_EU = True

# # COLLAPSE
# if DO_budget == True:

#     print('Making DO budget time series')

#     for i,station in enumerate(sta_dict):
            
#             # get interface depth from csv file
#             with open('interface_depths.csv', 'r') as f:
#                 for line in f:
#                     inlet, interface_depth = line.strip().split(',')
#                     interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
#             z_interface = float(interface_dict[station])

#             # initialize figure
#             plt.close('all')
#             fig, ax = plt.subplots(4,1,figsize = (12,9),sharex=True)
#             # format figure
#             plt.suptitle(station + ': DO Budget (Godin Filter)',size=14)
#             for axis in [ax[0],ax[1],ax[2],ax[3]]:
#                 # axis.plot([dates_local[0],dates_local[-1]],[0,0],color='k')
#                 axis.set_xlim([dates_local[0],dates_local[-1]])
#                 if axis == ax[3]:
#                     axis.set_ylabel(r'DO [mg L$^{-1}$]')
#                 else:
#                     axis.set_ylabel('DO transport\n' + r'[kmol O$_2$ s$^{-1}$]')
#                 axis.set_facecolor('#EEEEEE')
#                 axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#                 axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
#                 axis.tick_params(axis='x', labelrotation=30)
#                 for border in ['top','right','bottom','left']:
#                     axis.spines[border].set_visible(False)
#             ax[3].set_xlabel(year)
#             ax[0].text(0.02,0.9,'(a) Surface [shallower than {} m]'.format(-1*z_interface),
#                             ha='left', va='top', transform=ax[0].transAxes, fontsize=12)
#             ax[1].text(0.02,0.9,'(b) Bottom [deeper than {} m]'.format(-1*z_interface),
#                             ha='left', va='top', transform=ax[1].transAxes, fontsize=12)
#             ax[2].text(0.02,0.9,r'(c) Error [sum of surface and bottom vertical exchange]',
#                             ha='left', va='top', transform=ax[2].transAxes, fontsize=12)
#             ax[3].text(0.02,0.9,'(d) Average DO concentrations',
#                             ha='left', va='top', transform=ax[3].transAxes, fontsize=12)
        

#             # plot surface
#             # exchange flow
#             if residual == False:
#                 ax[0].plot(dates_local_daily,surfacelay_dict[station]['TEF Exchange Flow'],color=exchange_color,
#                         linewidth=1,label='TEF Exchange Flow')
#                 if show_EU:
#                     ax[0].plot(dates_local_daily,surfacelay_dict[station]['EU Exchange Flow'],color=exchange_color,
#                             linewidth=3,alpha=0.3,label='EU Exchange Flow')
#             # recirculation
#             else:
#                 ax[0].plot(dates_local_daily,surfacelay_dict[station]['TEF Recirculation'],color=exchange_color,
#                         linewidth=2,label='Recirculation TEF')
#                 if show_EU:
#                     ax[0].plot(dates_local_daily,surfacelay_dict[station]['EU Recirculation'],color=exchange_color,
#                             linewidth=3,alpha=0.3,label='Recirculation EU')
#             # rivers
#             ax[0].plot(dates_local_daily,surfacelay_dict[station]['TRAPS'],color=traps_color,
#                        linewidth=3,zorder=2,label='TRAPS')
#             # photosynthesis
#             ax[0].plot(dates_local_daily,surfacelay_dict[station]['Photosynthesis'],color=photo_color,
#                        linewidth=2,label='Photosynthesis')
#             # consumption
#             ax[0].plot(dates_local_daily,surfacelay_dict[station]['Bio Consumption'],color=cons_color,
#                        linewidth=2,linestyle=':',zorder=9,label='Bio Consumption')
#             # air-sea gas exchange
#             ax[0].plot(dates_local_daily,surfacelay_dict[station]['Air-Sea Transfer'],color=airsea_color,
#                        linewidth=1,zorder=8,label='Air-Sea Transfer')
#             # storage
#             ax[0].plot(dates_local_daily,surfacelay_dict[station]['Storage'],color=ddtDOV_color,
#                        linewidth=1,zorder=7,#alpha=0.6,
#                        label=r'$\frac{\mathrm{d}}{\mathrm{dt}}(\mathrm{DO}\cdot V)$')
#             # vertical exchange
#             if residual == False:
#                 ax[0].plot(dates_local_daily,surfacelay_dict[station]['TEF Vertical'],color=vertX_color,
#                        linewidth=1,label='TEF Vertical')
#                 if show_EU:
#                     ax[0].plot(dates_local_daily,surfacelay_dict[station]['EU Vertical'],color=vertX_color,
#                             linewidth=3,alpha=0.3,label='EU Vertical')
#             if residual == True:
#                 ax[0].legend(loc='lower center',ncol=6)
#             else:
#                 ax[0].legend(loc='lower center',ncol=5)
            
#             # plot deep
#             if residual == False:
#                 ax[1].plot(dates_local_daily,bottomlay_dict[station]['TEF Exchange Flow'],color=exchange_color,
#                         linewidth=1)
#                 if show_EU:
#                     ax[1].plot(dates_local_daily,bottomlay_dict[station]['EU Exchange Flow'],color=exchange_color,
#                             linewidth=3,alpha=0.3)
#             else:
#                 ax[1].plot(dates_local_daily,bottomlay_dict[station]['TEF Recirculation'],color=exchange_color,
#                         linewidth=2)
#                 if show_EU:
#                     ax[1].plot(dates_local_daily,bottomlay_dict[station]['EU Recirculation'],color=exchange_color,
#                             linewidth=3,alpha=0.3)
#             ax[1].plot(dates_local_daily,bottomlay_dict[station]['TRAPS'],color=traps_color,
#                        linewidth=3,zorder=2)
#             ax[1].plot(dates_local_daily,bottomlay_dict[station]['Photosynthesis'],color=photo_color,
#                        linewidth=2)
#             ax[1].plot(dates_local_daily,bottomlay_dict[station]['Bio Consumption'],color=cons_color,
#                        linewidth=2,linestyle=':',zorder=9)
#             ax[1].plot(dates_local_daily,bottomlay_dict[station]['Storage'],color=ddtDOV_color,
#                        linewidth=1,zorder=7)#,alpha=0.6)
#             if residual == False:
#                 ax[1].plot(dates_local_daily,bottomlay_dict[station]['TEF Vertical'],color=vertX_color,
#                         linewidth=1)
#                 if show_EU:
#                     ax[1].plot(dates_local_daily,bottomlay_dict[station]['EU Vertical'],color=vertX_color,
#                             linewidth=3,alpha=0.3)

#             # plot error
#             error_TEF = surfacelay_dict[station]['TEF Vertical']+bottomlay_dict[station]['TEF Vertical']
#             error_EU = surfacelay_dict[station]['EU Vertical']+bottomlay_dict[station]['EU Vertical']
#             ax[2].plot(dates_local_daily,error_TEF,color='crimson',
#                        linewidth=2,label='TEF')
#             if show_EU:
#                 ax[2].plot(dates_local_daily,error_EU,color='k',
#                         linewidth=1,linestyle='--',label='EU')
#                 ax[2].legend(loc='upper right')
                
#             for axis in [ax[0],ax[1],ax[2]]:
#                 if residual == False:
#                     ylimval = np.nanmax([np.abs(surfacelay_dict[station]['TEF Exchange Flow']),
#                                          np.abs(surfacelay_dict[station]['TEF Vertical']),
#                                          np.abs(surfacelay_dict[station]['Photosynthesis']),
#                                          np.abs(surfacelay_dict[station]['TRAPS']),
#                                          np.abs(surfacelay_dict[station]['Bio Consumption']),])*1.2
#                 else:
#                     ylimval = np.nanmax([np.abs(surfacelay_dict[station]['TEF Recirculation']),
#                                          np.abs(surfacelay_dict[station]['Photosynthesis']),
#                                          np.abs(surfacelay_dict[station]['TRAPS']),
#                                          np.abs(surfacelay_dict[station]['Bio Consumption']),])*1.2
#                 axis.set_ylim([-1*ylimval,ylimval])

# # plot DO concentrations -------------------------------------------------------------------------------

#             # plot DO
#             ax[3].plot(dates_local_daily,DOconcen_dict[station]['Surface Layer'],color='royalblue',
#                        linewidth=2,label='Avg. surface layer')
#             ax[3].plot(dates_local_daily,DOconcen_dict[station]['Deep Layer'],color='mediumorchid',
#                        linewidth=2,label='Avg. deep layer')
#             ax[3].plot(dates_local_daily,DOconcen_dict[station]['Qin DO'],color='black',
#                        linewidth=1,linestyle=':',label='Avg. Qin DO')
#             ax[3].plot(dates_local_daily,DOconcen_dict[station]['Bottom Sigma DO'],color='crimson',
#                        linewidth=1,linestyle='--',label=r'Avg. bottom $\sigma$-layer')
#             ax[3].plot(dates_local_daily,DOconcen_dict[station]['Minimum Bottom Layer DO'],color='crimson',
#                        linewidth=1,linestyle='-',label='Minimum Bottom Layer DO')
#             ax[3].legend(loc='lower left',ncol=2)


#             ax[3].set_ylim([0,12])

# # ---------------------------------- save figure --------------------------------------------

#             plt.subplots_adjust(hspace=0.06, bottom=0.06, top=0.94)
#             if residual == True:
#                 if show_EU == True:
#                     out_dir_budget = out_dir / 'recirculation_withEU'
#                 else:
#                     out_dir_budget = out_dir / 'recirculation'
#             else:
#                 if show_EU == True:
#                     out_dir_budget = out_dir / 'horizvert_withEU'
#                 else:
#                     out_dir_budget = out_dir / 'horizvert'

#             Lfun.make_dir(out_dir_budget)
#             plt.savefig(out_dir_budget / (station + '.png') )

# ##########################################################
# ##               Deep Budget Error Analysis             ##
# ##########################################################

# deep_budget_error = False

# # COLLAPSE
# if deep_budget_error == True:

#     print('Making deep DO budget time series with error')

#     for i,station in enumerate(sta_dict):
            
#             # get interface depth from csv file
#             with open('interface_depths.csv', 'r') as f:
#                 for line in f:
#                     inlet, interface_depth = line.strip().split(',')
#                     interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
#             z_interface = float(interface_dict[station])

#             # initialize figure
#             plt.close('all')
#             fig, ax = plt.subplots(1,1,figsize = (12,9))
#             # format figure
#             plt.suptitle(station + ': Entire inlet (surf+deep) DO Budget (10-day Hanning Window)',size=14)
#             for axis in [ax]:
#                 # axis.plot([dates_local[0],dates_local[-1]],[0,0],color='k')
#                 axis.set_xlim([dates_local[0],dates_local[-1]])
#                 axis.set_ylabel('DO transport (mg/L per day)')
#                 axis.set_facecolor('#EEEEEE')
#                 axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#                 axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
#                 axis.tick_params(axis='x', labelrotation=30)
#                 for border in ['top','right','bottom','left']:
#                     axis.spines[border].set_visible(False)
#             ax.set_xlabel(year)
#             # ax.text(0.02,0.9,'Deep Layer [deeper than {} m]'.format(-1*z_interface),
#             #                 ha='left', va='top', transform=ax.transAxes, fontsize=12)

#             # # plot surf + deep
#             # ax.plot(dates_local_daily,
#             #         zfun.lowpass(bottomlay_dict[station]['TEF Exchange Flow'].values,n=10)+
#             #         zfun.lowpass(surfacelay_dict[station]['TEF Exchange Flow'].values,n=10),
#             #         color=exchange_color,
#             #         linewidth=1, label='Exchange Flow')
#             # ax.plot(dates_local_daily,
#             #         zfun.lowpass(bottomlay_dict[station]['TRAPS'].values,n=10)+
#             #         zfun.lowpass(surfacelay_dict[station]['TRAPS'].values,n=10),
#             #         color=ddtDOV_color,linewidth=1,zorder=7, label='TRAPS')
#             # ax.plot(dates_local_daily,
#             #         zfun.lowpass(bottomlay_dict[station]['Photosynthesis'].values,n=10)+
#             #         zfun.lowpass(surfacelay_dict[station]['Photosynthesis'].values,n=10),
#             #         color=photo_color,
#             #            linewidth=2, label='Photosynthesis')
#             # ax.plot(dates_local_daily,
#             #         zfun.lowpass(bottomlay_dict[station]['Bio Consumption'].values,n=10)+
#             #         zfun.lowpass(surfacelay_dict[station]['Bio Consumption'].values,n=10),
#             #         color=cons_color,
#             #            linewidth=2,linestyle=':',zorder=9, label='Bio Consumption')
#             # ax.plot(dates_local_daily,
#             #         zfun.lowpass(bottomlay_dict[station]['Storage'].values,n=10)+
#             #         zfun.lowpass(surfacelay_dict[station]['Storage'].values,n=10),
#             #         color=traps_color,linewidth=3,zorder=2, label='Storage')#,alpha=0.6)
#             # ax.plot(dates_local_daily,
#             #         zfun.lowpass(surfacelay_dict[station]['Air-Sea Transfer'].values,n=10),
#             #         color=airsea_color,
#             #            linewidth=1,zorder=8,label='Air-Sea Transfer')

#             # plot surf + deep (volume normalized: units of mg/L per day)
#             # TEF
#             TEF = (bottomlay_dict[station]['TEF Exchange Flow'].values +
#                 surfacelay_dict[station]['TEF Exchange Flow'].values) / (
#                 dimensions_dict[station]['Inlet volume'].values) * (1000 * 32 * 60 * 60 * 24)
#             ax.plot(dates_local_daily,
#                     zfun.lowpass(TEF,n=10),
#                     color=exchange_color,
#                     linewidth=1, label='Exchange Flow')
#             # TRAPS
#             TRAPS = (bottomlay_dict[station]['TRAPS'].values +
#                 surfacelay_dict[station]['TRAPS'].values) / (
#                 dimensions_dict[station]['Inlet volume'].values) * (1000 * 32 * 60 * 60 * 24)
#             ax.plot(dates_local_daily,
#                     zfun.lowpass(TRAPS,n=10),
#                     color=ddtDOV_color,linewidth=1,zorder=7, label='TRAPS')
#             # Photosynthesis
#             photosynthesis = (bottomlay_dict[station]['Photosynthesis'].values +
#                 surfacelay_dict[station]['Photosynthesis'].values) / (
#                 dimensions_dict[station]['Inlet volume'].values) * (1000 * 32 * 60 * 60 * 24)
#             ax.plot(dates_local_daily,
#                     zfun.lowpass(photosynthesis,n=10),
#                     color=photo_color,
#                        linewidth=2, label='Photosynthesis')
#             # Consumption
#             consumption = (bottomlay_dict[station]['Bio Consumption'].values +
#                 surfacelay_dict[station]['Bio Consumption'].values) / (
#                 dimensions_dict[station]['Inlet volume'].values) * (1000 * 32 * 60 * 60 * 24)
#             ax.plot(dates_local_daily,
#                     zfun.lowpass(consumption,n=10),
#                     color=cons_color,
#                        linewidth=2,linestyle=':',zorder=9, label='Bio Consumption')
#             # Storage
#             storage = (bottomlay_dict[station]['Storage'].values +
#                 surfacelay_dict[station]['Storage'].values) / (
#                 dimensions_dict[station]['Inlet volume'].values) * (1000 * 32 * 60 * 60 * 24)
#             ax.plot(dates_local_daily,
#                     zfun.lowpass(storage,n=10),
#                     color=traps_color,linewidth=3,zorder=2, label='Storage')#,alpha=0.6)
#             # Air-sea
#             airsea = surfacelay_dict[station]['Air-Sea Transfer'].values / (
#                 dimensions_dict[station]['Inlet volume'].values) * (1000 * 32 * 60 * 60 * 24)
#             ax.plot(dates_local_daily,
#                     zfun.lowpass(airsea,n=10),
#                     color=airsea_color,
#                        linewidth=1,zorder=8,label='Air-Sea Transfer')
            
#             # # plot deep
#             # if residual == False:
#             #     ax.plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TEF Exchange Flow'].values,n=10),color=exchange_color,
#             #             linewidth=1, label='Exchange Flow')
#             # else:
#             #     ax.plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TEF Recirculation'].values,n=10),color=exchange_color,
#             #             linewidth=2, label='Exchange Flow & Vertical Transport')
#             # ax.plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TRAPS'].values,n=10),color=traps_color,
#             #            linewidth=3,zorder=2, label='WWTPs')
#             # ax.plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['Photosynthesis'].values,n=10),color=photo_color,
#             #            linewidth=2, label='Photosynthesis')
#             # ax.plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['Bio Consumption'].values,n=10),color=cons_color,
#             #            linewidth=2,linestyle=':',zorder=9, label='Bio Consumption')
#             # ax.plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['Storage'].values,n=10),color=ddtDOV_color,
#             #            linewidth=1,zorder=7, label='Storage')#,alpha=0.6)
#             # if residual == False:
#             #     ax.plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TEF Vertical'].values,n=10),color=vertX_color,
#             #             linewidth=1,label='Vertical Transport')

#             # plot error
#             error_TEF = (surfacelay_dict[station]['TEF Vertical']+bottomlay_dict[station]['TEF Vertical'])/ (
#                 dimensions_dict[station]['Inlet volume'].values) * (1000 * 32 * 60 * 60 * 24)
#             ax.plot(dates_local_daily,zfun.lowpass(error_TEF.values,n=10),color='crimson',
#                        linewidth=2,label='Error (Sum of deep&surf vertical transport)')
#             ax.legend(loc='lower center',ncol=2)
                
#             print(station)
#             QinDOin = (bottomlay_dict[station]['TEF Exchange Flow'].values/dimensions_dict[station]['Inlet volume'].values) * (1000 * 32 * 60 * 60 * 24)
#             print('     Error (ann avg.) = ' + str(round(np.nanmean(np.abs(error_TEF)),3)))
#             print('     QinDOin (ann avg.) = ' + str(round(np.nanmean(np.abs(QinDOin)),3)))
#             print('     perc of QinDOin (ann avg.) = ' + str(round((np.nanmean(np.abs(error_TEF))/np.nanmean(np.abs(QinDOin)))*100,2)) + '%')

# # ---------------------------------- save figure --------------------------------------------

#             plt.subplots_adjust(hspace=0.06, bottom=0.06, top=0.94)
#             # out_dir_budget = out_dir / 'deep_budget_error'
#             out_dir_budget = out_dir / 'deepANDsurf_budget_error'

#             Lfun.make_dir(out_dir_budget)
#             plt.savefig(out_dir_budget / (station + '.png') )

#     plt.close('all')


# ##########################################################
# ##               Mid-Jul to mid-Aug error               ##
# ##########################################################

# # mid July to mid August
# # minday = 194
# # maxday = 225
# minday = 0
# maxday = -1

# tef_error = np.array([])
# cons_error = np.array([])

# for i,station in enumerate(sta_dict):
#     # get error
#     error = np.nanmean(bottomlay_dict[station]['TEF Vertical'][minday:maxday] + surfacelay_dict[station]['TEF Vertical'][minday:maxday])
#     # get other terms
#     exchange_flow  = np.nanmean(bottomlay_dict[station]['TEF Exchange Flow'][minday:maxday])
#     vertical_trans = np.nanmean(bottomlay_dict[station]['TEF Vertical'][minday:maxday])
#     photosynthesis = np.nanmean(bottomlay_dict[station]['Photosynthesis'][minday:maxday])
#     consumption    = np.nanmean(bottomlay_dict[station]['Bio Consumption'][minday:maxday])

#     print('\n--------------------------------')
#     print(station)
#     print('EF       BC\n')

#     print('{}%    {}%'.format(
#         round(error/exchange_flow*100,1),
#         round(error/consumption*100,1)
#     ))

#     tef_error = np.concatenate((tef_error,[error/exchange_flow*100]))
#     cons_error = np.concatenate((cons_error,[error/consumption*100]))

# # average errors
# print('\nAVERAGE ERRORS (TEF, cons) [%]')
# print(np.nanmean(np.abs(tef_error)))
# print(np.nanmean(np.abs(cons_error)))

# # ##########################################################
# # ##         Lynch Cove example budget time series        ##
# # ##########################################################

# # LC_budget = True

# # # COLLAPSE
# # if LC_budget == True:
# #     plt.close('all')
# #     for i,station in enumerate(['lynchcove']):
            
# #             # get interface depth from csv file
# #             with open('interface_depths.csv', 'r') as f:
# #                 for line in f:
# #                     inlet, interface_depth = line.strip().split(',')
# #                     interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
# #             z_interface = float(interface_dict[station])

# #             # initialize figure
            
# #             fig, ax = plt.subplots(2,1,figsize=(9,8))
# #             # format figure
# #             # plt.suptitle(station + ': DO Budget (10-day Hanning Window)',size=14)
# #             ax[0].set_title('(a) 2017 Lynch Cove Deep DO Budget (10-day Hanning Window)',size=14, loc='left')
# #             ax[0].set_xlim([dates_local[0],dates_local[-25]])
# #             ax[0].set_ylabel('DO transport ' + r'[kmol O$_2$ s$^{-1}$]',size=12)
# #             # ax[0].set_facecolor('#EEEEEE')
# #             # ax[0].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
# #             ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
# #             ax[0].tick_params(axis='x', labelrotation=30, labelsize=12)
# #             ax[0].tick_params(axis='y', labelsize=12)
# #             loc = mdates.MonthLocator(interval=1)
# #             ax[0].xaxis.set_major_locator(loc)
# #             ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
# #             # ax[0].tick_params(axis='both', labelsize=12)
# #             # for border in ['top','right','bottom','left']:
# #             #     ax[0].spines[border].set_visible(False)
# #             # ax[0].set_xlabel(year, fontsize=14)
# #             ax[0].set_ylim([-0.2,0.2])
# #             # plot deep budget
# #             nwin = 10 # hanning window length
# #             ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['Storage'].values,n=nwin),color='k',
# #                        linewidth=1,label=r'$\frac{d}{dt}\int_V$DO dV') #,alpha=0.6)
# #             ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TEF Vertical'].values + surfacelay_dict[station]['TEF Vertical'].values,n=nwin),
# #                        color='darkorange', linewidth=1,label='Error')
# #             ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TEF Exchange Flow'].values,n=nwin),color='#0D4B91',
# #                     linewidth=2,label='Exchange Flow')
# #             ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TEF Vertical'].values,n=nwin),color='#99C5F7',
# #                     linewidth=2,label='Vertical Transport')
# #             # ax.plot(dates_local_daily,bottomlay_dict[station]['TRAPS'],color=traps_color,
# #             #            linewidth=3,zorder=2)
# #             ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['Photosynthesis'].values,n=nwin),color='#8F0445',
# #                        linewidth=2, label='Photosynthesis')
# #             ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['Bio Consumption'].values,n=nwin),color='#FCC2DD',
# #                        linewidth=2,label='Bio Consumption')
# #             ax[0].legend(loc='upper center',ncol=3, fontsize=12)
# #             # add drawdown period
# #             # mid Jul - mid Aug
# #             minday = 194
# #             maxday = 225
# #             ax[0].axvline(dates_local_daily[minday],0,12,color='grey')#,linestyle=':')
# #             ax[0].axvline(dates_local_daily[maxday],0,12,color='grey')#,linestyle=':')

# #             # Bar chart     
# #             # ax[1].set_facecolor('#EEEEEE')
# #             ax[1].tick_params(axis='x', labelrotation=30)
# #             ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='y')
# #             # ax[1].grid(True,color='w',linewidth=1,linestyle='-',axis='y')
# #             # for border in ['top','right','bottom','left']:
# #             #     ax[1].spines[border].set_visible(False)
# #             ax[1].tick_params(axis='y', labelsize=12)
# #             ax[1].set_xticklabels([])
# #             ax[1].set_ylabel('mg/L per day',fontsize=12)
# #             # ax[1].set_title(station + ' volume-averaged deep layer DO budget')
# #             ax[1].set_title('(b) Lynch Cove volume-averaged deep layer DO budget (drawdown period)',size=14, loc='left')
# #             ax[1].set_xlim([-0.3,1.05])
# #             ax[1].set_ylim([-0.3,0.3])
# #             width = 0.2
# #             # create bar chart
# #             for attribute, measurement in bottomlay_dict[station].items():
# #                 # skip variables we are not interested in
# #                 if attribute in ['EU Exchange Flow',
# #                                     'TRAPS',
# #                                     'EU Recirculation',
# #                                     'EU Vertical',
# #                                     'TEF Recirculation',
# #                                     'Photosynthesis & Consumption',
# #                                     'Volume',
# #                                     'Qin m3/s']:
# #                     continue
# #                 # assign colors
# #                 if attribute == 'TEF Exchange Flow':
# #                     color = '#0D4B91'
# #                     label = 'Exchange Flow'
# #                     pos = 0.15
# #                 if attribute == 'TEF Vertical':
# #                     color = '#99C5F7'
# #                     label = 'Vertical Transport'
# #                     pos = 0.35
# #                 if attribute == 'Photosynthesis':
# #                     color = '#8F0445'
# #                     label = attribute
# #                     pos = 0.65
# #                 if attribute == 'Bio Consumption':
# #                     color = '#FCC2DD'
# #                     label = attribute
# #                     pos = 0.85
# #                 if attribute == 'Storage':
# #                     color = 'black'
# #                     label = r'$\frac{d}{dt}$DO (net decrease)'
# #                     pos = -0.15

# #                 # calculate time average
# #                 time_avg = np.nanmean(measurement[minday:maxday])
# #                 # get volume average
# #                 avg = time_avg/(np.nanmean(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
# #                 # convert to mg/L per day
# #                 avg = avg * 1000 * 32 * 60 * 60 * 24

# #                 # calculate standard deviation
# #                 # get volume average
# #                 Vavg = measurement[minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday] # kmol O2 /s /m3
# #                 # get standard deviation
# #                 std = np.std(Vavg)
# #                 std = std * 1000 * 32 * 60 * 60 * 24

# #                 ax[1].bar(pos, avg, width, zorder=5, align='center', edgecolor=color,color=color, label=label)
# #                 # ax[1].errorbar(pos, avg, yerr=std, capsize=5, color='grey',zorder=6)
# #                 if avg < 0:
# #                     wiggle = 0.03
# #                 if avg > 0:
# #                     wiggle = -0.03
# #                 ax[1].text(pos, wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                         color='black',fontsize=12)
                
# #                 # # zero line
# #                 # ax[1].axhline(0,-0.3,1.05,color='dimgrey')
                
# #                 ax[1].legend(loc='upper right', fontsize=12, ncol=2)

# #                 plt.tight_layout()
# #                 plt.show()


# # ##########################################################
# # ##              Comparing deep DO budgets               ## 
# # ##########################################################

# # # June and July
# # minday = 150
# # maxday = 211

# # # bar width
# # width = 0.1

# # letters = ['a','b']

# # # mid July to mid August
# # minday = 194
# # maxday = 225

# # hypinlet_budget_comparison_v2 = True # ------------------------------------------- MONEY BAR CHART & t-tests
# # # COLLAPSE
# # if hypinlet_budget_comparison_v2 == True:

# #     multiplier_deep1 = 0
# #     multiplier_deep2 = 0

# #     fig, ax = plt.subplots(2,1,figsize = (11,9))

# #     # format grid
# #     for axis in [ax[0],ax[1]]:
# #         # axis.set_facecolor('#EEEEEE')
# #         axis.tick_params(axis='x', labelrotation=30)
# #         # axis.grid(True,color='w',linewidth=1,linestyle='-',axis='y')
# #         # for border in ['top','right','bottom','left']:
# #         #     axis.spines[border].set_visible(False)
# #         axis.grid(True,color='silver',linewidth=1,linestyle='--',axis='y')
# #         axis.tick_params(axis='y', labelsize=12)
# #         axis.set_xticklabels([])
# #         axis.set_ylabel('mg/L per day',fontsize=14)
# #         axis.set_xlim([-0.5,1.05])
# #     ax[1].set_ylim([-0.25,0.25])

# #     # create a new dictionary of results
# #     oxy_dict = {}
# #     hyp_dict = {}
# #     oxy_dict_std = {}
# #     hyp_dict_std = {}

# #     # part 1 with distinct physical and biological terms
# #     for station in sta_dict:
# #         for attribute, measurement in bottomlay_dict[station].items():
# #             # skip variables we are not interested in
# #             if attribute in ['EU Exchange Flow',
# #                              'TRAPS',
# #                              'EU Recirculation',
# #                              'EU Vertical',
# #                              'TEF Recirculation',
# #                              'Photosynthesis & Consumption',
# #                              'Volume',
# #                              'Qin m3/s']:
# #                 continue
# #             # calculate time average normalized by volume
# #             avg = np.nanmean(measurement[minday:maxday]/(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
# #             # convert to mg/L per day
# #             avg = avg * 1000 * 32 * 60 * 60 * 24

# #             # get standard deviations
# #             Vavg = measurement[minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday] # kmol O2 /s /m3
# #             # calculate standard deviation in mg/L per day
# #             std = np.std(Vavg) * 1000 * 32 * 60 * 60 * 24

# #             # save values in dictionary
# #             if station in ['penn','case','holmes','portsusan','lynchcove','dabob']:
# #                 if attribute in hyp_dict.keys():
# #                     hyp_dict[attribute].append(avg)
# #                     hyp_dict_std[attribute].append(std)
# #                 else:
# #                     hyp_dict[attribute] = [avg]
# #                     hyp_dict_std[attribute] = [std]
# #             else:
# #                 if attribute in oxy_dict.keys():
# #                     oxy_dict[attribute].append(avg)
# #                     oxy_dict_std[attribute].append(std)
# #                 else:
# #                     oxy_dict[attribute] = [avg]
# #                     oxy_dict_std[attribute] = [std]
# #     # t-test
# #     print('DISTINCT TERMS --------------------------------------')
# #     for attribute in oxy_dict:
# #         print('\n========================')
# #         print(attribute)
# #         a = oxy_dict[attribute]
# #         b = hyp_dict[attribute]
# #         ttest = ttest_ind(a, b, axis=0, equal_var=False)
# #         print(ttest)
# #     print('\n')
# #     for i,dict in enumerate([oxy_dict,hyp_dict]):
# #     # average all oxygenated and hypoxic inlet rate values, and calculate standard deviations
# #         if i ==0:
# #             shift = 0
# #         else:
# #             shift = 0.1
# #         for attribute, measurement in dict.items():
# #             # choose color
# #             if attribute == 'TEF Exchange Flow':
# #                 color = '#0D4B91'
# #                 label = 'Exchange Flow'
# #                 pos = 0.15
# #             if attribute == 'TEF Vertical':
# #                 color = '#99C5F7'
# #                 label = 'Vertical Transport'
# #                 pos = 0.35
# #             if attribute == 'Photosynthesis':
# #                 color = '#8F0445'
# #                 label = attribute
# #                 pos = 0.65
# #             if attribute == 'Bio Consumption':
# #                 color = '#FCC2DD'
# #                 label = attribute
# #                 pos = 0.85
# #             if attribute == 'Storage':
# #                 color = 'black'
# #                 label = r'$\frac{d}{dt}$DO (net decrease)'
# #                 pos = -0.25

# #             # calculate average and standard deviation
# #             avg = np.nanmean(measurement)
# #             # get error
# #             if dict == oxy_dict:
# #                 std_dict = oxy_dict_std
# #             elif dict == hyp_dict:
# #                 std_dict = hyp_dict_std
# #             # propagate error
# #             error = 1/len(std_dict.keys()) * np.sqrt(np.sum(np.square(std_dict[attribute])))

# #             if avg < 0:
# #                 wiggle = 0.35
# #                 ha = 'center'
# #             if avg > 0:
# #                 wiggle = -0.4
# #                 ha = 'center'
# #             # plot
# #             offset = width * multiplier_deep1
# #             if attribute == 'Storage':
# #                 hatchcolor = 'white'
# #             else:
# #                 hatchcolor = 'white'
# #             if i == 0:
# #                 # rects = ax[0].bar(offset, avg, width, zorder=5, edgecolor=hatchcolor,color=color, hatch='xx')
# #                 # rects = ax[0].bar(offset, avg, width, zorder=5, edgecolor=color,color='none')
# #                 rects = ax[0].bar(pos + shift, avg, width, zorder=5, align='center', edgecolor=hatchcolor,color=color, hatch='xx')
# #                 rects = ax[0].bar(pos + shift, avg, width, zorder=5, align='center', edgecolor=color,color='none')
# #                 # ax[0].errorbar(pos+shift, avg, yerr=error, capsize=5, color='red',zorder=6)
# #                 if attribute == 'Storage':
# #                     # ax[0].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                     #         color=color,fontsize=12, fontweight='bold')
# #                     ax[0].text(pos + shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
# #                             color=color,fontsize=12, fontweight='bold', rotation=45)
# #                 else:
# #                     # ax[0].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                     #         color='gray',fontsize=12)
# #                     ax[0].text(pos + shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
# #                             color='gray',fontsize=12, rotation=45)
# #                 multiplier_deep1 += 2
# #             elif i == 1:
# #                 offset = width * multiplier_deep2
# #                 # rects = ax[0].bar(offset+width, avg, width, zorder=5, edgecolor=color,color=color,label=label)
# #                 rects = ax[0].bar(pos + shift, avg, width, zorder=5, align='center', edgecolor=color,color=color, label=label)
# #                 # ax[0].errorbar(pos+shift, avg, yerr=error, capsize=5, color='red',zorder=6)
# #                 if attribute == 'Storage':
# #                     # ax[0].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                     #         color=color,fontsize=12, fontweight='bold')
# #                     ax[0].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
# #                             color=color,fontsize=12, fontweight='bold', rotation=45)
# #                 else:
# #                     # ax[0].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                     #         color='gray',fontsize=12)
# #                     ax[0].text(pos + shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
# #                             color='gray',fontsize=12, rotation=45)
# #                 multiplier_deep2 += 2
            
# #             ax[0].legend(loc='lower left', fontsize=12, ncol=1)

# #     # part 2 with combined physical and biological terms ---------------------------------

# #     multiplier_deep1 = 0
# #     multiplier_deep2 = 0

# #     # create a new dictionary of results
# #     oxy_dict = {}
# #     hyp_dict = {}
# #     oxy_dict_std = {}
# #     hyp_dict_std = {}

# #     for station in sta_dict:
# #         for attribute, measurement in bottomlay_dict[station].items():
# #             # skip variables we are not interested in
# #             if attribute in ['EU Exchange Flow',
# #                              'TEF Exchange Flow',
# #                              'TRAPS',
# #                              'EU Recirculation',
# #                              'EU Vertical',
# #                              'TEF Vertical',
# #                              'Photosynthesis',
# #                              'Bio Consumption',
# #                              'Volume',
# #                              'Qin m3/s']:
# #                 continue
# #              # calculate time average normalized by volume
# #             avg = np.nanmean(measurement[minday:maxday]/(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
# #             # convert to mg/L per day
# #             avg = avg * 1000 * 32 * 60 * 60 * 24

# #             # get standard deviations
# #             Vavg = measurement[minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday] # kmol O2 /s /m3
# #             # calculate standard deviation in mg/L per day
# #             std = np.std(Vavg) * 1000 * 32 * 60 * 60 * 24

# #             # save values in dictionary
# #             if station in ['penn','case','holmes','portsusan','lynchcove','dabob']:
# #                 if attribute in hyp_dict.keys():
# #                     hyp_dict[attribute].append(avg)
# #                     hyp_dict_std[attribute].append(std)
# #                 else:
# #                     hyp_dict[attribute] = [avg]
# #                     hyp_dict_std[attribute] = [std]
# #             else:
# #                 if attribute in oxy_dict.keys():
# #                     oxy_dict[attribute].append(avg)
# #                     oxy_dict_std[attribute].append(std)
# #                 else:
# #                     oxy_dict[attribute] = [avg]
# #                     oxy_dict_std[attribute] = [std]

# #     # t-test
# #     print('COMBINED TERMS --------------------------------------')
# #     width = 0.2
# #     for attribute in oxy_dict:
# #         print('\n========================')
# #         print(attribute)
# #         a = oxy_dict[attribute]
# #         b = hyp_dict[attribute]
# #         ttest = ttest_ind(a, b, axis=0, equal_var=False)
# #         # if attribute == 'Photosynthesis & Consumption':
# #         #     print('one-sided')
# #         #     # one-sided t-test 
# #         #     # alternative hypothesis: a < b (drawdown of oxygenated is more negative than drawdown of hypoxic)
# #         #     # the null hypothesis is thus that drawdown of hypoxic is more negative
# #         #     ttest = ttest_ind(a, b, axis=0, equal_var=False, alternative='less')
# #         # else:
# #         #     # two-sided t-test with null hypothesis that they are the same
# #         #     ttest = ttest_ind(a, b, axis=0, equal_var=False)
# #         print(ttest)

# #     print('\n')


# #     for i,dict in enumerate([oxy_dict,hyp_dict]):
# #     # average all oxygenated and hypoxic inlet rate values, and calculate standard deviations
# #         if i ==0:
# #             shift = 0
# #         else:
# #             shift = 0.2
# #         for attribute, measurement in dict.items():
# #             # choose color
# #             if attribute == 'TEF Recirculation':
# #                 color = '#488DDB'
# #                 label = 'Exchange Flow & Vertical Transport'
# #                 pos = 0.2
# #             if attribute == 'Photosynthesis & Consumption':
# #                 color = '#F069A8'
# #                 label = attribute + ' (drawdown)'
# #                 pos = 0.7
# #             if attribute == 'Storage':
# #                 color = 'black'
# #                 label = r'$\frac{d}{dt}$DO (net decrease)'
# #                 pos = -0.3
        
# #             # calculate average and standard deviation
# #             avg = np.nanmean(measurement)
# #             # get error
# #             if dict == oxy_dict:
# #                 std_dict = oxy_dict_std
# #             elif dict == hyp_dict:
# #                 std_dict = hyp_dict_std
# #             # propagate error
# #             error = 1/len(std_dict.keys()) * np.sqrt(np.sum(np.square(std_dict[attribute])))

# #             if avg < 0:
# #                 wiggle = 0.02
# #             if avg > 0:
# #                 wiggle = -0.02
# #             # plot
# #             offset = width * multiplier_deep1
# #             if attribute == 'Storage':
# #                 hatchcolor = 'white'
# #             else:
# #                 hatchcolor = 'white'
# #             if i == 0:
# #                 # rects = ax[1].bar(offset, avg, width, zorder=5, edgecolor=hatchcolor,color=color, hatch='xx')
# #                 # rects = ax[1].bar(offset, avg, width, zorder=5, edgecolor=color,color='none')
# #                 rects = ax[1].bar(pos + shift, avg, width, zorder=5, edgecolor=hatchcolor,color=color, hatch='xx')
# #                 rects = ax[1].bar(pos + shift, avg, width, zorder=5, edgecolor=color,color='none')
# #                 # ax[1].errorbar(pos+shift, avg, yerr=error, capsize=5, color='red',zorder=6)
# #                 if attribute == 'Storage':
# #                     # ax[1].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                     #         color=color,fontsize=12, fontweight='bold')
# #                     ax[1].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                             color=color,fontsize=12, fontweight='bold')
# #                 else:
# #                     # ax[1].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                     #         color='gray',fontsize=12)
# #                     ax[1].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                             color='gray',fontsize=12)
                    
# #                 multiplier_deep1 += 2
# #             elif i == 1:
# #                 offset = width * multiplier_deep2
# #                 # rects = ax[1].bar(offset+width, avg, width, zorder=5, edgecolor=color,color=color,label=label)
# #                 rects = ax[1].bar(pos + shift, avg, width, zorder=5, edgecolor=color,color=color,label=label)
# #                 # ax[1].errorbar(pos+shift, avg, yerr=error, capsize=5, color='red',zorder=6)
# #                 if attribute == 'Storage':
# #                     # ax[1].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                     #         color=color,fontsize=12, fontweight='bold')
# #                     ax[1].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                             color=color,fontsize=12, fontweight='bold')
# #                 else:
# #                     # ax[1].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                     #         color='gray',fontsize=12)
# #                     ax[1].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
# #                             color='gray',fontsize=12)
# #                 multiplier_deep2 += 2
            
# #             ax[1].legend(loc='lower left', fontsize=12)

# #     ax[0].text(0.98, 0.1, 'HATCHED: oxygenated (n={})\nSOLID: hypoxic (n={})      '.format(len(oxy_dict['Storage']),len(hyp_dict['Storage'])),
# #             color='black', verticalalignment='bottom', horizontalalignment='right',zorder=6,
# #             transform=ax[0].transAxes, fontsize=12)
# #     plt.suptitle('Volume-averaged deep layer budgets (mid-Jul to mid-Aug)',fontsize=16) 

# #     ax[0].set_title('(a) All budget rate terms', loc='left',fontsize=14)
# #     ax[1].set_title('(b) Combined physical and biological processes', loc='left',fontsize=14)

# #     plt.subplots_adjust(left=0.1, wspace=0.02, top=0.9, bottom=0.1, right=0.9)
# #     plt.show()


# ##########################################################
# ##         BAR CHARTS (Lynch Cove and 13 inlets)        ##
# ##########################################################

# budget_barchart = True

# # COLLAPSE
# if budget_barchart == True:
#     plt.close('all')
            
#     fig, ax = plt.subplots(4,1,figsize=(9.8,9.5))

#     # lynch cove example budget
#     station = 'lynchcove'

#     # format figure
#     # plt.suptitle(station + ': DO Budget (10-day Hanning Window)',size=14)
#     # ax[0].set_title('(a) 2017 Lynch Cove deep oxygen budget (10-day Hanning Window)',size=12, loc='left')
#     # ax[0].set_title('(a) Lynch Cove oxygen budget',size=12, loc='left', fontweight='bold')
#     ax[0].text(0.02, 0.88,'(a) Lynch Cove',fontsize=12, fontweight='bold',transform=ax[0].transAxes,)
#     ax[0].set_xlim([dates_local[0],dates_local[-25]])
#     ax[0].set_ylabel('DO transport ' + r'[kmol O$_2$ s$^{-1}$]',size=10)
#     # ax[0].set_facecolor('#EEEEEE')
#     # ax[0].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#     ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#     ax[0].tick_params(axis='x', labelrotation=30, labelsize=10)
#     ax[0].tick_params(axis='y', labelsize=10)
#     loc = mdates.MonthLocator(interval=1)
#     ax[0].xaxis.set_major_locator(loc)
#     ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
#     # ax[0].tick_params(axis='both', labelsize=12)
#     # for border in ['top','right','bottom','left']:
#     #     ax[0].spines[border].set_visible(False)
#     # ax[0].set_xlabel(year, fontsize=14)
#     ax[0].set_ylim([-0.23,0.18])
#     # plot deep budget
#     nwin = 10 # hanning window length
#     ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['Storage'].values,n=nwin),color='k',
#                 linewidth=1,label=r'$\frac{d}{dt}\int_V$DO dV') #,alpha=0.6)
#     ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TEF Vertical'].values + surfacelay_dict[station]['TEF Vertical'].values,n=nwin),
#                 color='darkorange', linewidth=1,label='Error')
#     ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TEF Exchange Flow'].values,n=nwin),color='#0D4B91',
#             linewidth=2,label='Exchange Flow')
#     ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['TEF Vertical'].values,n=nwin),color='#99C5F7',
#             linewidth=2,label='Vertical')
#     # ax.plot(dates_local_daily,bottomlay_dict[station]['TRAPS'],color=traps_color,
#     #            linewidth=3,zorder=2)
#     ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['Photosynthesis'].values,n=nwin),color='#8F0445',
#                 linewidth=2, label='Photosynthesis')
#     ax[0].plot(dates_local_daily,zfun.lowpass(bottomlay_dict[station]['Bio Consumption'].values,n=nwin),color='#FCC2DD',
#                 linewidth=2,label='Consumption')
#     ax[0].legend(loc='lower right',ncol=6, fontsize=9, handletextpad=0.15)
#     # add drawdown period
#     # mid Jul - mid Aug
#     minday = 194
#     maxday = 225
#     ax[0].axvline(dates_local_daily[minday],0,12,color='grey')#,linestyle=':')
#     ax[0].axvline(dates_local_daily[maxday],0,12,color='grey')#,linestyle=':')

#     # Bar chart     
#     # ax[1].set_facecolor('#EEEEEE')
#     ax[1].tick_params(axis='x', labelrotation=30)
#     # ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='y')
#     ax[1].axhline(y=0, xmin=-0.5, xmax=1.05,color='silver',linewidth=1,linestyle='--')
#     # ax[1].grid(True,color='w',linewidth=1,linestyle='-',axis='y')
#     # for border in ['top','right','bottom','left']:
#     #     ax[1].spines[border].set_visible(False)
#     ax[1].tick_params(axis='y', labelsize=10)
#     ax[1].set_xticklabels([])
#     ax[1].set_ylabel('mg/L per day',fontsize=10)
#     # ax[1].set_title(station + ' volume-averaged deep layer DO budget')
#     # ax[1].set_title('(b) Lynch Cove volume-averaged deep layer oxygen budget (drawdown period)',size=12, loc='left')
#     ax[1].text(0.02, 0.88,'(b) Lynch Cove',fontsize=12, fontweight='bold',transform=ax[1].transAxes,)
#     ax[1].set_xlim([-0.5,1.05])
#     ax[1].set_ylim([-0.25,0.25])
#     width = 0.2
#     # create bar chart
#     for attribute, measurement in bottomlay_dict[station].items():
#         # skip variables we are not interested in
#         if attribute in ['EU Exchange Flow',
#                             'TRAPS',
#                             'EU Recirculation',
#                             'EU Vertical',
#                             'TEF Recirculation',
#                             'Photosynthesis & Consumption',
#                             'Volume',
#                             'Qin m3/s']:
#             continue
#         # assign colors
#         if attribute == 'TEF Exchange Flow':
#             color = '#0D4B91'
#             label = 'Exchange Flow'
#             pos = 0.2
#         if attribute == 'TEF Vertical':
#             color = '#99C5F7'
#             label = 'Vertical'
#             pos = 0.4
#         if attribute == 'Photosynthesis':
#             color = '#8F0445'
#             label = attribute
#             pos = 0.7
#         if attribute == 'Bio Consumption':
#             color = '#FCC2DD'
#             label = 'Consumption'
#             pos = 0.9
#         if attribute == 'Storage':
#             color = 'black'
#             label = r'$\frac{d}{dt}$DO (net decrease)'
#             pos = -0.2

#         # calculate time average
#         time_avg = np.nanmean(measurement[minday:maxday])
#         # get volume average
#         avg = time_avg/(np.nanmean(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
#         # convert to mg/L per day
#         avg = avg * 1000 * 32 * 60 * 60 * 24

#         # calculate standard deviation
#         # get volume average
#         Vavg = measurement[minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday] # kmol O2 /s /m3
#         # get standard deviation
#         std = np.std(Vavg)
#         std = std * 1000 * 32 * 60 * 60 * 24

#         ax[1].bar(pos, avg, width, zorder=5, align='center', edgecolor=color,color=color, label=label)
#         # ax[1].errorbar(pos, avg, yerr=std, capsize=5, color='grey',zorder=6)
#         if avg < 0:
#             wiggle = 0.05
#         if avg > 0:
#             wiggle = -0.05
#         ax[1].text(pos, wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                 color='black',fontsize=10)
        
#         # # zero line
#         # ax[1].axhline(0,-0.3,1.05,color='dimgrey')
        
#         # ax[1].legend(loc='lower center', fontsize=9, ncol=5)
#         ax[1].legend(bbox_to_anchor=(0.5, -0.3), loc='lower center', fontsize=9, ncol=5, handletextpad=0.15)




# ##########################################################
# ##              Comparing deep DO budgets               ## 
# ##########################################################


#     # mid July to mid August
#     minday = 194
#     maxday = 225


#     width = 0.1

#     multiplier_deep1 = 0
#     multiplier_deep2 = 0

#     # format grid
#     for axis in [ax[2],ax[3]]:
#         # axis.set_facecolor('#EEEEEE')
#         axis.tick_params(axis='x', labelrotation=30)
#         # axis.grid(True,color='w',linewidth=1,linestyle='-',axis='y')
#         # for border in ['top','right','bottom','left']:
#         #     axis.spines[border].set_visible(False)
#         # axis.grid(True,color='silver',linewidth=1,linestyle='--',axis='y')
#         axis.axhline(y=0, xmin=-0.5, xmax=1.05,color='silver',linewidth=1,linestyle='--')
#         axis.tick_params(axis='y', labelsize=10)
#         axis.set_xticklabels([])
#         axis.set_ylabel('mg/L per day',fontsize=10)
#         axis.set_xlim([-0.5,1.05])
#     ax[2].set_ylim([-2.5,2.5])
#     ax[3].set_ylim([-0.35,0.25])

#     # create a new dictionary of results
#     oxy_dict = {}
#     hyp_dict = {}
#     oxy_dict_std = {}
#     hyp_dict_std = {}

#     # part 1 with distinct physical and biological terms
#     for station in sta_dict:
#         for attribute, measurement in bottomlay_dict[station].items():
#             # skip variables we are not interested in
#             if attribute in ['EU Exchange Flow',
#                              'TRAPS',
#                              'EU Recirculation',
#                              'EU Vertical',
#                              'TEF Recirculation',
#                              'Photosynthesis & Consumption',
#                              'Volume',
#                              'Qin m3/s']:
#                 continue
#             # calculate time average normalized by volume
#             avg = np.nanmean(measurement[minday:maxday]/(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
#             # convert to mg/L per day
#             avg = avg * 1000 * 32 * 60 * 60 * 24

#             # get standard deviations
#             Vavg = measurement[minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday] # kmol O2 /s /m3
#             # calculate standard deviation in mg/L per day
#             std = np.std(Vavg) * 1000 * 32 * 60 * 60 * 24

#             # save values in dictionary
#             if station in ['penn','case','holmes','portsusan','lynchcove','dabob']:
#                 if attribute in hyp_dict.keys():
#                     hyp_dict[attribute].append(avg)
#                     hyp_dict_std[attribute].append(std)
#                 else:
#                     hyp_dict[attribute] = [avg]
#                     hyp_dict_std[attribute] = [std]
#             else:
#                 if attribute in oxy_dict.keys():
#                     oxy_dict[attribute].append(avg)
#                     oxy_dict_std[attribute].append(std)
#                 else:
#                     oxy_dict[attribute] = [avg]
#                     oxy_dict_std[attribute] = [std]
#     # t-test
#     print('DISTINCT TERMS --------------------------------------')
#     for attribute in oxy_dict:
#         print('\n========================')
#         print(attribute)
#         a = oxy_dict[attribute]
#         b = hyp_dict[attribute]
#         ttest = ttest_ind(a, b, axis=0, equal_var=False)
#         print(ttest)
#     print('\n')
#     for i,dict in enumerate([oxy_dict,hyp_dict]):
#     # average all oxygenated and hypoxic inlet rate values, and calculate standard deviations
#         if i ==0:
#             shift = 0
#         else:
#             shift = 0.1
#         for attribute, measurement in dict.items():
#             # choose color
#             if attribute == 'TEF Exchange Flow':
#                 color = '#0D4B91'
#                 label = 'Exchange Flow'
#                 pos = 0.15
#             if attribute == 'TEF Vertical':
#                 color = '#99C5F7'
#                 label = 'Vertical'
#                 pos = 0.35
#             if attribute == 'Photosynthesis':
#                 color = '#8F0445'
#                 label = attribute
#                 pos = 0.65
#             if attribute == 'Bio Consumption':
#                 color = '#FCC2DD'
#                 label = 'Consumption'
#                 pos = 0.85
#             if attribute == 'Storage':
#                 color = 'black'
#                 label = r'$\frac{d}{dt}$DO (net decrease)'
#                 pos = -0.25

#             # calculate average and standard deviation
#             avg = np.nanmean(measurement)
#             # get error
#             if dict == oxy_dict:
#                 std_dict = oxy_dict_std
#             elif dict == hyp_dict:
#                 std_dict = hyp_dict_std
#             # propagate error
#             error = 1/len(std_dict.keys()) * np.sqrt(np.sum(np.square(std_dict[attribute])))

#             if avg < 0:
#                 wiggle = 0.65
#                 ha = 'center'
#             if avg > 0:
#                 wiggle = -0.7
#                 ha = 'center'
#             # plot
#             offset = width * multiplier_deep1
#             if attribute == 'Storage':
#                 hatchcolor = 'white'
#             else:
#                 hatchcolor = 'white'
#             if i == 0:
#                 # rects = ax[0].bar(offset, avg, width, zorder=5, edgecolor=hatchcolor,color=color, hatch='xx')
#                 # rects = ax[0].bar(offset, avg, width, zorder=5, edgecolor=color,color='none')
#                 rects = ax[2].bar(pos + shift, avg, width, zorder=5, align='center', edgecolor=hatchcolor,color=color, hatch='xx')
#                 rects = ax[2].bar(pos + shift, avg, width, zorder=5, align='center', edgecolor=color,color='none')
#                 # ax[0].errorbar(pos+shift, avg, yerr=error, capsize=5, color='red',zorder=6)
#                 if attribute == 'Storage':
#                     # ax[0].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                     #         color=color,fontsize=12, fontweight='bold')
#                     ax[2].text(pos + shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
#                             color=color,fontsize=10, fontweight='bold', rotation=45)
#                 else:
#                     # ax[0].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                     #         color='gray',fontsize=12)
#                     ax[2].text(pos + shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
#                             color='gray',fontsize=10, rotation=45)
#                 multiplier_deep1 += 2
#             elif i == 1:
#                 offset = width * multiplier_deep2
#                 # rects = ax[0].bar(offset+width, avg, width, zorder=5, edgecolor=color,color=color,label=label)
#                 rects = ax[2].bar(pos + shift, avg, width, zorder=5, align='center', edgecolor=color,color=color, label=label)
#                 # ax[0].errorbar(pos+shift, avg, yerr=error, capsize=5, color='red',zorder=6)
#                 if attribute == 'Storage':
#                     # ax[0].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                     #         color=color,fontsize=12, fontweight='bold')
#                     ax[2].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
#                             color=color,fontsize=10, fontweight='bold', rotation=45)
#                 else:
#                     # ax[0].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                     #         color='gray',fontsize=12)
#                     ax[2].text(pos + shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
#                             color='gray',fontsize=10, rotation=45)
#                 multiplier_deep2 += 2
            
#             # ax[2].legend(loc='lower center', fontsize=9, ncol=5)

#     # part 2 with combined physical and biological terms ---------------------------------

#     multiplier_deep1 = 0
#     multiplier_deep2 = 0

#     # create a new dictionary of results
#     oxy_dict = {}
#     hyp_dict = {}
#     oxy_dict_std = {}
#     hyp_dict_std = {}

#     for station in sta_dict:
#         for attribute, measurement in bottomlay_dict[station].items():
#             # skip variables we are not interested in
#             if attribute in ['EU Exchange Flow',
#                              'TEF Exchange Flow',
#                              'TRAPS',
#                              'EU Recirculation',
#                              'EU Vertical',
#                              'TEF Vertical',
#                              'Photosynthesis',
#                              'Bio Consumption',
#                              'Volume',
#                              'Qin m3/s']:
#                 continue
#              # calculate time average normalized by volume
#             avg = np.nanmean(measurement[minday:maxday]/(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
#             # convert to mg/L per day
#             avg = avg * 1000 * 32 * 60 * 60 * 24

#             # get standard deviations
#             Vavg = measurement[minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday] # kmol O2 /s /m3
#             # calculate standard deviation in mg/L per day
#             std = np.std(Vavg) * 1000 * 32 * 60 * 60 * 24

#             # save values in dictionary
#             if station in ['penn','case','holmes','portsusan','lynchcove','dabob']:
#                 if attribute in hyp_dict.keys():
#                     hyp_dict[attribute].append(avg)
#                     hyp_dict_std[attribute].append(std)
#                 else:
#                     hyp_dict[attribute] = [avg]
#                     hyp_dict_std[attribute] = [std]
#             else:
#                 if attribute in oxy_dict.keys():
#                     oxy_dict[attribute].append(avg)
#                     oxy_dict_std[attribute].append(std)
#                 else:
#                     oxy_dict[attribute] = [avg]
#                     oxy_dict_std[attribute] = [std]

#     # t-test
#     print('COMBINED TERMS --------------------------------------')
#     width = 0.2
#     for attribute in oxy_dict:
#         print('\n========================')
#         print(attribute)
#         a = oxy_dict[attribute]
#         b = hyp_dict[attribute]
#         ttest = ttest_ind(a, b, axis=0, equal_var=False)
#         # if attribute == 'Photosynthesis & Consumption':
#         #     print('one-sided')
#         #     # one-sided t-test 
#         #     # alternative hypothesis: a < b (drawdown of oxygenated is more negative than drawdown of hypoxic)
#         #     # the null hypothesis is thus that drawdown of hypoxic is more negative
#         #     ttest = ttest_ind(a, b, axis=0, equal_var=False, alternative='less')
#         # else:
#         #     # two-sided t-test with null hypothesis that they are the same
#         #     ttest = ttest_ind(a, b, axis=0, equal_var=False)
#         print(ttest)

#     print('\n')


#     for i,dict in enumerate([oxy_dict,hyp_dict]):
#     # average all oxygenated and hypoxic inlet rate values, and calculate standard deviations
#         if i ==0:
#             shift = 0
#         else:
#             shift = 0.2
#         for attribute, measurement in dict.items():
#             # choose color
#             if attribute == 'TEF Recirculation':
#                 color = '#488DDB'
#                 label = 'Exchange Flow & Vertical Transport'
#                 pos = 0.2
#             if attribute == 'Photosynthesis & Consumption':
#                 color = '#F069A8'
#                 label = attribute + ' (drawdown)'
#                 pos = 0.7
#             if attribute == 'Storage':
#                 color = 'black'
#                 label = r'$\frac{d}{dt}$DO'
#                 pos = -0.3
        
#             # calculate average and standard deviation
#             avg = np.nanmean(measurement)
#             # get error
#             if dict == oxy_dict:
#                 std_dict = oxy_dict_std
#             elif dict == hyp_dict:
#                 std_dict = hyp_dict_std
#             # propagate error
#             error = 1/len(std_dict.keys()) * np.sqrt(np.sum(np.square(std_dict[attribute])))

#             if avg < 0:
#                 wiggle = 0.04
#             if avg > 0:
#                 wiggle = -0.04
#             # plot
#             offset = width * multiplier_deep1
#             if attribute == 'Storage':
#                 hatchcolor = 'white'
#             else:
#                 hatchcolor = 'white'
#             if i == 0:
#                 # rects = ax[1].bar(offset, avg, width, zorder=5, edgecolor=hatchcolor,color=color, hatch='xx')
#                 # rects = ax[1].bar(offset, avg, width, zorder=5, edgecolor=color,color='none')
#                 rects = ax[3].bar(pos + shift, avg, width, zorder=5, edgecolor=hatchcolor,color=color, hatch='xx')
#                 rects = ax[3].bar(pos + shift, avg, width, zorder=5, edgecolor=color,color='none')
#                 # ax[1].errorbar(pos+shift, avg, yerr=error, capsize=5, color='red',zorder=6)
#                 if attribute == 'Storage':
#                     # ax[1].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                     #         color=color,fontsize=12, fontweight='bold')
#                     ax[3].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                             color=color,fontsize=10, fontweight='bold')
#                 else:
#                     # ax[1].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                     #         color='gray',fontsize=12)
#                     ax[3].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                             color='gray',fontsize=10)
                    
#                 multiplier_deep1 += 2
#             elif i == 1:
#                 offset = width * multiplier_deep2
#                 # rects = ax[1].bar(offset+width, avg, width, zorder=5, edgecolor=color,color=color,label=label)
#                 rects = ax[3].bar(pos + shift, avg, width, zorder=5, edgecolor=color,color=color,label=label)
#                 # ax[1].errorbar(pos+shift, avg, yerr=error, capsize=5, color='red',zorder=6)
#                 if attribute == 'Storage':
#                     # ax[1].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                     #         color=color,fontsize=12, fontweight='bold')
#                     ax[3].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                             color=color,fontsize=10, fontweight='bold')
#                 else:
#                     # ax[1].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                     #         color='gray',fontsize=12)
#                     ax[3].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
#                             color='gray',fontsize=10)
#                 multiplier_deep2 += 2
            
#             ax[3].legend(loc='lower center', fontsize=9, ncol=3, handletextpad=0.15)

#     ax[2].text(0.3, 0.12, 'HATCHED: oxygenated (n={})\nSOLID: hypoxic (n={})      '.format(len(oxy_dict['Storage']),len(hyp_dict['Storage'])),
#             color='black', verticalalignment='bottom', horizontalalignment='right',zorder=6,
#             transform=ax[2].transAxes, fontsize=9, fontweight='bold')
#     # plt.suptitle('Volume-averaged deep layer budgets (mid-Jul to mid-Aug)',fontsize=16) 

#     # ax[2].set_title('(c) All Inlet volume-averaged deep layer oxygen budgets (drawdown period)', loc='left',fontsize=12)
#     # ax[3].set_title('(d) All Inlet oxygen budgets; combined physical and biological terms', loc='left',fontsize=12)
#     ax[2].text(0.02, 0.88,'(c) All inlets',fontsize=12, fontweight='bold',transform=ax[2].transAxes,)
#     ax[3].text(0.02, 0.88,'(d) All inlets',fontsize=12, fontweight='bold',transform=ax[3].transAxes,)

#     plt.subplots_adjust(left=0.1, top=0.95, bottom=0.05, right=0.9, hspace=0.3)
#     # plt.tight_layout()
#     plt.show()


# #########################################################
# ##     biological vs. physical mid-jul to mid-aug      ##
# #########################################################

# # mid July to mid August
# minday = 194
# maxday = 225

# # initialize figure
# fig, ax = plt.subplots(1,1,figsize = (5,5))

# # format figure
# ax.set_title('Consumption vs. QinDOin\n(mid-Jul through mid-Aug)',size=14)
# ax.tick_params(axis='x', labelrotation=30)
# ax.grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
# ax.tick_params(axis='both', labelsize=12)
# ax.set_xlabel('Exchange Flow\n[mg/L per day]', fontsize=12)
# ax.set_ylabel('Bio Consumption\n[mg/L per day]', fontsize=12)

# # # add line with slope -1
# # ax.plot([0,0.45], [0,-0.45], color='grey', linestyle='-')

# physics_all = np.array([])
# biology_all = np.array([])

# for i,station in enumerate(sta_dict):
    
#     # get physical
#     # calculate time average normalized by volume
#     tef =  np.nanmean(bottomlay_dict[station]['TEF Exchange Flow'][minday:maxday]/
#                       (bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
#     vert =  np.nanmean(bottomlay_dict[station]['TEF Vertical'][minday:maxday]/
#                        (bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
#     # convert to mg/L per day
#     tef = tef * 1000 * 32 * 60 * 60 * 24
#     vert = vert * 1000 * 32 * 60 * 60 * 24
#     # sum
#     physics = tef + vert
#     physics = tef

#     # get biological
#     # calculate time average normalized by volume
#     photo =  np.nanmean(bottomlay_dict[station]['Photosynthesis'][minday:maxday]/
#                       (bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
#     cons =  np.nanmean(bottomlay_dict[station]['Bio Consumption'][minday:maxday]/
#                        (bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
#     # convert to mg/L per day
#     photo = photo * 1000 * 32 * 60 * 60 * 24
#     cons = cons * 1000 * 32 * 60 * 60 * 24
#     # sum
#     biology = photo + cons
#     biology = cons

#     # plot
#     ax.scatter(physics,biology,color='navy',alpha=0.3, s=60, zorder=5)
#     ax.set_ylim([-0.6,0])
#     ax.set_xlim([0,6])

#     # add to array
#     physics_all = np.concatenate((physics_all,[physics]))
#     biology_all = np.concatenate((biology_all,[biology]))

# # calculate correlation coefficient (Pearson)
# r,p = pearsonr(physics_all,biology_all)
# ax.text(0.4, 0.85, r'$r =$' + str(round(r,3)) + r'; $r^2 =$' + str(round(r**2,3)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=ax.transAxes, fontsize=12, fontweight='bold')
# ax.text(0.4, 0.79, r'$p =$' + '{:.1e}'.format(p) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=ax.transAxes, fontsize=12, fontweight='bold')
# # ax.text(0.4, 0.79, r'$p =$' + str(round(p,3)) ,color='black',
# #                         verticalalignment='bottom', horizontalalignment='left',
# #                         transform=ax.transAxes, fontsize=12, fontweight='bold')

# # ax.set_aspect('equal', adjustable='box')
# ax.xaxis.set_major_locator(plt.MaxNLocator(5))
# ax.yaxis.set_major_locator(plt.MaxNLocator(5))
# plt.tight_layout()
# plt.show()

# ##########################################################
# ##                   DO scatterplots                    ## 
# ##########################################################

# DO_analysis = True
# if DO_analysis == True:

#     # initialize arrays for plotting
#     intervals = 12
#     deep_lay_DO = np.zeros(len(sta_dict)*intervals)
#     bott_sig_DO = np.zeros(len(sta_dict)*intervals)
#     min_bott_sig_DO = np.zeros(len(sta_dict)*intervals)
#     annual_mean_DO = np.zeros(len(sta_dict)*intervals)
#     perc_hyp_vol = np.zeros(len(sta_dict)*intervals)
#     mean_DOin = np.zeros(len(sta_dict)*intervals)
#     mean_Tflush = np.zeros(len(sta_dict)*intervals)
#     mean_TEFin = np.zeros(len(sta_dict)*intervals)
#     mean_recirc = np.zeros(len(sta_dict)*intervals)
#     mean_cons = np.zeros(len(sta_dict)*intervals)
#     mean_depth = np.zeros(len(sta_dict)*intervals)
#     inlet_vol = np.zeros(len(sta_dict)*intervals)
#     aspect_ratio = np.zeros(len(sta_dict)*intervals)

#     annmean_DOin = np.zeros(len(sta_dict))
#     annmin_DOin = np.zeros(len(sta_dict))
#     colors_twentyone = []

#     colors = []

#     # get values for plotting and calculating r value
#     for i,station in enumerate(sta_dict):
#         # get interface depth from csv file
#         with open('interface_depths.csv', 'r') as f:
#             for line in f:
#                 inlet, interface_depth = line.strip().split(',')
#                 interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
#         z_interface = float(interface_dict[station])
#         fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_2014.01.01_2014.12.31.nc'
#         ds = xr.open_dataset(fn)
#         moor_depth = ds.h.values
        
#         for month in range(intervals):
#             if month == 0:
#                 minday = 0 #1
#                 maxday = 30 #32
#             elif month == 1:
#                 minday = 30# 32
#                 maxday = 58# 60
#             elif month == 2:
#                 minday = 58# 60
#                 maxday = 89# 91
#             elif month == 3:
#                 minday = 89# 91
#                 maxday = 119# 121
#             elif month == 4:
#                 minday = 119# 121
#                 maxday = 150# 152
#             elif month == 5:
#                 minday = 150# 152
#                 maxday = 180# 182
#             elif month == 6:
#                 minday = 180# 182
#                 maxday = 211# 213
#             elif month == 7:
#                 minday = 211# 213
#                 maxday = 242# 244
#             elif month == 8:
#                 minday = 242# 244
#                 maxday = 272# 274
#             elif month == 9:
#                 minday = 272# 274
#                 maxday = 303# 305
#             elif month == 10:
#                 minday = 303# 305
#                 maxday = 332# 335
#             elif month == 11:
#                 minday = 332# 335
#                 maxday = 363
#             # minday=month
#             # maxday=month+1
#             # save values
#             deep_lay_DO[i*intervals+month] =  np.nanmean(DOconcen_dict[station]['Deep Layer'][minday:maxday])
#             bott_sig_DO[i*intervals+month] =  np.nanmean(DOconcen_dict[station]['Bottom Sigma DO'][minday:maxday])
#             min_bott_sig_DO[i*intervals+month] =  np.nanmin(DOconcen_dict[station]['Minimum Bottom Layer DO'][minday:maxday])
#             annual_mean_DO[i*intervals+month] = np.nanmean(DOconcen_dict[station]['Deep Layer'])
#             perc_hyp_vol[i*intervals+month] = np.nanmean(DOconcen_dict[station]['percent hypoxic volume'][minday:maxday])
#             mean_DOin[i*intervals+month] = np.nanmean(DOconcen_dict[station]['Qin DO'][minday:maxday])
#             mean_Tflush[i*intervals+month] = np.nanmean(dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'][minday:maxday]) / (60*60*24)
#             mean_TEFin[i*intervals+month] = np.nanmean(bottomlay_dict[station]['TEF Exchange Flow'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]) * (
#                                     32 * 1000) * (60*60*24)
#             mean_recirc[i*intervals+month] = np.nanmean(bottomlay_dict[station]['TEF Recirculation'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]) * (
#                                     32 * 1000) * (60*60*24)
#             mean_cons[i*intervals+month] = np.nanmean(bottomlay_dict[station]['Bio Consumption'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]) * (
#                                     32 * 1000) * (60*60*24)
#             mean_depth[i*intervals+month] = dimensions_dict[station]['Mean depth'][0]
#             inlet_vol[i*intervals+month] = dimensions_dict[station]['Inlet volume'][0]
#             aspect_ratio[i*intervals+month] = dimensions_dict[station]['L/W aspect ratio'][0]
#             colors.append(basin_color_dict[basin_dict[station]])

#         annmean_DOin[i] = np.nanmean(DOconcen_dict[station]['Qin DO'])
#         annmin_DOin[i] = np.nanmin(DOconcen_dict[station]['Qin DO'])
#         colors_twentyone.append(basin_color_dict[basin_dict[station]])
        


#     # DOin vs. Tflush colored by percent hypoxic volume ============== 4PART MONEY PLOT SCATTER ================================
#     percent_hypoxic = True
#     if percent_hypoxic == True:
#         # initialize figure
#         fig, ax = plt.subplots(2,2,figsize = (10,9))
#         ax = ax.ravel()

#         # format figure
#         # ax[0].set_title('(a) ' + year + r' monthly mean DO$_{deep}$ vs. DO$_{in}$' + '\n' + r'colored by mean T$_{flush}$',
#         #                 size=14, loc='left')
#         ax[0].set_title('(a) All Inlets', size=14, loc='left', fontweight='bold')
#         # ax[0].set_title(year + r' monthly mean DO$_{deep}$ vs. DO$_{in}$',
#         #                 size=16, loc='left')
#         # format grid
#         # ax[0].set_facecolor('#EEEEEE')
#         ax[0].tick_params(axis='x', labelrotation=30)
#         # ax[0].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#         # for border in ['top','right','bottom','left']:
#         #     ax[0].spines[border].set_visible(False)
#         ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#         ax[0].tick_params(axis='both', labelsize=12)
#         ax[0].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
#         ax[0].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
#         # plot
#         cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
#         cmap_oxy = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
#         ax[0].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, c='k', alpha=0.5)
#         ax[0].plot([0,12],[0,12],color='gray')
#         # ax[0].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4, fontsize=10)
#         ax[0].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
#         cs = ax[0].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
#         # create colorbarlegend
#         cbar = fig.colorbar(cs)
#         cbar.ax.tick_params(labelsize=12)
#         cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
#         cbar.outline.set_visible(False)
#         ax[0].set_xlim([0,10])
#         ax[0].set_ylim([0,10])

#         # format figure
#         # ax[1].set_title('(b) ' + year + ' monthly mean '+r'DO$_{in}$ vs. T$_{flush}$'+'\ncolored by % hypoxic volume',
#         #                 size=14, loc='left')
#         ax[1].set_title('(b) All Inlets', size=14, loc='left', fontweight='bold')
#         # format grid
#         # ax[1].set_facecolor('#EEEEEE')
#         ax[1].tick_params(axis='x', labelrotation=30)
#         # ax[1].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#         # for border in ['top','right','bottom','left']:
#         #     ax[1].spines[border].set_visible(False)
#         ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#         ax[1].tick_params(axis='both', labelsize=12)
#         ax[1].set_xlabel(r'Monthly mean T$_{flush}$ [days]', fontsize=12)
#         ax[1].set_ylabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
#         ax[1].set_ylim([0,10])
#         ax[1].set_xlim([0,85])
#         # plot
#         cmap_hyp = plt.cm.get_cmap('gist_heat_r')
#         cs_DO = ax[1].scatter(mean_Tflush,mean_DOin,s=60,zorder=5,edgecolor='gray',c=perc_hyp_vol,cmap=cmap_hyp)
#         # create colorbarlegend
#         cbar = fig.colorbar(cs_DO)
#         cbar.ax.tick_params(labelsize=12)
#         cbar.ax.set_ylabel('Monthly mean % hypoxic volume', rotation=90, fontsize=12)
#         cbar.outline.set_visible(False)

#         # crescent bay
#         # ax[2].set_title('(c) Crescent Bay 2017 monthly mean \n' + r'DO$_{deep}$ vs. DO$_{in}$ colored by T$_{flush}$', loc='left', size=14)
#         ax[2].set_title('(c) Crescent Bay', size=14, loc='left', fontweight='bold')
#         # format grid
#         # ax[2].set_facecolor('#EEEEEE')
#         ax[2].tick_params(axis='x', labelrotation=30)
#         # ax[2].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#         # for border in ['top','right','bottom','left']:
#         #     ax[2].spines[border].set_visible(False)
#         ax[2].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#         ax[2].tick_params(axis='both', labelsize=12)
#         ax[2].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
#         ax[2].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
#         # plot
#         cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
#         cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
#         ax[2].plot([0,11],[0,11],color='dimgray')
#         ax[2].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
#         # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
#         ax[2].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, color='gray',alpha=0.5, edgecolor='none')
#         for i,station in enumerate(sta_dict):
#             if station == 'crescent':
#                 cs = ax[2].scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
#                                 s=150, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_tflush,
#                             linewidth=2, vmin=0, vmax=40)
#             else:
#                 continue
#         # create colorbarlegend
#         cbar = fig.colorbar(cs)
#         cbar.ax.tick_params(labelsize=12)
#         cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
#         cbar.outline.set_visible(False)
#         ax[2].set_xlim([0,11])
#         ax[2].set_ylim([0,11])

#         # lynch cove
#         # ax[3].set_title('(d) Lynch Cove 2017 monthly mean \n' + r'DO$_{deep}$ vs. DO$_{in}$ colored by T$_{flush}$', loc='left', size=14)
#         ax[3].set_title('(d) Lynch Cove', size=14, loc='left', fontweight='bold')
#         # format grid
#         # ax[3].set_facecolor('#EEEEEE')
#         ax[3].tick_params(axis='x', labelrotation=30)
#         # ax[3].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#         # for border in ['top','right','bottom','left']:
#         #     ax[3].spines[border].set_visible(False)
#         ax[3].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#         ax[3].tick_params(axis='both', labelsize=12)
#         ax[3].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
#         ax[3].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
#         # plot
#         cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
#         cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
#         ax[3].plot([0,11],[0,11],color='dimgray')
#         ax[3].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
#         # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
#         ax[3].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, color='gray',alpha=0.5, edgecolor='none')
#         for i,station in enumerate(sta_dict):
#             if station == 'lynchcove':
#                 cs = ax[3].scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
#                                 s=150, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_tflush,
#                             linewidth=2, vmin=0, vmax=40)
#             else:
#                 continue
#         # create colorbarlegend
#         cbar = fig.colorbar(cs)
#         cbar.ax.tick_params(labelsize=12)
#         cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
#         cbar.outline.set_visible(False)
#         ax[3].set_xlim([0,11])
#         ax[3].set_ylim([0,11])


#         plt.tight_layout()
#         # save figure
#         plt.show()

# ##########################################################
# ## MEAN DEEP DO vs % HYPOXIC VOLUME and  Deep DO time series ## 
# ##########################################################

# # mid Jul - mid Aug
# minday = 194
# maxday = 225

# # initialize figure
# fig, ax = plt.subplots(1,2,figsize = (12,5),gridspec_kw={'width_ratios': [1, 1.5]})


# # Deep DO vs. % hypoxic volume
# # ax[0].set_title('(a) ' + year + ' monthly mean \n     ' + r'% hypoxic volume vs. DO$_{deep}$',
# #                 size=14, loc='left')
# ax[0].set_title('(a) All inlets', size=14, loc='left', fontweight='bold')
# # format grid
# # ax[0].set_facecolor('#EEEEEE')
# ax[0].tick_params(axis='x', labelrotation=30)
# ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
# # ax[0].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
# # for border in ['top','right','bottom','left']:
# #     ax[0].spines[border].set_visible(False)
# ax[0].tick_params(axis='both', labelsize=12)
# ax[0].set_xlabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=14)
# ax[0].set_ylabel('Monthly mean % hypoxic volume', fontsize=14)
# # plot
# ax[0].scatter(deep_lay_DO,perc_hyp_vol,alpha=0.3,s=80,zorder=5,color='navy')
#         # color=colors)
# ax[0].set_xlim([0,10])
# ax[0].set_ylim([0,100])


# # Deep DO timeseries
# # ax[1].set_title('(b) ' + year + r' DO$_{deep}$ time series [mg/L]' +  '\n     (30-day Hanning Window)',
# #                 size=14, loc='left')
# ax[1].set_title('(b) All inlets', size=14, loc='left', fontweight='bold')
# # format grid
# # ax[1].set_facecolor('#EEEEEE')
# ax[1].tick_params(axis='x', labelrotation=30)
# loc = mdates.MonthLocator(interval=1)
# ax[1].xaxis.set_major_locator(loc)
# ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
# # ax[1].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
# # for border in ['top','right','bottom','left']:
# #     ax[1].spines[border].set_visible(False)
# ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
# ax[1].tick_params(axis='both', labelsize=12)
# # add drawdown period
# # ax[1].axvline(dates_local_daily[minday],0,12,color='pink')
# # ax[1].axvline(dates_local_daily[maxday],0,12,color='pink')
# ax[1].axvline(dates_local_daily[minday],0,12,color='grey')#,linestyle=':')
# ax[1].axvline(dates_local_daily[maxday],0,12,color='grey')#,linestyle=':')
# # loop through stations
# for i,station in enumerate(sta_dict):
#     # get average deep layer DO
#     deep_lay_DO_alltime = DOconcen_dict[station]['Deep Layer']
#     # 30-day hanning windwo
#     deep_lay_DO_alltime = zfun.lowpass(deep_lay_DO_alltime.values,n=30)
#     ax[1].plot(dates_local_daily,deep_lay_DO_alltime,linewidth=1,color='navy',alpha=0.5)

# # format labels
# ax[1].set_xlim([dates_local[0],dates_local[-2]])
# ax[1].set_ylim([0,10])
# ax[1].set_ylabel(r'DO$_{deep}$ [mg/L]',fontsize=14)
# plt.tight_layout()
# plt.show()

# ##########################################################
# ## net decrease (mid-July to mid-August) boxplots ## 
# ##########################################################

# # mid July to mid August
# minday = 194
# maxday = 225

# # initialize figure
# fig, ax = plt.subplots(1,1,figsize = (10,5))

# # format figure
# ax.set_title('d/dt(DO) (mid-Jul through mid-Aug)',size=14)
# ax.tick_params(axis='x', labelrotation=30)
# ax.grid(True,color='silver',linewidth=1,linestyle='--',axis='x')
# ax.tick_params(axis='both', labelsize=12)
# ax.set_ylabel('d/dt(DO) [mg/L per day]', fontsize=12)

# # # add line with slope -1
# # ax.plot([0,0.45], [0,-0.45], color='grey', linestyle='-')

# storage_all = []

# for i,station in enumerate(sta_dict):
    
#     # get daily net decrease rate
#     storage_daily =  bottomlay_dict[station]['Storage'][minday:maxday]/(bottomlay_dict[station]['Volume'][minday:maxday]) # kmol O2 /s /m3
    
#     # convert to mg/L per day
#     storage_daily = storage_daily.values * 1000 * 32 * 60 * 60 * 24

#     # add to array
#     storage_all.append(list(storage_daily))

# # create boxplot
# ax.axhline(y=0, xmin=-0.5, xmax=1.05,color='silver',linewidth=1,linestyle='--')
# bplot = plt.boxplot(storage_all, patch_artist=True, labels=sta_dict.keys(),
#             showmeans=True)

# for patch in bplot['boxes']:
#     patch.set_facecolor('honeydew')

# # condudct anova test
# # anova0 = f_oneway(storage_all[0],storage_all[1],storage_all[2],storage_all[3],storage_all[4],
# #          storage_all[5],storage_all[6],storage_all[7],storage_all[8],storage_all[9],
# #          storage_all[10],storage_all[11],storage_all[12])
# # anova1 = kruskal(storage_all[0],storage_all[1],storage_all[2],storage_all[3],storage_all[4],
# #          storage_all[5],storage_all[6],storage_all[7],storage_all[8],storage_all[9],
# #          storage_all[10],storage_all[11],storage_all[12])
# # anova2 = alexandergovern(storage_all[0],storage_all[1],storage_all[2],storage_all[3],storage_all[4],
# #          storage_all[5],storage_all[6],storage_all[7],storage_all[8],storage_all[9],
# #          storage_all[10],storage_all[11],storage_all[12])
# flat_storage = [x for xs in storage_all for x in xs]
# df = pd.DataFrame({'storage':flat_storage,
#                    'inlet': np.repeat(list(sta_dict.keys()), repeats=31)})
# # print(df)
# pingu = pg.welch_anova(data=df, dv='storage', between='inlet')


# print('\n===============================\n')
# # print(anova0)
# # print(anova1)
# # print(anova2)
# print(pingu)

# plt.tight_layout()
# plt.show()

# ##########################################################
# ##               Multiple regression (1)                ## 
# ##########################################################

# # create array of predictors

# input_array = np.array([mean_DOin, mean_Tflush, [1]*len(mean_DOin)]).T


# B,a,b,c = lstsq(input_array,deep_lay_DO)
# slope_DOin = B[0]
# slope_Tflush = B[1]
# intercept = B[2]

# print('\nMean deep layer DO [mg/L] = {}*DOin + {}*Tflush + {}'.format(
#     round(slope_DOin,2),round(slope_Tflush,2),round(intercept,2)))


# # calculate r^2 and p value
# r,p = pearsonr(mean_DOin,deep_lay_DO)
# print('\n===============================================')
# print('DO_deep dependence on DO_in')
# print('   r = {}'.format(r))
# print('   R^2 = {}'.format(r**2))
# print('   p = {}'.format(p))
# print('===============================================\n')

# # calculate r^2 and p value
# predicted_DOdeep = slope_DOin * mean_DOin + slope_Tflush * mean_Tflush + intercept
# r,p = pearsonr(deep_lay_DO,predicted_DOdeep)
# print('\n===============================================')
# print('DO_deep dependence on DO_in and T_flush')
# print('   r = {}'.format(r))
# print('   R^2 = {}'.format(r**2))
# print('   p = {}'.format(p))
# print('===============================================\n')
