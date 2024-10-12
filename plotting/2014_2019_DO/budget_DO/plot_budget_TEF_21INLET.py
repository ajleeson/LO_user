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

# group = '5Inlet'
# stations = ['lynchcove','penn','case','carr','budd'] # 5 inlets

# group = 'Whidbey'
# stations = ['similk','oak','crescent','penn','portsusan','holmes']

# group = 'HoodCanal'
# stations = ['lynchcove','dabob'] # Hood Canal

# group = 'MainBasin'
# stations = ['dyes','sinclair','elliot','quartermaster','commencement']

# group = 'SouthSound'
# stations = ['case','carr','hammersley','totten','eld','budd','henderson']

# group = 'Admiralty Inlet'
# stations = ['killsut']

# group = 'Deep Inlets'
# stations = ['dabob','carr','portsusan','commencement','elliot']

# group = 'Problem5Inlets'
# stations = ['lynchcove','dabob','case','carr','holmes'] # 5 inlets

# group = '15mOrDeeper'
# stations = ['case','lynchcove','carr','holmes','portsusan','elliot','commencement','dabob','penn','crescent']

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
    # # remove shallow inlets (< 10 m deep)
    # del sta_dict['hammersley']
    # del sta_dict['henderson']
    # del sta_dict['oak']
    # del sta_dict['totten']
    # del sta_dict['similk']
    # del sta_dict['budd']
    # del sta_dict['eld']
    # del sta_dict['killsut']
else:
    sta_dict = stations

    # set_difference = set(job_lists.get_sta_dict(jobname)) - set(stations)
    # sta_dict = list(set_difference)
    
# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
Lfun.make_dir(out_dir)

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]
# crop time vector (because we only have jan 2 - dec 30)
dates_no_crop = dates_local_daily
dates_local_daily = dates_local_daily[1:-1]

# create color dictionary for stations
basin_dict = {'similk': 'Whidbey',
              'oak': 'Whidbey',
              'crescent': 'Whidbey',
              'penn': 'Whidbey',
              'portsusan': 'Whidbey',
              'holmes': 'Whidbey',
              'dabob': 'Hood Canal',
              'lynchcove': 'Hood Canal',
              'dyes': 'Main Basin',
              'sinclair': 'Main Basin',
              'elliot': 'Main Basin',
              'quartermaster': 'Main Basin',
              'commencement': 'Main Basin',
              'case': 'South Sound',
              'carr': 'South Sound',
              'hammersley': 'South Sound', 
              'totten': 'South Sound', 
              'eld': 'South Sound', 
              'budd': 'South Sound', 
              'henderson': 'South Sound',
              'killsut': 'Admiralty Inlet'}
basin_color_dict = {'Whidbey': 'limegreen',
                    'Hood Canal': 'hotpink',
                    'Main Basin': 'deepskyblue',
                    'South Sound': 'blueviolet',
                    'Admiralty Inlet': 'black'}

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
        surfacelay_dict[station]['EU Exchange Flow'] = EU_surf
        surfacelay_dict[station]['TEF Exchange Flow'] = TEF_surf
        surfacelay_dict[station]['EU Recirculation'] = EU_surf + vertX_surf_EU
        surfacelay_dict[station]['TEF Recirculation'] = TEF_surf + vertX_surf_TEF
        surfacelay_dict[station]['TRAPS'] = traps_surf
        surfacelay_dict[station]['Photosynthesis'] = photo_surf
        surfacelay_dict[station]['Bio Consumption'] = cons_surf
        surfacelay_dict[station]['Air-Sea Transfer'] = airsea_surf
        surfacelay_dict[station]['EU Vertical'] = vertX_surf_EU
        surfacelay_dict[station]['TEF Vertical'] = vertX_surf_TEF
        surfacelay_dict[station]['Storage'] = ddtDOV_surf
        surfacelay_dict[station]['Volume'] = surf_V

        # bottom layer
        bottomlay_dict[station]['EU Exchange Flow'] = EU_deep
        bottomlay_dict[station]['TEF Exchange Flow'] = TEF_deep
        bottomlay_dict[station]['EU Recirculation'] = EU_deep + vertX_deep_EU
        bottomlay_dict[station]['TEF Recirculation'] = TEF_deep + vertX_deep_TEF
        bottomlay_dict[station]['EU Vertical'] = vertX_deep_EU
        bottomlay_dict[station]['TEF Vertical'] = vertX_deep_TEF
        bottomlay_dict[station]['TRAPS'] = traps_deep
        bottomlay_dict[station]['Photosynthesis'] = photo_deep
        bottomlay_dict[station]['Bio Consumption'] = cons_deep
        bottomlay_dict[station]['Photosynthesis & Consumption'] = photo_deep + cons_deep
        bottomlay_dict[station]['Storage'] = ddtDOV_deep
        bottomlay_dict[station]['Volume'] = deep_V
        bottomlay_dict[station]['Qin m3/s'] = Q_p.values 

# ------------------------- save DO concentrations in dataframe dict -----------------------------------

        # mg/L units

        # DO concentrations
        DOconcen_dict[station]['Surface Layer'] = o2_surf
        DOconcen_dict[station]['Deep Layer'] = o2_deep
        DOconcen_dict[station]['Bottom Sigma DO'] = bottom_oxygen
        DOconcen_dict[station]['Minimum Bottom Layer DO'] = oxygen_min
        DOconcen_dict[station]['Qin DO'] = DO_in
        DOconcen_dict[station]['percent hypoxic volume'] = hyp_vol/[inlet_vol] * 100

# ------------------------ save inlet dimensions in dataframe dict ----------------------------------------

        dimensions_dict[station]['Inlet volume'] = [inlet_vol] # m3
        dimensions_dict[station]['Mean depth'] = [mean_depth] # m
        dimensions_dict[station]['Mouth width'] = mouth_width
        dimensions_dict[station]['L/W aspect ratio'] = [est_aspect_ratio] # dimensionless

        # print keys
        if i == 0:
            print(list(surfacelay_dict[station].keys()))
            print(list(DOconcen_dict[station].keys()))
            print(list(dimensions_dict[station].keys()))

##########################################################
##              Plot DO budget of every inlet           ##
##########################################################

DO_budget = False

residual = True # recirculation (combined horizontal and vertical exchange)
show_EU = True

# COLLAPSE
if DO_budget == True:

    print('Making DO budget time series')

    for i,station in enumerate(sta_dict):
            
            # get interface depth from csv file
            with open('interface_depths.csv', 'r') as f:
                for line in f:
                    inlet, interface_depth = line.strip().split(',')
                    interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
            z_interface = float(interface_dict[station])

            # initialize figure
            plt.close('all')
            fig, ax = plt.subplots(4,1,figsize = (12,9),sharex=True)
            # format figure
            plt.suptitle(station + ': DO Budget (Godin Filter)',size=14)
            for axis in [ax[0],ax[1],ax[2],ax[3]]:
                # axis.plot([dates_local[0],dates_local[-1]],[0,0],color='k')
                axis.set_xlim([dates_local[0],dates_local[-1]])
                if axis == ax[3]:
                    axis.set_ylabel(r'DO [mg L$^{-1}$]')
                else:
                    axis.set_ylabel('DO transport\n' + r'[kmol O$_2$ s$^{-1}$]')
                axis.set_facecolor('#EEEEEE')
                axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
                axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
                axis.tick_params(axis='x', labelrotation=30)
                for border in ['top','right','bottom','left']:
                    axis.spines[border].set_visible(False)
            ax[3].set_xlabel(year)
            ax[0].text(0.02,0.9,'(a) Surface [shallower than {} m]'.format(-1*z_interface),
                            ha='left', va='top', transform=ax[0].transAxes, fontsize=12)
            ax[1].text(0.02,0.9,'(b) Bottom [deeper than {} m]'.format(-1*z_interface),
                            ha='left', va='top', transform=ax[1].transAxes, fontsize=12)
            ax[2].text(0.02,0.9,r'(c) Error [sum of surface and bottom vertical exchange]',
                            ha='left', va='top', transform=ax[2].transAxes, fontsize=12)
            ax[3].text(0.02,0.9,'(d) Average DO concentrations',
                            ha='left', va='top', transform=ax[3].transAxes, fontsize=12)
        

            # plot surface
            # exchange flow
            if residual == False:
                ax[0].plot(dates_local_daily,surfacelay_dict[station]['TEF Exchange Flow'],color=exchange_color,
                        linewidth=1,label='TEF Exchange Flow')
                if show_EU:
                    ax[0].plot(dates_local_daily,surfacelay_dict[station]['EU Exchange Flow'],color=exchange_color,
                            linewidth=3,alpha=0.3,label='EU Exchange Flow')
            # recirculation
            else:
                ax[0].plot(dates_local_daily,surfacelay_dict[station]['TEF Recirculation'],color=exchange_color,
                        linewidth=2,label='Recirculation TEF')
                if show_EU:
                    ax[0].plot(dates_local_daily,surfacelay_dict[station]['EU Recirculation'],color=exchange_color,
                            linewidth=3,alpha=0.3,label='Recirculation EU')
            # rivers
            ax[0].plot(dates_local_daily,surfacelay_dict[station]['TRAPS'],color=traps_color,
                       linewidth=3,zorder=2,label='TRAPS')
            # photosynthesis
            ax[0].plot(dates_local_daily,surfacelay_dict[station]['Photosynthesis'],color=photo_color,
                       linewidth=2,label='Photosynthesis')
            # consumption
            ax[0].plot(dates_local_daily,surfacelay_dict[station]['Bio Consumption'],color=cons_color,
                       linewidth=2,linestyle=':',zorder=9,label='Bio Consumption')
            # air-sea gas exchange
            ax[0].plot(dates_local_daily,surfacelay_dict[station]['Air-Sea Transfer'],color=airsea_color,
                       linewidth=1,zorder=8,label='Air-Sea Transfer')
            # storage
            ax[0].plot(dates_local_daily,surfacelay_dict[station]['Storage'],color=ddtDOV_color,
                       linewidth=1,zorder=7,#alpha=0.6,
                       label=r'$\frac{\mathrm{d}}{\mathrm{dt}}(\mathrm{DO}\cdot V)$')
            # vertical exchange
            if residual == False:
                ax[0].plot(dates_local_daily,surfacelay_dict[station]['TEF Vertical'],color=vertX_color,
                       linewidth=1,label='TEF Vertical')
                if show_EU:
                    ax[0].plot(dates_local_daily,surfacelay_dict[station]['EU Vertical'],color=vertX_color,
                            linewidth=3,alpha=0.3,label='EU Vertical')
            if residual == True:
                ax[0].legend(loc='lower center',ncol=6)
            else:
                ax[0].legend(loc='lower center',ncol=5)
            
            # plot deep
            if residual == False:
                ax[1].plot(dates_local_daily,bottomlay_dict[station]['TEF Exchange Flow'],color=exchange_color,
                        linewidth=1)
                if show_EU:
                    ax[1].plot(dates_local_daily,bottomlay_dict[station]['EU Exchange Flow'],color=exchange_color,
                            linewidth=3,alpha=0.3)
            else:
                ax[1].plot(dates_local_daily,bottomlay_dict[station]['TEF Recirculation'],color=exchange_color,
                        linewidth=2)
                if show_EU:
                    ax[1].plot(dates_local_daily,bottomlay_dict[station]['EU Recirculation'],color=exchange_color,
                            linewidth=3,alpha=0.3)
            ax[1].plot(dates_local_daily,bottomlay_dict[station]['TRAPS'],color=traps_color,
                       linewidth=3,zorder=2)
            ax[1].plot(dates_local_daily,bottomlay_dict[station]['Photosynthesis'],color=photo_color,
                       linewidth=2)
            ax[1].plot(dates_local_daily,bottomlay_dict[station]['Bio Consumption'],color=cons_color,
                       linewidth=2,linestyle=':',zorder=9)
            ax[1].plot(dates_local_daily,bottomlay_dict[station]['Storage'],color=ddtDOV_color,
                       linewidth=1,zorder=7)#,alpha=0.6)
            if residual == False:
                ax[1].plot(dates_local_daily,bottomlay_dict[station]['TEF Vertical'],color=vertX_color,
                        linewidth=1)
                if show_EU:
                    ax[1].plot(dates_local_daily,bottomlay_dict[station]['EU Vertical'],color=vertX_color,
                            linewidth=3,alpha=0.3)

            # plot error
            error_TEF = surfacelay_dict[station]['TEF Vertical']+bottomlay_dict[station]['TEF Vertical']
            error_EU = surfacelay_dict[station]['EU Vertical']+bottomlay_dict[station]['EU Vertical']
            ax[2].plot(dates_local_daily,error_TEF,color='crimson',
                       linewidth=2,label='TEF')
            if show_EU:
                ax[2].plot(dates_local_daily,error_EU,color='k',
                        linewidth=1,linestyle='--',label='EU')
                ax[2].legend(loc='upper right')
                
            for axis in [ax[0],ax[1],ax[2]]:
                if residual == False:
                    ylimval = np.nanmax([np.abs(surfacelay_dict[station]['TEF Exchange Flow']),
                                         np.abs(surfacelay_dict[station]['TEF Vertical']),
                                         np.abs(surfacelay_dict[station]['Photosynthesis']),
                                         np.abs(surfacelay_dict[station]['TRAPS']),
                                         np.abs(surfacelay_dict[station]['Bio Consumption']),])*1.2
                else:
                    ylimval = np.nanmax([np.abs(surfacelay_dict[station]['TEF Recirculation']),
                                         np.abs(surfacelay_dict[station]['Photosynthesis']),
                                         np.abs(surfacelay_dict[station]['TRAPS']),
                                         np.abs(surfacelay_dict[station]['Bio Consumption']),])*1.2
                axis.set_ylim([-1*ylimval,ylimval])

# plot DO concentrations -------------------------------------------------------------------------------

            # plot DO
            ax[3].plot(dates_local_daily,DOconcen_dict[station]['Surface Layer'],color='royalblue',
                       linewidth=2,label='Avg. surface layer')
            ax[3].plot(dates_local_daily,DOconcen_dict[station]['Deep Layer'],color='mediumorchid',
                       linewidth=2,label='Avg. deep layer')
            ax[3].plot(dates_local_daily,DOconcen_dict[station]['Qin DO'],color='black',
                       linewidth=1,linestyle=':',label='Avg. Qin DO')
            ax[3].plot(dates_local_daily,DOconcen_dict[station]['Bottom Sigma DO'],color='crimson',
                       linewidth=1,linestyle='--',label=r'Avg. bottom $\sigma$-layer')
            ax[3].plot(dates_local_daily,DOconcen_dict[station]['Minimum Bottom Layer DO'],color='crimson',
                       linewidth=1,linestyle='-',label='Minimum Bottom Layer DO')
            ax[3].legend(loc='lower left',ncol=2)


            ax[3].set_ylim([0,12])

# ---------------------------------- save figure --------------------------------------------

            plt.subplots_adjust(hspace=0.06, bottom=0.06, top=0.94)
            if residual == True:
                if show_EU == True:
                    out_dir_budget = out_dir / 'recirculation_withEU'
                else:
                    out_dir_budget = out_dir / 'recirculation'
            else:
                if show_EU == True:
                    out_dir_budget = out_dir / 'horizvert_withEU'
                else:
                    out_dir_budget = out_dir / 'horizvert'

            Lfun.make_dir(out_dir_budget)
            plt.savefig(out_dir_budget / (station + '.png') )

##########################################################
##          DO rate budget summary bar chart           ## 
##########################################################

budget_barchart = False

# COLLAPSE
if budget_barchart == True:

    seasons = ['Annual','Jan/Feb','Mar/Apr','May/Jun','Jul/Aug',
            'Sep/Oct','Nov/Dec']

    x = np.arange(len(seasons))

    # create list of station colors
    # station_color = ['hotpink','orange','gold','limegreen','deepskyblue'] # rainbow basic
    # station_color = ['#81667A','#8C8A93','#92B4A7','#B6CB9E','#D1F0B1']   # lilypad theme
    # station_color = ['#db5375','#E99F7D','#F4D35E','#B3C169','#6EB0C1']   # rainbow muted
    # station_color = ['#EF476F','darkorange','#FFD166','#B3C169','#6EB0C1']    # rainbow muted 2
    station_color = ['#F9627D','#62B6CB','#A8C256','#96031A','#957FEF','#476ad1','darkorange']    # pretty

    # initialize figure
    plt.close('all')
    fig, ax = plt.subplots(4,1,figsize = (13,8.6),sharex=True)

    # format figure
    plt.suptitle(group + ' summary of DO Budget Processes ('+year+')',size=16,fontweight='bold')
    # subtitles
    ax[0].set_title('(a) DO concentration',loc='left', fontsize=12, fontweight='bold')
    ax[1].set_title('(b) Surf Layer: average DO transport rate per unit volume',loc='left', fontsize=12, fontweight='bold')
    ax[2].set_title('(c) Deep Layer: average DO transport rate per unit volume',loc='left', fontsize=12, fontweight='bold')
    # ax[3].set_title('(d) Average ratio of sources to sinks',loc='left', fontsize=12, fontweight='bold')
    ax[3].set_title(r'(d) $\frac{d}{dt}DO$ per unit volume',loc='left', fontsize=12, fontweight='bold')
    # format grid
    for axis in [ax[0],ax[1],ax[2],ax[3]]:
        axis.set_facecolor('#EEEEEE')
        axis.grid(True,color='w',linewidth=1,linestyle='-',axis='y')
        for border in ['top','right','bottom','left']:
            axis.spines[border].set_visible(False)
        axis.tick_params(axis='y', labelsize=12)
    # label y-axis
    ax[0].set_ylabel('[mg/L]',fontsize=14)
    ax[1].set_ylabel(r'[mg/L day$^{-1}$]',fontsize=14)
    ax[2].set_ylabel(r'[mg/L day$^{-1}$]',fontsize=14) #(r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]',fontsize=14)
    ax[3].set_ylabel(r'[mg/L day$^{-1}$]',fontsize=14)

    width = 0.75/len(sta_dict)
    multiplier_conc = 0
    multiplier_rate_int = 0
    multiplier_rate_avg = 0
    multiplier_ratio = 0

    loop = 0
    # loop through stations
    for i,station in enumerate(sta_dict):
            
        # loop through different time intervals
        for j,season in enumerate(seasons):

            if season == 'Annual':
                minday = 0
                maxday = 363
            elif season == 'Jan/Feb':
                minday = 0
                maxday = 60
            elif season == 'Mar/Apr':
                minday = 60
                maxday = 121
            elif season == 'May/Jun':
                minday = 121
                maxday = 182
            elif season == 'Jul/Aug':
                minday = 182
                maxday = 244
            elif season == 'Sep/Oct':
                minday = 244
                maxday = 305
            elif season == 'Nov/Dec':
                minday = 305
                maxday = 363


            # DO concentrations --------------------------------------------------------
            # calculate average DO over time interval
            surf_lay_DO = np.nanmean(DOconcen_dict[station]['Surface Layer'][minday:maxday])
            deep_lay_DO = np.nanmean(DOconcen_dict[station]['Deep Layer'][minday:maxday])
            bott_DO = np.nanmean(DOconcen_dict[station]['Bottom Sigma DO'][minday:maxday])
            # calculate minimum DO in each season
            min_DO = np.nanmin(DOconcen_dict[station]['Minimum Bottom Layer DO'][minday:maxday])
            # get average inflowing DO concentration
            DO_in = np.nanmean(DOconcen_dict[station]['Qin DO'][minday:maxday])
            # plot
            offset = width * multiplier_conc
            if j == 0 and i==0:
                ax[0].scatter(j + offset, deep_lay_DO, zorder=5, s=120, color='white',edgecolor=station_color[i], label='Deep layer average')
                ax[0].scatter(j + offset, surf_lay_DO, zorder=5, s=60,color=station_color[i], label='Surface layer average')
                ax[0].scatter(j + offset, min_DO, zorder=5, s=30,marker='s',color=station_color[i],
                                    edgecolor='white', label='Absolute minimum recorded')
                ax[0].scatter(j + offset, DO_in, zorder=5, s=40, marker='+', color='black', label='Inflowing DO (TEF)')
            else:
                ax[0].scatter(j + offset, deep_lay_DO, zorder=5, s=120, color='white',edgecolor=station_color[i])
                ax[0].scatter(j + offset, surf_lay_DO, zorder=5, s=60, color=station_color[i])
                ax[0].scatter(j + offset, min_DO, zorder=5, s=30,marker='s',color=station_color[i], edgecolor='white')
                ax[0].scatter(j + offset, DO_in, zorder=5, s=40, marker='+', color='black')

            # surface layer transport terms (volume averaged)--------------------------------------------------------
            bottom_surf_pos = 0 # location of bottom of bar
            bottom_surf_neg = 0 # location of bottom of bar
            for attribute, measurement in surfacelay_dict[station].items():
                # skip variables we are not interested in
                if attribute in ['EU Exchange Flow', 'TEF Exchange Flow', 'EU Recirculation', 'EU Vertical', 'TEF Vertical', 'Volume', 'Storage']:
                    continue
                # calculate time average
                time_avg = np.nanmean(measurement[minday:maxday])
                # get volume average
                avg = time_avg/(np.nanmean(surfacelay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
                # convert to umol O2 /s /m3
                avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
                # convert to mg/L per day
                avg = avg * (32/1000/1000) * (60*60*24)
                # pick hatching pattern
                if attribute == 'Photosynthesis':
                    custom_color = 'white'
                elif attribute == 'TRAPS' or attribute == 'Air-Sea Transfer':
                    custom_color = 'black'
                else:
                    custom_color = station_color[i]
                if attribute == 'Bio Consumption':
                    hatch='//'
                    edgecolor = None
                    alpha=0.2
                elif attribute == 'Photosynthesis':
                    hatch='XX'
                    edgecolor = station_color[i]
                elif attribute == 'Air-Sea Transfer':
                    hatch = 'oo'
                    alpha = 0.4
                else:
                    hatch=None
                    edgecolor = None
                    alpha=1
                # decide where to start the bar, and increment where next bar should start
                if avg < 0:
                    bottom_surf = bottom_surf_neg
                    bottom_surf_neg += avg
                if avg > 0:
                    bottom_surf = bottom_surf_pos
                    bottom_surf_pos += avg
                # plot
                offset = width * multiplier_rate_int
                if loop == 0:
                    ax[1].bar(j + offset, avg, width, label=attribute, zorder=5, color=custom_color,
                        bottom=bottom_surf, hatch=hatch, edgecolor=edgecolor, alpha=alpha)
                else:
                    ax[1].bar(j + offset, avg, width, zorder=5, color=custom_color,
                        bottom=bottom_surf, hatch=hatch, edgecolor=edgecolor, alpha=alpha)
                    
            # deep layer transport terms (volume averaged)--------------------------------------------------------
            bottom_deep_pos = 0 # location of bottom of bar
            bottom_deep_neg = 0 # location of bottom of bar
            for attribute, measurement in bottomlay_dict[station].items():
                # skip variables we are not interested in
                if attribute in ['EU Exchange Flow', 'TEF Exchange Flow', 'EU Recirculation', 'EU Vertical', 'TEF Vertical', 'Volume', 'Storage','Qin m3/s']:
                    continue
                # calculate time average
                time_avg = np.nanmean(measurement[minday:maxday])
                # get volume average
                avg = time_avg/(np.nanmean(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
                # convert to umol O2 /s /m3
                avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
                # convert to mg/L per day
                avg = avg * (32/1000/1000) * (60*60*24)
                # pick hatching pattern
                if attribute == 'Photosynthesis':
                    custom_color = 'white'
                elif attribute == 'TRAPS':
                    custom_color = 'black'
                else:
                    custom_color = station_color[i]
                if attribute == 'Bio Consumption':
                    hatch='//'
                    edgecolor = None
                    alpha=0.2
                elif attribute == 'Photosynthesis':
                    hatch='XX'
                    edgecolor = station_color[i]
                else:
                    hatch=None
                    edgecolor = None
                    alpha=1
                # decide where to start the bar, and increment where next bar should start
                if avg < 0:
                    bottom_deep = bottom_deep_neg
                    bottom_deep_neg += avg
                if avg > 0:
                    bottom_deep = bottom_deep_pos
                    bottom_deep_pos += avg
                # plot
                offset = width * multiplier_rate_avg
                if loop == 0:
                    ax[2].bar(j + offset, avg, width, label=attribute, zorder=5, color=custom_color,
                        bottom=bottom_deep, hatch=hatch, edgecolor=edgecolor, alpha=alpha)
                else:
                    ax[2].bar(j + offset, avg, width, zorder=5, color=custom_color,
                        bottom=bottom_deep, hatch=hatch, edgecolor=edgecolor, alpha=alpha)

            # increment loop after first iteration
            loop = 1

            # storage terms --------------------------------------------------------
            # calculate d/dt(DO)
            storage_all_time_surf = surfacelay_dict[station]['Storage'][minday:maxday]/surfacelay_dict[station]['Volume'][minday:maxday]
            storage_all_time_deep = bottomlay_dict[station]['Storage'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]
            storage_surf = np.nanmean(storage_all_time_surf)
            storage_deep = np.nanmean(storage_all_time_deep)
            # convert from kmol to umol
            storage_surf = storage_surf * 1000 * 1000 * 1000 # umol O2 /s /m3
            storage_deep = storage_deep * 1000 * 1000 * 1000 # umol O2 /s /m3
            # convert to mg/L per day
            storage_surf = storage_surf * (32/1000/1000) * (60*60*24) 
            storage_surf = storage_deep * (32/1000/1000) * (60*60*24)
            # decide where to start next bar
            if storage_surf > 0 and storage_deep > 0:
                bottom_storage = storage_deep
            elif storage_surf < 0 and storage_deep < 0:
                bottom_storage = storage_deep
            else:
                bottom_storage = 0
            # plot
            offset = width * multiplier_ratio
            if j == 0 and i == 0:
                ax[3].scatter(j + offset, storage_deep, zorder=5, s=120, color='white',edgecolor=station_color[i], label='Deep layer average')
                ax[3].scatter(j + offset, storage_surf, zorder=5, s=60,color=station_color[i], label='Surface layer average')
            else:
                ax[3].scatter(j + offset, storage_deep, zorder=5, s=120, color='white',edgecolor=station_color[i])
                ax[3].scatter(j + offset, storage_surf, zorder=5, s=60,color=station_color[i])

        # advance multipliers
        multiplier_conc += 1
        multiplier_rate_int += 1
        multiplier_rate_avg += 1
        multiplier_ratio += 1
        
        # add legend for stations
        ax[3].text(0.2+0.11*i, 4.64,'•'+station,fontsize=12,fontweight='bold',ha='left', va='bottom', transform=ax[3].transAxes,
               color=station_color[i])

    # format figure --------------------------------------------------------
    # add legend
    ax[0].legend(loc='upper left', fontsize=10, ncol=4)
    ax[1].legend(loc='upper left', fontsize=10, ncol=5)
    ax[2].legend(loc='upper left', fontsize=10, ncol=4)
    ax[3].legend(loc='upper left', fontsize=10, ncol=2)
    # legend for entire plot (with stations labeled)
    # ax[3].legend(loc='upper center', ncol=5, fontsize=12,handletextpad=0.1,
    #             bbox_to_anchor=(0.64, 4.86), frameon=False)
    

    # restrict y axis
    ax[0].set_ylim([0,14])
    # ax[1].set_ylim([-1.2,1.2])
    # ax[2].set_ylim([-1.2,1.2])
    # ax[3].set_ylim([-0.6,0.6])

    # lines
    ax[0].axhline(2,0,5, linewidth=2, linestyle=':', color='k', zorder=4) # hypoxia
    ax[1].axhline(0,0,5, linewidth=2, linestyle=':', color='k', zorder=6) # zero
    ax[2].axhline(0,0,5, linewidth=2, linestyle=':', color='k', zorder=6) # zero
    ax[3].axhline(0,0,5, linewidth=2, linestyle=':', color='k', zorder=6) # one
    # vertical lines
    for addit in [0,1,2,3,4,5]:
        if addit == 0:
            color = 'gray'
        else:
            color = 'white'
        line_loc = 0.5*(1+(len(sta_dict)-1)*width)
        ax[0].axvline(line_loc+addit,0,12, linewidth=1, color=color)
        ax[1].axvline(line_loc+addit,-0.3,1, linewidth=1, color=color)
        ax[2].axvline(line_loc+addit,-0.3,1, linewidth=1, color=color)
        ax[3].axvline(line_loc+addit,0,5, linewidth=1, color=color)
    # add category labels on x-axis
    ax[3].set_xticks(line_loc-0.5+x, seasons, fontsize=14)

    # # add category labels on x-axis
    # ax[3].set_xticks(x*(0.25 + width*5) + 0.3, seasons, fontsize=14)

    plt.subplots_adjust(hspace=0.2, top=0.9, bottom=0.05, right=0.95)
    plt.show()

    plt.savefig(out_dir / (group + '_barchart_'+year+'.png') )


##########################################################
##              Comparing deep DO budgets               ## 
##########################################################

# June and July
minday = 150
maxday = 211

# bar width
width = 0.1

letters = ['a','b']

twoinlet_budget_comparison = False
# COLLAPSE
if twoinlet_budget_comparison == True:

    multiplier_deep = 0

    fig, ax = plt.subplots(1,2,figsize = (12,5),sharey=True)
    ax = ax.ravel()

    for i,station in enumerate(['crescent','lynchcove']):

        # format grid
        ax[i].set_facecolor('#EEEEEE')
        ax[i].tick_params(axis='x', labelrotation=30)
        ax[i].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax[i].spines[border].set_visible(False)
        ax[i].tick_params(axis='both', labelsize=12)

        for attribute, measurement in bottomlay_dict[station].items():

            # skip variables we are not interested in
            if attribute in ['EU Exchange Flow',
                             'TEF Exchange Flow',
                             'TRAPS',
                             'EU Recirculation',
                             'EU Vertical',
                             'TEF Vertical',
                             'Volume',
                             'Qin m3/s']:
                continue

             # calculate time average
            time_avg = np.nanmean(measurement[minday:maxday])
            # get volume average
            avg = time_avg/(np.nanmean(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
            # convert to umol O2 /s /m3
            avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
            # convert to mg/L per day
            avg = avg * (32/1000/1000) * (60*60*24)

            # choose color
            if attribute == 'TEF Recirculation':
                color = '#0C6291'
                label = 'Advection & Mixing'
            if attribute == 'Photosynthesis':
                color = '#A1CE37'
                label = attribute
            if attribute == 'Bio Consumption':
                color = '#A63446'
                label = attribute
            if attribute == 'Storage':
                color = '#000004'
                label = r'$\frac{d}{dt}$(DO)'

            # plot
            offset = width * multiplier_deep
            rects = ax[i].bar(i + offset, avg, width, zorder=5, color=color,label=label)
            multiplier_deep += 1

            # format
            ax[i].set_ylim([-0.3,0.3])
            ax[i].set_title('('+letters[i]+') '+ station, loc='left',fontsize=12)

            # add mean Aug/Sep DO
            augsep_DO = np.nanmean(DOconcen_dict[station]['Deep Layer'][211:272])
            ax[i].text(0.1, 0.1, 'Mean Aug & Sep DO: {} mg/L'.format(round(augsep_DO,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',zorder=6,
                        transform=ax[i].transAxes, fontsize=12, fontweight='bold')
            

    ax[0].legend(loc='upper center',ncol=2,fontsize=12)
    ax[0].set_ylabel('mg/L per day',fontsize=12)

    plt.suptitle('Volume-averaged deep layer budgets (Jun & Jul)',fontsize=14) 

    plt.subplots_adjust(left=0.1, wspace=0.02, top=0.85, bottom=0.1, right=0.9)
    plt.show()


# mid July to mid August
# June and July
minday = 194
maxday = 225

hypinlet_budget_comparison = False 
# COLLAPSE
if hypinlet_budget_comparison == True:

    multiplier_deep1 = 0
    multiplier_deep2 = 0

    fig, ax = plt.subplots(1,1,figsize = (8,5))
    # format grid
    ax.set_facecolor('#EEEEEE')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='w',linewidth=1,linestyle='-',axis='y')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_ylim([-2,1.4])
    ax.set_ylabel('mg/L per day',fontsize=12)

    # create a new dictionary of results
    oxy_dict = {}
    hyp_dict = {}

    for station in sta_dict:
        for attribute, measurement in bottomlay_dict[station].items():
            # skip variables we are not interested in
            if attribute in ['EU Exchange Flow',
                             'TEF Exchange Flow',
                             'TRAPS',
                             'EU Recirculation',
                             'EU Vertical',
                             'TEF Vertical',
                             'Volume',
                             'Photosynthesis & Consumption'
                             'Qin m3/s']:
                continue
             # calculate time average
            time_avg = np.nanmean(measurement[minday:maxday])
            # get volume average
            avg = time_avg/(np.nanmean(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
            # convert to umol O2 /s /m3
            avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
            # convert to mg/L per day
            avg = avg * (32/1000/1000) * (60*60*24)
            # save values in dictionary
            if station in ['penn','case','holmes','portsusan','lynchcove','dabob']:
                if attribute in hyp_dict.keys():
                    hyp_dict[attribute].append(avg)
                else:
                    hyp_dict[attribute] = [avg]
            else:
                if attribute in oxy_dict.keys():
                    oxy_dict[attribute].append(avg)
                else:
                    oxy_dict[attribute] = [avg]

    # t-test
    for attribute in oxy_dict:
        print('\n========================')
        print(attribute)
        a = oxy_dict[attribute]
        b = hyp_dict[attribute]
        ttest = ttest_ind(a, b, axis=0, equal_var=False)
        print(ttest)


    for i,dict in enumerate([oxy_dict,hyp_dict]):
    # average all oxygenated and hypoxic inlet rate values, and calculate standard deviations
        for attribute, measurement in dict.items():
            # choose color
            if attribute == 'TEF Recirculation':
                color = '#0C6291'
                label = 'Advection & Mixing'
            if attribute == 'Photosynthesis':
                color = '#A1CE37'
                label = attribute
            if attribute == 'Bio Consumption':
                color = '#A63446'
                label = attribute
            if attribute == 'Storage':
                color = 'gray'#'#000004'
                label = r'$\frac{d}{dt}$(DO)'
            # calculate average and standard deviation
            avg = np.nanmean(measurement)
            std = np.std(measurement)
            if avg < 0:
                wiggle = 0.12
            if avg > 0:
                wiggle = -0.12
            # plot
            offset = width * multiplier_deep1
            if i == 0:
                rects = ax.bar(offset, avg, width, zorder=5, edgecolor=color,color=color,label=label)
                ax.errorbar(offset, avg, yerr=std, capsize=5, color="k",zorder=6)
                ax.text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                        fontsize=12)
                multiplier_deep1 += 2
            elif i == 1:
                offset = width * multiplier_deep2
                rects = ax.bar(offset+width, avg, width, zorder=5, edgecolor=color,color='white', hatch='xx')
                ax.errorbar(offset+width, avg, yerr=std, capsize=5, color="k",zorder=6)
                ax.text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                        fontsize=12)
                multiplier_deep2 += 2
            
            ax.legend(loc='upper right',ncol=2,fontsize=12)

    ax.text(0.42, 0.05, 'SOLID: oxygenated (n={})\nHATCHED: hypoxic (n={})        '.format(len(oxy_dict['Storage']),len(hyp_dict['Storage'])),
            color='black', verticalalignment='bottom', horizontalalignment='right',zorder=6,
            transform=ax.transAxes, fontsize=12)
    plt.suptitle('Volume-averaged deep layer budgets (mid-Jul to mid-Aug)',fontsize=14) 

    plt.subplots_adjust(left=0.1, wspace=0.02, top=0.92, bottom=0.1, right=0.9)
    plt.show()


hypinlet_budget_comparison_v2 = True # ------------------------------------------- MONEY BAR CHART & t-tests
# COLLAPSE
if hypinlet_budget_comparison_v2 == True:

    multiplier_deep1 = 0
    multiplier_deep2 = 0

    fig, ax = plt.subplots(2,1,figsize = (9,9))

    # format grid
    for axis in [ax[0],ax[1]]:
        axis.set_facecolor('#EEEEEE')
        axis.tick_params(axis='x', labelrotation=30)
        axis.grid(True,color='w',linewidth=1,linestyle='-',axis='y')
        for border in ['top','right','bottom','left']:
            axis.spines[border].set_visible(False)
        axis.tick_params(axis='y', labelsize=12)
        axis.set_xticklabels([])
        axis.set_ylabel('mg/L per day',fontsize=12)
    ax[1].set_ylim([-1,1])

    # create a new dictionary of results
    oxy_dict = {}
    hyp_dict = {}

    # part 1 with distinct physical and biological terms
    for station in sta_dict:
        for attribute, measurement in bottomlay_dict[station].items():
            # skip variables we are not interested in
            if attribute in ['EU Exchange Flow',
                             'TRAPS',
                             'EU Recirculation',
                             'EU Vertical',
                             'TEF Recirculation',
                             'Photosynthesis & Consumption',
                             'Volume',
                             'Qin m3/s']:
                continue
             # calculate time average
            time_avg = np.nanmean(measurement[minday:maxday])
            # get volume average
            avg = time_avg/(np.nanmean(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
            # convert to umol O2 /s /m3
            avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
            # convert to mg/L per day
            avg = avg * (32/1000/1000) * (60*60*24)
            # save values in dictionary
            if station in ['penn','case','holmes','portsusan','lynchcove','dabob']:
                if attribute in hyp_dict.keys():
                    hyp_dict[attribute].append(avg)
                else:
                    hyp_dict[attribute] = [avg]
            else:
                if attribute in oxy_dict.keys():
                    oxy_dict[attribute].append(avg)
                else:
                    oxy_dict[attribute] = [avg]
    # t-test
    print('DISTINCT TERMS --------------------------------------')
    for attribute in oxy_dict:
        print('\n========================')
        print(attribute)
        a = oxy_dict[attribute]
        b = hyp_dict[attribute]
        ttest = ttest_ind(a, b, axis=0, equal_var=False)
        print(ttest)
    print('\n')
    for i,dict in enumerate([oxy_dict,hyp_dict]):
    # average all oxygenated and hypoxic inlet rate values, and calculate standard deviations
        for attribute, measurement in dict.items():
            # choose color
            if attribute == 'TEF Exchange Flow':
                color = '#a3bcf7'
                label = 'Exchange Flow'
            if attribute == 'TEF Vertical':
                color = '#b9a3f7'
                label = 'Vertical Transport'
            if attribute == 'Photosynthesis':
                color = '#88d7d8'
                label = attribute
            if attribute == 'Bio Consumption':
                color = 'darkgray'
                label = attribute
            if attribute == 'Storage':
                color = 'black'
                label = r'$\frac{d}{dt}$(DO)'
            # calculate average and standard deviation
            avg = np.nanmean(measurement)
            std = np.std(measurement)
            if avg < 0:
                wiggle = 0.4
            if avg > 0:
                wiggle = -0.6
            # plot
            offset = width * multiplier_deep1
            if attribute == 'Storage':
                hatchcolor = 'white'
            else:
                hatchcolor = 'black'
            if i == 0:
                rects = ax[0].bar(offset, avg, width, zorder=5, edgecolor=hatchcolor,color=color, hatch='xx')
                rects = ax[0].bar(offset, avg, width, zorder=5, edgecolor=color,color='none')
                # ax[1].errorbar(offset, avg, yerr=std, capsize=5, color="k",zorder=6)
                # ax[0].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                #         fontsize=12)
                if attribute == 'Storage':
                    ax[0].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color=color,fontsize=12, fontweight='bold')
                else:
                    ax[0].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color='gray',fontsize=12)
                multiplier_deep1 += 2
            elif i == 1:
                offset = width * multiplier_deep2
                rects = ax[0].bar(offset+width, avg, width, zorder=5, edgecolor=color,color=color,label=label)
                # ax[1].errorbar(offset+width, avg, yerr=std, capsize=5, color="k",zorder=6)
                # ax[0].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                #         fontsize=12)
                if attribute == 'Storage':
                    ax[0].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color=color,fontsize=12, fontweight='bold')
                else:
                    ax[0].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color='gray',fontsize=12)
                multiplier_deep2 += 2
            
            ax[0].legend(loc='upper right', fontsize=12, ncol=2)

    # part 2 with combined physical and biological terms ---------------------------------

    multiplier_deep1 = 0
    multiplier_deep2 = 0

    # create a new dictionary of results
    oxy_dict = {}
    hyp_dict = {}

    
    for station in sta_dict:
        for attribute, measurement in bottomlay_dict[station].items():
            # skip variables we are not interested in
            if attribute in ['EU Exchange Flow',
                             'TEF Exchange Flow',
                             'TRAPS',
                             'EU Recirculation',
                             'EU Vertical',
                             'TEF Vertical',
                             'Photosynthesis',
                             'Bio Consumption',
                             'Volume',
                             'Qin m3/s']:
                continue
             # calculate time average
            time_avg = np.nanmean(measurement[minday:maxday])
            # get volume average
            avg = time_avg/(np.nanmean(bottomlay_dict[station]['Volume'][minday:maxday])) # kmol O2 /s /m3
            # convert to umol O2 /s /m3
            avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
            # convert to mg/L per day
            avg = avg * (32/1000/1000) * (60*60*24)
            # save values in dictionary
            if station in ['penn','case','holmes','portsusan','lynchcove','dabob']:
                if attribute in hyp_dict.keys():
                    hyp_dict[attribute].append(avg)
                else:
                    hyp_dict[attribute] = [avg]
            else:
                if attribute in oxy_dict.keys():
                    oxy_dict[attribute].append(avg)
                else:
                    oxy_dict[attribute] = [avg]

    # t-test
    print('COMBINED TERMS --------------------------------------')
    for attribute in oxy_dict:
        print('\n========================')
        print(attribute)
        a = oxy_dict[attribute]
        b = hyp_dict[attribute]
        if attribute == 'Photosynthesis & Consumption':
            print('one-sided')
            # one-sided t-test 
            # alternative hypothesis: a < b (drawdown of oxygenated is more negative than drawdown of hypoxic)
            # the null hypothesis is thus that drawdown of hypoxic is more negative
            ttest = ttest_ind(a, b, axis=0, equal_var=False, alternative='less')
        else:
            # two-sided t-test with null hypothesis that they are the same
            ttest = ttest_ind(a, b, axis=0, equal_var=False)
        print(ttest)

    print('\n')


    for i,dict in enumerate([oxy_dict,hyp_dict]):
    # average all oxygenated and hypoxic inlet rate values, and calculate standard deviations
        for attribute, measurement in dict.items():
            # choose color
            if attribute == 'TEF Recirculation':
                color = '#b9a3f7'
                label = 'Advection & Mixing'
            # if attribute == 'Photosynthesis':
            #     color = '#A1CE37'
            #     label = attribute
            # if attribute == 'Bio Consumption':
            #     color = '#A63446'
            #     label = attribute
            if attribute == 'Photosynthesis & Consumption':
                color = '#88d7d8'
                label = attribute
            if attribute == 'Storage':
                color = 'black'
                label = r'$\frac{d}{dt}$(DO)'
            # calculate average and standard deviation
            avg = np.nanmean(measurement)
            std = np.std(measurement)
            if avg < 0:
                wiggle = 0.06
            if avg > 0:
                wiggle = -0.1
            # plot
            offset = width * multiplier_deep1
            if attribute == 'Storage':
                hatchcolor = 'white'
            else:
                hatchcolor = 'black'
            if i == 0:
                rects = ax[1].bar(offset, avg, width, zorder=5, edgecolor=hatchcolor,color=color, hatch='xx')
                rects = ax[1].bar(offset, avg, width, zorder=5, edgecolor=color,color='none')
                # ax[1].errorbar(offset, avg, yerr=std, capsize=5, color="k",zorder=6)
                # ax[1].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                #         color=color,fontsize=12)
                if attribute == 'Storage':
                    ax[1].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color=color,fontsize=12, fontweight='bold')
                else:
                    ax[1].text(offset, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color='gray',fontsize=12)
                    
                multiplier_deep1 += 2
            elif i == 1:
                offset = width * multiplier_deep2
                rects = ax[1].bar(offset+width, avg, width, zorder=5, edgecolor=color,color=color,label=label)
                # ax[1].errorbar(offset+width, avg, yerr=std, capsize=5, color="k",zorder=6)
                if attribute == 'Storage':
                    ax[1].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color=color,fontsize=12, fontweight='bold')
                else:
                    ax[1].text(offset+width, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color='gray',fontsize=12)
                multiplier_deep2 += 2
            
            ax[1].legend(loc='upper right', fontsize=12)

    ax[0].text(0.88, 0.1, 'HATCHED: oxygenated (n={})\nSOLID: hypoxic (n={})        '.format(len(oxy_dict['Storage']),len(hyp_dict['Storage'])),
            color='black', verticalalignment='bottom', horizontalalignment='right',zorder=6,
            transform=ax[0].transAxes, fontsize=12)
    plt.suptitle('Volume-averaged deep layer budgets (mid-Jul to mid-Aug)',fontsize=14) 

    ax[0].set_title('(a) All budget rate terms', loc='left',fontsize=12)
    ax[1].set_title('(b) Combined physical and biological processes', loc='left',fontsize=12)

    plt.subplots_adjust(left=0.1, wspace=0.02, top=0.92, bottom=0.1, right=0.9)
    plt.show()

# scatter plot of aug/sep min bottom sigma DO vs. dDO/dt during jun/jul
jjdDOdt_asbotDO = False
# COLLLAPSE
if jjdDOdt_asbotDO == True:
    # get values
    dDOdt = np.zeros(len(sta_dict))
    augsep_DO = np.zeros(len(sta_dict))
    for i,station in enumerate(sta_dict):
        # jun/jul d(DO)/dt -----------------
        time_avg = np.nanmean(bottomlay_dict[station]['Storage'][150:211])
        # get volume average
        avg = time_avg/(np.nanmean(bottomlay_dict[station]['Volume'][150:211])) # kmol O2 /s /m3
        # convert to umol O2 /s /m3
        avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
        # convert to mg/L per day
        dDOdt[i] = avg * (32/1000/1000) * (60*60*24)
        # aug/sep deep DO ----------------
        # augsep_DO[i] = np.nanmean(DOconcen_dict[station]['Minimum Bottom Layer DO'][211:272])
        # augsep_DO[i] = np.nanmean(DOconcen_dict[station]['Deep Layer'][211:272])
        augsep_DO[i] = np.nanmean(DOconcen_dict[station]['percent hypoxic volume'][211:272])

    # plot
    # initialize figure
    fig, ax = plt.subplots(1,1,figsize = (6,6))
    # format figure
    # plt.suptitle(r'Aug/Sep mean deep layer DO vs. Jun/Jul d(DO)/dt', size=14)
    plt.suptitle(r'Aug/Sep mean % hypoxic volume vs. Jun/Jul d(DO)/dt', size=14)
    # format grid
    ax.set_facecolor('#EEEEEE')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_xlabel('Jun/Jul d(DO)/dt [mg/L per day]', fontsize=12)
    # ax.set_ylabel(r'Aug/Sep mean deep layer DO [mg/L]', fontsize=12)
    ax.set_ylabel(r'Aug/Sep mean % hypoxic volume', fontsize=12)
    # plot
    ax.scatter(dDOdt,augsep_DO,alpha=0.5,s=100,zorder=5)
    ax.set_xlim([-0.06,0])
    ax.set_ylim([0,100])
    # calculate correlation coefficient (Pearson)
    r,p = pearsonr(dDOdt,augsep_DO)
    ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)) ,color='black',
                            verticalalignment='bottom', horizontalalignment='left',
                            transform=ax.transAxes, fontsize=12, fontweight='bold')
    ax.text(0.1, 0.79, r'$p =$' + str(round(p,5)) ,color='black',
                            verticalalignment='bottom', horizontalalignment='left',
                            transform=ax.transAxes, fontsize=12, fontweight='bold')
    for i,station in enumerate(sta_dict):
        ax.text(dDOdt[i],augsep_DO[i],station,horizontalalignment='right')
    # save figure
    plt.show()


# ##########################################################
# ##                Seasonal cycle analysis               ## 
# ##########################################################

# seasonal_cycle = True
# if seasonal_cycle == True:

#     # initialize figure
#     fig, ax = plt.subplots(1,1,figsize = (10,5))
#     # format figure
#     # plt.suptitle(year + r' DO$_{in}$-DO$_{deep}$ vs.' +
#     #               '\nExpected difference based on drawdown and flushing time',
#     #                 size=16)
#     plt.suptitle(year + r' (DO$_{in}$-DO$_{deep}$) vs. T$_{flush,ex}$ [60-day Hanning Window]',
#                     size=16)
#     # format grid
#     ax.set_facecolor('#EEEEEE')
#     ax.tick_params(axis='x', labelrotation=30)
#     ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#     for border in ['top','right','bottom','left']:
#         ax.spines[border].set_visible(False)
#     ax.tick_params(axis='both', labelsize=14)
#     ax.set_xlabel(r'T$_{flush,ex}$ [days]', fontsize=14)
#     ax.set_ylabel(r'DO$_{in}$-DO$_{deep}$ [mg/L]', fontsize=14)
#     # # colormap
#     # cmap_hyp = plt.cm.get_cmap('gist_heat_r')
#     # colors = plt.cm.gist_heat_r(np.linspace(0,1,80))
#     colors = plt.cm.tab20(np.linspace(0,1,20))
#     for i,station in enumerate(sta_dict):
#         # calculate flushing time
#         Tflush = dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'] / (60*60*24)
#         # calculate DOin - DOdeep
#         deltaDO = DOconcen_dict[station]['Qin DO'] - DOconcen_dict[station]['Deep Layer']

#         # get average DO consumption rate in mg/L per day
#         # get volume average
#         storage = bottomlay_dict[station]['Storage']
#         avg = storage/(np.nanmean(bottomlay_dict[station]['Volume'])) # kmol O2 /s /m3
#         # convert to umol O2 /s /m3
#         avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
#         # convert to mg/L per day
#         avg = avg * (32/1000/1000) * (60*60*24)
#         drawdown = avg

#         # get DO consumption
#         consumption = drawdown * Tflush # mg/L/day * days
#         drawdown = DOconcen_dict[station]['Qin DO'] - consumption

#         # filter
#         Tflush = zfun.lowpass(Tflush.values,n=60)
#         deltaDO = zfun.lowpass(deltaDO.values,n=60)
#         drawdown = zfun.lowpass(drawdown.values,n=30)
#         # remove padded nans
#         Tflush = Tflush[~np.isnan(Tflush)]
#         deltaDO = deltaDO[~np.isnan(deltaDO)]
#         drawdown = drawdown[~np.isnan(drawdown)]
#         # # calculate aug/sep % hypoxic volume
#         # augsep_hyp = np.nanmean(DOconcen_dict[station]['percent hypoxic volume'][211:272])
#         # plot
#         ax.plot(np.append(Tflush,Tflush[0]),np.append(deltaDO,deltaDO[0]),linewidth=3.5,
#                 color='black')
#         if station == 'henderson':
#             color = 'black'
#         else:
#             color = colors[i]
#         ax.plot(np.append(Tflush,Tflush[0]),np.append(deltaDO,deltaDO[0]),linewidth=2.5,
#                 color=color, label=station)#,
#         # ax.plot(np.append(drawdown,drawdown[0]),np.append(deltaDO,deltaDO[0]),linewidth=3,
#         #         color='k')
#         # ax.plot(np.append(drawdown,drawdown[0]),np.append(deltaDO,deltaDO[0]),linewidth=2)#,
#                 # color=colors[round(augsep_hyp)])#,
#                 # linewidth=1,color='k',alpha=0.5)
#         ax.set_xlim([0,80])
#     # cmap_oxy = plt.cm.get_cmap('rainbow_r')
#     # # cs_DO = ax.scatter(mean_Tflush,mean_DOin,s=45,zorder=5, edgecolor='k',c=deep_lay_DO, cmap=cmap_oxy)
#     # cs_DO = ax.scatter(mean_Tflush,mean_DOin,s=45,zorder=5, edgecolor='k',c=perc_hyp_vol, cmap=cmap_hyp)
#     # # create colorbarlegend
#     # cbar = fig.colorbar(cs_DO)
#     # cbar.ax.tick_params(labelsize=12)
#     # # cbar.ax.set_ylabel('Monthly mean deep layer DO [mg/L]', rotation=90, fontsize=14)
#     # cbar.ax.set_ylabel('Aug & Sep mean % hypoxic volume', rotation=90, fontsize=14)
#     # cbar.outline.set_visible(False)
#     # plt.xscale('log')
#     plt.legend(loc='lower right',ncol=3)
#     plt.tight_layout()
#     # save figure
#     plt.show()

##########################################################
##                   DO scatterplots                    ## 
##########################################################

DO_analysis = True
if DO_analysis == True:

    # initialize arrays for plotting
    intervals = 12
    deep_lay_DO = np.zeros(len(sta_dict)*intervals)
    bott_sig_DO = np.zeros(len(sta_dict)*intervals)
    min_bott_sig_DO = np.zeros(len(sta_dict)*intervals)
    annual_mean_DO = np.zeros(len(sta_dict)*intervals)
    perc_hyp_vol = np.zeros(len(sta_dict)*intervals)
    mean_DOin = np.zeros(len(sta_dict)*intervals)
    mean_Tflush = np.zeros(len(sta_dict)*intervals)
    mean_TEFin = np.zeros(len(sta_dict)*intervals)
    mean_recirc = np.zeros(len(sta_dict)*intervals)
    mean_cons = np.zeros(len(sta_dict)*intervals)
    mean_depth = np.zeros(len(sta_dict)*intervals)
    inlet_vol = np.zeros(len(sta_dict)*intervals)
    aspect_ratio = np.zeros(len(sta_dict)*intervals)

    annmean_DOin = np.zeros(len(sta_dict))
    annmin_DOin = np.zeros(len(sta_dict))
    colors_twentyone = []

    colors = []

    # get values for plotting and calculating r value
    for i,station in enumerate(sta_dict):
        # get interface depth from csv file
        with open('interface_depths.csv', 'r') as f:
            for line in f:
                inlet, interface_depth = line.strip().split(',')
                interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
        z_interface = float(interface_dict[station])
        fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_2014.01.01_2014.12.31.nc'
        ds = xr.open_dataset(fn)
        moor_depth = ds.h.values
        
        for month in range(intervals):
            if month == 0:
                minday = 1
                maxday = 32
            elif month == 1:
                minday = 32
                maxday = 60
            elif month == 2:
                minday = 60
                maxday = 91
            elif month == 3:
                minday = 91
                maxday = 121
            elif month == 4:
                minday = 121
                maxday = 152
            elif month == 5:
                minday = 152
                maxday = 182
            elif month == 6:
                minday = 182
                maxday = 213
            elif month == 7:
                minday = 213
                maxday = 244
            elif month == 8:
                minday = 244
                maxday = 274
            elif month == 9:
                minday = 274
                maxday = 305
            elif month == 10:
                minday = 305
                maxday = 335
            elif month == 11:
                minday = 335
                maxday = 363
            # minday=month
            # maxday=month+1
            # save values
            deep_lay_DO[i*intervals+month] =  np.nanmean(DOconcen_dict[station]['Deep Layer'][minday:maxday])
            bott_sig_DO[i*intervals+month] =  np.nanmean(DOconcen_dict[station]['Bottom Sigma DO'][minday:maxday])
            min_bott_sig_DO[i*intervals+month] =  np.nanmin(DOconcen_dict[station]['Minimum Bottom Layer DO'][minday:maxday])
            annual_mean_DO[i*intervals+month] = np.nanmean(DOconcen_dict[station]['Deep Layer'])
            perc_hyp_vol[i*intervals+month] = np.nanmean(DOconcen_dict[station]['percent hypoxic volume'][minday:maxday])
            mean_DOin[i*intervals+month] = np.nanmean(DOconcen_dict[station]['Qin DO'][minday:maxday])
            mean_Tflush[i*intervals+month] = np.nanmean(dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'][minday:maxday]) / (60*60*24)
            mean_TEFin[i*intervals+month] = np.nanmean(bottomlay_dict[station]['TEF Exchange Flow'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]) * (
                                    32 * 1000) * (60*60*24)
            # TESTING ------------- TESTING ------------- TESTING ------------- TESTING ------------- TESTING ------------- TESTING ------------- TESTING 
            # normalize by total inlet volume, rather than just bottom layer volume
            # mean_TEFin[i*12+month] = np.nanmean(bottomlay_dict[station]['TEF Exchange Flow'][minday:maxday]/dimensions_dict[station]['Inlet volume'][0]) * (
            #                         32 * 1000) * (60*60*24)
            # TESTING ------------- TESTING ------------- TESTING ------------- TESTING ------------- TESTING ------------- TESTING ------------- TESTING
            mean_recirc[i*intervals+month] = np.nanmean(bottomlay_dict[station]['TEF Recirculation'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]) * (
                                    32 * 1000) * (60*60*24)
            mean_cons[i*intervals+month] = np.nanmean(bottomlay_dict[station]['Bio Consumption'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]) * (
                                    32 * 1000) * (60*60*24)
            mean_depth[i*intervals+month] = dimensions_dict[station]['Mean depth'][0]
            inlet_vol[i*intervals+month] = dimensions_dict[station]['Inlet volume'][0]
            aspect_ratio[i*intervals+month] = dimensions_dict[station]['L/W aspect ratio'][0]
            colors.append(basin_color_dict[basin_dict[station]])

        annmean_DOin[i] = np.nanmean(DOconcen_dict[station]['Qin DO'])
        annmin_DOin[i] = np.nanmin(DOconcen_dict[station]['Qin DO'])
        colors_twentyone.append(basin_color_dict[basin_dict[station]])
        
        # # if annnual average
        # deep_lay_DO[i] =  np.nanmean(DOconcen_dict[station]['Deep Layer'])
        # mean_DOin[i] = np.nanmean(DOconcen_dict[station]['Qin DO'])
        # mean_Tflush[i] = np.nanmean(dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s']) / (60*60*24)
        # colors.append(basin_color_dict[basin_dict[station]])
        



    meanDOdeep_bottsigDO = False
    # COLLLAPSE
    if meanDOdeep_bottsigDO == True:
        # initialize figure
        # plt.close('all')
        fig, ax = plt.subplots(1,1,figsize = (5,5))
        # format figure
        plt.suptitle(year + ' monthly mean \n' + r'bottom $\sigma$ layer DO vs. deep layer DO',
                        size=14)
        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel('Monthly mean deep layer DO [mg/L]', fontsize=12)
        ax.set_ylabel(r'Monthly mean bottom $\sigma$ layer DO [mg/L]', fontsize=12)
        # plot
        ax.scatter(deep_lay_DO,bott_sig_DO,alpha=0.5,s=10,zorder=5,
                color=colors)
        ax.set_xlim([0,12])
        ax.set_ylim([0,12])
        # calculate correlation coefficient (Pearson)
        r,p = pearsonr(deep_lay_DO,bott_sig_DO)
        ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        ax.text(0.1, 0.79, r'$p =$' + str(round(p,10)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        # create legend
        whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='lower right')
        # save figure
        plt.show()
        # plt.savefig(out_dir / 'scatter_meanDO_bottsigDO.png' )

    meanDOdeep_perchypvol = True # ---------------------------  MEAN DEEP DO AND % HYPOXIC VOLUME
    # COLLLAPSE
    if meanDOdeep_perchypvol == True:
        # initialize figure
        # plt.close('all')
        fig, ax = plt.subplots(1,1,figsize = (5,5))
        # format figure
        plt.suptitle(year + ' monthly mean \n' + r'% hypoxic volume vs. DO$_{deep}$',
                        size=14)
        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='both', labelsize=12)
        ax.set_xlabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=14)
        ax.set_ylabel('Monthly mean % hypoxic volume', fontsize=14)
        # plot
        ax.scatter(deep_lay_DO,perc_hyp_vol,alpha=0.3,s=80,zorder=5,color='k')
                # color=colors)
        ax.set_xlim([0,12])
        ax.set_ylim([0,100])
        # calculate correlation coefficient (Pearson)
        # r,p = pearsonr(deep_lay_DO,perc_hyp_vol)
        # ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)) ,color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12)
        # ax.text(0.1, 0.79, r'$p =$' + str(round(p,10)) ,color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # # create legend
        # whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        # hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        # mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        # southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        # admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        # ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='upper right',fontsize=14)
        # save figure
        plt.tight_layout()
        plt.show()
        # plt.savefig(out_dir / 'scatter_meanDO_bottsigDO.png' )

    meanDOdeep_minbottsigDO = False
    # COLLLAPSE
    if meanDOdeep_minbottsigDO == True:
        # initialize figure
        fig, ax = plt.subplots(1,1,figsize = (5,5))
        # format figure
        plt.suptitle(year + ' monthly \n' + r'min bottom $\sigma$ layer DO vs. mean deep layer DO',
                        size=14)
        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel('Monthly mean deep layer DO [mg/L]', fontsize=12)
        ax.set_ylabel(r'Monthly minimum bottom $\sigma$ layer DO [mg/L]', fontsize=12)
        # plot
        ax.scatter(deep_lay_DO,min_bott_sig_DO,alpha=0.5,s=10,zorder=5,
                color=colors)
        ax.set_xlim([0,12])
        ax.set_ylim([0,12])
        # calculate correlation coefficient (Pearson)
        r,p = pearsonr(deep_lay_DO,min_bott_sig_DO)
        ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        ax.text(0.1, 0.79, r'$p =$' + str(round(p,10)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        # create legend
        whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='center left')
        # save figure
        plt.show()
        # plt.savefig(out_dir / 'scatter_meanDO_bottsigDO.png' )

    QinDO_meanDOdeep = False
    # COLLLAPSE
    if QinDO_meanDOdeep == True:
        # initialize figure
        plt.close('all')
        fig, ax = plt.subplots(1,1,figsize = (5,5))
        # format figure
        plt.suptitle(year + ' monthly mean \ndeep layer DO vs. inflowing TEF DO',
                        size=14)
        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel('Monthly mean inflowing TEF DO [mg/L]', fontsize=12)
        ax.set_ylabel('Monthly mean deep layer DO [mg/L]', fontsize=12)
        # plot
        ax.scatter(mean_DOin,deep_lay_DO,alpha=0.5,s=10,zorder=5,
                color=colors)
        ax.set_xlim([0,12])
        ax.set_ylim([0,12])
        # calculate correlation coefficient (Pearson)
        r,p = pearsonr(mean_DOin,deep_lay_DO)
        ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        ax.text(0.1, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        # create legend
        whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='lower right')
        # save figure
        plt.savefig(out_dir / 'scatter_TEFDOin_meanDO.png' )

    TEFDO_meanDOdeep = False
    # COLLLAPSE
    if TEFDO_meanDOdeep == True:
        # initialize figure
        plt.close('all')
        fig, ax = plt.subplots(1,1,figsize = (5,5))
        # format figure
        plt.suptitle(year + ' monthly mean \ndeep layer DO vs. TEF transport rate (vol-avg)',
                        size=14)
        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel('Monthly mean TEF transport rate [mg/L per day]', fontsize=12)
        ax.set_ylabel('Monthly mean deep layer DO [mg/L]', fontsize=12)
        # plot
        ax.scatter(mean_TEFin,deep_lay_DO,alpha=0.5,s=10,zorder=5,
                color=colors)
        # ax.set_xlim([0,12])
        ax.set_ylim([0,12])
        # calculate correlation coefficient (Pearson)
        r,p = pearsonr(mean_TEFin,deep_lay_DO)
        ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        ax.text(0.1, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        # create legend
        whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='lower right')
        plt.tight_layout()
        # save figure
        plt.savefig(out_dir / 'scatter_TEFin_meanDO.png' )

    recirc_meanDOdeep = False
    # COLLLAPSE
    if recirc_meanDOdeep == True:
        # initialize figure
        plt.close('all')
        fig, ax = plt.subplots(1,1,figsize = (5,5))
        # format figure
        plt.suptitle(year + ' monthly mean \ndeep layer DO vs. recirculation rate (vol-avg)',
                        size=14)
        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel('Monthly mean recirculation transport rate [mg/L per day]', fontsize=12)
        ax.set_ylabel('Monthly mean deep layer DO [mg/L]', fontsize=12)
        # plot
        ax.scatter(mean_recirc,deep_lay_DO,alpha=0.5,s=10,zorder=5,
                color=colors)
        # ax.set_xlim([0,12])
        ax.set_ylim([0,12])
        # calculate correlation coefficient (Pearson)
        r,p = pearsonr(mean_recirc,deep_lay_DO)
        ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        ax.text(0.1, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        # create legend
        whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='lower right')
        plt.tight_layout()
        # save figure
        plt.savefig(out_dir / 'scatter_recirc_meanDO.png' )

    cons_meanDOdeep = False
    # COLLLAPSE
    if cons_meanDOdeep == True:
        # initialize figure
        plt.close('all')
        fig, ax = plt.subplots(1,1,figsize = (5,5))
        # format figure
        plt.suptitle(year + ' monthly mean \ndeep layer DO vs. bio consumption rate (vol-avg)',
                        size=14)
        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel('Monthly mean bio consumption rate [mg/L per day]', fontsize=12)
        ax.set_ylabel('Monthly mean deep layer DO [mg/L]', fontsize=12)
        # plot
        ax.scatter(mean_cons,deep_lay_DO,alpha=0.5,s=10,zorder=5,
                color=colors)
        # ax.set_xlim([0,12])
        ax.set_ylim([0,12])
        # calculate correlation coefficient (Pearson)
        r,p = pearsonr(mean_cons,deep_lay_DO)
        ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        ax.text(0.1, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        # create legend
        whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='lower right')
        plt.tight_layout()
        # save figure
        plt.savefig(out_dir / 'scatter_cons_meanDO.png' )


    # DOin vs. Tflush colored by percent hypoxic volume ============== 2PART MONEY PLOT SCATTER ================================
    percent_hypoxic = True
    if percent_hypoxic == True:
        # initialize figure
        fig, ax = plt.subplots(1,2,figsize = (10,5))

        # format figure
        ax[0].set_title('(a) ' + year + r' monthly mean DO$_{deep}$ vs. DO$_{in}$' + '\n' + r'colored by mean T$_{flush}$',
                        size=16, loc='left')
        # format grid
        ax[0].set_facecolor('#EEEEEE')
        ax[0].tick_params(axis='x', labelrotation=30)
        ax[0].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax[0].spines[border].set_visible(False)
        ax[0].tick_params(axis='both', labelsize=14)
        ax[0].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=14)
        ax[0].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=14)
        # plot
        cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
        cmap_oxy = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
        ax[0].plot([0,12],[0,12],color='gray')
        ax[0].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4, fontsize=12)
        cs = ax[0].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
        # create colorbarlegend
        cbar = fig.colorbar(cs)
        cbar.ax.tick_params(labelsize=14)
        cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=14)
        cbar.outline.set_visible(False)
        ax[0].set_xlim([0,12])
        ax[0].set_ylim([0,12])


        # format figure
        ax[1].set_title('(b) ' + year + ' monthly mean '+r'DO$_{in}$ vs. T$_{flush}$'+'\ncolored by % hypoxic volume',
                        size=16, loc='left')
        # format grid
        ax[1].set_facecolor('#EEEEEE')
        ax[1].tick_params(axis='x', labelrotation=30)
        ax[1].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax[1].spines[border].set_visible(False)
        ax[1].tick_params(axis='both', labelsize=14)
        ax[1].set_xlabel(r'Monthly mean T$_{flush}$ [days]', fontsize=14)
        ax[1].set_ylabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=14)
        ax[1].set_ylim([0,11])
        ax[1].set_xlim([0,85])
        # plot
        cmap_hyp = plt.cm.get_cmap('gist_heat_r')
        cs_DO = ax[1].scatter(mean_Tflush,mean_DOin,s=60,zorder=5,edgecolor='gray',c=perc_hyp_vol,cmap=cmap_hyp)
        # create colorbarlegend
        cbar = fig.colorbar(cs_DO)
        cbar.ax.tick_params(labelsize=12)
        cbar.ax.set_ylabel('Monthly mean % hypoxic volume', rotation=90, fontsize=14)
        cbar.outline.set_visible(False)
        plt.tight_layout()
        # save figure
        plt.show()


    # difference between DOin and deep layer DO vs. Tflush
    deltaDO_Tflush = True
    if deltaDO_Tflush == True:
        # initialize figure
        fig, ax = plt.subplots(1,1,figsize = (6,6))
        # format figure
        plt.suptitle(year + ' monthly mean '+r' DO$_{in}$-DO$_{deep}$ vs. T$_{flush}$',
                        size=16)
        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='both', labelsize=14)
        ax.set_xlabel(r'Monthly mean T$_{flush}$ [days]', fontsize=14)
        ax.set_ylabel(r'Monthly mean DO$_{in}$-DO$_{deep}$ [mg/L]', fontsize=14)
        # plot
        ax.scatter(mean_Tflush,mean_DOin-deep_lay_DO,s=70,zorder=5,color='k',alpha=0.3)
        # get correlation
        r,p = pearsonr(mean_Tflush,mean_DOin-deep_lay_DO)
        ax.text(0.6, 0.15, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)),color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        ax.text(0.7, 0.09, r'$p =$' + str(round(p,8)) ,color='black',
                                verticalalignment='bottom', horizontalalignment='left',
                                transform=ax.transAxes, fontsize=12, fontweight='bold')
        plt.tight_layout()
        # save figure
        plt.show()

    # results!!!
    results_scatter = False
    # COLLAPSE
    if results_scatter == True:

        # # DOin vs. Tflush colored by mean deep layer DO
        # # initialize figure
        # fig, ax = plt.subplots(1,1,figsize = (5,5))
        # # format figure
        # plt.suptitle(year + ' monthly mean \ninflowing TEF DO vs. Tflush\ncolored by mean deep layer DO',
        #                 size=14)
        # # format grid
        # ax.set_facecolor('#EEEEEE')
        # ax.tick_params(axis='x', labelrotation=30)
        # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        # for border in ['top','right','bottom','left']:
        #     ax.spines[border].set_visible(False)
        # ax.tick_params(axis='y', labelsize=12)
        # ax.set_xlabel('Monthly mean Tflush = Vinlet/Qin [days]', fontsize=12)
        # ax.set_ylabel('Monthly mean inflowing TEF DO [mg/L]', fontsize=12)
        # # plot
        # cmap_oxy = plt.cm.get_cmap('rainbow_r')
        # cs_DO = ax.scatter(mean_Tflush,mean_DOin,s=30,zorder=5, c=deep_lay_DO, cmap=cmap_oxy)
        # # get correlation
        # r,p = pearsonr(mean_Tflush,mean_DOin)
        # ax.text(0.6, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)),color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # ax.text(0.7, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # # create colorbarlegend
        # cbar = fig.colorbar(cs_DO)
        # cbar.ax.tick_params(labelsize=12)
        # cbar.ax.set_ylabel('Monthly mean deep layer DO [mg/L]', rotation=90, fontsize=12)
        # cbar.outline.set_visible(False)
        # plt.tight_layout()
        # # save figure
        # plt.show()

        # # DOin vs. Tflush colored by mean deep layer DO (circles to show seasonal cycle) -------------------
        # # initialize figure
        # fig, ax = plt.subplots(1,1,figsize = (10,5))
        # # format figure
        # plt.suptitle(year + ' monthly mean'+r' TEF DO$_{in}$ vs. T$_{flush}$'+' colored by mean deep layer DO',
        #                 size=16)
        # # format grid
        # ax.set_facecolor('#EEEEEE')
        # ax.tick_params(axis='x', labelrotation=30)
        # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        # for border in ['top','right','bottom','left']:
        #     ax.spines[border].set_visible(False)
        # ax.tick_params(axis='both', labelsize=14)
        # ax.set_xlabel(r'Monthly mean T$_{flush}$ [days]', fontsize=14)
        # ax.set_ylabel(r'Monthly mean TEF DO$_{in}$ [mg/L]', fontsize=14)
        # for i,station in enumerate(sta_dict):
        #     # plot
        #     ax.plot(np.append(mean_Tflush[i*intervals:(i+1)*intervals],mean_Tflush[i*intervals]),
        #             np.append(mean_DOin[i*intervals:(i+1)*intervals],mean_DOin[i*intervals]),
        #             linewidth=1,color='k',alpha=0.5)
        # cmap_oxy = plt.cm.get_cmap('rainbow_r')
        # # cs_DO = ax.scatter(mean_Tflush,mean_DOin,s=45,zorder=5, edgecolor='k',c=deep_lay_DO, cmap=cmap_oxy)
        # cs_DO = ax.scatter(mean_Tflush,mean_DOin,s=45,zorder=5, edgecolor='k',c=perc_hyp_vol, cmap=cmap_hyp)
        # # create colorbarlegend
        # cbar = fig.colorbar(cs_DO)
        # cbar.ax.tick_params(labelsize=12)
        # # cbar.ax.set_ylabel('Monthly mean deep layer DO [mg/L]', rotation=90, fontsize=14)
        # cbar.ax.set_ylabel('Monthly mean % hypoxic volume', rotation=90, fontsize=14)
        # cbar.outline.set_visible(False)
        # # plt.xscale('log')
        # plt.tight_layout()
        # # save figure
        # plt.show()

        # # deltaDO vs. Tflush colored by mean deep layer DO (circles to show seasonal cycle) -------------------
        # # initialize figure
        # fig, ax = plt.subplots(1,1,figsize = (6,6))
        # # format figure
        # plt.suptitle(year + ' monthly mean'+r' DO$_{in}$-DO$_{deep}$ vs. T$_{flush}$'+' colored by mean deep layer DO',
        #                 size=16)
        # # format grid
        # ax.set_facecolor('#EEEEEE')
        # ax.tick_params(axis='x', labelrotation=30)
        # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        # for border in ['top','right','bottom','left']:
        #     ax.spines[border].set_visible(False)
        # ax.tick_params(axis='both', labelsize=14)
        # ax.set_xlabel(r'Monthly mean T$_{flush}$ [days]', fontsize=14)
        # ax.set_ylabel(r'Monthly mean DO$_{in}$-DO$_{deep}$ [mg/L]', fontsize=14)
        # for i,station in enumerate(sta_dict):
        #     # plot
        #     ax.plot(np.append(mean_Tflush[i*intervals:(i+1)*intervals],mean_Tflush[i*intervals]),
        #             np.append(mean_DOin[i*intervals:(i+1)*intervals],mean_DOin[i*intervals])-
        #             np.append(deep_lay_DO[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals]))#,
        #             # linewidth=1,color='k',alpha=0.5)
        # # cmap_oxy = plt.cm.get_cmap('rainbow_r')
        # # # cs_DO = ax.scatter(mean_Tflush,mean_DOin,s=45,zorder=5, edgecolor='k',c=deep_lay_DO, cmap=cmap_oxy)
        # # cs_DO = ax.scatter(mean_Tflush,mean_DOin-deep_lay_DO,s=45,zorder=5, edgecolor='k',c=perc_hyp_vol, cmap=cmap_hyp)
        # # # create colorbarlegend
        # # cbar = fig.colorbar(cs_DO)
        # # cbar.ax.tick_params(labelsize=12)
        # # # cbar.ax.set_ylabel('Monthly mean deep layer DO [mg/L]', rotation=90, fontsize=14)
        # # cbar.ax.set_ylabel('Monthly mean % hypoxic volume', rotation=90, fontsize=14)
        # # cbar.outline.set_visible(False)
        # # # plt.xscale('log')
        # plt.tight_layout()
        # # save figure
        # plt.show()

        # # bio consumption vs. recirculation colored by mean deep layer DO
        # # initialize figure
        # fig, ax = plt.subplots(1,1,figsize = (5,5))
        # # format figure
        # plt.suptitle(year + ' monthly mean \nbio consumption vs. recirculation\ncolored by mean deep layer DO',
        #                 size=14)
        # # format grid
        # ax.set_facecolor('#EEEEEE')
        # ax.tick_params(axis='x', labelrotation=30)
        # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        # for border in ['top','right','bottom','left']:
        #     ax.spines[border].set_visible(False)
        # ax.tick_params(axis='y', labelsize=12)
        # ax.set_xlabel('Monthly mean bio consumption [mg/L per day]', fontsize=12)
        # ax.set_ylabel('Monthly mean recirculation [mg/L per day]', fontsize=12)
        # # plot
        # cmap_oxy = plt.cm.get_cmap('rainbow_r')
        # cs_DO = ax.scatter(mean_cons,mean_recirc,s=30,zorder=5, c=deep_lay_DO, cmap=cmap_oxy)
        # # get correlation
        # r,p = pearsonr(mean_cons,mean_recirc)
        # ax.text(0.6, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)),color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # ax.text(0.7, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # # create colorbarlegend
        # cbar = fig.colorbar(cs_DO)
        # cbar.ax.tick_params(labelsize=12)
        # cbar.ax.set_ylabel('Monthly mean deep layer DO [mg/L]', rotation=90, fontsize=12)
        # cbar.outline.set_visible(False)
        # plt.tight_layout()
        # # save figure
        # plt.show()

        # # mean deep layer DO vs. DO in colored by Tflush (each inlet) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # for i,station in enumerate(sta_dict):
        #     plt.close('all')
        #     # initialize figure
        #     fig, ax = plt.subplots(1,1,figsize = (8,8))
        #     # format figure
        #     # plt.suptitle(year + ' monthly mean \ndeep layer DO vs. TEF ' + r'DO$_{in}$' + '\ncolored by mean flushing time',
        #     #                 size=18)
        #     plt.suptitle(year + ' monthly mean \ndeep layer DO vs. ' + r'DO$_{in}$', size=20)
        #     # format grid
        #     ax.set_facecolor('#EEEEEE')
        #     ax.tick_params(axis='x', labelrotation=30)
        #     ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        #     for border in ['top','right','bottom','left']:
        #         ax.spines[border].set_visible(False)
        #     ax.tick_params(axis='both', labelsize=16)
        #     ax.set_xlabel(r'Monthly mean TEF DO$_{in}$ [mg/L]', fontsize=18)
        #     ax.set_ylabel('Monthly mean deep layer DO [mg/L]', fontsize=18)
        #     # plot
        #     cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
        #     cmap_oxy = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
        #     ax.plot([0,11],[0,11],color='gray')
        #     ax.text(0.8,0.8,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4, fontsize=18)
        #     # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
        #     ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, color='gray',alpha=0.5, edgecolor='none')
        #     cs = ax.scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
        #                     s=150, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_oxy,
        #                     linewidth=2, vmin=0, vmax=80)
        #     ax.text(0.1, 0.85, station ,color='black',
        #                             verticalalignment='bottom', horizontalalignment='left',
        #                             transform=ax.transAxes, fontsize=16, fontweight='bold')
        #     # get correlation
        #     # r,p = pearsonr(mean_DOin,deep_lay_DO)
        #     # ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)) ,color='black',
        #     #                         verticalalignment='bottom', horizontalalignment='left',
        #     #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        #     # ax.text(0.1, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
        #     #                         verticalalignment='bottom', horizontalalignment='left',
        #     #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        #     # create colorbarlegend
        #     cbar = fig.colorbar(cs)
        #     cbar.ax.tick_params(labelsize=16)
        #     cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=18)
        #     cbar.outline.set_visible(False)
        #     ax.set_xlim([0,11])
        #     ax.set_ylim([0,11])
        #     plt.tight_layout()
        #     # save figure
        #     # plt.show()
        #     out_dir_DOinTflush = out_dir / 'DOin_Tflush_scatter'
        #     Lfun.make_dir(out_dir_DOinTflush)
        #     plt.savefig(out_dir_DOinTflush / (station + '.png') )

        # # mean deep layer DO vs. DO in colored by Tflush (only color the hypoxic inlets)
        # # initialize figure
        # fig, ax = plt.subplots(1,1,figsize = (8,8))
        # # format figure
        # # plt.suptitle(year + ' monthly mean \ndeep layer DO vs. TEF ' + r'DO$_{in}$' + '\ncolored by mean flushing time',
        # #                 size=18)
        # plt.suptitle(year + ' monthly mean \ndeep layer DO vs. ' + r'DO$_{in}$', size=20)
        # # format grid
        # ax.set_facecolor('#EEEEEE')
        # ax.tick_params(axis='x', labelrotation=30)
        # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        # for border in ['top','right','bottom','left']:
        #     ax.spines[border].set_visible(False)
        # ax.tick_params(axis='both', labelsize=16)
        # ax.set_xlabel(r'Monthly mean TEF DO$_{in}$ [mg/L]', fontsize=18)
        # ax.set_ylabel('Monthly mean deep layer DO [mg/L]', fontsize=18)
        # # plot
        # cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
        # cmap_oxy = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
        # ax.plot([0,11],[0,11],color='gray')
        # ax.text(0.8,0.8,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4, fontsize=18)
        # for i,station in enumerate(sta_dict):
        #     # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
        #     if station in ['case','penn','portsusan','holmes','lynchcove','dabob']:
        #         cs = ax.scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
        #                     s=50, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_oxy,
        #                     linewidth=2, vmin=0, vmax=80)
        #     # ax.text(0.1, 0.85, station ,color='black',
        #     #                         verticalalignment='bottom', horizontalalignment='left',
        #     #                         transform=ax.transAxes, fontsize=16, fontweight='bold')
        #     # get correlation
        #     # r,p = pearsonr(mean_DOin,deep_lay_DO)
        #     # ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)) ,color='black',
        #     #                         verticalalignment='bottom', horizontalalignment='left',
        #     #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        #     # ax.text(0.1, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
        #     #                         verticalalignment='bottom', horizontalalignment='left',
        #     #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        #     # create colorbarlegend
        # ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, color='gray',alpha=0.5, edgecolor='none')
        # ax.text(0.1, 0.79, 'hypoxic inlets colored' ,color='black',
        #         verticalalignment='bottom', horizontalalignment='left',
        #         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # cbar = fig.colorbar(cs)
        # cbar.ax.tick_params(labelsize=16)
        # cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=18)
        # cbar.outline.set_visible(False)
        # ax.set_xlim([0,11])
        # ax.set_ylim([0,11])
        # plt.tight_layout()
        # # save figure
        # plt.show()

        # mean deep layer DO vs. DO in colored by Tflush !!!!!!!!!!!!!!!!!!!!(MONEY PLOT)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # initialize figure
        fig, ax = plt.subplots(1,1,figsize = (8,8))
        # format figure
        plt.suptitle(year + ' monthly mean \ndeep layer DO vs. TEF ' + r'DO$_{in}$' + '\ncolored by mean flushing time',
                        size=16)
        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='both', labelsize=16)
        ax.set_xlabel(r'Monthly mean TEF DO$_{in}$ [mg/L]', fontsize=16)
        ax.set_ylabel('Monthly mean deep layer DO [mg/L]', fontsize=16)
        # plot
        cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
        cmap_oxy = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
        ax.plot([0,12],[0,12],color='gray')
        ax.text(0.8,0.8,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4, fontsize=14)
        cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
        # create colorbarlegend
        cbar = fig.colorbar(cs)
        cbar.ax.tick_params(labelsize=14)
        cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=16)
        cbar.outline.set_visible(False)
        ax.set_xlim([0,12])
        ax.set_ylim([0,12])
        plt.tight_layout()
        # save figure
        plt.show()

        # # mean deep layer DO vs. Tflush colored by DO in
        # # initialize figure
        # fig, ax = plt.subplots(1,1,figsize = (6,6))
        # # format figure
        # plt.suptitle(year + ' monthly mean \ndeep layer DO vs. mean flushing time\ncolored by inflowing TEF DO',
        #                 size=14)
        # # format grid
        # ax.set_facecolor('#EEEEEE')
        # ax.tick_params(axis='x', labelrotation=30)
        # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        # for border in ['top','right','bottom','left']:
        #     ax.spines[border].set_visible(False)
        # ax.tick_params(axis='y', labelsize=12)
        # ax.set_xlabel('Monthly mean Tflush = Vinlet/Qin [days]', fontsize=12)
        # ax.set_ylabel('Monthly mean deep layer DO [mg/L]', fontsize=12)
        # # plot
        # cmap_oxy = plt.cm.get_cmap('rainbow_r')
        # cs = ax.scatter(mean_Tflush,deep_lay_DO,s=20, zorder=5, c=mean_DOin, cmap=cmap_oxy)
        # # get correlation
        # r,p = pearsonr(mean_Tflush,deep_lay_DO)
        # ax.text(0.5, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)) ,color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # ax.text(0.5, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # # create colorbarlegend
        # cbar = fig.colorbar(cs)
        # cbar.ax.tick_params(labelsize=12)
        # cbar.ax.set_ylabel('Monthly mean inflowing TEF DO [mg/L]', rotation=90, fontsize=12)
        # cbar.outline.set_visible(False)
        # # ax.set_xlim([0,11])
        # ax.set_ylim([0,11])
        # plt.tight_layout()
        # # save figure
        # plt.show()


        # # min inflowing DO vs. annual mean inflowing DO colored by basin
        # # initialize figure
        # fig, ax = plt.subplots(1,1,figsize = (5,5))
        # # format figure
        # plt.suptitle(year + ' annual min DOin vs. annual mean DOin',
        #                 size=14)
        # # format grid
        # ax.set_facecolor('#EEEEEE')
        # ax.tick_params(axis='x', labelrotation=30)
        # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        # for border in ['top','right','bottom','left']:
        #     ax.spines[border].set_visible(False)
        # ax.tick_params(axis='y', labelsize=12)
        # ax.set_xlabel('Annual mean inflowing DO [mg/L]', fontsize=12)
        # ax.set_ylabel('Annual min inflowing DO [mg/L]', fontsize=12)
        # # plot
        # ax.scatter(annmean_DOin,annmin_DOin,alpha=0.5,s=60,zorder=5,
        #         color=colors_twentyone)
        # # ax.set_xlim([0,8])
        # # ax.set_ylim([0,8])
        # # calculate correlation coefficient (Pearson)
        # r,p = pearsonr(annmean_DOin,annmin_DOin)
        # ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) + r'; $r^2 =$' + str(round(r**2,2)) ,color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # ax.text(0.1, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
        #                         verticalalignment='bottom', horizontalalignment='left',
        #                         transform=ax.transAxes, fontsize=12, fontweight='bold')
        # # create legend
        # whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        # hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        # mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        # southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        # admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        # ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='lower right')
        # # save figure
        # plt.show()

##########################################################
##                Seasonal cycle analysis               ## 
##########################################################

seasonal_cycle = True
if seasonal_cycle == True:

    # initialize figure
    fig = plt.figure(figsize = (12,9))

    # crescent bay
    ax = fig.add_subplot(2, 2, 1)
    ax.set_title('(a) Crescent Bay 2017 monthly mean \n' + r'DO$_{deep}$ vs. DO$_{in}$ colored by T$_{flush}$', loc='left', size=14)
    # format grid
    ax.set_facecolor('#EEEEEE')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    ax.set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
    # plot
    cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
    cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
    ax.plot([0,11],[0,11],color='dimgray')
    ax.text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4, fontsize=10)
    # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
    ax.scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, color='gray',alpha=0.3, edgecolor='none')
    for i,station in enumerate(sta_dict):
        if station == 'crescent':
            cs = ax.scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
                            s=150, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_tflush,
                        linewidth=2, vmin=0, vmax=40)
        else:
            continue
    # create colorbarlegend
    cbar = fig.colorbar(cs)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
    cbar.outline.set_visible(False)
    ax.set_xlim([0,11])
    ax.set_ylim([0,11])

    # lynch cove
    ax = fig.add_subplot(2, 2, 3)
    ax.set_title('(b) Lynch Cove 2017 monthly mean \n' + r'DO$_{deep}$ vs. DO$_{in}$ colored by T$_{flush}$', loc='left', size=14)
    # format grid
    ax.set_facecolor('#EEEEEE')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    ax.set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
    # plot
    cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
    cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
    ax.plot([0,11],[0,11],color='dimgray')
    ax.text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4, fontsize=10)
    # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
    ax.scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, color='gray',alpha=0.3, edgecolor='none')
    for i,station in enumerate(sta_dict):
        if station == 'lynchcove':
            cs = ax.scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
                            s=150, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_tflush,
                        linewidth=2, vmin=0, vmax=40)
        else:
            continue
    # create colorbarlegend
    cbar = fig.colorbar(cs)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
    cbar.outline.set_visible(False)
    ax.set_xlim([0,11])
    ax.set_ylim([0,11])


    # spaghetti
    ax = fig.add_subplot(1, 2, 2)
    ax.set_title(r'(c) 2017 monthly mean $\Delta DO$ vs. T$_{flush}$',
                    size=14, loc='left')
    # format grid
    ax.set_facecolor('#EEEEEE')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_xlabel(r'Monthly mean T$_{flush}$ [days]', fontsize=12)
    ax.set_ylabel(r'Monthly mean $\Delta DO$ [mg/L]', fontsize=12)
    hypoxic_counter = 0
    oxygenated_counter = 0
    for i,station in enumerate(sta_dict):
        # plot
        zorder=4
        outline = 'black'
        if station in ['penn','holmes','portsusan','dabob','case']:
            if hypoxic_counter == 0:
                ax.scatter(mean_Tflush[i*intervals:(i+1)*intervals],
                       mean_DOin[i*intervals:(i+1)*intervals] - deep_lay_DO[i*intervals:(i+1)*intervals],
                       marker='o', s=100, alpha=0.5, zorder=5, c='hotpink', edgecolor='none', label='Hypoxic')
                hypoxic_counter += 1
            else:
                ax.scatter(mean_Tflush[i*intervals:(i+1)*intervals],
                    mean_DOin[i*intervals:(i+1)*intervals] - deep_lay_DO[i*intervals:(i+1)*intervals],
                    marker='o', s=100, alpha=0.5, zorder=5, c='hotpink', edgecolor='none')
        else:
            if oxygenated_counter == 0:
                ax.scatter(mean_Tflush[i*intervals:(i+1)*intervals],
                       mean_DOin[i*intervals:(i+1)*intervals] - deep_lay_DO[i*intervals:(i+1)*intervals],
                       marker='o', s=100, alpha=0.5, zorder=5, c='cornflowerblue', edgecolor='none', label='Oxygenated')
                oxygenated_counter += 1
            else:
                ax.scatter(mean_Tflush[i*intervals:(i+1)*intervals],
                       mean_DOin[i*intervals:(i+1)*intervals] - deep_lay_DO[i*intervals:(i+1)*intervals],
                       marker='o', s=100, alpha=0.5, zorder=5, c='cornflowerblue', edgecolor='none')
        if station == 'lynchcove':
            ax.scatter(mean_Tflush[i*intervals:(i+1)*intervals],
                       mean_DOin[i*intervals:(i+1)*intervals] - deep_lay_DO[i*intervals:(i+1)*intervals],
                       marker='s', s=150, zorder=6, c='deeppink', edgecolor='black', linewidth=2,
                       label='Lynch Cove')
        if station == 'crescent':
            ax.scatter(mean_Tflush[i*intervals:(i+1)*intervals],
                       mean_DOin[i*intervals:(i+1)*intervals] - deep_lay_DO[i*intervals:(i+1)*intervals],
                       marker='s', s=150, zorder=6, c='royalblue', edgecolor='black', linewidth=2,
                       label='Crescent Bay')
        
        ax.plot([0,85],[0,0],color='black',zorder=4,linewidth=0.5)
        ax.spines['bottom'].set_position('zero')
        ax.set_xlim([0,85])
    ax.legend(loc='lower right')

    plt.tight_layout()
    # save figure
    plt.show()

##########################################################
##                  Deep DO time series               ## 
##########################################################

# initialize figure
fig, ax = plt.subplots(1,1,figsize = (12,5))
# format figure
plt.suptitle(year + r' DO$_{deep}$ time series [mg/L] (30-day Hanning Window)',
                size=16)
# format grid
ax.set_facecolor('#EEEEEE')
ax.tick_params(axis='x', labelrotation=30)
loc = mdates.MonthLocator(interval=1)
ax.xaxis.set_major_locator(loc)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)
ax.tick_params(axis='both', labelsize=12)

# loop through stations
for i,station in enumerate(sta_dict):

    # get average deep layer DO
    deep_lay_DO = DOconcen_dict[station]['Deep Layer']

    # 30-day hanning windwo
    deep_lay_DO = zfun.lowpass(deep_lay_DO.values,n=30)

    # plot DO colored by flushing time
    # ax.plot(dates_local_daily,deep_lay_DO,linewidth=3.5,color='black', alpha=0.8)
    # ax.plot(dates_local_daily,deep_lay_DO,linewidth=1.7,color='silver')

    ax.plot(dates_local_daily,deep_lay_DO,linewidth=1,color='black',alpha=0.5)

# format labels
ax.set_xlim([dates_local[0],dates_local[-2]])
ax.set_ylim([0,12])
ax.set_ylabel(r'DO$_{deep}$ [mg/L]',fontsize=12)

plt.show()

##########################################################
##           Deep DO and Tflush time series             ## 
##########################################################

timeseries = False
if timeseries == True:
    for i,station in enumerate(sta_dict):
        plt.close('all')
        # initialize figure
        fig, ax = plt.subplots(1,1,figsize = (8,5.5))
        ax.set_title(station,fontsize=18)

        # format grid
        ax.set_facecolor('#EEEEEE')
        ax.tick_params(axis='x', labelrotation=30)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.tick_params(axis='both', labelsize=12)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        ax.set_xlim([dates_no_crop[0],dates_no_crop[-1]])
        # ax.set_ylim([0,12])
        ax.set_ylabel('Deep layer DO [mg/L]', fontsize=12)

        # get variables for plotting
        deep_layer_DO = DOconcen_dict[station]['Deep Layer']
        DOin = DOconcen_dict[station]['Qin DO']
        Tflush = dimensions_dict[station]['Inlet volume'][0] / bottomlay_dict[station]['Qin m3/s'] / (60*60*24)
        Qin = bottomlay_dict[station]['Qin m3/s']

        # filter so I can see what's going on
        DOin = zfun.lowpass(DOin.values,n=30)
        Qin = zfun.lowpass(Qin.values,n=30)
        Tflush = zfun.lowpass(Tflush.values,n=30)

        # get timeseries of monthly intervals
        days = [dates_local_daily[i] for i in [13,43,72,103,133,164,194,225,256,286,317,347]]

        # plot
        deepDO_color = 'black'
        ax.plot(dates_local_daily,deep_layer_DO,color=deepDO_color,linewidth=2,alpha=0.5,
                linestyle='-',label='Deep Layer DO')
        ax.scatter(days,deep_lay_DO[i*intervals:(i+1)*intervals],s=80,facecolor=deepDO_color,edgecolor='white',zorder=7)
        # ax.plot(dates_local_daily,DOin,color='deeppink',linewidth=2,linestyle='--',label=r'TEF DO$_{in}$ (30-day Hanning Window)')
        ax.legend(loc='upper left',fontsize=12)
        ax.set_ylim([0,13])
        ticks = 5

        # add Tflush
        Tflush_color = 'deeppink'
        ax2 = ax.twinx()
        for border in ['top','right','bottom','left']:
            ax2.spines[border].set_visible(False)
        ax2.plot(dates_local_daily,Tflush,linewidth=2,linestyle='-',color=Tflush_color,
                 alpha=0.5,label=r'T$_{flush}$ (30-day Hanning Window)')
        # add scatter points of monthly mean
        ax2.scatter(days,mean_Tflush[i*intervals:(i+1)*intervals],s=80,facecolor=Tflush_color,edgecolor='white',zorder=7)
        # ax2.plot(dates_local_daily,Qin,linewidth=2,linestyle='-',color='darkturquoise',label='Qin (30-day Hanning Window)')
        # ax2.set_ylim([0,np.nanmax(Tflush)*4])
        # ax2.set_ylim([0,np.nanmax(Qin)*4])
        ax2.set_ylabel(r'T$_{flush}$ [days]',color=Tflush_color,fontsize=12)
        ax2.legend(loc='upper right',fontsize=12)
        ymax = 1.2*np.nanmax(Tflush)
        # round to nearest multiple of 5
        ymax = int(5 * math.ceil(float(ymax)/5))
        ax2.set_ylim([0,ymax])
        ax2.set_yticks(np.arange(0, ymax+(ymax/ticks), ymax/ticks))
        ax2.tick_params(axis='y', labelcolor=Tflush_color, labelsize=12)

        out_dir_DOinTflush_ts = out_dir / 'DOin_Tflush_timeseries'
        Lfun.make_dir(out_dir_DOinTflush_ts)
        plt.savefig(out_dir_DOinTflush_ts / (station + '.png') )


    ##########################################################
    ##                 Scatter matrix                       ##
    ##########################################################

    scatter_matrix = False
    # COLLAPSE
    if scatter_matrix == True:


        fig, axes = plt.subplots(1,1,figsize = (8,8))

        # put all variables in pandas dataframe
        d = {'Deep layer\nDO [mg/L]': deep_lay_DO,
            'DOin\n[mg/L]': mean_DOin,
            'Tflush\n[days]': mean_Tflush,
            # 'TEF DO\n[mg/L/day]': mean_TEFin,
            'Recirculation\n[mg/L/day]': mean_recirc,
            'Consumption\n[mg/L/day]':  mean_cons,
            'Depth\n[m]': mean_depth,
            'Volume\n[m3]': inlet_vol,
            'Aspect ratio\n(L/W)': aspect_ratio
            }
        df = pd.DataFrame(data=d)

        # create colormap
        cmap_oxy = plt.cm.get_cmap('rainbow_r')

        pd.plotting.scatter_matrix(df[['Deep layer\nDO [mg/L]',
                                    'DOin\n[mg/L]',
                                    'Tflush\n[days]',
                                    # 'TEF DO\n[mg/L/day]',
                                    'Recirculation\n[mg/L/day]',
                                    'Consumption\n[mg/L/day]',
                                    'Depth\n[m]',
                                    'Volume\n[m3]',
                                    'Aspect ratio\n(L/W)']],
                                figsize=(15,15), range_padding=0.2,
                                ax = axes,
                                hist_kwds={'ec':'gray',
                                            'facecolor':'silver',
                                            'alpha':0.5,
                                            'bins':8},
                                c=deep_lay_DO, marker='.',s=30, #edgecolors='k'
                                alpha=1,cmap=cmap_oxy)

        plt.suptitle('Scatter matrix of monthly mean variables')

        # add colorbar
        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.91, 0.12, 0.01, 0.74])
        cbar = fig.colorbar(cs_DO, cax=cbar_ax)
        cbar.outline.set_visible(False)
        cbar.ax.set_ylabel('Monthly mean deep layer DO [mg/L]', rotation=90, fontsize=12)
        plt.subplots_adjust(left=0.05, top=0.95, bottom=0.07, right=0.9)

        plt.show()


    ##########################################################
    ##                Multiple regression (2)               ## 
    ##########################################################

    MLR = False
    if MLR == True:
        # split into shallower and deeper groups
        cutoff =  14.5 #14.5
        shallower = [1 if x <= cutoff else np.nan for x in mean_depth]
        deeper = [1 if x > cutoff else np.nan for x in mean_depth]
        # drop nans
        mean_DOin_shallower = [x for x in mean_DOin*shallower if str(x) != 'nan']
        mean_DOin_deeper = [x for x in mean_DOin*deeper if str(x) != 'nan']
        mean_Tflush_shallower = [x for x in mean_Tflush*shallower if str(x) != 'nan']
        mean_Tflush_deeper = [x for x in mean_Tflush*deeper if str(x) != 'nan']
        deep_lay_DO_shallower = [x for x in deep_lay_DO*shallower if str(x) != 'nan']
        deep_lay_DO_deeper = [x for x in deep_lay_DO*deeper if str(x) != 'nan']
        # create array of predictors
        input_array_shallow = np.array([mean_DOin_shallower, mean_Tflush_shallower, [1]*len(mean_DOin_shallower)]).T
        # input_array_shallow = np.array([mean_DOin_shallower, [0]*len(mean_DOin_shallower)]).T
        input_array_deep = np.array([mean_DOin_deeper, mean_Tflush_deeper, [1]*len(mean_DOin_deeper)]).T
        # calculate MLR
        B_shallow,_,_,_ = lstsq(input_array_shallow,deep_lay_DO_shallower)
        slope_DOin_shallow = B_shallow[0]
        intercept_shallow = B_shallow[1]
        # slope_Tflush_shallow = 0
        slope_Tflush_shallow = B_shallow[1]
        intercept_shallow = B_shallow[2]
        print(B_shallow)
        B_deep,_,_,_ = lstsq(input_array_deep,deep_lay_DO_deeper)
        slope_DOin_deep = B_deep[0]
        slope_Tflush_deep = B_deep[1]
        intercept_deep = B_deep[2]


        # # Time series -------------------------------------------------------------------------------------
        # # initialize figure
        # if year == '2017':
        #     # fig, ax = plt.subplots(7,3,figsize = (15,9),sharex=True,sharey=True)
        #     fig, ax = plt.subplots(5,4,figsize = (15,9),sharex=True,sharey=True)
        # elif year == '2014':
        #     fig, ax = plt.subplots(5,1,figsize = (12,5),sharex=True,sharey=True)
        # ax = ax.ravel()
        # equation_shallow = 'shallower than {} m: Deep DO [mg/L] = ({}) * DOin + ({} mg/L/day) *Tflush + ({} mg/L)'.format(str(cutoff),
        #                 round(slope_DOin_shallow,2),round(slope_Tflush_shallow,2),round(intercept_shallow,2))
        # equation_deep = 'deeper than {} m: Deep DO [mg/L] = ({}) * DOin + ({} mg/L/day) *Tflush + ({} mg/L)'.format(str(cutoff),
        #                 round(slope_DOin_deep,2),round(slope_Tflush_deep,2),round(intercept_deep,2))
        
        # # format figure
        # if year == '2017':
        #     plt.suptitle(year + ' deep layer DO time series\n' + equation_shallow +'\n' + equation_deep, size=10)
        # elif year == '2014':
        #     plt.suptitle(year + ' deep layer DO time series [mg/L]', size=14)


        # # loop through test stations
        # for i,station in enumerate(sta_dict):

        #     # format grid
        #     ax[i].set_facecolor('#EEEEEE')
        #     ax[i].tick_params(axis='x', labelrotation=30)
        #     ax[i].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        #     for border in ['top','right','bottom','left']:
        #         ax[i].spines[border].set_visible(False)
        #     ax[i].tick_params(axis='y', labelsize=12)
        #     ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        #     ax[i].set_xlim([dates_local_daily[0],dates_local_daily[-1]])
        #     ax[i].set_ylim([0,12])
        #     ax[i].text(0.02,0.9,station,ha='left', va='top', transform=ax[i].transAxes, fontsize=12,fontweight='bold')

        #     # get inflowing DO and flushing time
        #     DOin = DOconcen_dict[station]['Qin DO'] # mg/L
        #     Tflush = dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'] / (60*60*24) # days
        #     depth = dimensions_dict[station]['Mean depth'][0]
        #     aspect_ratio = dimensions_dict[station]['L/W aspect ratio'][0]
        #     deep_cons = bottomlay_dict[station]['Bio Consumption']

        #     # plot LiveOcean output
        #     ax[i].plot(dates_local_daily,DOconcen_dict[station]['Deep Layer'], color=basin_color_dict[basin_dict[station]],
        #             alpha=0.8,linewidth=3,label=(station + ' raw LO output'))
        #     # calculate and plot prediction
        #     if year == '2017':
        #         if depth <= cutoff:
        #             # prediction = slope_DOin_shallow*DOin + intercept_shallow
        #             prediction = slope_DOin_shallow*DOin + slope_Tflush_shallow*Tflush + intercept_shallow

        #         else:
        #             prediction = slope_DOin_deep*DOin + slope_Tflush_deep*Tflush + intercept_deep
        #     elif year == '2014':
        #         prediction = 0.9535885895820221*DOin + -0.0701628377978838*Tflush + 0.4763992129208265
        #     ax[i].plot(dates_local_daily,prediction,linestyle='-', color='k',#color=basin_color_dict[basin_dict[station]],
        #             linewidth=1,label=(station + ' regression prediction'))
            
        #     # add legend
        #     if i == 0:
        #         whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        #         hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        #         mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        #         southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        #         admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        #         ax[i].legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],
        #                     ncol=5,bbox_to_anchor=(3.38, 1.22))
        # plt.subplots_adjust(hspace=0.02,wspace=0.02,left=0.05, top=0.89, bottom=0.05, right=0.95)


        # # Residuals -------------------------------------------------------------------------------------
        # # initialize figure
        # if year == '2017':
        #     # fig, ax = plt.subplots(7,3,figsize = (15,9),sharex=True,sharey=True)
        #     fig, ax = plt.subplots(5,4,figsize = (15,9),sharex=True,sharey=True)
        # elif year == '2014':
        #     fig, ax = plt.subplots(5,1,figsize = (12,5),sharex=True,sharey=True)
        # ax = ax.ravel()

        # # format figure
        # if year == '2017':
        #     plt.suptitle(year + ' Residuals (LO-MLR)\n' + equation_shallow +'\n' + equation_deep, size=10)
        # elif year == '2014':
        #     plt.suptitle(year + ' Residuals (LO-MLR)', size=14)

        # # loop through test stations
        # for i,station in enumerate(sta_dict):

        #     # format grid
        #     ax[i].set_facecolor('#EEEEEE')
        #     ax[i].tick_params(axis='x', labelrotation=30)
        #     ax[i].grid(True,color='w',linewidth=1,linestyle='-',axis='x')
        #     for border in ['top','right','bottom','left']:
        #         ax[i].spines[border].set_visible(False)
        #     ax[i].tick_params(axis='y', labelsize=12)
        #     ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        #     ax[i].set_xlim([dates_local_daily[0],dates_local_daily[-1]])
        #     ax[i].set_ylim([-8,8])
        #     ax[i].axhline(0,0,1,color='gray')
        #     ax[i].text(0.02,0.9,station,ha='left', va='top', transform=ax[i].transAxes, fontsize=12,fontweight='bold')

        #     # get inflowing DO and flushing time
        #     DOin = DOconcen_dict[station]['Qin DO'] # mg/L
        #     Tflush = dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'] / (60*60*24) # days
        #     depth = dimensions_dict[station]['Mean depth'][0]
        #     aspect_ratio = dimensions_dict[station]['L/W aspect ratio'][0]
        #     deep_cons = bottomlay_dict[station]['Bio Consumption']

        #     # calculate and plot prediction
        #     if year == '2017':
        #         if depth <= cutoff:
        #             prediction = slope_DOin_shallow*DOin + slope_Tflush_shallow*Tflush + intercept_shallow
        #         else:
        #             prediction = slope_DOin_deep*DOin + slope_Tflush_deep*Tflush + intercept_deep
        #     elif year == '2014':
        #         prediction = 0.9535885895820221*DOin + -0.0701628377978838*Tflush + 0.4763992129208265
        #     # plot residuals
        #     ax[i].plot(dates_local_daily,DOconcen_dict[station]['Deep Layer']-prediction, color=basin_color_dict[basin_dict[station]],
        #             alpha=0.8,linewidth=1.5)
            
        #     # add legend
        #     if i == 0:
        #         whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        #         hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        #         mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        #         southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        #         admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        #         ax[i].legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],
        #                     ncol=5,bbox_to_anchor=(3.38, 1.22))
        # plt.subplots_adjust(hspace=0.02,wspace=0.02,left=0.05, top=0.89, bottom=0.05, right=0.95)

        # compare Aug/Sep mean deep layer DO -------------------------------------------------------------------------------------
        deep_lay_DO_summer_vals = np.zeros(len(sta_dict))
        prediction_summer_vals = np.zeros(len(sta_dict))
        fig, ax = plt.subplots(1,1,figsize = (8,8))
        # format grid
        ax.set_facecolor('#EEEEEE')
        # ax.tick_params(axis='x', labelrotation=30)
        ax.tick_params(axis='both', labelsize=16)
        ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        ax.set_xlabel('LiveOcean mean deep layer DO [mg/L]', fontsize=16)
        ax.set_ylabel('MLR predicted mean deep layer DO [mg/L]', fontsize=16)
        ax.set_ylim([0,11])
        ax.set_xlim([0,11])
        # ax.set_title('Aug/Sep mean deep layer DO: MLR prediction vs. LO\n'+ equation_shallow +'\n' + equation_deep, size=16)
        ax.set_title('Aug & Sep mean deep layer DO:\nMLR prediction vs. LiveOcean', size=16)
        # loop through stations
        for i,station in enumerate(sta_dict):
            # get data
            deep_lay_DO = DOconcen_dict[station]['Deep Layer']
            DOin = DOconcen_dict[station]['Qin DO'] # mg/L
            Tflush = dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'] / (60*60*24) # days
            depth = dimensions_dict[station]['Mean depth'][0]
            if depth <= cutoff:
                prediction = slope_DOin_shallow*DOin + slope_Tflush_shallow*Tflush + intercept_shallow
                color = 'teal'
            else:
                prediction = slope_DOin_deep*DOin + slope_Tflush_deep*Tflush + intercept_deep
                color = 'mediumorchid'
            # crop to Aug/Sep
            deep_lay_DO_summer = np.nanmean(deep_lay_DO[211:272])
            prediction_summer = np.nanmean(prediction[211:272])
            deep_lay_DO_summer_vals[i] = deep_lay_DO_summer
            prediction_summer_vals[i] = prediction_summer
            # plot
            ax.scatter(deep_lay_DO_summer,prediction_summer,s=120,alpha=0.5,zorder=5,
                    color=color)
                    # color=basin_color_dict[basin_dict[station]])
        # add unity line
        ax.plot([0,11],[0,11],color='gray')
        ax.text(8.8,8.8,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4,fontsize=16)
        # # create legend
        # whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
        # hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
        # mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
        # southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
        # admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
        # ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='lower right')
        # bias and rmse
        bias = np.nanmean(prediction_summer_vals-deep_lay_DO_summer_vals)
        rmse = np.sqrt(np.nanmean((prediction_summer_vals-deep_lay_DO_summer_vals)**2))
        ax.text(0.05, 0.85, 'bias: {}'.format(str(round(bias,2))),
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes, fontsize=16, color = 'k')
        ax.text(0.05, 0.8, 'rmse: {}'.format(str(round(rmse,2))),
            verticalalignment='top', horizontalalignment='left',
            transform=ax.transAxes, fontsize=16, color = 'k')

   
   
   
   
   
   
   
   
   
   
   
   
   
   
    # ##########################################################
    # ##               Multiple regression (1)                ## 
    # ##########################################################

    # # create array of predictors


    # input_array = np.array([mean_DOin, mean_Tflush, [1]*len(mean_DOin)]).T
    
    # # input_array = np.array([mean_DOin, [1]*len(mean_DOin)]).T
    # # input_array = np.array([mean_DOin, mean_Tflush, mean_depth, [1]*len(mean_DOin)]).T
    # # input_array = np.array([mean_DOin, mean_depth, aspect_ratio, mean_Tflush, mean_cons, [1]*len(mean_DOin)]).T
    # # calculate coefficients using multiple linear regression

    # B,a,b,c = lstsq(input_array,deep_lay_DO)
    # slope_DOin = B[0]
    # slope_Tflush = B[1]
    # intercept = B[2]
    # # slope_DOin = B[0]
    # # intercept = B[1]
    # # slope_DOin = B[0]
    # # slope_meandepth = B[1]
    # # slope_aspect = B[2]
    # # slope_Tflush = B[3]
    # # slope_cons = B[4]
    # # intercept = B[5]
    # # slope_DOin = B[0]
    # # slope_Tflush = B[1]
    # # slope_meandepth = B[2]
    # # intercept = B[3]

    # # print('\nMean deep layer DO [mg/L] = {}*DOin + {}*Tflush + {}'.format(
    # #     round(slope_DOin,2),round(slope_Tflush,2),round(intercept,2)))
    # # print('\nMean deep layer DO [mg/L] = {}*DOin + {}'.format(
    # #     round(slope_DOin,2),round(intercept,2)))
    # # print('\nMean deep layer DO [mg/L] = {}*DOin + {}*depth + {}'.format(
    # #     round(slope_DOin,2),round(slope_meandepth,2),round(intercept,2)))
    # # print('\nMean deep layer DO [mg/L] = {}*DOin + {}*depth + {}*aspect ratio + {}*Tflush + {}*consumption + {}'.format(
    # #     round(slope_DOin,2),round(slope_meandepth,2),round(slope_aspect,2),round(slope_Tflush,2),round(slope_cons,2),round(intercept,2)))


    # # Time series -------------------------------------------------------------------------------------
    # # initialize figure
    # if year == '2017':
    #     fig, ax = plt.subplots(7,3,figsize = (15,9),sharex=True,sharey=True)
    # elif year == '2014':
    #     fig, ax = plt.subplots(5,1,figsize = (12,5),sharex=True,sharey=True)
    # ax = ax.ravel()
    # equation = 'Mean deep layer DO [mg/L] = ({}) * DOin + ({} mg/L/day) *Tflush + ({} mg/L)'.format(
    #                 round(slope_DOin,2),round(slope_Tflush,2),round(intercept,2))
    # # equation = 'Mean deep layer DO [mg/L] = ({}) * DOin + ({} mg/L)'.format(
    # #                 round(slope_DOin,2),round(intercept,2))
    # # equation = 'Mean deep layer DO [mg/L] = ({}) * DOin + ({} mg/L/day) *Tflush + ({} mg/L/m) *depth + ({} mg/L)'.format(
    # #                 round(slope_DOin,2),round(slope_Tflush,2),round(slope_meandepth,2),round(intercept,2))
    
    # # format figure
    # if year == '2017':
    #     plt.suptitle(year + ' deep layer DO time series\n' + equation, size=14)
    # elif year == '2014':
    #     plt.suptitle(year + ' deep layer DO time series [mg/L]', size=14)
    # # plt.suptitle(year + ' deep layer DO time series' + 
    # #              '\nMean deep layer DO [mg/L] = ({}) * DOin + ({} mg/L/day) *Tflush + ({} mg/L)'.format(
    # #             slope_DOin,slope_Tflush,intercept), size=14)

    # # loop through test stations
    # for i,station in enumerate(sta_dict):

    #     # format grid
    #     ax[i].set_facecolor('#EEEEEE')
    #     ax[i].tick_params(axis='x', labelrotation=30)
    #     ax[i].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    #     for border in ['top','right','bottom','left']:
    #         ax[i].spines[border].set_visible(False)
    #     ax[i].tick_params(axis='y', labelsize=12)
    #     # ax[i].set_xlabel(year, fontsize=12)
    #     # ax[i].set_ylabel(r'Deep layer DO [mg/L]', fontsize=12)
    #     # ax.set_ylabel(r'Bottom $\sigma$-layer DO [mg/L]', fontsize=12)
    #     ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    #     ax[i].set_xlim([dates_local_daily[0],dates_local_daily[-1]])
    #     ax[i].set_ylim([0,12])
    #     ax[i].text(0.02,0.9,station,ha='left', va='top', transform=ax[i].transAxes, fontsize=12,fontweight='bold')

    #     # get inflowing DO and flushing time
    #     DOin = DOconcen_dict[station]['Qin DO'] # mg/L
    #     Tflush = dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'] / (60*60*24) # days
    #     depth = dimensions_dict[station]['Mean depth'][0]
    #     aspect_ratio = dimensions_dict[station]['L/W aspect ratio'][0]
    #     deep_cons = bottomlay_dict[station]['Bio Consumption']

    #     # plot LiveOcean output
    #     ax[i].plot(dates_local_daily,DOconcen_dict[station]['Deep Layer'], color=basin_color_dict[basin_dict[station]],
    #             alpha=0.8,linewidth=3,label=(station + ' raw LO output'))
    #     # calculate and plot prediction
    #     if year == '2017':
    #         prediction = slope_DOin*DOin + slope_Tflush*Tflush + intercept
    #         # prediction = slope_DOin*DOin + intercept
    #         # prediction = slope_DOin*DOin + slope_Tflush*Tflush + slope_meandepth*depth + intercept
    #     elif year == '2014':
    #         prediction = 0.9535885895820221*DOin + -0.0701628377978838*Tflush + 0.4763992129208265
    #     # prediction = slope_DOin*DOin + intercept
    #     # prediction = slope_DOin*DOin + slope_meandepth*depth + + slope_aspect*aspect_ratio + slope_Tflush*Tflush + slope_cons*deep_cons + intercept
    #     ax[i].plot(dates_local_daily,prediction,linestyle='-', color='k',#color=basin_color_dict[basin_dict[station]],
    #             linewidth=1,label=(station + ' regression prediction'))
        
    #     # add legend
    #     if i == 0:
    #         whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
    #         hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
    #         mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
    #         southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
    #         admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
    #         ax[i].legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],
    #                     ncol=5,bbox_to_anchor=(2.38, 1.4))
    # plt.subplots_adjust(hspace=0.02,wspace=0.02,left=0.05, top=0.89, bottom=0.05, right=0.95)


    # # Residuals -------------------------------------------------------------------------------------
    # # initialize figure
    # if year == '2017':
    #     fig, ax = plt.subplots(7,3,figsize = (15,9),sharex=True,sharey=True)
    # elif year == '2014':
    #     fig, ax = plt.subplots(5,1,figsize = (12,5),sharex=True,sharey=True)
    # ax = ax.ravel()

    # # format figure
    # if year == '2017':
    #     plt.suptitle(year + ' Residuals (LO-MLR)\n' + equation, size=14)
    # elif year == '2014':
    #     plt.suptitle(year + ' Residuals (LO-MLR)', size=14)
    # # plt.suptitle(year + ' deep layer DO time series' + 
    # #              '\nMean deep layer DO [mg/L] = ({}) * DOin + ({} mg/L/day) *Tflush + ({} mg/L)'.format(
    # #             slope_DOin,slope_Tflush,intercept), size=14)

    # # loop through test stations
    # for i,station in enumerate(sta_dict):

    #     # format grid
    #     ax[i].set_facecolor('#EEEEEE')
    #     ax[i].tick_params(axis='x', labelrotation=30)
    #     ax[i].grid(True,color='w',linewidth=1,linestyle='-',axis='x')
    #     for border in ['top','right','bottom','left']:
    #         ax[i].spines[border].set_visible(False)
    #     ax[i].tick_params(axis='y', labelsize=12)
    #     # ax[i].set_xlabel(year, fontsize=12)
    #     # ax[i].set_ylabel(r'Deep layer DO [mg/L]', fontsize=12)
    #     # ax.set_ylabel(r'Bottom $\sigma$-layer DO [mg/L]', fontsize=12)
    #     ax[i].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    #     ax[i].set_xlim([dates_local_daily[0],dates_local_daily[-1]])
    #     ax[i].set_ylim([-8,8])
    #     ax[i].axhline(0,0,1,color='gray')
    #     ax[i].text(0.02,0.9,station,ha='left', va='top', transform=ax[i].transAxes, fontsize=12,fontweight='bold')

    #     # get inflowing DO and flushing time
    #     DOin = DOconcen_dict[station]['Qin DO'] # mg/L
    #     Tflush = dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'] / (60*60*24) # days
    #     depth = dimensions_dict[station]['Mean depth'][0]
    #     aspect_ratio = dimensions_dict[station]['L/W aspect ratio'][0]
    #     deep_cons = bottomlay_dict[station]['Bio Consumption']

    #     # calculate and plot prediction
    #     if year == '2017':
    #         prediction = slope_DOin*DOin + slope_Tflush*Tflush + intercept
    #         # prediction = slope_DOin*DOin + intercept
    #         # prediction = slope_DOin*DOin + slope_Tflush*Tflush + slope_meandepth*depth + intercept
    #     elif year == '2014':
    #         prediction = 0.9535885895820221*DOin + -0.0701628377978838*Tflush + 0.4763992129208265
    #     # prediction = slope_DOin*DOin + intercept
    #     # prediction = slope_DOin*DOin + slope_meandepth*depth + + slope_aspect*aspect_ratio + slope_Tflush*Tflush + slope_cons*deep_cons + intercept
    #     # plot residuals
    #     ax[i].plot(dates_local_daily,DOconcen_dict[station]['Deep Layer']-prediction, color=basin_color_dict[basin_dict[station]],
    #             alpha=0.8,linewidth=1.5)
        
    #     # add legend
    #     if i == 0:
    #         whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
    #         hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
    #         mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
    #         southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
    #         admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
    #         ax[i].legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],
    #                     ncol=5,bbox_to_anchor=(2.38, 1.4))
    # plt.subplots_adjust(hspace=0.02,wspace=0.02,left=0.05, top=0.89, bottom=0.05, right=0.95)

    # # compare Aug/Sep mean deep layer DO -------------------------------------------------------------------------------------
    # deep_lay_DO_summer_vals = np.zeros(len(sta_dict))
    # prediction_summer_vals = np.zeros(len(sta_dict))
    # fig, ax = plt.subplots(1,1,figsize = (5,5))
    # # format grid
    # ax.set_facecolor('#EEEEEE')
    # ax.tick_params(axis='x', labelrotation=30)
    # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    # for border in ['top','right','bottom','left']:
    #     ax.spines[border].set_visible(False)
    # ax.tick_params(axis='y', labelsize=12)
    # ax.set_xlabel('LiveOcean mean deep layer DO [mg/L]', fontsize=12)
    # ax.set_ylabel('MLR predicted mean deep layer DO [mg/L]', fontsize=12)
    # ax.set_ylim([0,11])
    # ax.set_xlim([0,11])
    # ax.set_title('Aug/Sep mean deep layer DO: MLR prediction vs. LO\n'+equation)
    # # loop through stations
    # for i,station in enumerate(sta_dict):
    #     # get data
    #     deep_lay_DO = DOconcen_dict[station]['Deep Layer']
    #     DOin = DOconcen_dict[station]['Qin DO'] # mg/L
    #     Tflush = dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'] / (60*60*24) # days
    #     depth = dimensions_dict[station]['Mean depth'][0]
    #     prediction = slope_DOin*DOin + slope_Tflush*Tflush + intercept
    #     # prediction = slope_DOin*DOin + intercept
    #     # prediction = slope_DOin*DOin + slope_Tflush*Tflush + slope_meandepth*depth + intercept
    #     # crop to Aug/Sep
    #     deep_lay_DO_summer = np.nanmean(deep_lay_DO[211:272])
    #     prediction_summer = np.nanmean(prediction[211:272])
    #     deep_lay_DO_summer_vals[i] = deep_lay_DO_summer
    #     prediction_summer_vals[i] = prediction_summer
    #     # plot
    #     ax.scatter(deep_lay_DO_summer,prediction_summer,s=120,alpha=0.6,zorder=5,
    #             color=basin_color_dict[basin_dict[station]])
    # # add unity line
    # ax.plot([0,11],[0,11],color='gray')
    # ax.text(8.8,8.8,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4)
    # # create legend
    # whidbey = mpatches.Patch(color='limegreen', label='Whidbey Basin', alpha=0.5)
    # hoodcanal = mpatches.Patch(color='hotpink', label='Hood Canal', alpha=0.5)
    # mainbasin = mpatches.Patch(color='deepskyblue', label='Main Basin', alpha=0.5)
    # southsound = mpatches.Patch(color='blueviolet', label='South Sound', alpha=0.5)
    # admiralty = mpatches.Patch(color='black', label='Admiralty Inlet', alpha=0.5)
    # ax.legend(handles=[whidbey,hoodcanal,mainbasin,southsound,admiralty],loc='lower right')
    # # bias and rmse
    # bias = np.nanmean(prediction_summer_vals-deep_lay_DO_summer_vals)
    # rmse = np.sqrt(np.nanmean((prediction_summer_vals-deep_lay_DO_summer_vals)**2))
    # ax.text(0.05, 0.95, 'bias: {}'.format(str(round(bias,2))),
    #     verticalalignment='top', horizontalalignment='left',
    #     transform=ax.transAxes, fontsize=12, color = 'k')
    # ax.text(0.05, 0.9, 'rmse: {}'.format(str(round(rmse,2))),
    #     verticalalignment='top', horizontalalignment='left',
    #     transform=ax.transAxes, fontsize=12, color = 'k')
















# # Annual mean DO and summer DO ------------------------------------------------
# # initialize figure
# plt.close('all')
# fig, ax = plt.subplots(1,1,figsize = (5,5))
# # format figure
# plt.suptitle('Aug/Sep mean deep DO vs. annual mean inflowing DO', size=12)
# # format grid
# ax.set_facecolor('#EEEEEE')
# ax.tick_params(axis='x', labelrotation=30)
# ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
# for border in ['top','right','bottom','left']:
#     ax.spines[border].set_visible(False)
# ax.tick_params(axis='y', labelsize=12)
# ax.set_xlabel('Mean annual inflowing DO [mg/L]')
# ax.set_ylabel('Mean Aug/Sep deep layer DO [mg/L]')
# ax.set_xlim([0,11])
# ax.set_ylim([0,11])
# # plot
# ax.scatter(mean_DOin,deep_lay_DO,alpha=0.5,s=100,zorder=5,
#         color=colors)
# # calculate correlation coefficient (Pearson)
# r,p = pearsonr(mean_DOin,deep_lay_DO)
# ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=ax.transAxes, fontsize=12, fontweight='bold')
# ax.text(0.1, 0.79, r'$p =$' + str(round(p,3)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=ax.transAxes, fontsize=12, fontweight='bold')
# plt.show()

# # physical dimensions ------------------------------------------------
# # initialize figure
# plt.close('all')
# fig, ax = plt.subplots(1,1,figsize = (5,5))
# # format figure
# plt.suptitle('Aspect Ratio (L/W) vs. Mean Depth', size=12)
# # format grid
# ax.set_facecolor('#EEEEEE')
# ax.tick_params(axis='x', labelrotation=30)
# ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
# for border in ['top','right','bottom','left']:
#     ax.spines[border].set_visible(False)
# ax.tick_params(axis='y', labelsize=12)
# ax.set_xlabel('Mean Depth [m]')
# ax.set_ylabel('Aspect Ratio (L/W)')
# # plot
# ax.scatter(mean_depth,aspect_ratio,alpha=0.5,s=100,zorder=5,
#         color=colors)
# # calculate correlation coefficient (Pearson)
# r,p = pearsonr(mean_depth,aspect_ratio)
# ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=ax.transAxes, fontsize=12, fontweight='bold')
# ax.text(0.1, 0.79, r'$p =$' + str(round(p,3)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=ax.transAxes, fontsize=12, fontweight='bold')
# plt.show()

# ##########################################################
# ##                  Bottom DO time series               ## 
# ##########################################################

# # initialize figure
# # plt.close('all')
# fig, ax = plt.subplots(2,1,figsize = (10,8),sharex=True)
# # format figure
# plt.suptitle(year + ' deep layer DO time series [mg/L]\ncolored by annual mean flushing time',
#                 size=16)
# # format grid
# for axis in [ax[0],ax[1]]:
#     axis.set_facecolor('#EEEEEE')
#     axis.tick_params(axis='x', labelrotation=30)
#     axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
#     axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#     for border in ['top','right','bottom','left']:
#         axis.spines[border].set_visible(False)
#     axis.tick_params(axis='y', labelsize=12)

# ax[0].set_title('(a) Average deep-layer DO',
#                 fontsize=12,loc='left')
# ax[1].set_title('(b) Average deep-layer DO with annual mean subtracted',
#                 fontsize=12,loc='left')

# n = 47#67#106#67
# # colors = plt.cm.bwr_r(np.linspace(0,1,n))
# # colors = plt.cm.rainbow_r(np.linspace(0,1,n))
# # colors = plt.cm.coolwarm_r(np.linspace(0,1,n))
# # colors = plt.cm.nipy_spectral(np.linspace(0,1,n))
# colors = plt.cm.magma_r(np.linspace(0,1,n))
# # colors = plt.cm.cubehelix_r(np.linspace(0,1,n))
        
# # loop through stations
# for i,station in enumerate(sta_dict):

#     # get average deep layer DO
#     deep_lay_DO = DOconcen_dict[station]['Deep Layer']

#     # get the Minimum Bottom Layer DO
#     mindeepDO = np.nanmin(DOconcen_dict[station]['Deep Layer'])
#     # print(round(mindeepDO))

#     # get annual mean deep layer DO
#     ann_mean_DO = np.nanmean(DOconcen_dict[station]['Deep Layer'])

#     # subtract annual mean
#     deep_lay_no_mean = deep_lay_DO - ann_mean_DO

#     # get mean depth
#     mean_depth = dimensions_dict[station]['Mean depth'].values[0]

#     # get volume
#     inlet_vol = dimensions_dict[station]['Inlet volume'].values[0]

#     # get aspect ratio
#     aspect_ratio = dimensions_dict[station]['L/W aspect ratio'].values[0]

#     # get mean annual flushing time
#     Tflush = np.nanmean(inlet_vol/bottomlay_dict[station]['Qin m3/s'] / (60*60*24))

#     # # 30 day hanning window
#     # deep_lay_DO = zfun.lowpass(deep_lay_DO.values,f='hanning',n=30)
#     # deep_lay_no_mean = zfun.lowpass(deep_lay_no_mean.values,f='hanning',n=90)

#     # plot DO colored by mean depth
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8)
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=colors[round(mean_depth)-1])
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8)
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=colors[round(mean_depth)-1])

#     # plot DO colored by flushing time
#     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8)
#     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=colors[round(Tflush)-1])
#     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8)
#     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=colors[round(Tflush)-1])

#     # # plot DO colored by inlet volume
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8)
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=colors[round(inlet_vol/10000000)-1])
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8)
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=colors[round(inlet_vol/10000000)-1])

#     # # plot DO colored by aspect ratio
#     # # if station == 'dyes':
#     # #     continue # dyes inlet aspect ratio is large, and likely an overestimate, given geometry
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8)
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=colors[round(aspect_ratio)-1])
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8)
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=colors[round(aspect_ratio)-1])

#     # plot raw deep layer DO colored by basin
#     # if basin_dict[station] == 'Hood Canal':
#         # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8,zorder=24)
#         # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=basin_color_dict[basin_dict[station]],zorder=25)#colors[round(mindeepDO)-1])
#         # # plot deep layer DO with annual mean subtracted out colored by min deep layer DO
#         # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8, zorder=24)
#         # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=basin_color_dict[basin_dict[station]],zorder=25)#colors[round(mindeepDO)-1])

#     # if station == 'dabob':
#     #     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8,zorder=24)
#     #     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=basin_color_dict[basin_dict[station]],zorder=25)#colors[round(mindeepDO)-1])
#     #     # plot deep layer DO with annual mean subtracted out colored by min deep layer DO
#     #     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8, zorder=24)
#     #     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=basin_color_dict[basin_dict[station]],zorder=25)#colors[round(mindeepDO)-1])
#     #     in_dir = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
#     #     bulk = xr.open_dataset(in_dir)
#     #     tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
#     #     Qin =  tef_df['q_p'] # Qin [m3/s]
#     #     ax2=ax[1].twinx()
#     #     ax2.plot(dates_local_daily,Qin,linewidth=2,color='black',zorder=26, label='Qin')
#     #     ax2.set_ylim(0,8000)
#     #     ax2.set_ylabel('TEF inflow [m3/s]')
#     #     ax2.legend(loc='upper right', ncol=2)
#     #     for border in ['top','right','bottom','left']:
#     #         ax2.spines[border].set_visible(False)

#         # if station == 'penn':
#         #     # add Skagit river hydrograph
#         #     ds = xr.open_dataset('2017rivflow.nc')
#         #     flow_snohomish = ds.sel(riv='skagit')['transport'].values
#         #     flow_skagit = ds.sel(riv='snohomish')['transport'].values
#         #     ax2=ax[1].twinx()
#         #     ax2.plot(dates_local_daily,flow_skagit[1:-1],color='black',zorder=30,linestyle=':', label='Skagit')
#         #     ax2.plot(dates_local_daily,flow_snohomish[1:-1],color='black',zorder=30, label='Snohomish')
#         #     ax2.set_ylim(0,6000)
#         #     ax2.set_ylabel('Hydrograph [m3/s]')
#         #     ax2.legend(loc='upper right', ncol=2)
#         #     for border in ['top','right','bottom','left']:
#         #         ax2.spines[border].set_visible(False)

#     # else:
#     #     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color=basin_color_dict[basin_dict[station]], alpha=0.3)
#     #     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color=basin_color_dict[basin_dict[station]], alpha=0.3)

#     # format labels
#     ax[0].set_xlim([dates_local[0],dates_local[-1]])
#     # ax.set_ylim([-5,5])

# # ax.legend(loc='upper right', ncol=5)

# plt.show()