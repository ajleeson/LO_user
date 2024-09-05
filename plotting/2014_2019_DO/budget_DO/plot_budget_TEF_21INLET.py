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
import math
import matplotlib.patches as patches
import csv
import cmocean
from scipy.stats import pearsonr
import matplotlib.pylab as plt
import gsw
import pickle
import get_two_layer

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

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
else:
    sta_dict = stations
    
# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
Lfun.make_dir(out_dir)

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]
# crop time vector (because we only have jan 2 - dec 30)
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
        ds = xr.open_dataset('../../../../LO_output/extract/'+gtagex+
                            '/tef2/extractions_'+startdate+'_'+enddate+
                            '/'+station+'.nc')
        # calculate average bottom oxygen throughout section
        bottom_oxygen_unfiltered = np.nanmean(ds['oxygen'].values[:,0,:],axis=1) * 32/1000 # mmol/m3 to mg/L (time,z,p)
        bottom_oxygen = zfun.lowpass(bottom_oxygen_unfiltered, f='godin')[36:-34:24]
        # calculate minimum DO throughout section
        oxygen_min_unfiltered = np.nanmin(np.nanmin(ds['oxygen'].values,axis=1),axis=1) * 32/1000
        oxygen_min = zfun.lowpass(oxygen_min_unfiltered, f='godin')[36:-34:24]

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
        bottomlay_dict[station]['TRAPS'] = traps_deep
        bottomlay_dict[station]['Photosynthesis'] = photo_deep
        bottomlay_dict[station]['Bio Consumption'] = cons_deep
        bottomlay_dict[station]['EU Vertical'] = vertX_deep_EU
        bottomlay_dict[station]['TEF Vertical'] = vertX_deep_TEF
        bottomlay_dict[station]['Storage'] = ddtDOV_deep
        bottomlay_dict[station]['Volume'] = deep_V

# ------------------------- save DO concentrations in dataframe dict -----------------------------------

        # mg/L units

        # DO concentrations
        DOconcen_dict[station]['Surface Layer'] = o2_surf
        DOconcen_dict[station]['Deep Layer'] = o2_deep
        DOconcen_dict[station]['Bottom Sigma DO'] = bottom_oxygen
        DOconcen_dict[station]['Minimum Deep Layer DO'] = oxygen_min
        DOconcen_dict[station]['Qin DO'] = DO_in

# ------------------------ save inlet dimensions in dataframe dict ----------------------------------------

        dimensions_dict[station]['Inlet volume'] = [inlet_vol] # m3
        dimensions_dict[station]['Mean depth'] = [mean_depth] # m
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

residual = False # recirculation (combined horizontal and vertical exchange)
show_EU = True

# COLLAPSE
if DO_budget == True:

    print('Making DO budget time series')

    for i,station in enumerate(sta_dict):
            
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
            ax[3].plot(dates_local_daily,DOconcen_dict[station]['Minimum Deep Layer DO'],color='crimson',
                       linewidth=1,linestyle='-',label='Minimum deep layer DO')
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
            min_DO = np.nanmin(DOconcen_dict[station]['Minimum Deep Layer DO'][minday:maxday])
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
                if attribute in ['EU Exchange Flow', 'TEF Exchange Flow', 'EU Recirculation', 'EU Vertical', 'TEF Vertical', 'Volume', 'Storage']:
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
##                   DO scatterplots                    ## 
##########################################################

# initialize arrays for plotting
deep_lay_DO = np.zeros(len(sta_dict))
bott_sig_DO = np.zeros(len(sta_dict))
annual_mean_DO = np.zeros(len(sta_dict))
mean_DOin = np.zeros(len(sta_dict))
mean_depth = np.zeros(len(sta_dict))
inlet_vol = np.zeros(len(sta_dict))
aspect_ratio = np.zeros(len(sta_dict))
colors = []

# get values for plotting and calculating r value
for i,station in enumerate(sta_dict):
    # get minday and maxday for aug 1 sep 30:
    minday = 212
    maxday = 243
    # save values
    deep_lay_DO[i] =  np.nanmean(DOconcen_dict[station]['Deep Layer'][minday:maxday])
    bott_sig_DO[i] =  np.nanmean(DOconcen_dict[station]['Bottom Sigma DO'][minday:maxday])
    annual_mean_DO[i] = np.nanmean(DOconcen_dict[station]['Deep Layer'])
    mean_DOin[i] = np.nanmean(DOconcen_dict[station]['Qin DO'])
    mean_depth[i] = dimensions_dict[station]['Mean depth'][0]
    inlet_vol[i] = dimensions_dict[station]['Inlet volume'][0]
    aspect_ratio[i] = dimensions_dict[station]['L/W aspect ratio'][0]
    colors.append(basin_color_dict[basin_dict[station]])
    
meanDO_bottsigDO = False
# COLLLAPSE
if meanDO_bottsigDO == True:
    # initialize figure
    plt.close('all')
    fig, ax = plt.subplots(1,1,figsize = (5,5))
    # format figure
    plt.suptitle(year + ' Aug/Sep \n' + r'Mean bottom $\sigma$ layer DO vs. mean deep layer DO',
                    size=12)
    # format grid
    ax.set_facecolor('#EEEEEE')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    ax.tick_params(axis='y', labelsize=12)
    ax.set_xlabel('Mean deep layer DO [mg/L]')
    ax.set_ylabel(r'Mean bottom $\sigma$ layer DO [mg/L]')
    # plot
    ax.plot([0,11],[0,11],color='gray')
    ax.scatter(deep_lay_DO,bott_sig_DO,alpha=0.5,s=100,zorder=5,
            color=colors)
    ax.set_xlim([0,11])
    ax.set_ylim([0,11])
    # calculate correlation coefficient (Pearson)
    r,p = pearsonr(deep_lay_DO,bott_sig_DO)
    ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
                            verticalalignment='bottom', horizontalalignment='left',
                            transform=ax.transAxes, fontsize=12, fontweight='bold')
    ax.text(0.1, 0.79, r'$p =$' + str(round(p,8)) ,color='black',
                            verticalalignment='bottom', horizontalalignment='left',
                            transform=ax.transAxes, fontsize=12, fontweight='bold')
    plt.savefig(out_dir / 'scatter_meanDO_bottsigDO.png' )


# Annual mean DO and summer DO ------------------------------------------------
# initialize figure
plt.close('all')
fig, ax = plt.subplots(1,1,figsize = (5,5))
# format figure
plt.suptitle('Aug/Sep mean deep DO vs. annual mean inflowing DO', size=12)
# format grid
ax.set_facecolor('#EEEEEE')
ax.tick_params(axis='x', labelrotation=30)
ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)
ax.tick_params(axis='y', labelsize=12)
ax.set_xlabel('Mean annual inflowing DO [mg/L]')
ax.set_ylabel('Mean Aug/Sep deep layer DO [mg/L]')
ax.set_xlim([0,11])
ax.set_ylim([0,11])
# plot
ax.scatter(mean_DOin,deep_lay_DO,alpha=0.5,s=100,zorder=5,
        color=colors)
# calculate correlation coefficient (Pearson)
r,p = pearsonr(mean_DOin,deep_lay_DO)
ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=12, fontweight='bold')
ax.text(0.1, 0.79, r'$p =$' + str(round(p,3)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=12, fontweight='bold')
plt.show()

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
# plt.close('all')
# fig, ax = plt.subplots(2,1,figsize = (10,8),sharex=True)
# # format figure
# plt.suptitle(year + ' deep layer DO time series [mg/L]\ncolored by aspect ration (L/W)',
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

# n = 67#106#67
# # colors = plt.cm.bwr_r(np.linspace(0,1,n))
# # colors = plt.cm.rainbow_r(np.linspace(0,1,n))
# # colors = plt.cm.coolwarm_r(np.linspace(0,1,n))
# # colors = plt.cm.nipy_spectral(np.linspace(0,1,n))
# colors = plt.cm.magma_r(np.linspace(0,1,n))
        
# # loop through stations
# for i,station in enumerate(sta_dict):

#     # get average deep layer DO
#     deep_lay_DO = DOconcen_dict[station]['Deep Layer']

#     # get the minimum deep layer DO
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

#     # # 30 day hanning window
#     # deep_lay_DO = zfun.lowpass(deep_lay_DO.values,f='hanning',n=30)
#     # deep_lay_no_mean = zfun.lowpass(deep_lay_no_mean.values,f='hanning',n=90)

#     # plot DO colored by mean depth
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8)
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=colors[round(mean_depth)-1])
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8)
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=colors[round(mean_depth)-1])

#     # # plot DO colored by inlet volume
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8)
#     # ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=colors[round(inlet_vol/10000000)-1])
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8)
#     # ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=colors[round(inlet_vol/10000000)-1])

#     # plot DO colored by aspect ratio
#     # if station == 'dyes':
#     #     continue # dyes inlet aspect ratio is large, and likely an overestimate, given geometry
#     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8)
#     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=colors[round(aspect_ratio)-1])
#     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8)
#     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=colors[round(aspect_ratio)-1])

#     # # plot raw deep layer DO colored by basin
#     # if basin_dict[station] == 'Admiralty Inlet':
#     #     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color='black', alpha=0.8,zorder=24)
#     #     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=2,color=basin_color_dict[basin_dict[station]],zorder=25)#colors[round(mindeepDO)-1])
#     #     # plot deep layer DO with annual mean subtracted out colored by min deep layer DO
#     #     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color='black', alpha=0.8, zorder=24)
#     #     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=2,color=basin_color_dict[basin_dict[station]],zorder=25)#colors[round(mindeepDO)-1])
#     # else:
#     #     ax[0].plot(dates_local_daily,deep_lay_DO,linewidth=3,color=basin_color_dict[basin_dict[station]], alpha=0.3)
#     #     ax[1].plot(dates_local_daily,deep_lay_no_mean,linewidth=3,color=basin_color_dict[basin_dict[station]], alpha=0.3)

#     # format labels
#     ax[0].set_xlim([dates_local[0],dates_local[-1]])
#     # ax.set_ylim([-5,5])

# # ax.legend(loc='upper right', ncol=5)

# plt.show()


# # # # ##########################################################
# # # # ##             Time series of storage term              ## ######################################################################
# # # # ##########################################################

# # # # initialize figure
# # # fig, ax = plt.subplots(5,1,figsize = (12,8),sharey=True,sharex=True)
# # # axis = ax.ravel()

# # # # loop through stations
# # # for i,station in enumerate(stations):

# # #     # format figure
# # #     plt.suptitle('Major Bottom Layer Rates\n' + r'[mg/L day$^{-1}$]',
# # #                  size=16)
# # #     # format grid
# # #     axis[i].set_facecolor('#EEEEEE')
# # #     axis[i].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
# # #     axis[i].tick_params(axis='x', labelrotation=30)
# # #     axis[i].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
# # #     for border in ['top','right','bottom','left']:
# # #         axis[i].spines[border].set_visible(False)
# # #     axis[i].tick_params(axis='y', labelsize=12)
# # #     axis[i].text(0.02,0.05,station,fontsize=12,fontweight='bold',ha='left', va='bottom', transform=axis[i].transAxes)

# # #     # get volume of bottom layer
# # #     fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
# # #     df_V = pd.read_pickle(fn)
# # #     # Godin filter already applied earlier in workflow
# # #     deep_V = df_V['deep [m3]'].values[1:-1]

# # #     # convert from kmol/s/m3 to  to mg/L per day
# # #     # calculate d/dt(DO)
# # #     storage = storage_deep_dict[station]['Storage']/deep_V * 1000 * (32) * (60*60*24)
# # #     # calculate major source terms
# # #     photo = bottomlay_dict[station]['Photosynthesis']/deep_V * 1000 * (32) * (60*60*24)
# # #     recirc = bottomlay_dict[station]['Recirculation']/deep_V * 1000 * (32) * (60*60*24)
# # #     wwtps = bottomlay_dict[station]['TRAPS']/deep_V * 1000 * (32) * (60*60*24)
# # #     # calculate major sink terms
# # #     cons = bottomlay_dict[station]['Bio Consumption']/deep_V * 1000 * (32) * (60*60*24)

# # #     # plot
# # #     axis[i].plot(dates_local_daily[1:-1],recirc,color='darkorange',linewidth=2,label='Recirculation')
# # #     axis[i].plot(dates_local_daily[1:-1],cons,color='deeppink',linewidth=2,label='Bio Consumption')
# # #     axis[i].plot(dates_local_daily[1:-1],photo,color='lightseagreen',linewidth=2,label='Photosynthesis')
# # #     axis[i].plot(dates_local_daily[1:-1],wwtps,color='mediumorchid',alpha=0.5,linewidth=3,label='WWTP inflow')
# # #     axis[i].plot(dates_local_daily[1:-1],storage,color='navy',linewidth=0.6,label='Storage')

# # #     # # get average bottom layer DO
# # #     # deep_lay_DO =DOconcen_dict[station]['Deep Layer']
# # #     # # bott_DO = DOconcen_dict[station]['Bottom DO']
# # #     # ax2 = axis[i].twinx()
# # #     # ax2.plot(dates_local_daily[1:-1],deep_lay_DO,color='limegreen',linewidth=2)
# # #     # # format grid
# # #     # ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
# # #     # ax2.tick_params(axis='x', labelrotation=30)
# # #     # for border in ['top','right','bottom','left']:
# # #     #     ax2.spines[border].set_visible(False)
# # #     # ax2.tick_params(axis='y', labelsize=12, labelcolor='limegreen')
# # #     # ax2.set_ylim([0,12])
# # #     # ax2.set_yticks(np.arange(0, 12, 3))

# # #     # format labels
# # #     axis[i].set_xlim([dates_local[0],dates_local[-1]])
# # #     axis[i].set_ylim([-1,1])

# # #     # legend
# # #     if i == 0:
# # #         axis[i].legend(loc='upper left',ncol=5)

# # # plt.subplots_adjust(hspace=0.1, top=0.9, bottom=0.05, right=0.95)
# # # plt.show()

# # ##########################################################
# # ##               Flushing time time series              ## 
# # ##########################################################

# # initialize figure
# fig, ax = plt.subplots(1,1,figsize = (10,6))
# # format figure
# plt.suptitle('Bottom layer flushing time time series [days]\n' +
#              r'$T_{flush} = V_{layer}/Q_{in,TEF}$; with 30-day Hanning filter',
#                 size=16)
# # format grid
# ax.set_facecolor('#EEEEEE')
# ax.tick_params(axis='x', labelrotation=30)
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
# ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
# for border in ['top','right','bottom','left']:
#     ax.spines[border].set_visible(False)
# ax.tick_params(axis='y', labelsize=12)
        
# # loop through stations
# for i,station in enumerate(sta_dict):

#     # get Qin to bottom layer (TEF)
#     in_dir = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
#     bulk = xr.open_dataset(in_dir)
#     tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
#     Qin = tef_df['q_p'].values # Qin [m3/s]

#     # calculate fushing time in days
#     Tflush = bottomlay_dict[station]['Volume']/Qin / (60*60*24) # convert seconds to days

#     # month average
#     Tflush_filtered = zfun.lowpass( Tflush.values,f='hanning',n=30)

#     # plot
#     # ax.plot(dates_local_daily[1:-1],Tflush,color=station_color[i],linewidth=2,label=station)
#     ax.plot(dates_local_daily,Tflush_filtered,linewidth=2,label=station)

#     # format labels
#     ax.set_xlim([dates_local[0],dates_local[-1]])
#     # ax.set_ylim([-5,5])

# ax.legend(loc='upper right')

# plt.show()