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
import math
import csv
import cmocean
import matplotlib.pylab as plt
import gsw
import pickle
import get_two_layer

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

residual = False
show_EU = False

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
jobname = 'twentyoneinlets'
startdate = '2014.01.01'
enddate = '2014.12.31'
enddate_hrly = '2015.01.01 00:00:00'
year = '2014' # for making a date label

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
Lfun.make_dir(out_dir)

# create time_vecotr
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]

print('\n')

##########################################################
##              Plot DO budget of every inlet           ##
##########################################################

stations = ['lynchcove','penn','case','carr','budd']

# create dictionaries with interface depths
interface_dict = dict()

# create dictionaries to save dataframe of budget terms for each station
bottomlay_dict = {'lynchcove': pd.DataFrame(),
                'penn': pd.DataFrame(),
                'budd': pd.DataFrame(),
                'carr': pd.DataFrame(),
                'case': pd.DataFrame()}
surfacelay_dict = {'lynchcove': pd.DataFrame(),
                'penn': pd.DataFrame(),
                'budd': pd.DataFrame(),
                'carr': pd.DataFrame(),
                'case': pd.DataFrame()}

for i,station in enumerate(stations): # enumerate(sta_dict):
        # print status
        print('({}/{}) Working on {}...'.format(i+1,len(sta_dict),station))

        # get interface depth from csv file
        with open('interface_depths.csv', 'r') as f:
            for line in f:
                inlet, interface_depth = line.strip().split(',')
                interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
        z_interface = float(interface_dict[station])

##########################################################
##                      One layer                       ##
##########################################################
        if math.isnan(z_interface):
            # one layer
            print('One layer....make code to deal with this case')
        
##########################################################
##                    Two layers                        ##
##########################################################
        else:

            # initialize figure
            plt.close('all')
            fig, ax = plt.subplots(3,1,figsize = (10,8),sharex=True)
            # format figure
            plt.suptitle(station + ': DO Budget (Godin Filter)',size=14)
            for axis in [ax[0],ax[1],ax[2]]:
                # axis.plot([dates_local[0],dates_local[-1]],[0,0],color='k')
                axis.set_xlim([dates_local[0],dates_local[-1]])
                axis.set_ylabel(r'DO transport [kmol O$_2$ s$^{-1}$]')
                axis.set_facecolor('#EEEEEE')
                axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
                axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
                axis.tick_params(axis='x', labelrotation=30)
                for border in ['top','right','bottom','left']:
                    axis.spines[border].set_visible(False)
            ax[2].set_xlabel(year)
            ax[0].text(0.02,0.9,'(a) Surface [shallower than {} m]'.format(-1*z_interface),
                            ha='left', va='top', transform=ax[0].transAxes, fontsize=12)
            ax[1].text(0.02,0.9,'(b) Bottom [deeper than {} m]'.format(-1*z_interface),
                            ha='left', va='top', transform=ax[1].transAxes, fontsize=12)
            ax[2].text(0.02,0.9,r'(c) Error [sum of surface and bottom vertical exchange]',
                            ha='left', va='top', transform=ax[2].transAxes, fontsize=12)
            
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
            in_dir = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'tef2' / 'bulk_2014.01.01_2014.12.31' / (station + '.nc')
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
            months = ['2014.01.01_2014.01.31',
                      '2014.02.01_2014.02.28',
                      '2014.03.01_2014.03.31',
                      '2014.04.01_2014.04.30',
                      '2014.05.01_2014.05.31',
                      '2014.06.01_2014.06.30',
                      '2014.07.01_2014.07.31',
                      '2014.08.01_2014.08.31',
                      '2014.09.01_2014.09.30',
                      '2014.10.01_2014.10.31',
                      '2014.11.01_2014.11.30',
                      '2014.12.01_2014.12.31',]
            
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
            # 10-day hanning window filter (10 days = 240 hours)
            # traps_surf = zfun.lowpass(rivers_surf_unfiltered, f='hanning', n=240)
            # traps_deep = zfun.lowpass(wwtps_deep_unfiltered, f='hanning', n=240)
            # Godin filter
            traps_surf = zfun.lowpass(rivers_surf_unfiltered, f='godin')[36:-34:24]
            traps_deep =  zfun.lowpass(wwtps_deep_unfiltered, f='godin')[36:-34:24]
            traps_color = 'black'

# ------------------------------- get vertical exchange ----------------------------------------

            # pick color
            vertX_color = 'mediumorchid'

            # calculate error term, and ascribe that the vertical exchange
            vertX_surf_TEF = ddtDOV_surf - (TEF_surf # ddtDOV_surf[1:179] - (TEF_surf[0:178]
                                        + photo_surf#[1:179]
                                        + cons_surf#[1:179]
                                        + airsea_surf#[1:179]
                                        + traps_surf)#[1:179])
            vertX_deep_TEF = ddtDOV_deep - (TEF_deep # ddtDOV_deep[1:179] - (TEF_deep[0:178]
                                        + photo_deep#[1:179]
                                        + cons_deep#[1:179]
                                        + traps_deep)#[1:179])
            
            # calculate error term, and ascribe that the vertical exchange
            vertX_surf_EU = ddtDOV_surf - (EU_surf # ddtDOV_surf[1:179] - (EU_surf[0:178]
                                        + photo_surf#[1:179]
                                        + cons_surf#[1:179]
                                        + airsea_surf#[1:179]
                                        + traps_surf)#[1:179])
            vertX_deep_EU = ddtDOV_deep - (EU_deep # ddtDOV_deep[1:179] - (EU_deep[0:178]
                                        + photo_deep#[1:179]
                                        + cons_deep#[1:179]
                                        + traps_deep)#[1:179])
        

# ---------------------------------- plot and save --------------------------------------------
            # plot surface
            ax[0].plot(dates_local_daily[1:-1],TEF_surf,color=exchange_color,
                       linewidth=1,label='TEF Exchange Flow')
            if show_EU:
                ax[0].plot(dates_local_daily[1:-1],EU_surf,color=exchange_color,
                        linewidth=3,alpha=0.2,label='EU Exchange Flow')
            ax[0].plot(dates_local_daily[1:-1],traps_surf,color=traps_color,
                       linewidth=3,zorder=2,label='TRAPS')
            ax[0].plot(dates_local_daily[1:-1],photo_surf,color=photo_color,
                       linewidth=2,label='Photosynthesis')
            ax[0].plot(dates_local_daily[1:-1],cons_surf,color=cons_color,
                       linewidth=2,linestyle=':',zorder=9,label='Bio Consumption')
            ax[0].plot(dates_local_daily[1:-1],airsea_surf,color=airsea_color,
                       linewidth=1,zorder=8,label='Air-Sea Transfer')
            ax[0].plot(dates_local_daily[1:-1],ddtDOV_surf,color=ddtDOV_color,
                       linewidth=1,zorder=7,#alpha=0.6,
                       label=r'$\frac{\mathrm{d}}{\mathrm{dt}}(\mathrm{DO}\cdot V)$')
            ax[0].plot(dates_local_daily[1:-1],vertX_surf_TEF,color=vertX_color,
                       linewidth=1,label='TEF Vertical')
            if show_EU:
                ax[0].plot(dates_local_daily[1:-1],vertX_surf_EU,color=vertX_color,
                        linewidth=3,alpha=0.2,label='EU Vertical')
            ax[0].legend(loc='lower center',ncol=4)
            
            # plot deep
            ax[1].plot(dates_local_daily[1:-1],TEF_deep,color=exchange_color,
                       linewidth=1)
            if show_EU:
                ax[1].plot(dates_local_daily[1:-1],EU_deep,color=exchange_color,
                        linewidth=3,alpha=0.2)
            ax[1].plot(dates_local_daily[1:-1],traps_deep,color=traps_color,
                       linewidth=3,zorder=2)
            ax[1].plot(dates_local_daily[1:-1],photo_deep,color=photo_color,
                       linewidth=2)
            ax[1].plot(dates_local_daily[1:-1],cons_deep,color=cons_color,
                       linewidth=2,linestyle=':',zorder=9)
            ax[1].plot(dates_local_daily[1:-1],ddtDOV_deep,color=ddtDOV_color,
                       linewidth=1,zorder=7)#,alpha=0.6)
            ax[1].plot(dates_local_daily[1:-1],vertX_deep_TEF,color=vertX_color,
                       linewidth=1)
            if show_EU:
                ax[1].plot(dates_local_daily[1:-1],vertX_deep_EU,color=vertX_color,
                        linewidth=3,alpha=0.2)

            # plot error
            error_TEF = vertX_surf_TEF+vertX_deep_TEF
            error_EU = vertX_surf_EU+vertX_deep_EU
            ax[2].plot(dates_local_daily[1:-1],error_TEF,color='crimson',
                       linewidth=2,label='TEF')
            if show_EU:
                ax[2].plot(dates_local_daily[1:-1],error_EU,color='k',
                        linewidth=1,linestyle='--',label='EU')
                ax[2].legend(loc='best')
            for axis in [ax[0],ax[1],ax[2]]:
                ylimval = np.nanmax(np.abs(TEF_surf))*1.1
                axis.set_ylim([-1*ylimval,ylimval])
                # axis.set_ylim([-0.15,0.15])

            plt.subplots_adjust(hspace = 0.02, top=0.93)
            plt.savefig(out_dir / (station+'.png'))

# ------------------------- save data in dataframe dict -----------------------------------
            # Note: everything is in units of kmol O2 /s

            # surface layer
            if residual == False:
                if show_EU:
                    surfacelay_dict[station]['EU Exchange Flow'] = EU_surf
                else:
                    surfacelay_dict[station]['TEF Exchange Flow'] = TEF_surf
            else:
                if show_EU:
                    surfacelay_dict[station]['Residual'] = EU_surf + vertX_surf_EU
                else:
                    surfacelay_dict[station]['Residual'] = TEF_surf + vertX_surf_TEF
            surfacelay_dict[station]['TRAPS'] = traps_surf
            surfacelay_dict[station]['Photosynthesis'] = photo_surf
            surfacelay_dict[station]['Bio Consumption'] = cons_surf
            surfacelay_dict[station]['Air-Sea Transfer'] = airsea_surf
            if residual == False:
                if show_EU:
                    surfacelay_dict[station]['EU Vertical'] = vertX_surf_EU
                else:
                    surfacelay_dict[station]['TEF Vertical'] = vertX_surf_TEF
            surfacelay_dict[station]['Storage'] = ddtDOV_surf

            # bottom layer
            if residual == False:
                if show_EU:
                    bottomlay_dict[station]['EU Exchange Flow'] = EU_deep
                else:
                    bottomlay_dict[station]['TEF Exchange Flow'] = TEF_deep
            else:
                if show_EU:
                    bottomlay_dict[station]['Residual'] = EU_deep + vertX_deep_EU
                else:
                    bottomlay_dict[station]['Residual'] = TEF_deep + vertX_deep_TEF
            bottomlay_dict[station]['TRAPS'] = traps_deep
            bottomlay_dict[station]['Photosynthesis'] = photo_deep
            bottomlay_dict[station]['Bio Consumption'] = cons_deep
            if residual == False:
                if show_EU:
                    bottomlay_dict[station]['EU Vertical'] = vertX_deep_EU
                else:
                    bottomlay_dict[station]['TEF Vertical'] = vertX_deep_TEF
            bottomlay_dict[station]['Storage'] = ddtDOV_deep


##########################################################
##                Create bar charts                     ##
##########################################################

# annual average summer bottom DO
stations_w_DO = ['Lynch Cove\n0.0 mg/L','Penn Cove\n1.8 mg/L',
                 'Case Inlet\n1.9 mg/L','Carr Inlet\n4.0 mg/L',
                 'Budd Inlet\n5.8 mg/L']

# initialize figure
plt.close('all')
fig, ax = plt.subplots(5,2,figsize = (13,10), sharey=True, sharex='col')
# format figure
if show_EU:
    exchange = 'Eulerian'
else:
    exchange = 'TEF'
plt.suptitle('Volume-averaged DO transport rates ({})'.format(exchange),size=14)
for axis in [ax[0,0],ax[0,1],ax[1,0],ax[1,1],ax[2,0],ax[2,1],ax[3,0],ax[3,1],ax[4,0],ax[4,1]]:
    axis.set_facecolor('#EEEEEE')
    axis.grid(True,color='w',linewidth=1,linestyle='-',axis='y')
    for border in ['top','right','bottom','left']:
        axis.spines[border].set_visible(False)
# label y-axis
ax[0,0].set_ylabel('Vol-avg DO transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
ax[1,0].set_ylabel('Vol-avg DO transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
ax[2,0].set_ylabel('Vol-avg DO transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
ax[3,0].set_ylabel('Vol-avg DO transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
ax[4,0].set_ylabel('Vol-avg DO transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
# add surface and bottom layer label
ax[0,0].set_title('Surface Layer',fontweight='bold')
ax[0,1].set_title('Bottom Layer',fontweight='bold')
# add subpanel labels
ax[0,0].text(0.02,0.93,'(a) Annual Average'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[0,0].transAxes, fontsize=12)
ax[0,1].text(0.02,0.93,'(b) Annual Average'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[0,1].transAxes, fontsize=12)
ax[1,0].text(0.02,0.93,'(c) Winter Average\n    Jan/Feb/Mar'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[1,0].transAxes, fontsize=12)
ax[1,1].text(0.02,0.93,'(d) Winter Average\n    Jan/Feb/Mar'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[1,1].transAxes, fontsize=12)
ax[2,0].text(0.02,0.93,'(e) Spring Average\n    Apr/May/Jun'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[2,0].transAxes, fontsize=12)
ax[2,1].text(0.02,0.93,'(f) Spring Average\n    Apr/May/Jun'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[2,1].transAxes, fontsize=12)
ax[3,0].text(0.02,0.93,'(g) Summer Average\n    Jul/Aug/Sep'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[3,0].transAxes, fontsize=12)
ax[3,1].text(0.02,0.93,'(h) Summer Average\n    Jul/Aug/Sep'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[3,1].transAxes, fontsize=12)
ax[4,0].text(0.02,0.93,'(i) Fall Average\n    Oct/Nov/Dec'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[4,0].transAxes, fontsize=12)
ax[4,1].text(0.02,0.93,'(j) Fall Average\n    Oct/Nov/Dec'.format(-1*z_interface),zorder=6,
                ha='left', va='top', transform=ax[4,1].transAxes, fontsize=12)

# Annual averages
x = np.arange(len(stations))  # the label locations
print(x)
width = 0.2  # the width of the bars

seasons = ['Annual','Winter','Spring','Summer','Fall']

# loop through different time intervals
for j,season in enumerate(seasons):
    if season == 'Annual':
        minday = 0
        maxday = 363
    elif season == 'Winter':
        minday = 0
        maxday = 91
    elif season == 'Spring':
        minday = 91
        maxday = 182
    elif season == 'Summer':
        minday = 182
        maxday = 274
    elif season == 'Fall':
        minday = 274
        maxday = 363

    multiplier_surf = 0
    multiplier_deep = 0

    for i,station in enumerate(stations):
        
        # get volume of two layers
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
        df_V = pd.read_pickle(fn)
        # Godin filter already applied earlier in workflow
        surf_V = df_V['surface [m3]'].values[1:-1]
        deep_V = df_V['deep [m3]'].values[1:-1]

        print(station)
        for attribute, measurement in surfacelay_dict[station].items():
            # calculate annual average
            time_avg = np.nanmean(measurement[minday:maxday])
            # get volume average
            avg = time_avg/(np.nanmean(surf_V[minday:maxday])) # kmol O2 /s /m3
            # convert to umol
            avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
            # plot
            offset = width * multiplier_surf
            if attribute == 'Residual':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = exchange_color, label=attribute)
            if attribute == 'TEF Exchange Flow':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = exchange_color, label=attribute)
            if attribute == 'EU Exchange Flow':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = exchange_color, label=attribute)
            if attribute == 'TRAPS':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = traps_color, label=attribute)
            if attribute == 'Photosynthesis':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = photo_color, label=attribute)
            if attribute == 'Bio Consumption':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = cons_color, label=attribute)
            if attribute == 'Air-Sea Transfer':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = airsea_color, label=attribute)
            if attribute == 'TEF Vertical':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = vertX_color, label=attribute)
            if attribute == 'EU Vertical':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = vertX_color, label=attribute)
            if attribute == 'Storage':
                rects = ax[j,0].bar(i + offset, avg, width, zorder=5,
                                    color = ddtDOV_color, label=attribute)
            multiplier_surf += 1

        for attribute, measurement in bottomlay_dict[station].items():
            # calculate annual average
            time_avg = np.nanmean(measurement[minday:maxday])
            # get volume average
            avg = time_avg/(np.nanmean(deep_V[minday:maxday]))
            # convert to umol
            avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
            # plot
            offset = width * multiplier_deep
            if attribute == 'Residual':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = exchange_color, label=attribute)
            if attribute == 'TEF Exchange Flow':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = exchange_color, label=attribute)
            if attribute == 'EU Exchange Flow':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = exchange_color, label=attribute)
            if attribute == 'TRAPS':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = traps_color, label=attribute)
            if attribute == 'Photosynthesis':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = photo_color, label=attribute)
            if attribute == 'Bio Consumption':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = cons_color, label=attribute)
            if attribute == 'Air-Sea Transfer':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = airsea_color, label=attribute)
            if attribute == 'TEF Vertical':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = vertX_color, label=attribute)
            if attribute == 'EU Vertical':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = vertX_color, label=attribute)
            if attribute == 'Storage':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = ddtDOV_color, label=attribute)
            multiplier_deep += 1

        if i == 0 and j == 0:
            ax[j,0].legend(loc='lower right', ncols=2, fontsize=9)

# label x axis
if residual == False:
    ax[0,0].set_xticks(x*(1 + width*7) + 0.3, stations_w_DO)
    ax[0,1].set_xticks(x*(1 + width*6) + 0.3, stations_w_DO)
else:
    ax[0,0].set_xticks(x*(1 + width*6) + 0.3, stations_w_DO)
    ax[0,1].set_xticks(x*(1 + width*5) + 0.3, stations_w_DO)

# format figure
plt.subplots_adjust(wspace=0.02, hspace=0.04, top=0.92)
# save figure
out_dir_barcharts = out_dir / 'bar_charts'
Lfun.make_dir(out_dir_barcharts)
if residual == True:
    if show_EU:
        plt.savefig(out_dir_barcharts / 'barchart_residual_EU.png')
    else:
        plt.savefig(out_dir_barcharts / 'barchart_residual_TEF.png')
else:
    if show_EU:
        plt.savefig(out_dir_barcharts / 'barchart_total_EU.png')
    else:
        plt.savefig(out_dir_barcharts / 'barchart_total_TEF.png')


# ##########################################################
# ##             DO usage ratio time series               ##
# ##########################################################

# # initialize figure
# plt.close('all')
# fig, ax = plt.subplots(2,1,figsize = (10,6),sharex=True)
# # format figure
# plt.suptitle('Ratio of Respiration to DO sources',size=14)
# for axis in [ax[0],ax[1]]:
#     axis.set_xlim([dates_local[0],dates_local[-1]])
#     axis.set_ylabel('Ratio [nondimensional]')
#     axis.set_facecolor('#EEEEEE')
#     axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#     axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
#     axis.tick_params(axis='x', labelrotation=30)
#     for border in ['top','right','bottom','left']:
#         axis.spines[border].set_visible(False)
# ax[1].set_xlabel(year)
# ax[0].text(0.02,0.9,'(a) Surface Layer'.format(-1*z_interface),
#                 ha='left', va='top', transform=ax[0].transAxes, fontsize=12)
# ax[1].text(0.02,0.9,'(b) Bottom Layer'.format(-1*z_interface),
#                 ha='left', va='top', transform=ax[1].transAxes, fontsize=12)

# for station in stations:
#     # calculate surface ratio
#     respiration = zfun.lowpass(surfacelay_dict[station]['Bio Consumption'].values,n=60) * -1
#     if residual == True:
#         sources = (zfun.lowpass(surfacelay_dict[station]['TRAPS'].values,n=60) +
#                    zfun.lowpass(surfacelay_dict[station]['Photosynthesis'].values,n=60) +
#                    zfun.lowpass(surfacelay_dict[station]['Air-Sea Transfer'].values,n=60) +
#                    zfun.lowpass(surfacelay_dict[station]['Residual'].values,n=60)*-1 )
#     # else:
#     #     if show_EU == True:
#     #         sources = (surfacelay_dict[station]['TRAPS'] +
#     #                    surfacelay_dict[station]['Photosynthesis'] + 
#     #                    surfacelay_dict[station]['Air-Sea Transfer'] +
#     #                    surfacelay_dict[station]['EU Exchange Flow'] +
#     #                    surfacelay_dict[station]['EU Vertical'])
#     #     else:
#     #         sources = (surfacelay_dict[station]['TRAPS'] +
#     #                    surfacelay_dict[station]['Photosynthesis'] + 
#     #                    surfacelay_dict[station]['Air-Sea Transfer'] +
#     #                    surfacelay_dict[station]['TEF Exchange Flow'] +
#     #                    surfacelay_dict[station]['TEF Vertical'])
#     ratio_surf = respiration / sources

#     # calculate bottom ratio
#     respiration = zfun.lowpass(bottomlay_dict[station]['Bio Consumption'].values,n=90) * -1
#     if residual == True:
#         sources = (zfun.lowpass(bottomlay_dict[station]['TRAPS'].values,n=90) +
#                    zfun.lowpass(bottomlay_dict[station]['Photosynthesis'].values,n=90) +
#                    zfun.lowpass(bottomlay_dict[station]['Residual'].values,n=90))
#     # else:
#     #     if show_EU == True:
#     #         sources = (bottomlay_dict[station]['TRAPS'] +
#     #                    bottomlay_dict[station]['Photosynthesis'] + 
#     #                    bottomlay_dict[station]['EU Exchange Flow'] +
#     #                    bottomlay_dict[station]['EU Vertical'])
#     #     else:
#     #         sources = (bottomlay_dict[station]['TRAPS'] +
#     #                    bottomlay_dict[station]['Photosynthesis'] + 
#     #                    bottomlay_dict[station]['TEF Exchange Flow'] +
#     #                    bottomlay_dict[station]['TEF Vertical'])
            
#     ratio_deep = respiration / sources

#     # plot
#     ax[0].plot(dates_daily[1:-1],ratio_surf,label=station)
#     ax[1].plot(dates_daily[1:-1],ratio_deep,label=station)
#     ax[0].legend(loc='best')
    