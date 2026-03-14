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
import matplotlib.patches as patches
import csv
import cmocean
import matplotlib.pylab as plt
import gsw
import pickle
import get_two_layer

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

residual = True # recirculation
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
storage_deep_dict = {'lynchcove': pd.DataFrame(),
                'penn': pd.DataFrame(),
                'budd': pd.DataFrame(),
                'carr': pd.DataFrame(),
                'case': pd.DataFrame()}
DOconcen_dict = {'lynchcove': pd.DataFrame(),
                'penn': pd.DataFrame(),
                'budd': pd.DataFrame(),
                'carr': pd.DataFrame(),
                'case': pd.DataFrame()}

for i,station in enumerate(stations): # enumerate(sta_dict):

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
            fig, ax = plt.subplots(4,1,figsize = (12,9),sharex=True)
            # format figure
            plt.suptitle(station + ': DO Budget (Godin Filter)',size=14)
            for axis in [ax[0],ax[1],ax[2],ax[3]]:
                # axis.plot([dates_local[0],dates_local[-1]],[0,0],color='k')
                axis.set_xlim([dates_local[0],dates_local[-1]])
                if axis == ax[3]:
                    axis.set_ylabel(r'DO [mg L$^{-1}$]')
                else:
                    axis.set_ylabel(r'DO transport [kmol O$_2$ s$^{-1}$]')
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
        

# ---------------------------------- plot budgets --------------------------------------------
            # plot surface
            if residual == False:
                ax[0].plot(dates_local_daily[1:-1],TEF_surf,color=exchange_color,
                        linewidth=1,label='TEF Exchange Flow')
                if show_EU:
                    ax[0].plot(dates_local_daily[1:-1],EU_surf,color=exchange_color,
                            linewidth=3,alpha=0.2,label='EU Exchange Flow')
            else:
                ax[0].plot(dates_local_daily[1:-1],TEF_surf + vertX_surf_TEF,color=exchange_color,
                        linewidth=2,label='Recirculation TEF')
                if show_EU:
                    ax[0].plot(dates_local_daily[1:-1],EU_surf + vertX_surf_EU,color=exchange_color,
                            linewidth=3,alpha=0.2,label='Recirculation EU')
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
            if residual == False:
                ax[0].plot(dates_local_daily[1:-1],vertX_surf_TEF,color=vertX_color,
                       linewidth=1,label='TEF Vertical')
                if show_EU:
                    ax[0].plot(dates_local_daily[1:-1],vertX_surf_EU,color=vertX_color,
                            linewidth=3,alpha=0.2,label='EU Vertical')
            ax[0].legend(loc='lower center',ncol=6)
            
            # plot deep
            if residual == False:
                ax[1].plot(dates_local_daily[1:-1],TEF_deep,color=exchange_color,
                        linewidth=1)
                if show_EU:
                    ax[1].plot(dates_local_daily[1:-1],EU_deep,color=exchange_color,
                            linewidth=3,alpha=0.5)
            else:
                ax[1].plot(dates_local_daily[1:-1],TEF_deep + vertX_deep_TEF,color=exchange_color,
                        linewidth=2)
                if show_EU:
                    ax[1].plot(dates_local_daily[1:-1],EU_deep + vertX_deep_EU,color=exchange_color,
                            linewidth=3,alpha=0.2)
            ax[1].plot(dates_local_daily[1:-1],traps_deep,color=traps_color,
                       linewidth=3,zorder=2)
            ax[1].plot(dates_local_daily[1:-1],photo_deep,color=photo_color,
                       linewidth=2)
            ax[1].plot(dates_local_daily[1:-1],cons_deep,color=cons_color,
                       linewidth=2,linestyle=':',zorder=9)
            ax[1].plot(dates_local_daily[1:-1],ddtDOV_deep,color=ddtDOV_color,
                       linewidth=1,zorder=7)#,alpha=0.6)
            if residual == False:
                ax[1].plot(dates_local_daily[1:-1],vertX_deep_TEF,color=vertX_color,
                        linewidth=1)
                if show_EU:
                    ax[1].plot(dates_local_daily[1:-1],vertX_deep_EU,color=vertX_color,
                            linewidth=3,alpha=0.5)

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
                if residual == False:
                    ylimval = np.nanmax(np.abs(TEF_surf))*1.1
                else:
                    ylimval = np.nanmax(np.abs(TEF_surf + vertX_surf_TEF))*1.2
                axis.set_ylim([-1*ylimval,ylimval])

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

            # plot DO
            ax[3].plot(dates_local_daily[1:-1],o2_surf,color='royalblue',
                       linewidth=2,label='Avg. Surface')
            ax[3].plot(dates_local_daily[1:-1],o2_deep,color='mediumorchid',
                       linewidth=2,label='Avg. Deep')
            ax[3].plot(dates_local_daily[1:-1],DO_in,color='black',
                       linewidth=1,linestyle=':',label='Avg. Qin')
            ax[3].plot(dates_local_daily[1:-1],bottom_oxygen,color='crimson',
                       linewidth=1,linestyle='--',label='Avg. bottom DO')
            ax[3].plot(dates_local_daily[1:-1],oxygen_min,color='crimson',
                       linewidth=1,linestyle='-',label='Minimum DO')
            ax[3].legend(loc='lower left',ncol=2)


            ax[3].set_ylim([0,12])


# ---------------------------------- save figure --------------------------------------------

            plt.subplots_adjust(hspace=0.06, bottom=0.06, top=0.94)
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
                    surfacelay_dict[station]['Recirculation'] = EU_surf + vertX_surf_EU
                else:
                    surfacelay_dict[station]['Recirculation'] = TEF_surf + vertX_surf_TEF
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
                    bottomlay_dict[station]['Recirculation'] = EU_deep + vertX_deep_EU
                else:
                    bottomlay_dict[station]['Recirculation'] = TEF_deep + vertX_deep_TEF
            bottomlay_dict[station]['TRAPS'] = traps_deep
            bottomlay_dict[station]['Photosynthesis'] = photo_deep
            bottomlay_dict[station]['Bio Consumption'] = cons_deep
            if residual == False:
                if show_EU:
                    bottomlay_dict[station]['EU Vertical'] = vertX_deep_EU
                else:
                    bottomlay_dict[station]['TEF Vertical'] = vertX_deep_TEF
            bottomlay_dict[station]['Storage'] = ddtDOV_deep
            storage_deep_dict[station]['Storage'] = ddtDOV_deep

# ------------------------- save DO concentrations in dataframe dict -----------------------------------

            # mg/L units

            # DO concentrations
            DOconcen_dict[station]['Surface Layer'] = o2_surf
            DOconcen_dict[station]['Deep Layer'] = o2_deep
            DOconcen_dict[station]['Bottom DO'] = bottom_oxygen
            DOconcen_dict[station]['Minimum DO'] = oxygen_min
            DOconcen_dict[station]['Qin DO'] = DO_in

##########################################################
##           Create seasonal bar charts                 ##
##########################################################

# annual average summer bottom DO
stations_looped = ['Lynch\nCove','Penn\nCove',
                 'Case\nInlet','Carr\nInlet',
                 'Budd\nInlet']

# initialize figure
plt.close('all')
fig, ax = plt.subplots(5,2,figsize = (11,8), sharey=True,sharex='col')
# format figure
if show_EU:
    exchange = 'Eulerian'
else:
    exchange = 'TEF'
plt.suptitle('Volume-averaged DO transport rates ({}) and DO concentrations'.format(exchange),size=14)
for axis in [ax[0,0],ax[0,1],ax[1,0],ax[1,1],ax[2,0],ax[2,1],ax[3,0],ax[3,1],ax[4,0],ax[4,1]]:
    axis.set_facecolor('#EEEEEE')
    axis.grid(True,color='w',linewidth=1,linestyle='-',axis='y')
    for border in ['top','right','bottom','left']:
        axis.spines[border].set_visible(False)
# label y-axis
ax[0,0].set_ylabel('Transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
ax[1,0].set_ylabel('Transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
ax[2,0].set_ylabel('Transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
ax[3,0].set_ylabel('Transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
ax[4,0].set_ylabel('Transport\n'+ r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]')
# ax[0,2].set_ylabel(r'DO [mg L$^{-1}$]')
# ax[1,2].set_ylabel(r'DO [mg L$^{-1}$]')
# ax[2,2].set_ylabel(r'DO [mg L$^{-1}$]')
# ax[3,2].set_ylabel(r'DO [mg L$^{-1}$]')
# ax[4,2].set_ylabel(r'DO [mg L$^{-1}$]')
# add surface and bottom layer label
ax[0,0].set_title('Surface Layer',fontweight='bold')
ax[0,1].set_title('Bottom Layer',fontweight='bold')
# ax[0,2].set_title('DO Concentrations',fontweight='bold')
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
width = 0.2  # the width of the bars

seasons = ['Annual','Winter','Spring','Summer','Fall']

# loop through different time intervals
for j,season in enumerate(seasons):

    # print('========================='+season+'==========================')

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
    multiplier_DO = 0

    for i,station in enumerate(stations):

        # print('\n' + station)
        
        # get volume of two layers
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
        df_V = pd.read_pickle(fn)
        # Godin filter already applied earlier in workflow
        surf_V = df_V['surface [m3]'].values[1:-1]
        deep_V = df_V['deep [m3]'].values[1:-1]

        # Surface layer -----------------------------------------------
        for attribute, measurement in surfacelay_dict[station].items():
            # calculate annual average
            time_avg = np.nanmean(measurement[minday:maxday])
            # get volume average
            avg = time_avg/(np.nanmean(surf_V[minday:maxday])) # kmol O2 /s /m3
            # avg = time_avg # skip volume averaging
            # convert to umol
            avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
            # plot
            offset = width * multiplier_surf
            if attribute == 'Recirculation':
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

        # Bottom layer ----------------------------------------------
        for attribute, measurement in bottomlay_dict[station].items():
            # calculate annual average
            time_avg = np.nanmean(measurement[minday:maxday])
            # get volume average
            avg = time_avg/(np.nanmean(deep_V[minday:maxday]))
            # avg = time_avg # skip volume averaging
            # convert to umol
            avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3
            # plot
            offset = width * multiplier_deep
            if attribute == 'Recirculation':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = exchange_color, label=attribute)
                # print('Recirculation: {}'.format(round(avg,3)))
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
                # print('Photosynthesis: {}'.format(round(avg,3)))
            if attribute == 'Bio Consumption':
                rects = ax[j,1].bar(i + offset, avg, width, zorder=5,
                                    color = cons_color, label=attribute)
                # print('Bio Consumption: {}'.format(round(avg,3)))
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

        # # DO concentrations ----------------------------------------
        # for attribute, measurement in DOconcen_dict[station].items():
        #     # calculate annual average
        #     avg = np.nanmean(measurement[minday:maxday])
        #     # plot
        #     offset = width * multiplier_DO
        #     if attribute == 'Surface Layer':
        #         rects = ax[j,2].bar(i + offset, avg, width, zorder=5,
        #                             color='deepskyblue', alpha=0.7, label=attribute)
        #     if attribute == 'Deep Layer':
        #         rects = ax[j,2].bar(i + offset, avg, width, zorder=5,
        #                             color = 'blueviolet', alpha=0.6, label=attribute)
        #     if attribute == 'Bottom DO':
        #         rects = ax[j,2].bar(i + offset, avg, width, zorder=5,
        #                             color = 'deeppink', alpha=0.7, label=attribute)
        #     if attribute == 'Qin DO':
        #         rects = ax[j,2].bar(i + offset, avg, width, zorder=5,
        #                             color = 'black', alpha=0.5, label=attribute)
        #     multiplier_DO += 1

        if i == 0 and j == 0:
            ax[j,0].legend(loc='lower right', ncols=2, fontsize=9)
            # ax[j,2].legend(loc='upper right', fontsize=9)

# label x axis
if residual == False:
    ax[0,0].set_xticks(x*(1 + width*7) + 0.3, stations_looped)
    ax[0,1].set_xticks(x*(1 + width*6) + 0.3, stations_looped)
else:
    ax[0,0].set_xticks(x*(1 + width*6) + 0.3, stations_looped)
    ax[0,1].set_xticks(x*(1 + width*5) + 0.3, stations_looped)
# ax[0,2].set_xticks(x*(1 + width*4) + 0.3, stations_looped)

# format figure
plt.subplots_adjust(wspace=0.02, hspace=0.1, top=0.92, right=0.98)
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

##########################################################
##             Ratio of sinks to sources                ## ######################################################################
##########################################################

stations = ['lynchcove','penn','case','carr','budd']

# create list of station colors
# station_color = ['hotpink','orange','gold','limegreen','deepskyblue'] # rainbow basic
# station_color = ['#81667A','#8C8A93','#92B4A7','#B6CB9E','#D1F0B1']   # lilypad theme
# station_color = ['#db5375','#E99F7D','#F4D35E','#B3C169','#6EB0C1']   # rainbow muted
# station_color = ['#EF476F','darkorange','#FFD166','#B3C169','#6EB0C1']    # rainbow muted 2
station_color = ['#F9627D','#62B6CB','#A8C256','#96031A','#957FEF']    # pretty

# create list of hatching
hatch_dict = {'Recirculation': 'OO',
              'TRAPS': None,
              'Photosynthesis': '**',
              'Bio Consumption': 'XX'}

# initialize figure
plt.close('all')
fig, ax = plt.subplots(4,1,figsize = (12,8.6),sharex=True)

# format figure
plt.suptitle('Bottom Layer DO Processes (2014)',size=16,fontweight='bold')
# subtitles
ax[0].set_title('(a) DO concentration',loc='left', fontsize=12, fontweight='bold')
ax[1].set_title('(b) Average DO transport rate (volume-integrated)',loc='left', fontsize=12, fontweight='bold')
ax[2].set_title('(c) Average DO transport rate per unit volume',loc='left', fontsize=12, fontweight='bold')
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
ax[1].set_ylabel(r'[kmol O$_2$ s$^{-1}$]',fontsize=14)
ax[2].set_ylabel(r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]',fontsize=14)
ax[3].set_ylabel(r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]',fontsize=14)

width = 0.15
multiplier_conc = 0
multiplier_rate_int = 0
multiplier_rate_avg = 0
multiplier_ratio = 0

# ratio of sources to sinks --------------------------------------------------------
# loop through stations
for i,station in enumerate(stations):
        
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

        # calculate average DO over time interval, and calculate standard deviation
        deep_lay_DO = np.nanmean(DOconcen_dict[station]['Deep Layer'][minday:maxday])
        deep_lay_std_dev = np.std(DOconcen_dict[station]['Deep Layer'][minday:maxday])
        bott_DO = np.nanmean(DOconcen_dict[station]['Bottom DO'][minday:maxday])
        bott_std_dev = np.std(DOconcen_dict[station]['Bottom DO'][minday:maxday])
        # calculate minimum DO in each season
        min_DO = np.nanmin(DOconcen_dict[station]['Minimum DO'][minday:maxday])


        # plot
        offset = width * multiplier_conc
        # if j == 0:
        #     rects = ax[0].bar(j + offset, bott_DO, width, zorder=5, color = station_color[i], label=station)
        # else:
        #     rects = ax[0].bar(j + offset, bott_DO, width, zorder=5, color = station_color[i])

        if j == 0 and i==0:
            rects = ax[0].scatter(j + offset, deep_lay_DO, zorder=5, s=120, color=station_color[i], label='Layer average')
            # rects = ax[0].scatter(j + offset, bott_DO, zorder=5, s=30,marker='s',color='white',
            #                        edgecolor='k', label='Bottom s-level')
            rects = ax[0].scatter(j + offset, min_DO, zorder=5, s=30,marker='s',color=station_color[i],
                                   edgecolor='white', label='Absolute minimum')
        else:
            rects = ax[0].scatter(j + offset, deep_lay_DO, zorder=5, s=120, color=station_color[i])
            # rects = ax[0].scatter(j + offset, bott_DO, zorder=5, s=30,marker='s',color='white', edgecolor='k')
            rects = ax[0].scatter(j + offset, min_DO, zorder=5, s=30,marker='s',color=station_color[i], edgecolor='white')
        # # add standard deviation
        # ax[0].errorbar(j + offset, deep_lay_DO, yerr=deep_lay_std_dev, color=station_color[i])
        # ax[0].errorbar(j + offset, bott_DO, yerr=bott_std_dev, color='black')

    multiplier_conc += 1

# transport terms (volume integrated)--------------------------------------------------------

loop = 0
# loop through stations
for i,station in enumerate(stations):

    # remove storage term from dict
    del bottomlay_dict[station]['Storage']
        
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
        
        bottom = np.zeros(5)
        for attribute, measurement in bottomlay_dict[station].items():

            # calculate time average
            time_avg = np.nanmean(measurement[minday:maxday])
            avg = time_avg # skip volume avg (kmol O2 /s)

            if attribute == 'Photosynthesis':
                custom_color = 'black'
            elif attribute == 'TRAPS':
                custom_color = 'white'
            else:
                custom_color = station_color[i]
            if attribute == 'Recirculation':
                hatch='//'
                edgecolor = None
                alpha=0.05
            else:
                hatch=None
                edgecolor = None
                alpha=1

            # plot
            offset = width * multiplier_rate_int
            if loop == 0 and attribute != 'TRAPS':
                ax[1].bar(j + offset, avg, width, label=attribute, zorder=5, color=custom_color,
                      bottom=bottom, hatch=hatch, edgecolor=edgecolor, alpha=alpha)
            else:
                ax[1].bar(j + offset, avg, width, zorder=5, color=custom_color,
                      bottom=bottom, hatch=hatch, edgecolor=edgecolor, alpha=alpha)
            
            if attribute == 'Photosynthesis':
                bottom = 0
            else:
                bottom += avg

        loop = 1

    multiplier_rate_int += 1

# transport terms (volume averaged)--------------------------------------------------------

loop = 0
# loop through stations
for i,station in enumerate(stations):

    # get volume of bottom layer
    fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
    df_V = pd.read_pickle(fn)
    # Godin filter already applied earlier in workflow
    deep_V = df_V['deep [m3]'].values[1:-1]
        
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
        
        bottom = np.zeros(5)
        for attribute, measurement in bottomlay_dict[station].items():

            # calculate time average
            time_avg = np.nanmean(measurement[minday:maxday])
            # get volume average
            avg = time_avg/(np.nanmean(deep_V[minday:maxday])) # kmol O2 /s /m3
            # convert to umol
            avg = avg * 1000 * 1000 * 1000 # umol O2 /s /m3

            if attribute == 'Photosynthesis':
                custom_color = 'black'
            elif attribute == 'TRAPS':
                custom_color = 'white'
            else:
                custom_color = station_color[i]
            if attribute == 'Recirculation':
                hatch='//'
                edgecolor = None
                alpha=0.05
            else:
                hatch=None
                edgecolor = None
                alpha=1

            # plot
            offset = width * multiplier_rate_avg
            if loop == 0 and attribute != 'TRAPS':
                ax[2].bar(j + offset, avg, width, label=attribute, zorder=5, color=custom_color,
                      bottom=bottom, hatch=hatch, edgecolor=edgecolor, alpha=alpha)
            else:
                ax[2].bar(j + offset, avg, width, zorder=5, color=custom_color,
                      bottom=bottom, hatch=hatch, edgecolor=edgecolor, alpha=alpha)
            
            if attribute == 'Photosynthesis':
                bottom = 0
            else:
                bottom += avg

        loop = 1

    multiplier_rate_avg += 1

# ratio of sources to sinks --------------------------------------------------------
# loop through stations
for i,station in enumerate(stations):

    # get volume of bottom layer
    fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
    df_V = pd.read_pickle(fn)
    # Godin filter already applied earlier in workflow
    deep_V = df_V['deep [m3]'].values[1:-1]
        
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

        # # calculate average ratio of sources to sinks over time interval
        # photo = np.nanmean(bottomlay_dict[station]['Photosynthesis'][minday:maxday])
        # recirc = np.nanmean(bottomlay_dict[station]['Recirculation'][minday:maxday])
        # wwtps = np.nanmean(bottomlay_dict[station]['TRAPS'][minday:maxday])
        # sources = photo + recirc + wwtps
        # sinks = np.nanmean(bottomlay_dict[station]['Bio Consumption'][minday:maxday] * -1)
        # ratio = sources/sinks
        # # ratio = sinks/sources

        # calculate average ratio of sources to sinks over time interval per unit volume
        storage_all_time = storage_deep_dict[station]['Storage'][minday:maxday]/deep_V[minday:maxday]
        storage = np.nanmean(storage_all_time)
        # convert from kmol to umol
        storage = storage * 1000 * 1000 * 1000 # umol O2 /s /m3

        # get min and max 
        upper_err = (np.nanmax(storage_all_time) * 1000 * 1000 * 1000 )- storage
        lower_err = storage - (np.nanmin(storage_all_time) * 1000 * 1000 * 1000 )
        y_error = np.reshape(np.array([lower_err,upper_err]),(2,1))

        # plot
        offset = width * multiplier_ratio
        # if j == 0:
        #     rects = ax[3].scatter(j + offset, storage, s=120, marker='o', zorder=5, color = station_color[i], label=station)
        # else:
        #     rects = ax[3].scatter(j + offset, storage, s=120, marker='o', zorder=5, color = station_color[i])
        if j == 0:
            rects = ax[3].bar(j + offset, storage, width, zorder=5, color = station_color[i], label=station)
        else:
            rects = ax[3].bar(j + offset, storage, width, zorder=5, color = station_color[i])
        # add error
        # ax[3].errorbar(j + offset, storage, yerr=y_error, color='black')

        # if storage > 0:
        #     ax[3].text(j + offset,storage+upper_err,str(round(storage,3)),zorder=6,
        #         ha='left', va='bottom', fontsize=10, rotation=75, color='blue')
        # else:
        #     ax[3].text(j + offset,storage+upper_err,str(round(storage,3)),zorder=6,
        #         ha='left', va='bottom', fontsize=10, rotation=75, color='red')

    multiplier_ratio += 1

# add legend
ax[0].legend(loc='upper right', fontsize=10)
ax[1].legend(loc='lower right', fontsize=10)
ax[3].legend(loc='upper center', ncol=5, fontsize=12,handletextpad=0.1,
             bbox_to_anchor=(0.64, 4.86), frameon=False)

# restrict y axis
ax[0].set_ylim([0,12])
# ax[3].set_ylim([-0.5,0.5])

# add description
ax[1].text(0.8,0.94,'*WWTP input is negligible',zorder=6,
                ha='left', va='top', transform=ax[1].transAxes, fontsize=10)
# ax[3].text(0.01,0.94,
#            r'ratio = $\frac{\mathrm{Sources}}{\mathrm{Sinks}}$ = $\frac{\mathrm{Photo + Recirc + WWTPs}}{\mathrm{Cons.}}$',
#            zorder=6, ha='left', va='top', transform=ax[3].transAxes, fontsize=14)
# ax[3].text(0.3,0.69,'ratio > 1 = gaining oxygen\nratio < 1 = losing oxygen',zorder=6,
#                 ha='left', va='top', transform=ax[3].transAxes, fontsize=12)

# lines
ax[0].axhline(2,0,5, linewidth=2, linestyle=':', color='k', zorder=4) # hypoxia
ax[1].axhline(0,0,5, linewidth=2, linestyle=':', color='k', zorder=6) # zero
ax[2].axhline(0,0,5, linewidth=2, linestyle=':', color='k', zorder=6) # zero
ax[3].axhline(0,0,5, linewidth=2, linestyle=':', color='k', zorder=6) # one
# vertical lines
for addit in [0,1,2,3]:
    ax[0].axvline(0.8+addit,0,12, linewidth=1, color='white')
    ax[1].axvline(0.8+addit,-0.3,1, linewidth=1, color='white')
    ax[2].axvline(0.8+addit,-0.3,1, linewidth=1, color='white')
    ax[3].axvline(0.8+addit,0,5, linewidth=1, color='white')

# add category labels on x-axis
ax[3].set_xticks(x*(0.25 + width*5) + 0.3, seasons, fontsize=14)

plt.subplots_adjust(hspace=0.2, top=0.9, bottom=0.05, right=0.95)
plt.show()

# # ##########################################################
# # ##             Time series of storage term              ## ######################################################################
# # ##########################################################

# # initialize figure
# fig, ax = plt.subplots(5,1,figsize = (12,8),sharey=True,sharex=True)
# axis = ax.ravel()

# # loop through stations
# for i,station in enumerate(stations):

#     # format figure
#     plt.suptitle('Major Bottom Layer Rates\n' + r'[$\mu$mol O$_2$ s$^{-1}$ m$^{-3}$]',
#                  size=16)
#     # format grid
#     axis[i].set_facecolor('#EEEEEE')
#     axis[i].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
#     axis[i].tick_params(axis='x', labelrotation=30)
#     axis[i].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#     for border in ['top','right','bottom','left']:
#         axis[i].spines[border].set_visible(False)
#     axis[i].tick_params(axis='y', labelsize=12)
#     axis[i].text(0.02,0.05,station,fontsize=12,fontweight='bold',ha='left', va='bottom', transform=axis[i].transAxes)

#     # get volume of bottom layer
#     fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
#     df_V = pd.read_pickle(fn)
#     # Godin filter already applied earlier in workflow
#     deep_V = df_V['deep [m3]'].values[1:-1]

#     # calculate d/dt(DO)
#     storage = storage_deep_dict[station]['Storage']/deep_V * 1000 * 1000 * 1000
#     # calculate major source terms
#     photo = bottomlay_dict[station]['Photosynthesis']/deep_V * 1000 * 1000 * 1000
#     recirc = bottomlay_dict[station]['Recirculation']/deep_V * 1000 * 1000 * 1000
#     # calculate major sink terms
#     cons = bottomlay_dict[station]['Bio Consumption']/deep_V * 1000 * 1000 * 1000

#     # plot
#     axis[i].plot(dates_local_daily[1:-1],recirc,color='darkorange',linewidth=2,label='Recirculation')
#     axis[i].plot(dates_local_daily[1:-1],cons,color='deeppink',linewidth=2,label='Bio Consumption')
#     axis[i].plot(dates_local_daily[1:-1],photo,color='lightseagreen',linewidth=2,label='Photosynthesis')
#     axis[i].plot(dates_local_daily[1:-1],storage,color='navy',linewidth=0.6,label='Storage')

#     # get average bottom layer DO
#     deep_lay_DO =DOconcen_dict[station]['Deep Layer']
#     # bott_DO = DOconcen_dict[station]['Bottom DO']
#     ax2 = axis[i].twinx()
#     ax2.plot(dates_local_daily[1:-1],deep_lay_DO,color='limegreen',linewidth=2)
#     # format grid
#     ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
#     ax2.tick_params(axis='x', labelrotation=30)
#     for border in ['top','right','bottom','left']:
#         ax2.spines[border].set_visible(False)
#     ax2.tick_params(axis='y', labelsize=12, labelcolor='limegreen')
#     ax2.set_ylim([0,12])
#     ax2.set_yticks(np.arange(0, 12, 3))

#     # format labels
#     axis[i].set_xlim([dates_local[0],dates_local[-1]])
#     axis[i].set_ylim([-0.5,0.5])

#     # legend
#     if i == 0:
#         axis[i].legend(loc='upper left',ncol=4)

# plt.subplots_adjust(hspace=0.1, top=0.9, bottom=0.05, right=0.95)
# plt.show()

# ##########################################################
# ##                  Bottom DO time series               ## ######################################################################
# ##########################################################

# initialize figure
fig, ax = plt.subplots(1,1,figsize = (10,6))
# format figure
plt.suptitle('Bottom layer DO time series [mg/L]\nwith annual mean subtracted and 30-day Hanning filter applied',
                size=16)
# format grid
ax.set_facecolor('#EEEEEE')
ax.tick_params(axis='x', labelrotation=30)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)
ax.tick_params(axis='y', labelsize=12)
        
# loop through stations
for i,station in enumerate(stations):

    # get average bottom layer DO
    deep_lay_DO = DOconcen_dict[station]['Deep Layer']

    # subtract annual mean
    deep_lay_DO_noDCoffset = deep_lay_DO - np.nanmean(DOconcen_dict[station]['Deep Layer'])
    # bott_DO = DOconcen_dict[station]['Bottom DO']
    deep_lay_DO_noDCoffset_filtered = zfun.lowpass(deep_lay_DO_noDCoffset.values,f='hanning',n=30)

    # ax.plot(dates_local_daily[1:-1],deep_lay_DO,color=station_color[i],linewidth=2, label=station)
    # ax.plot(dates_local_daily[1:-1],deep_lay_DO_noDCoffset,color=station_color[i],linewidth=3,alpha=0.4)
    ax.plot(dates_local_daily[1:-1],deep_lay_DO_noDCoffset_filtered,color=station_color[i],linewidth=2,label=station)

    # # boxplot of spread of annual DO with annual mean subtracted
    # dot = dict(markerfacecolor=station_color[i], marker='.',alpha=0.2,markeredgecolor='none')
    # bp = ax.boxplot(deep_lay_DO_noDCoffset,positions=[i],patch_artist = True, flierprops=dot)
    # plt.setp(bp['medians'], color='k')
    # for patch in bp['boxes']:
    #     patch.set(facecolor=station_color[i], alpha=0.7)

    # format labels
    ax.set_xlim([dates_local[0],dates_local[-1]])
    ax.set_ylim([-5,5])

ax.legend(loc='upper right')

plt.show()