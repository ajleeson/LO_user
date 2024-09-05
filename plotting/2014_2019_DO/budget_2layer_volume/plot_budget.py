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
from importlib import reload
reload(get_two_layer)

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
jobname = 'twentyoneinlets'
startdate = '2017.01.01'
enddate = '2017.12.31'
enddate_hrly = '2018.01.01 00:00:00'
year = '2017' # for making a date label

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
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_'+startdate+'_'+enddate) / '2layer_figures'
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

stations = ['lynchcove','penn','budd','case','carr']
# create dictionaries with interface depths
interface_dict = dict()

for i,station in enumerate(sta_dict): # enumerate(stations): 
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
            plt.suptitle(station + ': Volume Budget (Godin filter)',size=14)
            for axis in [ax[0],ax[1],ax[2]]:
                # axis.plot([dates_local[0],dates_local[-1]],[0,0],linewidth=0.5,color='k')
                axis.set_xlim([dates_local[0],dates_local[-1]])
                # axis.set_ylim([-1200,1200])
                axis.set_ylabel(r'Q [m$^3$ s$^{-1}$]')
                axis.set_facecolor('#EEEEEE')
                axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
                axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
                axis.tick_params(axis='x', labelrotation=30)
                for border in ['top','right','bottom','left']:
                    axis.spines[border].set_visible(False)
            ax[2].set_xlabel(year)
            ax[0].set_title('(a) Surface [shallower than {} m]'.format(-1*z_interface),
                            loc='left')
            ax[1].set_title('(b) Bottom [deeper than {} m]'.format(-1*z_interface),
                            loc='left')
            ax[2].set_title(r'(c) Error [sum of surface and bottom vertical exchange]',
                            loc='left')

# # --------------------------- get exchange flow terms ----------------------------------------
            fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layerEU_exchange_flow' / (station + '.p')
            df_exchange = pd.read_pickle(fn)
            exchange_surf_unfiltered = df_exchange['surface [m3/s]'].values 
            exchange_deep_unfiltered = df_exchange['deep [m3/s]'].values
            # 10-day hanning window filter (10 days = 240 hours)
            # exchange_surf = zfun.lowpass(exchange_surf_unfiltered, f='hanning', n=240) 
            # exchange_deep = zfun.lowpass(exchange_deep_unfiltered, f='hanning', n=240)
            # Godin filter
            EU_surf = zfun.lowpass(exchange_surf_unfiltered, f='godin')[36:-34:24]
            EU_deep = zfun.lowpass(exchange_deep_unfiltered, f='godin')[36:-34:24]
            EU_color = 'royalblue'

# --------------------------- get TEF exchange flow terms ----------------------------------------
            in_dir = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station+ '.nc')
            bulk = xr.open_dataset(in_dir)
            tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
            Q_p = tef_df['q_p'] # Qin
            Q_m = tef_df['q_m'] # Qout
            TEF_surf = Q_m.values
            TEF_deep = Q_p.values
            TEF_color = 'blue'

# ------------------------------- get rivers and WWTPs ----------------------------------------
            fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_traps' / (station + '.p')
            df_traps = pd.read_pickle(fn)
            rivers_surf_unfiltered = df_traps['surface [m3/s]'].values
            wwtps_deep_unfiltered = df_traps['deep [m3/s]'].values
            # 10-day hanning window filter (10 days = 240 hours)
            # traps_surf = zfun.lowpass(rivers_surf_unfiltered, f='hanning', n=240)
            # traps_deep = zfun.lowpass(wwtps_deep_unfiltered, f='hanning', n=240)
            # Pad with zeros if no rivers or WWTPs
            if rivers_surf_unfiltered.size == 0:
                 rivers_surf_unfiltered = np.zeros(8761)
            if wwtps_deep_unfiltered.size == 0:
                 wwtps_deep_unfiltered = np.zeros(8761)
            # Godin filter
            traps_surf = zfun.lowpass(rivers_surf_unfiltered, f='godin')[36:-34:24]
            traps_deep = zfun.lowpass(wwtps_deep_unfiltered, f='godin')[36:-34:24]
            # traps_color = 'turquoise'
            traps_color = 'black'

# ------------------------------- get volume storage ----------------------------------------
            fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
            df_V = pd.read_pickle(fn)
            # Godin filter already applied earlier in workflow
            surf_V = df_V['surface [m3]'].values
            deep_V = df_V['deep [m3]'].values
            print(surf_V[0]+deep_V[0])
            # get dVdt
            dVdt_surf = np.diff(surf_V)* (1/24) * (1/60) * (1/60) # m3/day to m3/s
            dVdt_deep = np.diff(deep_V)* (1/24) * (1/60) * (1/60) # m3/day to m3/s
            dVdt_color = 'deeppink'

# ------------------------------- get vertical exchange ----------------------------------------

            # pick color
            vertXEU_color = 'mediumorchid'
            vertX_color = 'rebeccapurple'

            # calculate error term, and acribe that the the vertical exchange
            vertXTEF_surf = dVdt_surf[1::] - (TEF_surf + traps_surf)
            vertXTEF_deep = dVdt_deep[1::] - (TEF_deep + traps_deep)
            vertXEU_surf = dVdt_surf[1::] - (EU_surf + traps_surf)
            vertXEU_deep = dVdt_deep[1::] - (EU_deep + traps_deep)
            

# ---------------------------------- plot and save --------------------------------------------
            # plot surface
            ax[0].plot(dates_local_daily[1:-1],EU_surf,color=EU_color,
                       linewidth=2,alpha=0.5,label='EU Exchange Flow')
            ax[0].plot(dates_local_daily[1:-1],TEF_surf,color=TEF_color,
                       linewidth=1,linestyle='--',label='TEF Exchange Flow')
            ax[0].plot(dates_local_daily[0:-1],dVdt_surf,color=dVdt_color,
                       linewidth=1.5,zorder=5,
                       label=r'$\frac{\mathrm{dV}}{\mathrm{dt}}$')
            ax[0].plot(dates_local_daily[1:-1],traps_surf,color=traps_color,
                       linewidth=3,zorder=4,label='TRAPS')
            ax[0].plot(dates_local_daily[1:-1],vertXEU_surf,color=vertXEU_color,
                       linewidth=2,alpha=0.5,label='EU Vertical')
            ax[0].plot(dates_local_daily[1:-1],vertXTEF_surf,color=vertX_color,
                       linewidth=1,linestyle='--',label='TEF Vertical')
            ax[0].legend(loc='best',ncol=3)
            
            # plot deep
            ax[1].plot(dates_local_daily[1:-1],EU_deep,color=EU_color,
                       linewidth=2,alpha=0.5)
            ax[1].plot(dates_local_daily[1:-1],TEF_deep,color=TEF_color,
                       linewidth=1,linestyle='--')
            ax[1].plot(dates_local_daily[0:-1],dVdt_deep,color=dVdt_color,
                       linewidth=1.5,zorder=5)
            ax[1].plot(dates_local_daily[1:-1],traps_deep,color=traps_color,
                       linewidth=3,zorder=4)
            ax[1].plot(dates_local_daily[1:-1],vertXEU_deep,color=vertXEU_color,
                       linewidth=2,alpha=0.5)
            ax[1].plot(dates_local_daily[1:-1],vertXTEF_deep,color=vertX_color,
                       linewidth=1,linestyle='--')

            # plot error
            error_TEF = vertXTEF_surf + vertXTEF_deep
            error_EU = vertXEU_surf + vertXEU_deep
            ax[2].plot(dates_local_daily[1:-1],error_EU,color='olivedrab',
                       linewidth=3,alpha=0.5,label='EU Error')
            ax[2].plot(dates_local_daily[1:-1],error_TEF,color='darkgreen',
                       linewidth=1,linestyle='--',label='TEF Error')
            ax[2].legend(loc='best')

            for axis in [ax[0],ax[1],ax[2]]:
                ylimval_TEF = np.nanmax(np.abs(TEF_surf))*1.1
                ylimval_EU = np.nanmax(np.abs(EU_surf))*1.1
                ylimval = np.nanmax([ylimval_TEF,ylimval_EU])
                axis.set_ylim([-1*ylimval,ylimval])


            plt.savefig(out_dir / (station+'.png'))