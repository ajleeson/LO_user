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
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('budget_'+startdate+'_'+enddate) / 'figures'
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

for i,station in enumerate(['lynchcove']): # enumerate(sta_dict):
        # print status
        print('({}/{}) Working on {}...'.format(i+1,len(sta_dict),station))

        # get interface depth from csv file
        with open('interface_depths.csv', 'r') as f:
            for line in f:
                inlet, z_interface = line.strip().split(',')
        z_interface = float(z_interface)

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
            fig, ax = plt.subplots(3,1,figsize = (12,7),sharex=True)
            # format figure
            plt.suptitle(station + ': DO Budget (Godin filter)',size=14)
            for axis in [ax[0],ax[1],ax[2]]:
                axis.plot([dates_local[0],dates_local[-1]],[0,0],linewidth=0.5,color='k')
                axis.set_xlim([dates_local[0],dates_local[-1]])
                axis.set_ylim([-1200,1200])
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
            ax[2].set_title(r'(c) Error = $\frac{\mathrm{dV}}{\mathrm{dt}}$ - (Exchange + TRAPS)',
                            loc='left')

# # --------------------------- get exchange flow terms ----------------------------------------
#             fn = Ldir['LOo'] / 'pugetsound_DO' / ('budget_' + startdate + '_' + enddate) / 'DO_exchange_flow' / (station + '.p')
#             df_exchange = pd.read_pickle(fn)
#             exchange_surf_unfiltered = df_exchange['surface [kmol/s]'].values 
#             exchange_deep_unfiltered = df_exchange['deep [kmol/s]'].values
#             # 10-day hanning window filter (10 days = 240 hours)
#             exchange_surf = zfun.lowpass(exchange_surf_unfiltered, f='hanning', n=240) 
#             exchange_deep = zfun.lowpass(exchange_deep_unfiltered, f='hanning', n=240)
            exchange_color = 'blue'

# --------------------------- get TEF exchange flow terms ----------------------------------------
            # TODO: need to multiply by DO
            in_dir = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'tef2' / 'bulk_2014.01.01_2014.12.31' / 'lynchcove.nc'
            bulk = xr.open_dataset(in_dir)
            tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
            Q_p = tef_df['q_p'] # Qout
            Q_m = tef_df['q_m'] # Qin
            TEF_surf = Q_m.values
            TEF_deep = Q_p.values

# # ---------------------------------- get BGC terms --------------------------------------------
#             bgc_dir = Ldir['LOo'] / 'pugetsound_DO' / ('budget_' + startdate + '_' + enddate) / 'DO_bgc' / station
#             # get months
#             months = ['2014.01.01_2014.01.31',
#                       '2014.02.01_2014.02.28',
#                       '2014.03.01_2014.03.31',
#                       '2014.04.01_2014.04.30',
#                       '2014.05.01_2014.05.31',
#                       '2014.06.01_2014.06.30',
#                       '2014.07.01_2014.07.31',
#                       '2014.08.01_2014.08.31',
#                       '2014.09.01_2014.09.30',
#                       '2014.10.01_2014.10.31',
#                       '2014.11.01_2014.11.30',
#                       '2014.12.01_2014.12.31',]
            
#             # initialize arrays to save values
#             photo_surf_unfiltered = []
#             photo_deep_unfiltered = []
#             photo_color = 'forestgreen'
#             cons_surf_unfiltered = [] # nitrification, respiration
#             cons_deep_unfiltered = [] # nitrification, respiration, sediment oxygen demand
#             cons_color = 'black'
#             airsea_surf_unfiltered = []
#             airsea_color = 'deeppink'
#             o2vol_surf_unfiltered = []
#             o2vol_deep_unfiltered = []
#             ddtDOV_color = 'rebeccapurple'

#             # combine all months
#             for month in months:
#                 fn = 'O2_bgc_shallow_deep_' + month + '.nc'
#                 ds = xr.open_dataset(bgc_dir/fn, decode_times=False)
#                 # conversion factor to go from mmol O2/hr to kmol O2/s
#                 conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
#                 # get photosynthesis
#                 photo_surf_unfiltered = np.concatenate((photo_surf_unfiltered, ds['Oxy_pro_sum_shallow'].values * conv)) # kmol/s
#                 photo_deep_unfiltered = np.concatenate((photo_deep_unfiltered, ds['Oxy_pro_sum_deep'].values * conv)) # kmol/s
#                 # get consumption
#                 surf_cons_terms = ds['Oxy_nitri_sum_shallow'].values + ds['Oxy_remi_sum_shallow'].values
#                 deep_cons_terms = ds['Oxy_nitri_sum_deep'].values + ds['Oxy_remi_sum_deep'].values + ds['Oxy_sed_sum2'].values
#                 cons_surf_unfiltered = np.concatenate((cons_surf_unfiltered, surf_cons_terms * conv * -1)) # kmol/s; multiply by -1 b/c loss term
#                 cons_deep_unfiltered = np.concatenate((cons_deep_unfiltered, deep_cons_terms * conv * -1)) # kmol/s; multiply by -1 b/c loss term
#                 # get air-sea gas exchange
#                 airsea_surf_unfiltered = np.concatenate((airsea_surf_unfiltered, ds['Oxy_air_flux_sum'].values * conv)) # kmol/s
#                 # get (DO*V)
#                 o2vol_surf_unfiltered = np.concatenate((o2vol_surf_unfiltered, ds['Oxy_vol_sum_shallow'].values)) # mmol
#                 o2vol_deep_unfiltered = np.concatenate((o2vol_deep_unfiltered, ds['Oxy_vol_sum_deep'].values)) # mmol

#             # take time derivative of (DO*V) to get d/dt (DO*V)
#             ddtDOV_surf_unfiltered = np.diff(o2vol_surf_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
#             ddtDOV_deep_unfiltered = np.diff(o2vol_deep_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s

#             # apply 10-day hanning window
#             photo_surf = zfun.lowpass(photo_surf_unfiltered, f='hanning', n=240)
#             photo_deep = zfun.lowpass(photo_deep_unfiltered, f='hanning', n=240)
#             cons_surf = zfun.lowpass(cons_surf_unfiltered, f='hanning', n=240)
#             cons_deep = zfun.lowpass(cons_deep_unfiltered, f='hanning', n=240)
#             airsea_surf = zfun.lowpass(airsea_surf_unfiltered, f='hanning', n=240)
#             ddtDOV_surf = zfun.lowpass(ddtDOV_surf_unfiltered, f='hanning', n=240)
#             ddtDOV_deep = zfun.lowpass(ddtDOV_deep_unfiltered, f='hanning', n=240)

# ------------------------------- get rivers and WWTPs ----------------------------------------
            fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_traps' / (station + '.p')
            df_traps = pd.read_pickle(fn)
            print(df_traps)
            rivers_surf_unfiltered = df_traps['surface [m3/s]'].values
            wwtps_deep_unfiltered = df_traps['deep [m3/s]'].values
            # 10-day hanning window filter (10 days = 240 hours)
            # traps_surf = zfun.lowpass(rivers_surf_unfiltered, f='hanning', n=240)
            # traps_deep = zfun.lowpass(wwtps_deep_unfiltered, f='hanning', n=240)
            # Godin filter
            traps_surf = zfun.lowpass(rivers_surf_unfiltered, f='godin')[36:-34:24]
            traps_deep = zfun.lowpass(wwtps_deep_unfiltered, f='godin')[36:-34:24]
            traps_color = 'turquoise'

# # ------------------------------- get vertical exchange ----------------------------------------

#             # pick color
#             vertX_color = 'sandybrown'

#             # calculate raw error term, and acribe that the the vertical exchange
#             vertX_surf_unfiltered = ddtDOV_surf_unfiltered - (exchange_surf_unfiltered[0:-2]
#                                                               + photo_surf_unfiltered[0:-1]
#                                                               + cons_surf_unfiltered[0:-1]
#                                                               + airsea_surf_unfiltered[0:-1]
#                                                               + traps_surf[0:-1])
#             vertX_deep_unfiltered = ddtDOV_deep_unfiltered - (exchange_deep_unfiltered[0:-2]
#                                                               + photo_deep_unfiltered[0:-1]
#                                                               + cons_deep_unfiltered[0:-1]
#                                                               + traps_deep[0:-1])
            
#             # apply 10-day hanning window
#             vertX_surf = zfun.lowpass(vertX_surf_unfiltered, f='hanning', n=240)
#             vertX_deep = zfun.lowpass(vertX_deep_unfiltered, f='hanning', n=240)

# ---------------------------------- plot and save --------------------------------------------
            # plot surface
            # ax[0].plot(dates_local,exchange_surf,color=exchange_color,linewidth=1,linestyle='--',label='EU Exchange Flow')
            ax[0].plot(dates_local_daily[1:-1],TEF_surf,color=exchange_color,linewidth=1,linestyle='--',label='EU Exchange Flow')
            # ax[0].plot(dates_local[0:-1],photo_surf,color=photo_color,linewidth=2,label='Photosynthesis')
            # ax[0].plot(dates_local[0:-1],cons_surf,color=cons_color,linewidth=2,linestyle=':',label='Bio Consumption')
            # ax[0].plot(dates_local[0:-1],airsea_surf,color=airsea_color,linewidth=1,label='Air-Sea Transfer')
            # ax[0].plot(dates_local[0:-2],ddtDOV_surf,color=ddtDOV_color,linewidth=2,alpha=0.6,
            #            label=r'$\frac{\mathrm{d}}{\mathrm{dt}}(\mathrm{DO}\cdot V)$')
            ax[0].plot(dates_local_daily[1:-1],traps_surf,color=traps_color,linewidth=3,alpha=0.6,label='TRAPS')
            # ax[0].plot(dates_local[0:-2],vertX_surf,color=vertX_color,linewidth=2,label='Vertical Exchange')
            ax[0].legend(loc='best',ncol=3)
            
            # plot deep
            # ax[1].plot(dates_local,exchange_deep,color=exchange_color,linewidth=1,linestyle='--')
            ax[1].plot(dates_local_daily[1:-1],TEF_deep,color=exchange_color,linewidth=1,linestyle='--')
            # ax[1].plot(dates_local[0:-1],photo_deep,color=photo_color,linewidth=2,)
            # ax[1].plot(dates_local[0:-1],cons_deep,color=cons_color,linewidth=2,linestyle=':')
            # ax[1].plot(dates_local[0:-2],ddtDOV_deep,color=ddtDOV_color,linewidth=2,alpha=0.6)
            ax[1].plot(dates_local_daily[1:-1],traps_deep,color=traps_color,linewidth=3,alpha=0.6)
            # ax[1].plot(dates_local[0:-2],vertX_deep,color=vertX_color,linewidth=2)

            # plot error
            # ax.plot(dates_local[1:-1],vertX_surf+vertX_deep)