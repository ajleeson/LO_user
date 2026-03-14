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

print('\n')

##########################################################
##              Plot DO budget of every inlet           ##
##########################################################

for i,station in enumerate(['lynchcove']): # enumerate(sta_dict):
        # print status
        print('({}/{}) Working on {}...'.format(i+1,len(sta_dict),station))

        # initialize figure
        plt.close('all')
        fig, ax = plt.subplots(2,1,figsize = (13,6),sharex=True)
        # format figure
        plt.suptitle(station + ': Eulerian Volume Budget (10-day hanning window)',size=14)
        for axis in [ax[0],ax[1]]:
            # axis.plot([dates_local[0],dates_local[-1]],[0,0],color='k')
            axis.set_xlim([dates_local[0],dates_local[-1]])
            axis.set_ylim([-60,60])
            axis.set_ylabel(r'Q [$m^3\ s^{-1}$]')
            axis.set_facecolor('#EEEEEE')
            axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
            axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
            axis.tick_params(axis='x', labelrotation=30)
            for border in ['top','right','bottom','left']:
                axis.spines[border].set_visible(False)
        ax[1].set_xlabel(year)
        ax[0].set_title('(a) Volume Budget',loc='left')
        ax[1].set_title(r'(b) Error = $\frac{\mathrm{dV}}{\mathrm{dt}}$ - (Exchange + TRAPS)', loc='left')

# --------------------------- get exchange flow terms ----------------------------------------
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / 'EU_exchange_flow' / (station + '.p')
        df_exchange = pd.read_pickle(fn)
        exchange_unfiltered = df_exchange['total [m3/s]'].values 
        # 10-day hanning window filter (10 days = 240 hours)
        exchange = zfun.lowpass(exchange_unfiltered, f='hanning', n=240)
        exchange_color = 'blue'

# --------------------------- get storage terms ----------------------------------------
        fn = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'tef2' / 'segments_2014.01.01_2014.12.31_cas7_c21_traps00.nc'
        ds_storage = xr.open_dataset(fn)
        volume = ds_storage.loc[dict(seg='lynchcove_p')]['volume'].values # m3
        dVdt_unfiltered = np.diff(volume) * (1/60) * (1/60) # m3/hr * (1/3600) = m3/s
        # 10-day hanning window filter (10 days = 240 hours)
        dVdt = zfun.lowpass(dVdt_unfiltered, f='hanning', n=240)
        dVdt_color = 'rebeccapurple'

# ------------------------------- get rivers and WWTPs ----------------------------------------
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / 'traps' / (station + '.p')
        df_traps = pd.read_pickle(fn)
        traps_unfiltered = df_traps['total [m3/s]'].values
        # 10-day hanning window filter (10 days = 240 hours)
        traps = zfun.lowpass(traps_unfiltered, f='hanning', n=240)
        traps_color = 'turquoise'

# ---------------------------------- plot and save --------------------------------------------
        # plot budget
        ax[0].plot([dates_local[0],dates_local[-1]],[0,0],'k-',linewidth=0.5)
        ax[0].plot(dates_local,exchange,color=exchange_color,linewidth=1,linestyle='--',label='Eulerian Exchange Flow')
        ax[0].plot(dates_local[0:-1],traps,color=traps_color,linewidth=3,alpha=0.6,label='TRAPS')
        ax[0].plot(dates_local[0:-1],dVdt,color=dVdt_color,linewidth=2,alpha=0.6,
                    label=r'Storage = $\frac{\mathrm{dV}}{\mathrm{dt}}$')
        ax[0].legend(loc='best')
        
        # plot error
        print(np.nanmax(traps))
        error = dVdt - (exchange[0:-1] + traps)
        ax[1].plot([dates_local[0],dates_local[-1]],[0,0],'k-',linewidth=0.5)
        ax[1].plot(dates_local[0:-1],error)
