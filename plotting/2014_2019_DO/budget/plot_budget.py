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

        # get interface depth from csv file
        with open('interface_depths.csv', 'r') as f:
            for line in f:
                inlet, z_interface = line.strip().split(',')
        z_interface = float(z_interface)

        # one or two layers?
        if math.isnan(z_interface):
            # one layer
            print('One layer....make code to deal with this case')
        
        # two layers -----------------------------------------------
        else:

            # initialize figure
            plt.close('all')
            fig, ax = plt.subplots(2,1,figsize = (13,6),sharex=True)
            # format figure
            plt.suptitle(station + ': DO Budget (10-day hanning window)',size=14)
            for axis in [ax[0],ax[1]]:
                axis.plot([dates_local[0],dates_local[-1]],[0,0],color='k')
                axis.set_xlim([dates_local[0],dates_local[-1]])
                axis.set_ylim([-120,120])
                axis.set_ylabel(r'DO transport [$mol \ O_2 \ s^{-1}$]')
                axis.set_facecolor('#EEEEEE')
                axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
                axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
                axis.tick_params(axis='x', labelrotation=30)
                for border in ['top','right','bottom','left']:
                    axis.spines[border].set_visible(False)
            ax[1].set_xlabel(year)
            ax[0].set_title('(a) Surface [shallower than {} m]'.format(-1*z_interface),
                            loc='left')
            ax[1].set_title('(b) Bottom [deeper than {} m]'.format(-1*z_interface),
                            loc='left')

            # get exchange flow terms
            fn = Ldir['LOo'] / 'pugetsound_DO' / ('budget_' + startdate + '_' + enddate) / 'DO_exchange_flow' / (station + '.p')
            df_exchange = pd.read_pickle(fn)
            exchange_surf = df_exchange['surface [mol/s]']
            exchange_deep = df_exchange['deep [mol/s]']
            exchange_color = 'mediumorchid'

            # plot surface
            ax[0].plot(dates_local,exchange_surf,color=exchange_color,label='Exchange Flow')
            ax[0].legend(loc='upper right')
            # plot deep
            ax[1].plot(dates_local,exchange_deep,color=exchange_color)
