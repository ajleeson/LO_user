"""
Generate depth vs. time property plots using mooring extraction data. 
Used to look at 21 inlets in Puget Sound, and compare to Ecology monitoring stations, if available
"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import gsw
import pinfo
import pickle
from importlib import reload
reload(pinfo)

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

vn_list = ['oxygen']
rows = len(vn_list)

# figure settings
fs = 14 # figure font size
ls = 14 # label size
ts = 14 # title size

# basin color list
color_list = ['limegreen',  # whidbey
            'hotpink',    # hood canal
            'deepskyblue',  # main basin
            'blueviolet',       # south sound
            'black']        # admiralty inlet

# inlets in different basins
whidbey = ['similk','oak','crescent','penn','portsusan','holmes']
hoodcanal = ['dabob','lynchcove']
mainbasin = ['dyes','sinclair','elliot','quartermaster','commencement']
southsound = ['case','carr','hammersley','totten','eld','budd','henderson']
admiraltysill = ['killsut']

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
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures' / 'twentyone' / 'timeseries'
Lfun.make_dir(out_dir)

##########################################################
##                      Plotting                        ##
##########################################################

# create time vector
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]


# Initialize Figure
plt.close('all')
fig, ax = plt.subplots(1,1,figsize = (15,5))

# Loop through all of the mooring stations
for i,station in enumerate(sta_dict): # enumerate(['commencement']): #

    # get data
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)

    # get color
    if station in whidbey:
        color = color_list[0]
    elif station in hoodcanal:
        color = color_list[1]
    elif station in mainbasin:
        color = color_list[2]
    elif station in southsound:
        color = color_list[3]
    else:
        color = color_list[4]

    bott_ox_raw = ds['oxygen'].values[:,0] # 0 = bottom layer
    val = bott_ox_raw *  pinfo.fac_dict['oxygen'] # convert to mg/L

    val = zfun.lowpass(val-np.nanmean(val),f='hanning',n=30)

    # plot
    ax.plot(dates_local,val,linewidth=1,color=color,alpha=0.6)

    # plot time and location of max bott_DO
    max_DO = np.max(val)
    max_i = np.argmax(val)
    max_time = dates_local[max_i]
    # ax.scatter(max_time,max_DO,color=color,s=70,edgecolor='k',zorder=3)

# FORMATTING ---------------------------------------------------------
ax.set_ylabel(r'DO$_{bot}$ [mg $L^{-1}$]',fontsize=ls+2)
# ax.set_ylim([0,14])
# ax.set_yticks(np.arange(0, 16, 2))
ax.set_xlim([dates_local[0],dates_local[-1]])
ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax.tick_params(axis='x', labelrotation=30, labelsize=ls)
ax.set_xlabel(year, fontsize = fs)
ax.set_title(r'DO$_{bot}$ time series in Puget Sound terminal inlets',fontsize=ts+2)

# add basin label
ax.text(0.97, 0.91, 'Whidbey Basin',color=color_list[0],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.97, 0.87, 'Hood Canal',color=color_list[1],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.97, 0.83, 'Main Basin',color=color_list[2],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.97, 0.79, 'South Sound',color=color_list[3],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.97, 0.75, 'Admiralty Inlet',color=color_list[4],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')

# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.subplots_adjust(bottom=0.16)
plt.savefig(out_dir / '21inlet_DO_timeseries.png')
