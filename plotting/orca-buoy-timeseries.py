"""
Plots orca timeseries.

Written to compare 2017 observations to other years

From the terminal: python orca-buoy-timeseries.py

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import pandas as pd
import cmocean
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# font sizes
fs_label = 12
fs_header = 12
fs_title = 14

# variable
vn = 'oxygen' # salt, NO3, oxygen
layer = 'bottom'


# Get correct variable information
if vn == 'oxygen':
    orca_vn = 'oxy'
    lim0 = 0
    lim1 = 12
    title = 'DO [mg/L]'
    var = 'DO'
elif vn == 'salt':
    orca_vn = 'sal'
    lim0 = 28
    lim1 = 32
    title = 'salinity'
    var = 'salt'
elif vn == 'NO3':
    orca_vn = 'nitrate'
    lim0 = 0
    lim1 = 60
    title = 'NO3 [uM]'
    var = 'NO3'

# get correct surface/bottom layer
if layer =='bottom':
    model_layer = 0
    orca_layer = -1
    layer_name = 'Bottom'
elif layer =='surface':
    model_layer = -1
    orca_layer = 0
    layer_name = 'Surface'

# Initialize figure
plt.close('all')
fs = 10
pfun.start_plot(fs=fs, figsize=(10,7))
fig = plt.figure()

###############################################################
##                  Add orca timeseries                      ##
###############################################################


# get orca data
orca_nc = ['PW','CI','HP','TW']
# fig numbers
fig_no = [1,2,3,4]#[3,4,7,8,11,12]
# label
labels = ['(a) ','(b) ','(c) ','(d) ']
# titles
titles = ['PW - Point Wells', 'CI - Carr Inlet', 'HP - Hoodsport', 'TW - Twanoh']#['TW - Twanoh','PW - Point Wells','NB - North Buoy', 'HP - Hoodsport', 'DB - Dabob Bay', 'CI - Carr Inlet']
# years
# years = [2017,2018,2019,2020,2021]
years = [2021,2020,2019,2018,2017]
colors = ['grey','grey','grey','grey','deeppink']

yrday = pd.date_range(start ='1/1/2020', end ='12/31/2020', freq ='D')

# loop through all of the stations
for i,orca in enumerate(orca_nc):
    # create subplot
    ax = fig.add_subplot(4,1,fig_no[i])
    
    # add title
    ax.text(0.05,0.8,labels[i] + titles[i],transform=ax.transAxes,fontsize=fs_header,fontweight='bold')
    # format y labels
    ax.set_ylim([lim0,lim1])
    plt.yticks(fontsize=fs_label)
    # format x labels
    ax.set_xlim([yrday[0],yrday[-1]])
    ax.xaxis.set_major_formatter(DateFormatter('%b'))
    ax.xaxis.set_major_locator(MonthLocator())
    ax.grid(True,color='w',linewidth=2)
    if i == 3:
        plt.tick_params(axis='x',rotation=30)
        plt.xticks(fontsize=fs_label)
    else:
        ax.set_xticklabels([])
    # format figure color
    ax.set_facecolor('#EEEEEE')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)

    # get orca observations
    ds_orca = xr.open_dataset('/home/aleeson/LO_data/obs/ORCA/LO_orca_moor/datasets/'+orca+'_ds.nc')
    orca_time = ds_orca.time.values
    orca_vals = ds_orca[orca_vn].values[:,orca_layer]

    # convert dataset to dataframe
    df = pd.DataFrame() 
    df['time'] = orca_time
    df[orca_vn] = orca_vals
    # get dates
    df['year'] = pd.DatetimeIndex(df['time']).year
    df['day'] = pd.DatetimeIndex(df['time']).day_of_year

    # plot observations
    for y,year in enumerate(years):
        values_to_plot = df[orca_vn].loc[df['year'] == year].values
        values_to_plot = values_to_plot.tolist()
        # skip leap years
        if np.mod(year,4) != 0:
            # pad Feb 29th with nan
            values_to_plot = values_to_plot[0:60] + [np.nan] + values_to_plot[60::]
        ax.plot(yrday,values_to_plot,label=year,marker='o',linestyle='none',alpha=0.5, color = colors[y])

    # add legend
    if i == 0:
        ax.legend(loc='upper right',fontsize=fs_label,ncol=1)

    # add legend
    ax.text(yrday[110],9.8, r'$z_{mooring}$' +' = {} m'.format(round(ds_orca.depth.values[orca_layer],1)),fontsize=fs_label)
 
# Generate plot
plt.subplots_adjust(wspace=0, hspace=0.1)
plt.suptitle('ORCA mooring ' + layer_name + ' ' + title,fontweight='bold',fontsize=16)
plt.tight_layout
plt.subplots_adjust(left=0.07, bottom=0.06, right=0.95, top=0.91, wspace=0, hspace=0.15)
plt.savefig(out_dir / ('orca_'+layer+'_'+var+'_timeseries.png'))

###############################################################
##                     Hoodsport Only                        ##
###############################################################
plt.close('all')

orca_layer = -42
orca = 'HP'

pfun.start_plot(fs=fs, figsize=(6,4))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

# add title
ax.text(yrday[10],11,'Hoodsport',fontsize=fs_header,fontweight='bold')
# format y labels
ax.set_ylim([lim0,lim1])
plt.yticks(fontsize=fs_label)
# format x labels
ax.set_xlim([yrday[0],yrday[-1]])
ax.xaxis.set_major_formatter(DateFormatter('%b'))
ax.xaxis.set_major_locator(MonthLocator())
ax.grid(True,color='w',linewidth=2)
if i == 3:
    plt.tick_params(axis='x',rotation=30)
    plt.xticks(fontsize=fs_label)
else:
    ax.set_xticklabels([])
# format figure color
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

# get orca observations
ds_orca = xr.open_dataset('/home/aleeson/LO_data/obs/ORCA/LO_orca_moor/datasets/'+orca+'_ds.nc')
orca_time = ds_orca.time.values
orca_vals = ds_orca[orca_vn].values[:,orca_layer]

# convert dataset to dataframe
df = pd.DataFrame() 
df['time'] = orca_time
df[orca_vn] = orca_vals
# get dates
df['year'] = pd.DatetimeIndex(df['time']).year
df['day'] = pd.DatetimeIndex(df['time']).day_of_year

# plot observations
for y,year in enumerate(years):
    values_to_plot = df[orca_vn].loc[df['year'] == year].values
    values_to_plot = values_to_plot.tolist()
    # skip leap years
    if np.mod(year,4) != 0:
        # pad Feb 29th with nan
        values_to_plot = values_to_plot[0:60] + [np.nan] + values_to_plot[60::]
    ax.plot(yrday,values_to_plot,label=year,marker='o',linestyle='none',alpha=0.5, color = colors[y])

# add legend
ax.legend(loc='upper right',fontsize=fs_label,ncol=1)
# add depth
ax.text(yrday[110],11, r'$z_{mooring}$' +' = {} m'.format(round(ds_orca.depth.values[orca_layer],1)),fontsize=fs_label)

plt.title('ORCA mooring '+ title,fontweight='bold',fontsize=14)

plt.show()