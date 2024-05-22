
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

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots' / 'numod_project'
Lfun.make_dir(out_dir)

## Define tags to loop through
tags = ['amtc','amtk','ahtc','ahtk']

colors = ['red','red','navy','navy']
styles = ['-',':','-',':']

location = 'near_river'

# font sizes
fs_label = 12
fs_header = 12
fs_title = 14
pfun.start_plot(fs=10, figsize=(8,3))
fig,ax = plt.subplots(1,1)

# Loop through all model runs, get data, and plot
for i,tag in enumerate(tags):

    # get executable
    if 'c' in tag:
        ex = '_xatc'
    else:
        ex = '_xatk'

    # get data
    ds = xr.open_dataset(Ldir['LOo'] / 'extract' / ('hcal_' + tag + ex) / 'moor' / 'hcal' / (location + '_2020.01.01_2020.01.15.nc'))
    
    # get salinity along transect
    s_top = ds.salt[:,-1].values
    s_bot = ds.salt[:,0].values

    # calcualte stratification
    S = s_bot.reshape(-1) - s_top.reshape(-1)

    # create time vector
    dates_local = [pfun.get_dt_local(pd.Timestamp(x)) for x in ds.ocean_time.values]

    # plot stratification along transect
    plt.plot(dates_local,S,label=tag,alpha=0.5,linewidth=2,color=colors[i],linestyle=styles[i])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b-%d"))
    plt.legend()

    # format figure
    # plt.ylim([10,17.5])
    ax.set_xlim([dates_local[0],dates_local[-1]])
    plt.title(r'$s_{bot}-s_{top}$' + ' time-series')
    plt.ylabel(r'Stratification [g kg$^{-1}$]')
    plt.xlabel('Date')
    # format figure
    ax.grid(visible=True, color='w')
    # format background color
    ax.set_facecolor('#EEEEEE')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)

plt.savefig(out_dir / 'stratification_time_series.png')

# plot sea surface height ---------------------------------------------


plt.close('all')

# font sizes
fs_label = 12
fs_header = 12
fs_title = 14
pfun.start_plot(fs=10, figsize=(8,4))
fig,ax = plt.subplots(1,1)

# plot stratification along transect
plt.plot(dates_local,ds.zeta.values,alpha=0.5,linewidth=2,color='k')
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b-%d"))

# format figure
# plt.ylim([10,17.5])
ax.set_xlim([dates_local[0],dates_local[-1]])
plt.title('Sea surface height time-series')
plt.ylabel('SSH [m]')
plt.xlabel('Date')
# format figure
ax.grid(visible=True, color='w')
# format background color
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig(out_dir / 'SSH_time_series.png')