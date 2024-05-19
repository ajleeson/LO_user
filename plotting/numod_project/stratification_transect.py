
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
out_dir = Ldir['LOo'] / 'AL_custom_plots' / 'numod_project'
Lfun.make_dir(out_dir)

## Define tags to loop through
tags = ['amtc','amtk','ahtc','ahtk']

colors = ['red','red','navy','navy']
styles = ['-',':','-',':']

date = '2020.01.01'
hn = '25'

# font sizes
fs_label = 12
fs_header = 12
fs_title = 14
pfun.start_plot(fs=10, figsize=(8,5))
fig,ax = plt.subplots(1,1)

# Loop through all model runs, get data, and plot
for i,tag in enumerate(tags):

    # get executable
    if 'c' in tag:
        ex = '_xatc'
    else:
        ex = '_xatk'

    # get data
    ds = xr.open_dataset(Ldir['roms_out'] / ('hcal_' + tag + ex) / ('f' + date) / ('ocean_his_00' + hn + '.nc'))

    start_ind = 23

    # get transect lats
    lat = ds.lat_rho.values[start_ind:222,95]

    # get transect distance
    x,y = zfun.ll2xy(0, lat, 0, lat[0])
    
    # get salinity along transect
    s_top = ds.salt[:,-1,start_ind:222,95].values
    s_bot = ds.salt[:,0,start_ind:222,95].values

    # calcualte stratification
    S = s_bot.reshape(-1) - s_top.reshape(-1)

    # plot stratification along transect
    plt.plot(y/1000,S,label=tag,alpha=0.5,linewidth=2,color=colors[i],linestyle=styles[i])
    plt.legend()

    # format figure
    plt.ylim([0,30])
    plt.xlim([0,100])
    plt.title(r'$s_{bot}-s_{top}$' + '\n' + date + ' hour ' + hn)
    plt.ylabel(r'Stratification [g kg$^{-1}$]')
    plt.xlabel('Distance [km]')
    # format figure
    ax.grid(visible=True, color='w')
    # format background color
    ax.set_facecolor('#EEEEEE')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)

plt.savefig(out_dir / (date + '_hr' + hn + '_stratification.png'))