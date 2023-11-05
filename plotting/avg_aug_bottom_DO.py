"""
Plot average August bottom DO in Puget Sound based on box extraction results

From ipython: run avg_aug_bottom_DO.py

"""

###################################################################
##                       import packages                         ##  
###################################################################      

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

###################################################################
##      load output folder, grid data, box extraction data       ##  
################################################################### 

gtagex = 'cas7_trapsV00_meV00'

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# Get grid data
G = zrfun.get_basic_info('/home/aleeson/LO_data/grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset('/home/aleeson/LO_data/grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values
px, py = pfun.get_plon_plat(G['lon_rho'],G['lat_rho'])

# layer, text, axis limits
vn = 'oxygen'
slev = 0
stext = 'Bottom'
vmin = 0
vmax = 10
# scale variable
scale =  pinfo.fac_dict[vn]

# lon/lat limits (Puget Sound)
xmin = -123.2
xmax = -122.1
ymin = 46.93
ymax = 48.45

# get box extraction data
fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / 'prelimDO_2017.08.01_2017.08.31.nc'
ds = xr.open_dataset(fp)

###################################################################
##                  Plotting and saving figure                   ##  
################################################################### 

# Initialize figure
fs = 10
pfun.start_plot(fs=fs, figsize=(8,11))
fig = plt.figure()

# create custom colormap
# newcmap = cmocean.tools.crop(cmocean.cm.amp_r, vmin = 0, vmax = 5, pivot= 2)
newcmap = cmocean.cm.oxy
# newcmap.set_bad('#ffffff',1.) # background color

# plot values
ax = fig.add_subplot(1,1,1)
v = ds[vn][:,slev,:,:].values * scale
# average over the month
v = np.nanmean(v,axis=0)
# only save values in that are hypoxic (so >2 mg/L is set to arbitrarily high value)
# v[v > 2] = 1000
cs = ax.pcolormesh(ds.coords['lon_rho'],ds.coords['lat_rho'],v,
                    vmin=vmin, vmax=10, cmap=newcmap)
                    # vmin=vmin, vmax=2, cmap=newcmap)
cbar = fig.colorbar(cs)
cbar.ax.tick_params(labelsize=14)
cbar.outline.set_visible(False)
# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.axis('off')
# pfun.add_coast(ax, color='k')
pfun.dar(ax)
ax.set_title('Average August Bottom DO [mg/L]', fontsize=16, fontweight='bold')

# add 10 km bar
lat0 = 47
lon0 = -122.4
lat1 = lat0
lon1 = -122.27
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=6)
ax.text((lon0+lon1)/2,lat0+0.02,'{} km'.format(x_dist_km),color='k',
        horizontalalignment='center', fontsize=15)

# Generate plot
plt.tight_layout
plt.savefig(out_dir / ('avg_aug_bottom_DO.png'))