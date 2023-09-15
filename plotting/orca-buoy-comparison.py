"""
Plots a comparison of bottom DO for two identical grids with different oxygen value.

This is a custom function for a particular experiment, but code can be adapted for other use cases in the future.

From the terminal: python orca-buoy-comparison.py

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

# Get grid data
G = zrfun.get_basic_info('/home/aleeson/LO_data/grids/cas6/grid.nc', only_G=True)
grid_ds = xr.open_dataset('/home/aleeson/LO_data/grids/cas6/grid.nc')
z = -grid_ds.h.values
mask_rho = np.transpose(grid_ds.mask_rho.values)
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
X = lon[0,:] # grid cell X values
Y = lat[:,0] # grid cell Y values
plon, plat = pfun.get_plon_plat(lon,lat)
# make a version of z with nans where masked
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan

# Initialize figure
plt.close('all')
fs = 10
pfun.start_plot(fs=fs, figsize=(18,10))
fig = plt.figure()

# Create bathymetry plot --------------------------------------------------------------
ax = fig.add_subplot(1,2,1)
newcmap = cmocean.tools.crop_by_percent(cmocean.cm.deep_r, 10, which='max', N=None)
cs = ax.pcolormesh(plon, plat, zm, vmin=-300, vmax=0, cmap=newcmap)
cbar = plt.colorbar(cs,ax=ax, location='left')
cbar.ax.tick_params(labelsize=14, rotation=30)
cbar.ax.set_ylabel('Depth [m]', fontsize=14)
pfun.dar(ax)
# pfun.add_coast(ax)
# Set axis limits
ax.set_xlim([-123.3,-122.1])
ax.set_ylim([46.93,48.45])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.axis('off')
# add title
ax.set_title('Orca Buoy Locations',fontsize=16)
# add 50 km bar
lat0 = 46.94
lon0 = -123.05
lat1 = lat0
lon1 = -122.91825
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=5)
ax.text(lon0,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=14)
# plot buoy locations
buoy_names = ['CI','PW','NB','DB','HP','TW']
buoy_lon = [-122.7300,-122.3972,-122.6270,-122.8029,-123.1126,-123.0083]
buoy_lat = [47.2800,47.7612,47.9073,47.8034,47.4218,47.3750]
ax.scatter(buoy_lon,buoy_lat,s=100,color='deeppink')
for i,name in enumerate(buoy_names):
    ax.text(buoy_lon[i]-0.04,buoy_lat[i]+0.02,name,fontsize='16',fontweight='bold',color='black')

# Loop through and add Orca Buoy data and model output ----------------------------------------

# get orca data
orca_nc = ['TW','PW','NB','HP','DB','CI']
# fig numbers
fig_no = [3,4,7,8,11,12]
# titles
titles = ['TW - Twanoh','PW - Point Wells','NB - North Buoy', 'HP - Hoodsport', 'DB - Dabob Bay', 'CI - Carr Inlet']

# loop through all of the stations
for i,orca in enumerate(orca_nc):
    # create subplot
    ax = fig.add_subplot(3,4,fig_no[i])
    # Add GRC case and create timevector
    ds_baseline = xr.open_dataset('/home/aleeson/LO_output/extract/cas6_traps2_x2b/moor/orca/'+orca+'_2017.01.01_2017.12.31.nc')
    DO_baseline = ds_baseline['oxygen'][:,0]*pinfo.fac_dict['oxygen'] # bottom DO
    year_time = ds_baseline['ocean_time']
    ax.plot(year_time,DO_baseline,linestyle='-',color='mediumturquoise',linewidth=2, label='GRC results (traps2_x2b)')
    # Add model evaluation results
    ds_WWTPDIN = xr.open_dataset('/home/aleeson/LO_output/extract/cas7_trapsV00_meV00/moor/orca/'+orca+'_2017.01.01_2017.08.15.nc')
    DO_WWTPDIN = ds_WWTPDIN['oxygen'][:,0]*pinfo.fac_dict['oxygen'] # bottom DO
    ax.plot(year_time[0:227],DO_WWTPDIN,linestyle='--',color='k',linewidth=1,label='New results (trapsV00_meV00)')
    # add title
    ax.set_title(titles[i] + ' Bottom DO [mg/L]',fontsize=16)
    # format y labels
    ax.set_ylim([0,12])
    if np.mod(i,2) != 0:
        ax.set_yticklabels([])
    else:
        plt.yticks(fontsize=12)
    # format x labels
    ax.set_xlim([year_time[0] - np.timedelta64(1,'D'),year_time[-1]])
    ax.xaxis.set_major_locator(MonthLocator())
    ax.grid(True,color='w',linewidth=2)
    if i < 4:
        ax.set_xticklabels([])
    else:
        plt.tick_params(axis='x',rotation=30)
        ax.xaxis.set_major_formatter(DateFormatter('%b'))
        plt.xticks(fontsize=12)
    # format figure color
    ax.set_facecolor('#EEEEEE')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    # add observations
    ds_orca = xr.open_dataset('/home/aleeson/LO_data/obs/ORCA/LO_orca_moor/datasets/'+orca+'_ds.nc')
    orca_time = ds_orca.time.values
    orca_do = ds_orca.oxy.values[:,-1]
    ax.plot(orca_time,orca_do,'o',color='darkmagenta',markersize=6,alpha=0.3,label='Observations*')
    # add legend
    if i == 0:
        ax.legend(loc='upper right',fontsize=12)
    # add water depth
    ax.text(year_time[15],1.2,'buoy depth = {} m'.format(round(ds_orca.depth.values[-1],1)),fontsize=12)
    ax.text(year_time[15],0.5,'model depth = {} m'.format(round(-1*np.mean(ds_baseline['z_rho'][:,0].values),1)),fontsize=12)


# Generate plot
plt.suptitle('LiveOcean and ORCA buoy bottom DO comparison',fontweight='bold',fontsize=20)
plt.tight_layout
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.90, wspace=0.05, hspace=0.3)
plt.savefig('orca_buoy_plot.png')