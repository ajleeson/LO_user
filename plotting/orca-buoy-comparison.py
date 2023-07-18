"""
Plots a comparison of bottom DO for two identical grids with different oxygen value.

This is a custom function for a particular experiment, but code can be adapted for other use cases in the future.

From the terminal: python money_plot.py

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

# Loop through and add Orca Buoy data and model output ----------------------------------------

# get orca data
orca_nc = ['TW_ds.nc','PW_ds.nc','NB_ds.nc','HP_ds.nc','DB_ds.nc','CI_ds.nc']
# fig numbers
fig_no = [3,4,7,8,11,12]
# titles
titles = ['TW - Twanoh','PW - Point Wells','NB - North Buoy', 'HP - Hoodsport', 'DB - Dabob Bay', 'CI - Carr Inlet']

# read in data model data at penn cove
ds_penn_withDIN = xr.open_dataset('/home/aleeson/LO_output/extract/cas6_traps2_x2b/moor/pennlynch/PennCove_2017.01.01_2017.12.30.nc')
DO_penn_withDIN = ds_penn_withDIN['oxygen'][:,0]*pinfo.fac_dict['oxygen'] # bottom DO
ds_penn_base = xr.open_dataset('/home/aleeson/LO_output/extract/cas6_traps3_x2b/moor/pennlynch/PennCove_2017.01.01_2017.12.31.nc')
DO_penn_base = ds_penn_base['oxygen'][:,0]*pinfo.fac_dict['oxygen'] # bottom DO
# get time vectors
time_withDIN = ds_penn_withDIN['ocean_time']
time_base = ds_penn_base['ocean_time']

# loop through all of the stations
for i,orca in enumerate(orca_nc):
    # create subplot
    ax = fig.add_subplot(3,4,fig_no[i])
    # add title
    ax.set_title(titles[i] + ' Bottom DO [mg/L]',fontsize=16)
    # format y labels
    ax.set_ylim([0,12])
    if np.mod(i,2) != 0:
        ax.set_yticklabels([])
    else:
        plt.yticks(fontsize=12)
    # format x labels
    ax.set_xlim([time_base[0] - np.timedelta64(1,'D'),time_base[-1]])
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
    # ax.plot(time_base,DO_lynch_base,linestyle='-',color='mediumturquoise',linewidth=2, label='Baseline [no WW discharge]')
    # ax.plot(time_withDIN,DO_lynch_withDIN,linestyle='--',color='k',linewidth=1,label='With WW discharge')
    # add observations
    ds_orca = xr.open_dataset('/home/aleeson/LO_data/obs/ORCA/LO_orca_moor/datasets/'+orca)
    orca_time = ds_orca.time.values
    orca_do = ds_orca.oxy.values[:,-1]
    ax.plot(orca_time,orca_do,'o',color='darkmagenta',markersize=6,alpha=0.3,label='Observations*')
    # add legend
    if i == 0:
        ax.legend(loc='upper right',fontsize=12)


# # Lynch Cove
# ds_lynch_withDIN = xr.open_dataset('/home/aleeson/LO_output/extract/cas6_traps2_x2b/moor/pennlynch/LynchCove_2017.01.01_2017.12.30.nc')
# DO_lynch_withDIN = ds_lynch_withDIN['oxygen'][:,0]*pinfo.fac_dict['oxygen'] # bottom DO
# ds_lynch_base = xr.open_dataset('/home/aleeson/LO_output/extract/cas6_traps3_x2b/moor/pennlynch/LynchCove_2017.01.01_2017.12.31.nc')
# DO_lynch_base = ds_lynch_base['oxygen'][:,0]*pinfo.fac_dict['oxygen'] # bottom DO
# ax = fig.add_subplot(3,4,3)
# ax.set_title('TW - Twanoh Bottom DO [mg/L]',fontsize=16)
# ax.set_ylim([0,12])
# ax.set_xticklabels([])
# ax.set_xlim([time_base[0] - np.timedelta64(1,'D'),time_base[-1]])
# ax.plot(time_base,DO_lynch_base,linestyle='-',color='mediumturquoise',linewidth=2, label='Baseline [no WW discharge]')
# ax.plot(time_withDIN,DO_lynch_withDIN,linestyle='--',color='k',linewidth=1,label='With WW discharge')
# # Add orca buoy data
# ds_orca = xr.open_dataset('/home/aleeson/LO_data/obs/ORCA/LO_orca_moor/datasets/TW_ds.nc')
# orca_time = ds_orca.time.values
# orca_do = ds_orca.oxy.values[:,-1]
# ax.plot(orca_time,orca_do,'o',color='darkmagenta',markersize=6,alpha=0.5,label='Observations*')
# ax.legend(loc='upper right',fontsize=12)
# ax.set_facecolor('#EEEEEE')
# for border in ['top','right','bottom','left']:
#     ax.spines[border].set_visible(False)
# ax.tick_params(axis='x', colors='w')
# ax.xaxis.set_major_locator(MonthLocator())
# ax.grid(True,color='w',linewidth=2)
# plt.yticks(fontsize=14)
# # add water depth
# ax.text(time_base[15],0.5,'depth = {} m'.format(round(-1*np.mean(ds_lynch_base['z_rho'][:,0].values),1)),fontsize=12)



# Generate plot
plt.tight_layout
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.05, hspace=0.3)
plt.savefig('orca_buoy_plot.png')
# plt.show()