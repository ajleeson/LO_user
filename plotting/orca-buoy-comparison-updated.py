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

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# font sizes
fs_label = 15
fs_header = 15
fs_title = 18

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

# Get grid data
G = zrfun.get_basic_info('/home/aleeson/LO_data/grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset('/home/aleeson/LO_data/grids/cas7/grid.nc')
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
pfun.start_plot(fs=fs, figsize=(18,11))
fig = plt.figure()

# Create bathymetry plot --------------------------------------------------------------
ax = fig.add_subplot(1,2,1)
newcmap = cmocean.tools.crop_by_percent(cmocean.cm.deep_r, 10, which='max', N=None)
# newcmap = cmocean.cm.deep_r
newcmap.set_bad('#EEEEEE',1.) # background color
cs = ax.pcolormesh(plon, plat, zm, vmin=-300, vmax=0, cmap=newcmap)
cbar = plt.colorbar(cs,ax=ax, location='left')
cbar.ax.tick_params(labelsize=fs_label, rotation=30)
cbar.ax.set_ylabel('Depth [m]', fontsize=fs_label)
cbar.outline.set_visible(False)
pfun.dar(ax)
pfun.add_coast(ax,color='gray')
# Set axis limits
ax.set_xlim([-123.29,-122.1])
ax.set_ylim([46.93,48.35])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.axis('off')
# add title
ax.set_title('Orca Buoy Locations',fontsize=fs_header,fontweight='bold')
# add 10 km bar
lat0 = 47
lon0 = -122.33
lat1 = lat0
lon1 = -122.2
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=5)
ax.text((lon0+lon1)/2,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=fs_label,
        horizontalalignment='center')
# plot buoy locations
buoy_names = ['CI','PW','HP','TW'] #['CI','PW','NB','DB','HP','TW']
buoy_lon = [-122.7300,-122.3972,-123.1126,-123.0083] #[-122.7300,-122.3972,-122.6270,-122.8029,-123.1126,-123.0083]
buoy_lat = [47.2800,47.7612,47.4218,47.3750]#[47.2800,47.7612,47.9073,47.8034,47.4218,47.3750]
ax.scatter(buoy_lon,buoy_lat,s=100,color='violet',edgecolors='k')
for i,name in enumerate(buoy_names):
    # ax.text(buoy_lon[i]-0.04,buoy_lat[i]+0.02,name,fontsize='16',fontweight='bold',color='black')
    txt = ax.text(buoy_lon[i]-0.04,buoy_lat[i]+0.02,name,fontsize=fs_title,fontweight='bold',color='black')
    txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])

# Loop through and add Orca Buoy data and model output ----------------------------------------

# get orca data
orca_nc = ['PW','CI','HP','TW']#['TW','PW','NB','HP','DB','CI']
# fig numbers
fig_no = [2,4,6,8]#[3,4,7,8,11,12]
# label
labels = ['(a) ','(b) ','(c) ','(d) ']
# titles
titles = ['PW - Point Wells', 'CI - Carr Inlet', 'HP - Hoodsport', 'TW - Twanoh']#['TW - Twanoh','PW - Point Wells','NB - North Buoy', 'HP - Hoodsport', 'DB - Dabob Bay', 'CI - Carr Inlet']

# loop through all of the stations
for i,orca in enumerate(orca_nc):
    # create subplot
    ax = fig.add_subplot(4,2,fig_no[i])#fig.add_subplot(3,4,fig_no[i])
    
    # Add GRC case and create timevector
    ds_baseline = xr.open_dataset('/home/aleeson/LO_output/extract/cas6_traps2_x2b/moor/orca/'+orca+'_2017.01.01_2017.12.31.nc')
    DO_baseline = ds_baseline[vn][:,model_layer]*pinfo.fac_dict[vn] # bottom DO
    year_time = ds_baseline['ocean_time']

    # # Add prior model evaluation results
    # ds_WWTPDIN0 = xr.open_dataset('/home/aleeson/LO_output/extract/cas7_trapsV00_meV00_AugVFCinis/moor/orca/'+orca+'_2017.01.01_2017.08.15.nc')
    # DO_WWTPDIN0 = ds_WWTPDIN0[vn][:,model_layer]*pinfo.fac_dict[vn] # bottom DO
    
    # Add model evaluation results
    ds_WWTPDIN = xr.open_dataset('/home/aleeson/LO_output/extract/cas7_trapsV00_meV00/moor/orca/'+orca+'_2017.01.01_2017.12.31.nc')
    DO_WWTPDIN = ds_WWTPDIN[vn][:,model_layer]*pinfo.fac_dict[vn] # bottom DO
    
    # add title
    ax.set_title(labels[i] + titles[i],fontsize=fs_header,fontweight='bold')
    # format y labels
    ax.set_ylim([lim0,lim1])
    plt.yticks(fontsize=fs_label)
    ax.set_ylabel(layer_name + ' ' + title, fontsize=fs_label)
    # format x labels
    ax.set_xlim([year_time[0] - np.timedelta64(1,'D'),year_time[-1]])
    ax.xaxis.set_major_locator(MonthLocator())
    ax.grid(True,color='w',linewidth=2)
    if i == 3:
        plt.tick_params(axis='x',rotation=30)
        ax.xaxis.set_major_formatter(DateFormatter('%b'))
        plt.xticks(fontsize=fs_label)
    else:
        ax.set_xticklabels([])
    # format figure color
    ax.set_facecolor('#EEEEEE')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    # add observations
    ds_orca = xr.open_dataset('/home/aleeson/LO_data/obs/ORCA/LO_orca_moor/datasets/'+orca+'_ds.nc')
    orca_time = ds_orca.time.values
    orca_do = ds_orca[orca_vn].values[:,orca_layer]
    ax.plot(orca_time,orca_do,'o',color='mediumorchid',markersize=8,alpha=0.3,label='Observations*')

    # model output
    ax.plot(year_time,DO_baseline,linestyle='-',color='darkturquoise',linewidth=4, label='Current LiveOcean')
    # ax.plot(year_time[0:227],DO_WWTPDIN0,linestyle='-',color='k',linewidth=0.5,label='Prior results (trapsV00_meV00-VFC ICs)')
    ax.plot(year_time,DO_WWTPDIN,linestyle='-',color='k',linewidth=2,label='Updated model')

    # add legend
    if i == 0:
        ax.legend(loc='best',fontsize=fs_label,ncol=3)
    # add water depth
    ax.text(year_time[15],2.4, r'$z_{buoy}$' +' = {} m'.format(round(ds_orca.depth.values[orca_layer],1)),fontsize=fs_label)
    ax.text(year_time[15],0.6, r'$z_{model}$' +' = {} m'.format(round(-1*np.mean(ds_baseline['z_rho'][:,model_layer].values),1)),fontsize=fs_label)


# Generate plot
plt.subplots_adjust(wspace=0, hspace=0)
plt.suptitle('LiveOcean and ORCA buoy ' + layer_name + ' ' + var + ' comparison',fontweight='bold',fontsize=20)
plt.tight_layout
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.90, wspace=0.05, hspace=0.3)
plt.savefig(out_dir / ('orca_buoy_plot_updated.png'))