"""
Code to draw plot regional bathymetry
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import xarray as xr
from matplotlib.ticker import MaxNLocator
import cmocean
import cmcrameri.cm as cmc
import pandas as pd
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()


# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
z = -ds.h.values
mask_rho = np.transpose(ds.mask_rho.values)
lon = ds.lon_rho.values
lat = ds.lat_rho.values
X = lon[0,:] # grid cell X values
Y = lat[:,0] # grid cell Y values
plon, plat = pfun.get_plon_plat(lon,lat)
# make a version of z with nans where masked
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
depth = -zm

# Create bathymetry plot --------------------------------------------------------------

# Initialize figure
plt.close('all')
fs = 10
pfun.start_plot(fs=fs, figsize=(12.3,8.5))
fig = plt.figure()
plt.subplots_adjust(wspace=0, hspace=0)
# create colormap
# newcmap = cmocean.tools.crop_by_percent(cmocean.cm.deep, 20, which='min')
# newcmap = cmc.tokyo_r
newcmap = cmc.batlow_r
newcmap.set_bad(color='white')

# Model Domain ----------------------------------------------------------
ax0 = fig.add_subplot(1,3,1)
# cs = ax0.pcolormesh(plon, plat, zm, vmin=-5, vmax=0, cmap=newcmap)
cs = ax0.pcolormesh(plon, plat, depth, vmin=0, vmax=800, cmap=newcmap)
cbar_ax = fig.add_axes([0.117, 0.14, 0.515, 0.02])
cbar = plt.colorbar(cs, cax=cbar_ax, orientation='horizontal')
cbar.ax.tick_params(labelsize=18)
cbar.set_label('Depth [m]', fontsize=18)
cbar.outline.set_visible(False)
pfun.dar(ax0)
pfun.add_coast(ax0, color='silver')
ax0.set_xlim(-130, -121.8)
ax0.set_ylim(42, 52)

# plt.xticks(rotation=30, color='gray')
plt.xticks(rotation=30,horizontalalignment='right',fontsize=18)
plt.yticks(fontsize=18, rotation=30)
ax0.set_ylabel('Latitude',  fontsize=14)
ax0.set_xlabel('Longitude', fontsize=14)
ax0.xaxis.set_major_locator(MaxNLocator(integer=True))
ax0.yaxis.set_major_locator(MaxNLocator(integer=True))
ax0.tick_params(axis='both', labelrotation=45, labelsize=12)
# add title
# ax0.set_title('(a) Model Domain',fontsize=20, loc='left', fontweight='bold')#,color='#EEEEEE')
ax0.text(0.06, 0.95, '(a)', transform=ax0.transAxes, fontsize=16, fontweight='bold')
# add distance bar
lat0 = 42.6
lon0 = -123.6
lat1 = lat0
lon1 = -122.988875
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax0.plot([lon0,lon1],[lat0,lat1],color='black',linewidth=1, marker='|')
ax0.text((lon0+lon1)/2,lat0+0.15,'{} km'.format(x_dist_km),color='k',fontsize=12,
         horizontalalignment='center')

# draw box around study domain
bordercolor = 'deeppink'#'#EEEEEE'
ax0.add_patch(Rectangle((-126, 45.5), 4, 5,
             edgecolor = bordercolor, facecolor='none', lw=1))

# Study Domain ----------------------------------------------------------
ax1 = fig.add_subplot(1,3,2)

cs = ax1.pcolormesh(plon, plat, depth, vmin=0, vmax=800, cmap=newcmap)

# format figure
pfun.dar(ax1)

# Set axis limits
ax1.set_xlim([-126,-122])
ax1.set_ylim([45.5,50.5])
ax1.set_yticklabels([])
ax1.set_xticklabels([])
# add title
# ax1.set_title('(b) Study Domain',fontsize=20,loc='left', fontweight='bold')# color='#EEEEEE')
ax1.text(0.06, 0.95, '(b)', transform=ax1.transAxes, fontsize=16, fontweight='bold')
pfun.add_coast(ax1, color='silver')

# draw box around Puget Sound
bordercolor = 'deeppink'#'#EEEEEE'
ax1.add_patch(Rectangle((-123.3, 46.95), 1.2, 1.52,
             edgecolor = bordercolor, facecolor='none', lw=1))

# add distance bar
lat0 = 45.8
lon0 = -123.63425
lat1 = lat0
lon1 = -122.988875
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax1.plot([lon0,lon1],[lat0,lat1],color='black',linewidth=1, marker='|')
ax1.text((lon0+lon1)/2,lat0+0.1,'{} km'.format(x_dist_km),color='k',fontsize=12,
         horizontalalignment='center')

# Puget Sound ----------------------------------------------------------
ax2 = fig.add_subplot(1,3,3)

# get the grid data
ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
z = -ds.h.values
mask_rho = np.transpose(ds.mask_rho.values)
lon = ds.lon_rho.values
lat = ds.lat_rho.values
X = lon[0,:] # grid cell X values
Y = lat[:,0] # grid cell Y values
plon, plat = pfun.get_plon_plat(lon,lat)
# make a version of z with nans where masked
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
zm[np.transpose(mask_rho) != 0] = -1

# Create map
cs = ax2.pcolormesh(plon, plat, depth, vmin=0, vmax=250, cmap=newcmap)
cbar_ax = fig.add_axes([0.65, 0.14, 0.24, 0.02])
cbar = plt.colorbar(cs, cax=cbar_ax, orientation='horizontal')
cbar.ax.tick_params(labelsize=18)
cbar.set_label('Depth [m]', fontsize=18)
cbar.outline.set_visible(False)
pfun.add_coast(ax2, color='silver')

# format figure
pfun.dar(ax2)
# Set axis limits
ax2.set_xlim(-123.3, -122.1)
ax2.set_ylim(46.95, 48.47)
ax2.set_yticklabels([])
ax2.set_xticklabels([])
# add title
# ax2.set_title('(c) Puget Sound',fontsize=20,loc='left', fontweight='bold')# color='#EEEEEE')
ax2.text(0.06, 0.95, '(c)', transform=ax2.transAxes, fontsize=16, fontweight='bold')

# add distance bar
lat0 = 47.05
lon0 = -122.38208
lat1 = lat0
lon1 = -122.25
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax2.plot([lon0,lon1],[lat0,lat1],color='black',linewidth=1, marker='|')
ax2.text((lon0+lon1)/2,lat0+0.02,'{} km'.format(x_dist_km),color='k',fontsize=12,
         horizontalalignment='center')


# plt.tight_layout()
plt.subplots_adjust(wspace = 0, bottom=0.26)
# plt.savefig(out_dir / ('model_bathy.png'))#,transparent='True')
