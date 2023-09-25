"""
Code to draw plot regional bathymetry
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import xarray as xr
import cmocean
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

# define grid indices to look at
j1 = 570
j2 = 1200
i1 = 250
i2 = 652

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../LO_data/grids/cas7/grid.nc')
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

# Create bathymetry plot --------------------------------------------------------------

# Initialize figure
plt.close('all')
fs = 10
pfun.start_plot(fs=fs, figsize=(12,7))
fig = plt.figure()
plt.subplots_adjust(wspace=0, hspace=0)
# create colormap
# newcmap = cmocean.tools.crop_by_percent(cmocean.cm.deep_r, 20, which='max', N=None)
newcmap = cmocean.cm.deep_r
newcmap.set_bad('silver',1.) # background color

# Salish Sea ----------------------------------------------------------
ax0 = fig.add_subplot(1,2,1)
cs = ax0.pcolormesh(plon, plat, zm, vmin=-600, vmax=0, cmap=newcmap)
cbar = plt.colorbar(cs,ax=ax0, location='left')
cbar.ax.tick_params(labelsize=14, rotation=30)
cbar.ax.set_ylabel('Depth [m]', fontsize=14)
cbar.outline.set_visible(False)
# format figure
pfun.dar(ax0)
pfun.add_coast(ax0, color='gray')
for border in ['top','right','bottom','left']:
        ax0.spines[border].set_visible(False)
# Set axis limits
ax0.set_xlim(X[i1],-122)#X[i2]) # Salish Sea
ax0.set_ylim(Y[j1],Y[j2]) # Salish Sea
plt.xticks(rotation=30, color='gray')
plt.yticks(color='gray')
# add title
ax0.set_title('Salish Sea',fontsize=16)
# add 10 km bar
lat0 = 47
lon0 = -124.63175
lat1 = lat0
lon1 = -124.36975
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax0.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=5)
ax0.text((lon0+lon1)/2,lat0+0.05,'{} km'.format(x_dist_km),color='k',fontsize=12,
         horizontalalignment='center')
# # draw box around Puget Sound
# ax0.add_patch(Rectangle((-123.2, 46.93), 1.1, 1.52,
#              edgecolor = 'blueviolet', facecolor = 'none', lw=2))

# Puget Sound ----------------------------------------------------------
ax1 = fig.add_subplot(1,2,2)
cs = ax1.pcolormesh(plon, plat, zm, vmin=-250, vmax=0, cmap=newcmap)
# # add isobath
# cl = ax1.contour(lon, lat, zm, [-100], colors='k',
#                  linestyles='solid', linewidths=1)
# ax0.clabel(cl, inline=True, fontsize=10)
# add colorbar
cbar = plt.colorbar(cs,ax=ax1, location='right')
cbar.ax.tick_params(labelsize=14, rotation=30)
cbar.ax.set_ylabel('Depth [m]', fontsize=14)
cbar.outline.set_visible(False)
# format figure
pfun.dar(ax1)
pfun.add_coast(ax1, color='gray')
for border in ['top','right','bottom','left']:
        ax1.spines[border].set_visible(False)
# Set axis limits
ax1.set_xlim([-123.3,-122.1])
ax1.set_ylim([46.93,48.45])
plt.xticks(rotation=30, color='gray')
plt.yticks(color='gray')
# add title
ax1.set_title('Puget Sound',fontsize=16)
# add 10 km bar
lat0 = 47
lon0 = -123.2
lat1 = lat0
lon1 = -123.06825
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax1.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=5)
ax1.text((lon0+lon1)/2,lat0+0.015,'{} km'.format(x_dist_km),color='k',fontsize=12,
         horizontalalignment='center')

plt.savefig('model_bathy.png')\
