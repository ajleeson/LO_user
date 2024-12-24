"""
Code to draw plot regional bathymetry
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import xarray as xr
from matplotlib.ticker import MaxNLocator
import cmocean
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

# define grid indices to look at
j1 = 570
j2 = 1220
i1 = 250
i2 = 652

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

background = 'white'#'whitesmoke' #'dimgray'

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
pfun.start_plot(fs=fs, figsize=(11,7))
fig = plt.figure()
plt.subplots_adjust(wspace=0, hspace=0)
# create colormap
# newcmap = cmocean.tools.crop_by_percent(cmocean.cm.deep_r, 20, which='max', N=None)
# newcmap = cmocean.cm.deep_r
# newcmap = cmocean.tools.crop_by_percent(cmocean.cm.deep_r, 15, which='max')
newcmap = cmocean.tools.crop_by_percent(cmocean.cm.thermal, 20, which='max')
newcmap.set_bad(background,1.) # background color

# zm[np.transpose(mask_rho) != 0] = -1
# newcmap = plt.get_cmap('Blues_r')
# newcmap.set_bad(background,1.)

# Salish Sea ----------------------------------------------------------
ax0 = fig.add_subplot(1,2,1)
# cs = ax0.pcolormesh(plon, plat, zm, vmin=-5, vmax=0, cmap=newcmap)
cs = ax0.pcolormesh(plon, plat, zm, vmin=-250, vmax=0, cmap=newcmap)
# cbar = plt.colorbar(cs,ax=ax0, location='left', pad=0.1)#pad=0.05)
# cbar.ax.tick_params(labelsize=11)#, color='#EEEEEE')#, rotation=30)
# cbar.ax.set_ylabel('Depth [m]', fontsize=11)#, color='#EEEEEE')
# cbar.outline.set_visible(False)
# cbar_yticks = plt.getp(cbar.ax.axes, 'yticklabels')
# plt.setp(cbar_yticks)#, color='#EEEEEE')
# format figure
pfun.dar(ax0)
# pfun.add_coast(ax0, color='gray')
# for border in ['top','right','bottom','left']:
#         ax0.spines[border].set_visible(False)
# Set axis limits
ax0.set_xlim(X[i1],-122)#X[i2]) # Salish Sea
ax0.set_ylim(Y[j1],Y[j2]) # Salish Sea

print('Lon={},{}'.format(X[i1],-122))
print('Lat={},{}'.format(Y[j1],Y[j2]))

# plt.xticks(rotation=30, color='gray')
plt.xticks(rotation=30,horizontalalignment='right',fontsize=12)
plt.yticks(fontsize=12)
ax0.set_ylabel('Latitude', fontsize=12)
ax0.set_xlabel('Longitude', fontsize=12)
ax0.xaxis.set_major_locator(MaxNLocator(integer=True))
ax0.yaxis.set_major_locator(MaxNLocator(integer=True))
# plt.yticks(color='gray')
# ax0.set_yticklabels([])
# ax0.set_xticklabels([])
# add title
ax0.set_title('(a) Salish Sea',fontsize=14, loc='left', fontweight='bold')#,color='#EEEEEE')
# add 10 km bar
lat0 = 47
lon0 = -124.63175
lat1 = lat0
lon1 = -124.36975
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax0.plot([lon0,lon1],[lat0,lat1],color='white',linewidth=5)
ax0.text((lon0+lon1)/2,lat0+0.05,'{} km'.format(x_dist_km),color='w',fontsize=12,
         horizontalalignment='center')
# draw box around Puget Sound
bordercolor = 'k'#'#EEEEEE'
ax0.add_patch(Rectangle((-123.3, 46.93), 1.2, 1.52,
             edgecolor = bordercolor, facecolor='none', lw=1))

# add major cities
# Seattle
ax0.scatter([-122.3328],[47.6061],s=[250],color='pink',
            marker='*',edgecolors='darkred')
ax0.text(-122.3328 + 0.1,47.6061,'Seattle',color='darkred', rotation=90,
         horizontalalignment='left',verticalalignment='center', size=12)
# Vancouver
ax0.scatter([-123.1207],[49.2827],s=[250],color='pink',
            marker='*',edgecolors='darkred')
ax0.text(-123.1207 + 0.1,49.2827,'Vancouver',color='darkred', rotation=0,
         horizontalalignment='left',verticalalignment='center', size=12)

# add major water bodies
ax0.text(-124.937095,47.782238,'Pacific Ocean',color='k', rotation=-75,
         horizontalalignment='left',verticalalignment='center', size=12,fontweight='bold')

# Puget Sound ----------------------------------------------------------
ax1 = fig.add_subplot(1,2,2)
# zm = z.copy()
# zm[np.transpose(mask_rho) == 0] = np.nan
# zm[np.transpose(mask_rho) != 0] = -1
# newcmap = plt.get_cmap('Blues_r')
# newcmap.set_bad(background,1.)
# cs = ax1.pcolormesh(plon, plat, zm, vmin=-5, vmax=0, cmap=newcmap)
cs = ax1.pcolormesh(plon, plat, zm, vmin=-250, vmax=0, cmap=newcmap)
cbar = plt.colorbar(cs,ax=ax1, location='right', pad=0.05)
cbar.ax.tick_params(labelsize=11)#, color='#EEEEEE')#, rotation=30)
cbar.ax.set_ylabel('Depth [m]', fontsize=11)#, color='#EEEEEE')
cbar.outline.set_visible(False)
cbar_yticks = plt.getp(cbar.ax.axes, 'yticklabels')
plt.setp(cbar_yticks)#, color='#EEEEEE')
# format figure
pfun.dar(ax1)
# pfun.add_coast(ax1, color='gray')
# for border in ['top','right','bottom','left']:
#         ax1.spines[border].set_visible(False)
#         # ax1.spines[border].set_color(bordercolor)
#         # ax1.spines[border].set_linewidth(2)
# Set axis limits
ax1.set_xlim([-123.3,-122.1])
ax1.set_ylim([46.93,48.45])
ax1.set_yticklabels([])
ax1.set_xticklabels([])
# plt.xticks(rotation=30,horizontalalignment='right',color='gray')
# ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
# ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
# plt.xticks(rotation=30, color='gray')
# plt.yticks(color='gray')
# add title
ax1.set_title('(b) Puget Sound',fontsize=14,loc='left', fontweight='bold')# color='#EEEEEE')
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

plt.subplots_adjust(wspace = -0.3)
plt.savefig(out_dir / ('model_bathy.png'))#,transparent='True')
