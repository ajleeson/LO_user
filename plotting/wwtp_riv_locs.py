"""
Code to draw plot regional bathymetry
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import xarray as xr
import cmocean
import pandas as pd
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

background = 'silver'

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

# Get wwtp and riv locations -------------------------------------------------------

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
zm[np.transpose(mask_rho) != 0] = -1

# get river and wwtp locations
df_LOriv = pd.read_csv('../../LO_data/grids/cas7/river_info.csv')
df_triv = pd.read_csv('../../LO_data/grids/cas7/triv_info.csv')
df_wwtp = pd.read_csv('../../LO_data/grids/cas7/wwtp_info.csv')
df_riv = pd.concat([df_LOriv,df_triv])

# get lat/lon coordinates for rivers and wwtps
lon_riv = [X[int(col)] for col in df_riv['col_py']]
lat_riv = [Y[int(row)] for row in df_riv['row_py']]
lon_wwtp = [X[int(col)] for col in df_wwtp['col_py']]
lat_wwtp = [Y[int(row)] for row in df_wwtp['row_py']]


# Create bathymetry plot --------------------------------------------------------------

# Initialize figure
plt.close('all')
fs = 10
pfun.start_plot(fs=fs, figsize=(10,7))
fig = plt.figure()
plt.subplots_adjust(wspace=0, hspace=0)
# create colormap
newcmap = cmocean.tools.crop_by_percent(cmocean.cm.deep_r, 10, which='max')
newcmap.set_bad(background,1.) # background color

zm[np.transpose(mask_rho) != 0] = -1
newcmap = plt.get_cmap('Blues_r')
newcmap.set_bad(background,1.)

# Salish Sea ----------------------------------------------------------
ax0 = fig.add_subplot(1,2,1)
cs = ax0.pcolormesh(plon, plat, zm, vmin=-10, vmax=0, cmap=newcmap)

# format figure
pfun.dar(ax0)
# pfun.add_coast(ax0, color='gray')
for border in ['top','right','bottom','left']:
        ax0.spines[border].set_visible(False)
# Set axis limits
ax0.set_xlim(X[i1],-122)#X[i2]) # Salish Sea
ax0.set_ylim(Y[j1],Y[j2]) # Salish Sea
ax0.set_xlabel('Longitude',color='gray', fontsize=12)
ax0.set_ylabel('Latitude', color='gray', fontsize=12)

print('Lon={},{}'.format(X[i1],-122))
print('Lat={},{}'.format(Y[j1],Y[j2]))

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
# draw box around Puget Sound
bordercolor = 'black'
ax0.add_patch(Rectangle((-123.2, 46.93), 1.1, 1.52,
             edgecolor = bordercolor, facecolor='none', lw=1))

# add river and wwtp locations
size = 20
ax0.scatter(lon_riv,lat_riv,color='royalblue',marker='d',alpha=0.9,s=size,label='Rivers')
ax0.scatter(lon_wwtp,lat_wwtp,color='deeppink',alpha=0.9,s=size,label='WWTPs')
ax0.legend(loc='upper right',fontsize=14, borderpad=0.3, handletextpad=0.02)


# Puget Sound ----------------------------------------------------------
ax1 = fig.add_subplot(1,2,2)
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
zm[np.transpose(mask_rho) != 0] = -1
newcmap = plt.get_cmap('Blues_r')
newcmap.set_bad(background,1.)
cs = ax1.pcolormesh(plon, plat, zm, vmin=-10, vmax=0, cmap=newcmap)

# add river and wwtp locations
size = 20
ax1.scatter(lon_riv,lat_riv,color='royalblue',marker='d',alpha=0.7,s=size,label='Rivers')
ax1.scatter(lon_wwtp,lat_wwtp,color='deeppink',alpha=0.7,s=size,label='WWTPs')

# format figure
pfun.dar(ax1)
# pfun.add_coast(ax1, color='gray')
for border in ['top','right','bottom','left']:
        # ax1.spines[border].set_visible(False)
        ax1.spines[border].set_color(bordercolor)
        ax1.spines[border].set_linewidth(2)
# Set axis limits
ax1.set_xlim([-123.3,-122.1])
ax1.set_ylim([46.93,48.45])
ax1.set_yticklabels([])
ax1.set_xticklabels([])
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

plt.subplots_adjust(hspace = 0.01)
plt.savefig(out_dir / ('wwtp_riv_locs.png'))
