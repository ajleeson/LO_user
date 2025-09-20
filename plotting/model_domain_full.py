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
Ldir = Lfun.Lstart()

# define grid indices to look at
j1 = 570
j2 = 1220
i1 = 250
i2 = 652

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

background = 'dimgray'

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../LO_data/grids/cas7/grid.nc')
z = ds.h.values
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
pfun.start_plot(fs=fs, figsize=(10,7))
fig = plt.figure()
plt.subplots_adjust(wspace=0, hspace=0)
# create colormap
# newcmap = cmocean.tools.crop_by_percent(cmocean.cm.deep_r, 20, which='max', N=None)
# newcmap = cmocean.cm.deep_r
newcmap = cmocean.tools.crop_by_percent(cmocean.cm.deep, 5, which='max')
# newcmap = cmocean.tools.crop_by_percent(cmocean.cm.rain_r, 10, which='max')
newcmap.set_bad(background,1.) # background color

# zm[np.transpose(mask_rho) != 0] = -1
# newcmap = plt.get_cmap('Blues_r')
# newcmap.set_bad(background,1.)

# Salish Sea ----------------------------------------------------------
ax0 = fig.add_subplot(1,2,1)
# cs = ax0.pcolormesh(plon, plat, zm, vmin=-5, vmax=0, cmap=newcmap)
cs = ax0.pcolormesh(plon, plat, zm, vmin=0, vmax=4000, cmap=newcmap)
# cs = ax0.pcolormesh(plon, plat, zm, vmin=-4000, vmax=0, cmap='gist_stern_r')
cbar = plt.colorbar(cs,ax=ax0, location='right', pad=0.05)
cbar.ax.tick_params(labelsize=11)#, color='#EEEEEE')#, rotation=30)
cbar.ax.set_ylabel('Depth [m]', fontsize=11)#, color='#EEEEEE')
cbar.outline.set_visible(False)
cbar_yticks = plt.getp(cbar.ax.axes, 'yticklabels')
# plt.setp(cbar_yticks, color='#EEEEEE')
# format figure
pfun.dar(ax0)
# pfun.add_coast(ax0, color='gray')
for border in ['top','right','bottom','left']:
        ax0.spines[border].set_visible(False)
# Set axis limits
# ax0.set_xlim(X[i1],-122)#X[i2]) # Salish Sea
# ax0.set_ylim(Y[j1],Y[j2]) # Salish Sea

ax0.set_ylabel('Latitude')
ax0.set_xlabel('Longitude')

print('Lon={},{}'.format(X[i1],-122))
print('Lat={},{}'.format(Y[j1],Y[j2]))

plt.xticks(rotation=30)#, color='#EEEEEE')
# plt.yticks(color='#EEEEEE')

# draw box around Puget Sound
bordercolor = 'violet'
ax0.add_patch(Rectangle((-123.2, 46.93), 1.1, 1.52,
             edgecolor = bordercolor, facecolor='none', lw=1.5))


plt.savefig(out_dir / ('model_domain_full.png'),transparent='True')
