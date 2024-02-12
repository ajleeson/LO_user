"""
Code to draw boz around Puget Sound in map
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import xarray as xr
import cmocean
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun

# define grid indices to look at
j1 = 570
j2 = 1170
i1 = 220
i2 = 652

# lon/lat limits (Puget Sound)
xmin = -123.2
xmax = -122.1
ymin = 46.93
ymax = 48.45

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
zm[np.transpose(mask_rho) != 0] = -1

# Create map
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1,1,1)
ax.add_patch(Rectangle((X[i1], Y[j1]), -121.4-X[i1],Y[j2]-Y[j1], facecolor='#EEEEEE'))
plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-6, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
# Add mooring location
plt.scatter(-122.813348, 47.667696, s=150, c='deeppink', marker='*')#, edgecolors='k')
# format
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
pfun.dar(ax)

# # draw box around Puget Sound
# ax.add_patch(Rectangle((-123.2, 46.93), 1.1, 1.52,
#              edgecolor = 'white', facecolor = 'none', lw=1.5))

# remove border
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig('hood_canal_mooring.png',transparent=True)
