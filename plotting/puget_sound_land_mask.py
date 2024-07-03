"""
Just Puget Sound land mask and water
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import xarray as xr
import cmocean
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun


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
fig = plt.figure(figsize=(8,15))
ax = fig.add_subplot(1,1,1)
ax.add_patch(Rectangle((X[0], Y[0]), X[-1]-X[0],Y[-1]-Y[0], facecolor='dimgray'))
plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-6, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
# format
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
ax.set_xlim(-123.29, -122.1) # Puget Sound DO
ax.set_ylim(46.95, 48.93) # Puget Sound DO
pfun.dar(ax)

# remove border
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig('puget_sound_landmask.png',transparent=True)
