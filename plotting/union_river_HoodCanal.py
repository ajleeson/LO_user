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
Ldir = Lfun.Lstart()

# define grid indices to look at
j1 = 590
j2 = 950
i1 = 475
i2 = 650

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc')
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
fig = plt.figure(figsize=(5,7))
ax = fig.add_subplot(1,1,1)
ax.add_patch(Rectangle((X[i1], Y[j1]), -X[i2]-X[i1],Y[j2]-Y[j1], facecolor='#EEEEEE'))
plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-6, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
# format
plt.xticks(rotation=30, ha='right')
ax.set_xlim(X[i1],X[i2])#X[i2]) # Salish Sea
ax.set_ylim(Y[j1],Y[j2]) # Salish Sea
pfun.dar(ax)

# remove border
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig(out_dir / ('Union_river.png'))
