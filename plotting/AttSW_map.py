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

lat_LB1 = 47.02
lat_UB1 = 50.29
lon_LB1 = -123.89

lon_LB2 = -125.31
lon_UB2 = lon_LB1
lat_LB2 = 49.13
lat_UB2 = 51.02

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc') # Change to cas7 later!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
fig = plt.figure(figsize=(8,10))
ax = fig.add_subplot(1,1,1)
ax.add_patch(Rectangle((X[0], Y[0]), X[-1]-X[0],Y[-1]-Y[0], facecolor='dimgray'))
plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-6, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
# format
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
ax.set_ylim(45.5,Y[-1])
pfun.dar(ax)

# draw box around Salish Sea
ax.add_patch(Rectangle((lon_LB1, lat_LB1), X[-1]-lon_LB1, lat_UB1-lat_LB1,
             edgecolor = 'white', facecolor = 'pink', alpha=0.5, lw=1.5))
ax.add_patch(Rectangle((lon_LB2, lat_LB2), lon_UB2-lon_LB2, lat_UB2-lat_LB2,
             edgecolor = 'white', facecolor = 'pink', alpha=0.5, lw=1.5))

# remove border
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

# label AttSW values
ax.text(0.10,0.35,'AttSW =\n' + r'0.05 m$^{-1}$',transform=ax.transAxes,size=16,color='teal')
ax.text(0.80,0.75,'AttSW =\n' + r'0.15 m$^{-1}$',transform=ax.transAxes,size=16,color='pink')

# # add test points
# ax.scatter([X[50],X[595],X[300]], [Y[756],Y[756],Y[1130]],s=75,color='white',edgecolor='k')
# ax.text(X[45],Y[680],'Coast\ntest point')
# ax.text(X[420],Y[680],'Main Basin\ntest point')
# ax.text(X[250],Y[1150],'North SoG\ntest point')

plt.savefig(out_dir / ('AttSW_map.png'))

# print('longitude: {}'.format(X[300]))
# print('latitude: {}'.format(Y[1130]))
