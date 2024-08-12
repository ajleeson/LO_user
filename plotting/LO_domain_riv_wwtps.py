"""
Code to draw boz around Puget Sound in map
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from matplotlib.patches import Rectangle
import cmocean
import pandas as pd
from matplotlib.ticker import MaxNLocator
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun


# Get LiveOcean grid info --------------------------------------------------
plt.close('all')

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

# Create map
fig = plt.figure(figsize=(6,8),facecolor='none')
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-6, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
# add river and wwtp locations
size = 20
ax.scatter(lon_riv,lat_riv,color='royalblue',marker='d',alpha=0.7,s=size,label='Rivers')
ax.scatter(lon_wwtp,lat_wwtp,color='tomato',alpha=0.7,s=size,label='WWTPs')
# format
ax.set_ylabel('Latitude',color='#EEEEEE', fontsize=14)
ax.set_xlabel('Longitude',color='#EEEEEE', fontsize=14)
ax.set_facecolor('#EEEEEE')
ax.legend(loc='lower left',fontsize=14)
ax.set_ylim([46.5,50])
ax.set_xlim([-125,-122])
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.tick_params(axis='x')#, labelcolor='#EEEEEE', labelsize=14,color='#EEEEEE')
ax.tick_params(axis='y')#, labelcolor='#EEEEEE', labelsize=14,color='#EEEEEE')
pfun.dar(ax)

# draw box around Puget Sound
ps_color = 'darkmagenta'
ax.add_patch(Rectangle((-123.2, 46.93), 1.1, 1.52,
             edgecolor = ps_color, facecolor = 'none', lw=1.5))
ax.text(0.8, 0.1, 'Puget Sound', color=ps_color,
            verticalalignment='center', horizontalalignment='center',
            transform=ax.transAxes, fontsize=14)


# remove border
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig('riv_wwtp_locs.png')#,transparent=True)
