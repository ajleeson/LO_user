"""
Plots a comparison of bottom DO for two identical grids with different oxygen value.

This is a custom function for a particular experiment, but code can be adapted for other use cases in the future.

From the terminal: python bottom_DO_nutrient_comparison.py

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

# ----------------------------------------------------------------------

# # Provide information about models to compare
# # basecase (has nutrients)
# bc_gtagex = 'cas6_traps2_x2b'
# # test condition (no nutrients)
# c1_gtagex = 'cas6_traps3_x2b'
# gtagexes = [bc_gtagex, c1_gtagex]

# hr = '0025'
# date = '2017.03.07'

#'Puget Sound','Whidbey Basin','North King County','Lynch Cove'
region = 'North King County'

# # Variables to compare
# vn_list = ['oxygen']#['NO3','phytoplankton','oxygen']

# Show WWTP locations?
WWTP_loc = True

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# get plotting limits based on region
if region == 'North King County':
    xmin = -122.55
    xmax = -122.35
    ymin = 47.55
    ymax = 47.75

# Get grid data
G = zrfun.get_basic_info('/home/aleeson/LO_data/grids/cas6/grid.nc', only_G=True)
grid_ds = xr.open_dataset('/home/aleeson/LO_data/grids/cas6/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values
px, py = pfun.get_plon_plat(G['lon_rho'],G['lat_rho'])
z = -grid_ds.h.values
mask_rho = np.transpose(grid_ds.mask_rho.values)
# get indices of min and max
imin_x = find_nearest(px[0,:],xmin)
imax_x = find_nearest(px[0,:],xmax)
imin_y = find_nearest(py[:,0],ymin)
imax_y = find_nearest(py[:,0],ymax)

# Initialize figure
fs = 10
pfun.start_plot(fs=fs, figsize=(8,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
# make a version of z with nans where masked
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
zm[np.transpose(mask_rho) != 0] = -1
ax.pcolormesh(px, py, zm, edgecolor='powderblue', linewidth=0.5, vmin=-1000, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
pfun.add_coast(ax)
pfun.dar(ax)
ax.set_title('Locations of WWTPs as Tiny Rivers')
ax.set_ylabel('Lat')
ax.set_xlabel('Lon')

# add locations as wwtps
wwtplabel=True
loc_df = pd.read_csv('../../LO_data/grids/cas6/wwtp_info.csv')
for rn in loc_df.index:
    ii = int(loc_df.loc[rn,'col_py'])
    jj = int(loc_df.loc[rn,'row_py'])
    if wwtplabel==True:
        ax.scatter(lon[jj,ii], lat[jj,ii],marker='o', s=60, color='#AEDC3C', edgecolor = 'k',label='Original WWTP Location')
        wwtplabel=False
    else:
        ax.scatter(lon[jj,ii], lat[jj,ii],marker='o', s=60, color='#AEDC3C', edgecolor = 'k')

# ax.scatter(lon_wwtps,lat_wwtps,color='lawngreen', edgecolors='green', alpha=0.7)

# add locations as tiny rivers
loc_df = pd.read_csv('../../LO_data/grids/cas6/wwtp_as_triv_info.csv')
rivlabel=True
for rn in loc_df.index:
# u or v grids.
    ii = int(loc_df.loc[rn,'col_py'])
    jj = int(loc_df.loc[rn,'row_py'])
    uv = loc_df.loc[rn,'uv']
    isign = loc_df.loc[rn,'isign']
    idir = loc_df.loc[rn,'idir']
    if uv == 'u' and isign == 1:
        if rivlabel==True:
            ax.scatter(lon[jj,ii+1], lat[jj,ii+1],marker='*', s=60, color='#7148BC', label='Moved Location')
            rivlabel=False
        else:
            ax.scatter(lon[jj,ii+1], lat[jj,ii+1],marker='*', s=60, color='#7148BC')
    if uv == 'u' and isign == -1:
        ax.scatter(lon[jj,ii], lat[jj,ii],marker='*', s=60, color='#7148BC')

    if uv == 'v' and isign == 1:
        ax.scatter(lon[jj+1,ii], lat[jj+1,ii],marker='*', s=60, color='#7148BC')
                
    if uv == 'v' and isign == -1:
        ax.scatter(lon[jj,ii], lat[jj,ii],marker='*', s=60, color='#7148BC')

plt.legend(loc='upper right')
                        

# Generate plot
plt.show()
