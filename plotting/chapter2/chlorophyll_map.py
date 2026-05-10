"""
Plot surface chlorophyll
"""

###################################################################
##                       import packages                         ##  
###################################################################      

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import matplotlib.colors
import cmocean
import numpy as np
import xarray as xr
import cmcrameri.cm as cmc
import csv
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()
Ldir = Lfun.Lstart()

plt.close('all')

###################################################################
##                          User Inputs                          ##  
################################################################### 

date = '2017.02.12'


###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

gtagex_loading = 'cas7_t1_x11b'#'cas7_t1_x11ab'
gtagex_noloading = 'cas7_t1noDIN_x11b'#'cas7_t1_x11ab'

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values

# lon/lat limits (Puget Sound)
xmin = -123.2
xmax = -122.1
ymin = 46.93
ymax = 48.45

fp = Ldir['roms_out'] / gtagex_loading / ('f'+date) / 'ocean_his_0001.nc'
ds_loading = xr.open_dataset(fp)

fp = Ldir['roms_out'] / gtagex_noloading / ('f'+date) / 'ocean_his_0001.nc'
ds_noloading = xr.open_dataset(fp)

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values

px, py = pfun.get_plon_plat(lon,lat)

###################################################################
## Determine flux of nutrients to surface due to vertical mixing ##  
################################################################### 

# get NH4
phytoplankton_loading = ds_loading['phytoplankton'].values[0,-1,:,:] # [mmol/m3]
phytoplankton_noloading = ds_noloading['phytoplankton'].values[0,-1,:,:] # [mmol/m3]


# ###################################################################
# ##                  Plotting and saving figure                   ##  
# ################################################################### 

# Initialize figure
fig,ax = plt.subplots(1,2, figsize=(12,9))

# plot no-loading nutrient flux
cs = ax[0].pcolormesh(px,py,phytoplankton_noloading,cmap=cmocean.cm.algae, vmin=-0, vmax=3)
cbar = fig.colorbar(cs)
cbar.ax.tick_params(labelsize=12)
cbar.outline.set_visible(False)
# format figure
ax[0].set_xlim([xmin,xmax])
ax[0].set_ylim([ymin,ymax])
ax[0].set_yticklabels([])
ax[0].set_xticklabels([])
ax[0].axis('off')
pfun.add_coast(ax[0], color='silver')
pfun.dar(ax[0])
ax[0].set_title('No-loading\n' + r'surface phytoplankton [mmol m$^{-3}$]' + '\n' + date + ' ocean_his_0001',
            fontsize=12, fontweight='bold')



# plot loading minus no-loading nutrient flux
cs = ax[1].pcolormesh(px,py,phytoplankton_noloading - phytoplankton_loading,cmap=cmc.vik, vmin=-0.5, vmax=0.5)
# cs = ax[1].pcolormesh(px,py,F_photic_loading - F_photic_noloading,cmap=cmc.vik,vmin=-1e-5, vmax=1e-5)
cbar = fig.colorbar(cs)
cbar.ax.tick_params(labelsize=12)
cbar.outline.set_visible(False)
# format figure
ax[1].set_xlim([xmin,xmax])
ax[1].set_ylim([ymin,ymax])
ax[1].set_yticklabels([])
ax[1].set_xticklabels([])
ax[1].axis('off')
pfun.add_coast(ax[1], color='silver')
pfun.dar(ax[1])
ax[1].set_title('Loading - No-loading\n' + r'surface phytoplankton [mmol m$^{-3}$]' + '\n' + date + ' ocean_his_0001',
            fontsize=12, fontweight='bold')


# Generate plot
plt.tight_layout
plt.show()