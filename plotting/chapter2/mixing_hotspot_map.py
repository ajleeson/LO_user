"""
Plot vertical flux of NH4
"""

###################################################################
##                       import packages                         ##  
###################################################################      

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import matplotlib.colors
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
NH4_loading = ds_loading['NH4'].values[0,:,:,:] # [mmol/m3]
NH4_noloading = ds_noloading['NH4'].values[0,:,:,:] # [mmol/m3]
# Get z_rho and z_w
# get S for the whole grid
Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
reader = csv.DictReader(open(Sfp))
S_dict = {}
for row in reader:
    S_dict[row['ITEMS']] = row['VALUES']
S = zrfun.get_S(S_dict)
h = ds_loading['h'].values # height of water column
zeta = ds_loading['zeta'].values
z_rho, z_w = zrfun.get_z(h, zeta[0, :, :], S) # [m from surface, negative down]

# get eddy diffusivity (note: this is for salt!!! DO WE NEED A DIFFERENT VALUE FOR NUTRIENTS????)
Aks_loading = ds_loading['AKs'].values[0,:,:,:] # [m2/s]
Aks_noloading = ds_loading['AKs'].values[0,:,:,:] # [m2/s]
# land mask
mask_rho = ds_loading['mask_rho'].values
# apply land mask
Aks_loading = np.where(mask_rho==1, Aks_loading, np.nan) # [m2/s]
Aks_noloading = np.where(mask_rho==1, Aks_loading, np.nan) # [m2/s]

# Initialize d/dz(NH4) array with same shape as eddy diffusivity (31, 1302, 663)
dNH4_dzw_loading = np.zeros_like(Aks_loading) # pad with zeros
dNH4_dzw_noloading = np.zeros_like(Aks_noloading) # pad with zeros

# Get d/dz(NH4) at the vertical midpoints (z_w) using central differences
dNH4_dzw_loading[1:-1, :, :] = np.diff(NH4_loading, axis=0) / np.diff(z_rho, axis=0)
dNH4_dzw_noloading[1:-1, :, :] = np.diff(NH4_noloading, axis=0) / np.diff(z_rho, axis=0)
# Note that at surface and boundary, d/dz(NH4) is zero

# Get flux of NH4 due to turbulent diffusion F = -Aks * d/dz(NH4)
F_loading = -Aks_loading * dNH4_dzw_loading # [mmol m^-2 s^-1]
F_noloading = -Aks_noloading * dNH4_dzw_noloading # [mmol m^-2 s^-1]

# Get turbulent diffusion neares to bottom of photic zone
photic_zone_depth = -10
# Get index of z_w that is closest to photic zone depth at each horizontal location
k_idx = np.argmin(np.abs(z_w - photic_zone_depth), axis=0) 
# Convert to 3D so we can apply it to eddy diffusivity, which is 3D
k_idx_3d = np.expand_dims(k_idx, axis=0)
# Get F value at the photic zone depth
F_photic_loading = np.take_along_axis(F_loading, k_idx_3d, axis=0).squeeze(axis=0)
F_photic_noloading = np.take_along_axis(F_noloading, k_idx_3d, axis=0).squeeze(axis=0)
# Mask out areas where the physical bottom is shallower than the photic zone depth
F_photic_loading = np.where(ds_loading['h'].values >= np.abs(photic_zone_depth), F_photic_loading, np.nan)
F_photic_noloading = np.where(ds_noloading['h'].values >= np.abs(photic_zone_depth), F_photic_noloading, np.nan)

# ###################################################################
# ##                  Plotting and saving figure                   ##  
# ################################################################### 

# Initialize figure
fig,ax = plt.subplots(1,2, figsize=(12,9))

# plot no-loading nutrient flux
cs = ax[0].pcolormesh(px,py,F_photic_noloading,cmap=cmc.bam,
                   norm=matplotlib.colors.SymLogNorm(linthresh=1e-7, linscale=1,
                                              vmin=-1e-4, vmax=1e-4, base=10))
# cs = ax[0].pcolormesh(px,py,F_photic_noloading,cmap=cmc.bam,vmin=-1e-5, vmax=1e-5)
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
ax[0].set_title('No-loading\n' + r'NH4 vertical flux [mmol m$^{-2}$ s$^{-1}$]' + '\n' + date + ' ocean_his_0001',
            fontsize=12, fontweight='bold')



# plot loading minus no-loading nutrient flux
cs = ax[1].pcolormesh(px,py,F_photic_loading - F_photic_noloading,cmap=cmc.vik,
                   norm=matplotlib.colors.SymLogNorm(linthresh=1e-7, linscale=1,
                                              vmin=-1e-5, vmax=1e-5, base=10))
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
ax[1].set_title('Loading - No-loading\n' + r'NH4 vertical flux [mmol m$^{-2}$ s$^{-1}$]' + '\n' + date + ' ocean_his_0001',
            fontsize=12, fontweight='bold')


# Generate plot
plt.tight_layout
plt.show()