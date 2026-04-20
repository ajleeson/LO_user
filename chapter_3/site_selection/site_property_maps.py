"""
Plot model domain properties for OAE site selection
"""

###################################################################
##                       import packages                         ##  
###################################################################      

import copy
import numpy as np
import xarray as xr
import cmocean
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import cmcrameri.cm as cmc

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

date = '2020.05.22'
date_formatted = '2020-05-22'

###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

gtagex = 'cas7_t1d_x11ad'

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
xmin = -126
xmax = -122
ymin = 45.5
ymax = 50.5

# get model output
fp = Ldir['roms_out'] / gtagex / ('f'+date) / 'ocean_avg_0001.nc'
ds_dye = xr.open_dataset(fp)

# get surface mixed layer depth
fp = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'sml_plus' / 'threshold_p125' / 'LO_domain_2020.01.01_2020.12.31' / 'LO_domain_sml_plus_2020.01.01_2020.12.31.nc'
ds_sml = xr.open_dataset(fp)
print(ds_sml['zSML'])
ds_sml_may22 = ds_sml.sel(ocean_time=date_formatted)

###################################################################
##                       Plot surface dye                        ##  
################################################################### 

px, py = pfun.get_plon_plat(lon,lat)
vmin = 0
vmax = 0.5

# surface dye concentration
surf_dye_01 = ds_dye['dye_01'][0,-1,:,:].values

# Initialize figure
fig,ax = plt.subplots(1,1, figsize=(7,9))

# plot values
cmap = copy.copy(plt.get_cmap('plasma'))
cmap.set_bad(color='white', alpha=1.0)
cmap = cmc.acton
cs = ax.pcolormesh(px,py,surf_dye_01,vmin=vmin, vmax=vmax, cmap=cmap)

# Add Puget Sound Inset
# [x0, y0, width, height]
axins = ax.inset_axes([0.71, 0.012, 0.45, 0.6])
# plot values in inset
axins.pcolormesh(px, py, surf_dye_01, vmin=vmin, vmax=vmax, cmap=cmap)
# Puget Sound limits
axins.set_xlim(-123.2, -122.1)
axins.set_ylim(46.95, 48.4)
axins.tick_params(left=False, bottom=False)
# format
axins.set_xticklabels([])
axins.set_yticklabels([])
for spine in axins.spines.values():
    spine.set_edgecolor('grey')
    spine.set_linewidth(2)

# add colorbar
cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
cbar.ax.tick_params(labelsize=18, rotation=30)
cbar.outline.set_visible(False)

# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.tick_params(left=False, bottom=False)
# pfun.add_coast(ax, color='k')
pfun.dar(ax)
pfun.dar(axins)
ax.set_title('Surface dye\n'+ r'concentration [kg m$^3$]', fontsize=20,
             loc='Left', fontweight='bold')
# ax.text()

# Generate plot
plt.tight_layout
plt.subplots_adjust(bottom=0.001, top=0.9)
plt.show()

###################################################################
##                Plot surface mixed layer depth                 ##  
################################################################### 

px, py = pfun.get_plon_plat(lon,lat)
vmin = 0
vmax = 35

# surface dye concentration
sml = ds_sml_may22['SML_thickness'][0,:,:].values

# Initialize figure
fig,ax = plt.subplots(1,1, figsize=(7,9))

# plot values
# cmap = copy.copy(plt.get_cmap('viridis_r'))
# cmap.set_bad(color='white', alpha=1.0)
cmap = cmc.batlow
cs = ax.pcolormesh(px,py,sml,vmin=vmin, vmax=vmax, cmap=cmap)

# Add Puget Sound Inset
# [x0, y0, width, height]
axins = ax.inset_axes([0.71, 0.012, 0.45, 0.6])
# plot values in inset
axins.pcolormesh(px, py, sml, vmin=vmin, vmax=vmax, cmap=cmap)
# Puget Sound limits
axins.set_xlim(-123.2, -122.1)
axins.set_ylim(46.95, 48.4)
axins.tick_params(left=False, bottom=False)
# format
axins.set_xticklabels([])
axins.set_yticklabels([])
for spine in axins.spines.values():
    spine.set_edgecolor('grey')
    spine.set_linewidth(2)

# add colorbar
cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
cbar.ax.tick_params(labelsize=18, rotation=30)
cbar.outline.set_visible(False)

# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.tick_params(left=False, bottom=False)
# pfun.add_coast(ax, color='k')
pfun.dar(ax)
pfun.dar(axins)
ax.set_title('Surface mixed layer\n'+ r'depth [m]', fontsize=20,
             loc='Left', fontweight='bold')
# ax.text()

# Generate plot
plt.tight_layout
plt.subplots_adjust(bottom=0.001, top=0.9)
plt.show()

###################################################################
##               Plot SALINITY in surface mixed layer            ##  
################################################################### 

px, py = pfun.get_plon_plat(lon,lat)
vmin = 25
vmax = 33

# surface dye concentration
SA = ds_sml_may22['SA'][0,:,:].values

# Initialize figure
fig,ax = plt.subplots(1,1, figsize=(7,9))

# plot values
# cmap = copy.copy(plt.get_cmap('viridis_r'))
# cmap.set_bad(color='white', alpha=1.0)
cmap = cmc.imola
cs = ax.pcolormesh(px,py,SA,vmin=vmin, vmax=vmax, cmap=cmap)

# Add Puget Sound Inset
# [x0, y0, width, height]
axins = ax.inset_axes([0.71, 0.012, 0.45, 0.6])
# plot values in inset
axins.pcolormesh(px, py, SA, vmin=vmin, vmax=vmax, cmap=cmap)
# Puget Sound limits
axins.set_xlim(-123.2, -122.1)
axins.set_ylim(46.95, 48.4)
axins.tick_params(left=False, bottom=False)
# format
axins.set_xticklabels([])
axins.set_yticklabels([])
for spine in axins.spines.values():
    spine.set_edgecolor('grey')
    spine.set_linewidth(2)

# add colorbar
cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
cbar.ax.tick_params(labelsize=18, rotation=30)
cbar.outline.set_visible(False)

# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.tick_params(left=False, bottom=False)
# pfun.add_coast(ax, color='k')
pfun.dar(ax)
pfun.dar(axins)
ax.set_title('Surface mixed layer\n'+ r'abs. salinity [g kg$^{-1}$]', fontsize=20,
             loc='Left', fontweight='bold')
# ax.text()

# Generate plot
plt.tight_layout
plt.subplots_adjust(bottom=0.001, top=0.9)
plt.show()


###################################################################
##             Plot TEMPERATURE in surface mixed layer           ##  
################################################################### 

px, py = pfun.get_plon_plat(lon,lat)
vmin = 8
vmax = 16

# surface dye concentration
CT = ds_sml_may22['CT'][0,:,:].values

# Initialize figure
fig,ax = plt.subplots(1,1, figsize=(7,9))

# plot values
# cmap = copy.copy(plt.get_cmap('viridis_r'))
# cmap.set_bad(color='white', alpha=1.0)
cmap = cmc.lajolla
cs = ax.pcolormesh(px,py,CT,vmin=vmin, vmax=vmax, cmap=cmap)

# Add Puget Sound Inset
# [x0, y0, width, height]
axins = ax.inset_axes([0.71, 0.012, 0.45, 0.6])
# plot values in inset
axins.pcolormesh(px, py, CT, vmin=vmin, vmax=vmax, cmap=cmap)
# Puget Sound limits
axins.set_xlim(-123.2, -122.1)
axins.set_ylim(46.95, 48.4)
axins.tick_params(left=False, bottom=False)
# format
axins.set_xticklabels([])
axins.set_yticklabels([])
for spine in axins.spines.values():
    spine.set_edgecolor('grey')
    spine.set_linewidth(2)

# add colorbar
cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
cbar.ax.tick_params(labelsize=18, rotation=30)
cbar.outline.set_visible(False)

# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.tick_params(left=False, bottom=False)
# pfun.add_coast(ax, color='k')
pfun.dar(ax)
pfun.dar(axins)
ax.set_title('Surface mixed layer\n'+ r'cons. temp [$\degree$C]', fontsize=20,
             loc='Left', fontweight='bold')
# ax.text()

# Generate plot
plt.tight_layout
plt.subplots_adjust(bottom=0.001, top=0.9)
plt.show()

###################################################################
##             Plot ALKALINITY in surface mixed layer           ##  
################################################################### 

px, py = pfun.get_plon_plat(lon,lat)
vmin = 1000
vmax = 2500

# surface dye concentration
alk = ds_sml_may22['ALK'][0,:,:].values

# Initialize figure
fig,ax = plt.subplots(1,1, figsize=(7,9))

# plot values
# cmap = copy.copy(plt.get_cmap('viridis_r'))
# cmap.set_bad(color='white', alpha=1.0)
cmap = cmc.devon
cs = ax.pcolormesh(px,py,alk,vmin=vmin, vmax=vmax, cmap=cmap)

# Add Puget Sound Inset
# [x0, y0, width, height]
axins = ax.inset_axes([0.71, 0.012, 0.45, 0.6])
# plot values in inset
axins.pcolormesh(px, py, alk, vmin=vmin, vmax=vmax, cmap=cmap)
# Puget Sound limits
axins.set_xlim(-123.2, -122.1)
axins.set_ylim(46.95, 48.4)
axins.tick_params(left=False, bottom=False)
# format
axins.set_xticklabels([])
axins.set_yticklabels([])
for spine in axins.spines.values():
    spine.set_edgecolor('grey')
    spine.set_linewidth(2)

# add colorbar
cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
cbar.ax.tick_params(labelsize=18, rotation=30)
cbar.outline.set_visible(False)

# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.tick_params(left=False, bottom=False)
# pfun.add_coast(ax, color='k')
pfun.dar(ax)
pfun.dar(axins)
ax.set_title('Surface mixed layer\n'+ r'alkalinity [???]', fontsize=20,
             loc='Left', fontweight='bold')
# ax.text()

# Generate plot
plt.tight_layout
plt.subplots_adjust(bottom=0.001, top=0.9)
plt.show()

