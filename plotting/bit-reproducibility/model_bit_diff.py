"""
Plot difference between surface/bottom values of specified state variable.
Calculates difference between two different runs
(Written to compare long model1 to N-less run)

From ipython: run model_bit_diff

"""

###################################################################
##                       import packages                         ##  
###################################################################      

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import pandas as pd
import cmocean
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

vns = ['NH4'] #['NH4','NO3','oxygen','u','v','w','salt','temp'] # u, v, w, NO3, NH4 oxygen
date = '2012.10.07'

# filetype = 'ocean_avg_0001.nc'
filetype = 'ocean_his_0002.nc'

###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

model1 = 'cas7_newtraps00_debugx11ab'
# model2 = 'cas7_newtraps01_debugx11ab'
model2 = 'cas7_newtraps00noWP_debugx11ab'
# model1 = 'cas7_newtraps_x11ab'
# model2 = 'cas7_newtrapsnoN_x11ab'

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values

# get model output
fp_model1 = Ldir['roms_out'] / model1 / ('f'+date) / filetype
fp_model2 = Ldir['roms_out'] / model2 / ('f'+date) / filetype
ds_model1 = xr.open_dataset(fp_model1)
ds_model2 = xr.open_dataset(fp_model2)


###################################################################
##                      Binary differences                       ##  
################################################################### 

for vn in vns:

    # Get data, and get rid of ocean_time dim (because this is at a single time)
    v1 = ds_model1[vn].squeeze()
    v2 = ds_model2[vn].squeeze()

    # Identify vertical and horizontal dims
    if 's_rho' in v1.dims:
        vert_dim = 's_rho'
    elif 's_w' in v1.dims:
        vert_dim = 's_w'
    else:
        raise ValueError(f"No vertical dimension found for {vn}")

    if ('eta_rho' in v1.dims) and ('xi_rho' in v1.dims):
        h_dims = ('eta_rho', 'xi_rho')
        lon = ds_model1['lon_rho']
        lat = ds_model1['lat_rho']
    elif ('eta_u' in v1.dims) and ('xi_u' in v1.dims):
        h_dims = ('eta_u', 'xi_u')
        lon = ds_model1['lon_u']
        lat = ds_model1['lat_u']
    elif ('eta_v' in v1.dims) and ('xi_v' in v1.dims):
        h_dims = ('eta_v', 'xi_v')
        lon = ds_model1['lon_v']
        lat = ds_model1['lat_v']
    elif ('eta_psi' in v1.dims) and ('xi_psi' in v1.dims):
        h_dims = ('eta_psi', 'xi_psi')
        lon = ds_model1['lon_psi']
        lat = ds_model1['lat_psi']
    else:
        raise ValueError(f"Unknown grid type for variable '{vn}'.")

    # Compute strict difference (True where different, False where equal)
    diff_mask = (v1 != v2) | (v1.isnull() != v2.isnull())
    diff_2d = diff_mask.any(dim=vert_dim)

    # Mask out locations where both are NaN at all depths
    both_nan = v1.isnull() & v2.isnull()
    diff_2d = diff_2d.where(~both_nan.any(dim=vert_dim))

    # Convert to numeric: 0 = same, 1 = diff, NaN = both missing
    plot_data = diff_2d.astype(float)

    # Set up colormap: black = 0, lightblue = 1, white = NaN
    cmap = mcolors.ListedColormap(['paleturquoise', 'black'])
    bounds = [-0.5, 0.5, 1.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    # Plotting -------------------------------------------------- 

    # Initialize figure
    fig, ax = plt.subplots(1,1, figsize=(10, 8))

    # plot
    plt.pcolormesh(lon, lat, plot_data, cmap=cmap, norm=norm, shading='auto')
    # cbar = plt.colorbar(map, ticks=[0, 1])
    # cbar.ax.set_yticklabels(['No differences', 'Differences'],fontsize=12)

    # add West Point location
    wwtp_fn = Ldir['data'] / 'trapsD01' / 'processed_data'/ 'wwtp_data_wasielewski_etal_2024.nc'
    wwtp_data_ds = xr.open_dataset(wwtp_fn)
    WP_lat = wwtp_data_ds.lat.sel(source = wwtp_data_ds.name=='King County West Point WWTP').values
    WP_lon = wwtp_data_ds.lon.sel(source = wwtp_data_ds.name=='King County West Point WWTP').values
    ax.scatter(WP_lon,WP_lat,s=80, facecolors='none', edgecolors='deeppink')

    # format figure
    ax.text(0.04, 0.95, 'Differences', color='black', fontweight='bold', fontsize=12,
            transform=ax.transAxes)
    ax.text(0.04, 0.92, 'No Differences', color='turquoise', fontweight='bold', fontsize=12,
            transform=ax.transAxes)
    plt.suptitle('Locations where {} differs between runs at any s-level'.format(vn),fontsize=14,fontweight='bold')
    ax.set_title('{} and {} ({})'.format(model1,model2, filetype))
    plt.xlabel('Lon', fontsize=12)
    plt.ylabel('Lat', fontsize=12)
    pfun.dar(ax)

    # # Salish Sea
    # ax.set_ylim([46.5,50])
    # ax.set_xlim([-126.5,-122])

    # West Point
    ax.set_ylim([47.4,48])
    ax.set_xlim([-122.9,-122.1])

    plt.tight_layout()
    plt.show()

###################################################################
##                    Colormap differences                       ##  
################################################################### 

for vn in vns:

    # Get data, and get rid of ocean_time dim (because this is at a single time)
    v1 = ds_model1[vn].squeeze()
    v2 = ds_model2[vn].squeeze()

    if ('eta_rho' in v1.dims) and ('xi_rho' in v1.dims):
        h_dims = ('eta_rho', 'xi_rho')
        lon = ds_model1['lon_rho']
        lat = ds_model1['lat_rho']
    elif ('eta_u' in v1.dims) and ('xi_u' in v1.dims):
        h_dims = ('eta_u', 'xi_u')
        lon = ds_model1['lon_u']
        lat = ds_model1['lat_u']
    elif ('eta_v' in v1.dims) and ('xi_v' in v1.dims):
        h_dims = ('eta_v', 'xi_v')
        lon = ds_model1['lon_v']
        lat = ds_model1['lat_v']
    elif ('eta_psi' in v1.dims) and ('xi_psi' in v1.dims):
        h_dims = ('eta_psi', 'xi_psi')
        lon = ds_model1['lon_psi']
        lat = ds_model1['lat_psi']
    else:
        raise ValueError(f"Unknown grid type for variable '{vn}'.")

    # set bounds
    if vn in ['NO3','NH4']:
        vn_name = vn
        vmin = -0.001
        vmax =  0.001
    elif vn in ['u','v','w']:
        vn_name = vn
        vmin = -0.00001#-0.01
        vmax =  0.00001#0.01
    elif vn == 'oxygen':
        vn_name = vn
        vmin = -0.001
        vmax =  0.001
    elif vn == 'salt':
        vn_name = vn
        vmin = -0.00001
        vmax =  0.00001
    elif vn == 'temp':
        vn_name = vn
        vmin = -0.00001
        vmax =  0.00001
    else:
        print('vmin and vmax not provided for '+ vn)

    # scale variable
    scale =  pinfo.fac_dict[vn_name]

    # Get model1 data
    surf_vn_model1 = ds_model1[vn][0,-1,:,:].values
    bott_vn_model1 = ds_model1[vn][0,0,:,:].values
    # Get model2 data
    surf_vn_model2 = ds_model2[vn][0,-1,:,:].values
    bott_vn_model2 = ds_model2[vn][0,0,:,:].values
    # Get difference
    surf_diff = (surf_vn_model1 - surf_vn_model2) * scale
    bott_diff = (bott_vn_model1 - bott_vn_model2) * scale

    # Initialize figure
    fig = plt.figure(figsize=(16,8)) # 15,11 for Puget sound and 18,8 for Salish Sea
    plt.tight_layout()

    subplotnums = [121,122]
    stexts = ['Surface','Bottom']
    values = [surf_diff,bott_diff]
    # values = [surf_vn_model1,bott_vn_model1]

    newcmap = cmocean.tools.crop_by_percent(cmocean.cm.balance_r, 20, which='both', N=None)

    # loop through all of the plots we need to make
    for i,stext in enumerate(stexts):

        # add water/land
        ax = fig.add_subplot(subplotnums[i])

        # plot values
        cs = ax.pcolormesh(lon, lat,values[i],vmin=vmin, vmax=vmax, cmap=newcmap)
        # cs = ax.pcolormesh(lon, lat,values[i],vmin=0, vmax=0.005)#, cmap=newcmap)
        cbar = fig.colorbar(cs)
        cbar.ax.tick_params(labelsize=14)
        cbar.outline.set_visible(False)
        # format figure
        # ax.set_xlim([xmin,xmax])
        # ax.set_ylim([ymin,ymax])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.axis('off')
        ax.set_ylim([47.4,48])
        ax.set_xlim([-122.9,-122.1])
        ax.scatter(WP_lon,WP_lat,s=80, facecolors='none', edgecolors='deeppink')
        # pfun.add_coast(ax, color='k')
        pfun.dar(ax)
        ax.set_title(vn + ' difference at ' + stext + pinfo.units_dict[vn_name], fontsize=16)
        fig.suptitle('{} minus {}\n'.format(model1,model2) + date + ' {}'.format(filetype),
                    fontsize=14, fontweight='bold')
        

###################################################################
##                    Property-property plot                     ##  
################################################################### 

for vn in vns:

    # West Point position
    lat0, lon0 = WP_lat, WP_lon

    # Compute residual and squeeze ocean_time
    res = ds_model1[vn].squeeze() - ds_model2[vn].squeeze()  # shape (s_rho, eta_rho, xi_rho)

    # Load lat/lon on rho grid (shape (eta_rho, xi_rho))
    lat = ds_model1['lat_rho'].values
    lon = ds_model1['lon_rho'].values

    # get residual data: (s_rho, eta_rho, xi_rho)
    res_data = res.values  # (30, 1302, 663)

    # Broadcast lat/lon to match res vertical dimension for flattening
    # lat and lon are 2D (eta_rho, xi_rho)
    # We want 3D (s_rho, eta_rho, xi_rho) by repeating lat/lon vertically
    lat3d = np.broadcast_to(lat, res_data.shape)
    lon3d = np.broadcast_to(lon, res_data.shape)

    # Flatten all arrays to 1D for scatter plot
    flat_res = res_data.flatten()
    flat_lat = lat3d.flatten()
    flat_lon = lon3d.flatten()

    # Mask to keep only finite, nonzero residuals
    valid = np.isfinite(flat_res) & (flat_res != 0)
    flat_res = flat_res[valid]
    flat_lat = flat_lat[valid]
    flat_lon = flat_lon[valid]

    # Get distance from West Point
    x, y = zfun.ll2xy(flat_lon, flat_lat, lon0, lat0)
    distance = np.sqrt(x**2 + y**2) / 1000  # convert to km

    # convert to proper units
    scale =  pinfo.fac_dict[vn]
    residuals = flat_res * scale
    # get positive and negative residual values, and take absolute value
    pos_residuals = flat_res[flat_res > 0] * scale
    pos_distances = distance[flat_res > 0] * scale
    neg_residuals = flat_res[flat_res < 0] * scale * -1
    neg_distances = distance[flat_res < 0] * scale

    # # get mean residual value
    # mean_res = np.nanmean(flat_res) * scale

    # Plot
    fig, axes = plt.subplots(2,1, figsize=(8, 8),sharex=True)
    ax = axes.ravel()

    # Residual
    ax[0].scatter(distance, residuals, s=5, alpha=0.3,color='black',zorder=5)
    # format figure
    ax[0].set_ylabel(f'{vn} residual {pinfo.units_dict[vn]}\n(West Point - no West Point)',fontsize=12)
    ax[0].grid(True,color='gainsboro')


    # Absolute value and log transform
    ax[1].scatter(pos_distances, pos_residuals, s=5, alpha=0.3,color='royalblue',
                  label='Positive Residual',zorder=5)
    ax[1].scatter(neg_distances, neg_residuals, s=5, alpha=0.3,color='crimson',
                  label='Negative Residual',zorder=5)
    # format figure
    ax[1].set_yscale('log')
    ax[1].set_xlim([0,25])
    ax[1].set_xlabel(f'Distance from West Point [km]',fontsize=12)
    ax[1].set_ylabel(f'{vn} residual {pinfo.units_dict[vn]}\n(West Point - no West Point)',fontsize=12)
    ax[1].grid(True,color='gainsboro')
    ax[1].legend(loc='upper right',fontsize=12)


    plt.suptitle(f'Residuals of {vn} vs. distance from West Point',fontsize=14)
    plt.tight_layout()
    plt.show()