"""
Helper script used to generate a roms history file for cas7
using a roms history files from cas6.

An intended use-case is to generate initial conditions for cas7
based on history files from a cas6 run.

From ipython:
run make_cas7his_from_cas6his.py -oldg cas6 -newg cas7
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Ellipse
import cmocean
import argparse
import sys
from pathlib import Path
from importlib import reload
import Ofun
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun, zrfun

reload(Ofun)

#################################################################################
#                              Argument parsing                                 #
#################################################################################

# read arguments
parser = argparse.ArgumentParser()
parser.add_argument('-oldg', '--oldgridname', type=str) # original grid
parser.add_argument('-newg', '--newgridname', type=str) # new grid
args = parser.parse_args()

# make sure grid names were provided
argsd = args.__dict__
if argsd['oldgridname'] == None or argsd['oldgridname'] == None:
    print('*** Missing required argument to forcing_argfun.intro(): gridnames')
    sys.exit()

Ldir = Lfun.Lstart()

#################################################################################
#                               Get grid data                                   #
#################################################################################

# get the old grid data
oldgrid_fn = Ldir['data'] / 'grids' / args.oldgridname / 'grid.nc'
grid_ds_old = xr.open_dataset(oldgrid_fn)
z_old = -grid_ds_old.h.values
# transpose mask to work with history file
# mask_rho_old = np.transpose(grid_ds_old.mask_rho.values)
mask_rho_old = grid_ds_old.mask_rho.values
mask_u_old = grid_ds_old.mask_u.values
mask_v_old = grid_ds_old.mask_v.values
lonr_old = grid_ds_old.lon_rho.values
latr_old = grid_ds_old.lat_rho.values
lonu_old = grid_ds_old.lon_u.values
latu_old = grid_ds_old.lat_u.values
lonv_old = grid_ds_old.lon_v.values
latv_old = grid_ds_old.lat_v.values
X_old = lonr_old[0,:] # grid cell X values
Y_old = latr_old[:,0] # grid cell Y values
plon, plat = pfun.get_plon_plat(lonr_old,latr_old)

# get the new grid data
newgrid_fn = Ldir['data'] / 'grids' / args.newgridname / 'grid.nc'
grid_ds_new = xr.open_dataset(newgrid_fn)
z_new = -grid_ds_new.h.values
# transpose mask to work with history file
# mask_rho_new = np.transpose(grid_ds_new.mask_rho.values)
mask_rho_new = grid_ds_new.mask_rho.values
mask_u_new = grid_ds_new.mask_u.values
mask_v_new = grid_ds_new.mask_v.values
lonr_new = grid_ds_new.lon_rho.values
latr_new = grid_ds_new.lat_rho.values
lonu_new = grid_ds_new.lon_u.values
latu_new = grid_ds_new.lat_u.values
lonv_new = grid_ds_new.lon_v.values
latv_new = grid_ds_new.lat_v.values
X_new = lonr_new[0,:] # grid cell X values
Y_new = latr_new[:,0] # grid cell Y values
plon, plat = pfun.get_plon_plat(lonr_new,latr_new)

#################################################################################
#                                  Plot grid                                    #
#################################################################################

plt.close('all')
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1,1,1)
mask_sum = mask_rho_old+mask_rho_new
cmap = colors.ListedColormap(['white', 'red', 'powderblue'])
plt.pcolormesh(plon, plat,mask_sum,vmin=0,vmax=2,cmap=cmap)
pfun.dar(ax)
plt.title('Grid Differences (Red)',fontsize=14)
# Puget Sound
ax.set_xlim([-123.3,-122.1])
ax.set_ylim([47,48.5])
plt.show()

#################################################################################
#                            Get history file data                              #
#################################################################################

# get old history file
# in_dir0 = Path('/pgdat1/auroral/LO_roms')
in_dir0 = Ldir['roms_out']
in_fn = in_dir0 / 'cas6_v00_uu0mb' / 'f2016.12.31' / 'ocean_his_0025.nc'

# create spot to save new history file
out_dir = Ldir['roms_out'] / 'cas7_trapsV00_meV00' / 'f2016.12.31'
Lfun.make_dir(out_dir)
out_fn = out_dir / 'ocean_his_0025.nc'

# create place to save figures
fig_dir =  Ldir['LOo'] / 'testing' / 'ocn_ini'

# get old grid history file in dataset form
ds_old = xr.open_dataset(in_fn, decode_times=False)

# get variable names
vars = list(ds_old.keys())

#################################################################################
#                   New history file is a copy of the old one                   #
#################################################################################

ds_new = ds_old.copy()

#################################################################################
#                           Get the new grid indices                           #
#################################################################################

# get the new mask minus the old mask
mask_diff_rho = mask_rho_new-mask_rho_old
mask_diff_u = mask_u_new-mask_u_old
mask_diff_v = mask_v_new-mask_v_old

# where mask_diff = 1 are locations of new water cells
# get the indices of these new water cells
i_rho = np.where(mask_diff_rho == 1)[0]
j_rho = np.where(mask_diff_rho == 1)[1]
i_u = np.where(mask_diff_u == 1)[0]
j_u = np.where(mask_diff_u == 1)[1]
i_v = np.where(mask_diff_v == 1)[0]
j_v = np.where(mask_diff_v == 1)[1]

#################################################################################
#                            Plot state variables                               #
#################################################################################

# generate directory to save figures
fig_dir = fig_dir / 'figures'
Lfun.make_dir(fig_dir)

# Plot all state variables
for ii,var in enumerate(vars):

    # get dimension of variables
    dim = ds_new[var].dims
    unit = ds_new[var].units

    # Figure out which grid the variable is on (rho, u, or v)
    if 'eta_rho' in dim:
        lat = latr_new
        lon = lonr_new
        i_inds = i_rho
        j_inds = j_rho
    elif 'eta_u' in dim:
        lat = latu_new
        lon = lonu_new
        i_inds = i_u
        j_inds = j_u
    elif 'eta_v' in dim:
        lat = latv_new
        lon = lonv_new
        i_inds = i_v
        j_inds = j_v
    else:
        print('Grid unclear')
        pass

    # Figure out if variable is 2D or 3D
    if 's_rho' in dim:
        for s in range(30):
            ds_new[var].values[0,s,:,:] = Ofun.extrap_nearest_to_masked(lon, lat,
                                        ds_old[var].values[0,s,:,:], i_inds, j_inds)
        depths = ['Surface','Bottom']
    else:
        ds_new[var].values[0,:,:] = Ofun.extrap_nearest_to_masked(lon, lat,
                                        ds_old[var].values[0,:,:], i_inds, j_inds)
        depths = ['Depth-uniform']

    # get lon and lat for plotting
    plon, plat = pfun.get_plon_plat(lon,lat)

    # loop through different depth layers
    for depth in depths:

        # Create figure
        fig = plt.figure(figsize=(12,8))

        if depth == 'Surface':
            val_old = ds_old[var].values[0,-1,:,:]
            val_new = ds_new[var].values[0,-1,:,:]
        elif depth == 'Bottom':
            val_old = ds_old[var].values[0,0,:,:]
            val_new = ds_new[var].values[0,0,:,:]
        elif depth == 'Depth-uniform':
            val_old = ds_old[var].values[0,:,:]
            val_new = ds_new[var].values[0,:,:]


        # get bounds
        vmin = np.nanmin(ds_new[var].values)
        vmax = np.nanmax(ds_new[var].values)

        if vmin == vmax:
            vmin = -0.1
            vmax = 0.1
        if var in ['u','ubar','v','vbar']:
            cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)
        else:
            cmap = cmocean.cm.thermal

        # Old history file #########################

        # Agate Pass (Old)
        ax1 = fig.add_subplot(2,2,1)
        ax1.pcolormesh(plon, plat,val_old, vmin=vmin, vmax=vmax, cmap=cmap)
        pfun.dar(ax1)
        ax1.set_title('Agate Pass ({})'.format(args.oldgridname))
        ax1.set_xlim([-122.7,-122.3])
        ax1.set_ylim([47.6,47.8])
        plt.xticks([])
        ax1.add_patch(Ellipse((-122.56, 47.71), 0.06,0.09, angle=135, linewidth = 2,
                            edgecolor='k', facecolor='none'))
        # Swinomish Channel (Old)
        ax2 = fig.add_subplot(2,2,2)
        ax2.pcolormesh(plon, plat,val_old, vmin=vmin, vmax=vmax, cmap=cmap)
        pfun.dar(ax2)
        ax2.set_title('Swinomish Channel ({})'.format(args.oldgridname))
        ax2.set_xlim([-122.8,-122.4])
        ax2.set_ylim([48.3,48.5])
        plt.xticks([])
        ax2.add_patch(Ellipse((-122.51, 48.41), 0.06,0.13, angle=0, linewidth = 2,
                            edgecolor='k', facecolor='none'))

        # Old history file #########################

        # Agate Pass (New)
        ax3 = fig.add_subplot(2,2,3)
        ax3.pcolormesh(plon, plat,val_new, vmin=vmin, vmax=vmax, cmap=cmap)
        pfun.dar(ax3)
        ax3.set_title('Agate Pass ({})'.format(args.newgridname))
        ax3.set_xlim([-122.7,-122.3])
        ax3.set_ylim([47.6,47.8])
        plt.xticks(rotation=30)
        ax3.add_patch(Ellipse((-122.56, 47.71), 0.06,0.09, angle=135, linewidth = 2,
                            edgecolor='k', facecolor='none'))
        # Swinomish Channel (New)
        ax4 = fig.add_subplot(2,2,4)
        im = ax4.pcolormesh(plon, plat,val_new, vmin=vmin, vmax=vmax, cmap=cmap)
        pfun.dar(ax4)
        ax4.set_title('Swinomish Channel ({})'.format(args.newgridname))
        ax4.set_xlim([-122.8,-122.4])
        ax4.set_ylim([48.3,48.5])
        plt.xticks(rotation=30)
        ax4.add_patch(Ellipse((-122.51, 48.41), 0.06,0.13, angle=0, linewidth = 2,
                            edgecolor='k', facecolor='none'))

        # add colorbar
        fig.subplots_adjust(right=0.9)
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        fig.colorbar(im, cax=cbar_ax)

        plt.suptitle('{} {} [{}]'.format(depth, var,unit),fontsize=14)
        # Save figure
        figname = var + '_' + depth + '.png'
        save_path = fig_dir / figname
        fig.savefig(save_path)
        plt.close('all')

print('Plotting done')

#################################################################################
#                       Check and save new history file                         #
#################################################################################

# count the number of nans before and after
for var in vars:

    # get dimension of variables
    dim = ds_new[var].dims

    # Figure out which grid the variable is on (rho, u, or v)
    if 'eta_rho' in dim:
        lat = latr_new
        lon = lonr_new
        i_inds = i_rho
        j_inds = j_rho
    elif 'eta_u' in dim:
        lat = latu_new
        lon = lonu_new
        i_inds = i_u
        j_inds = j_u
    elif 'eta_v' in dim:
        lat = latv_new
        lon = lonv_new
        i_inds = i_v
        j_inds = j_v
    else:
        print('Grid unclear')
        pass

    # Figure out if variable is 2D or 3D
    if 's_rho' in dim:
        val_old = ds_old[var].values[0,:,i_inds,j_inds]
        val_new = ds_new[var].values[0,:,i_inds,j_inds]
        dimension = '3D'
    else:
        val_old = ds_old[var].values[0,i_inds,j_inds]
        val_new = ds_new[var].values[0,i_inds,j_inds]
        dimension = '2D'

    # print and count the number of nans in the new grid cells
    num_nans_old = np.count_nonzero(np.isnan(val_old))
    num_nans_new = np.count_nonzero(np.isnan(val_new))
    print('\n=======================')
    print('{} ({})'.format(var,dimension))
    print('Old cells: {} Nans'.format(num_nans_old))
    print('New cells: {} Nans'.format(num_nans_new))

# Save the new history file
# and save to NetCDF
Enc_dict = {vn:zrfun.enc_dict for vn in ds_new.data_vars}
ds_new.to_netcdf(out_fn, encoding=Enc_dict)