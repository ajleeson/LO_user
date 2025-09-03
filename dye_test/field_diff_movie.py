"""
Plot difference movie between surface/bottom values of specified state variable.
Calculates difference between two different runs
(Written to compare long model1 to N-less run)

From ipython: run field_diff_movie
Figures saved in LO_output/AL_custom_plots/noise_WestPoint/[vn]_diff.mp3

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

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys

Ldir = Lfun.Lstart()

plt.close('all')

###################################################################
##                          User Inputs                          ##  
################################################################### 

# USER OPTIONS ----------------------------------------------------

d0= '2012.10.07'
d1 = '2012.10.07'

list_type = 'hourly' #'snapshot', 'daily', 'hourly ', 'allhours'
his_num = 25# 11
timestep_interval = 60

# ----------------------------------------------------------------

# gtagex of files to difference
Ldir_WWTP   = Lfun.Lstart(gridname='cas7', tag='exdye2', ex_name='x11exdye2')
# Ldir_noWWTP = Lfun.Lstart(gridname='cas7', tag='exdye2', ex_name='x11exdye2')
Ldir_noWWTP = Lfun.Lstart(gridname='cas7', tag='exdye2duplicate', ex_name='x11exdye2')

# get list of history files to plot (and skip ocean_his_0025 from previous day)
fn_list_WWTP   = Lfun.get_fn_list(list_type, Ldir_WWTP,
    d0, d1, his_num)[0:his_num:]
fn_list_noWWTP = Lfun.get_fn_list(list_type, Ldir_noWWTP,
    d0, d1, his_num)[0:his_num:]

fn_list_WWTP[0]   = Ldir['roms_out'] / 'cas7_exdye2_x11exdye2' / 'f2012.10.07' / 'ocean_his_0001.nc'
# fn_list_noWWTP[0] = Ldir['roms_out'] / 'cas7_exdye2_x11exdye2' / 'f2012.10.07' / 'ocean_his_0001.nc'
fn_list_noWWTP[0] = Ldir['roms_out'] / 'cas7_exdye2duplicate_x11exdye2' / 'f2012.10.07' / 'ocean_his_0001.nc'


# vns = ['NH4'] #['NH4','NO3','oxygen','u','v','w','salt','temp'] # u, v, w, NO3, NH4 oxygen
vn = 'dye_01' # one variable at a time
date = '2012.10.07'

# filetype = 'ocean_avg_0001.nc'
# filetype = 'ocean_his_0002.nc'

# PLOTTING
# outdir0 = Ldir['LOo'] / 'AL_custom_plots' / 'expdecay_dye2' / 'same_run'
outdir0 = Ldir['LOo'] / 'AL_custom_plots' / 'expdecay_dye2' / 'two_different_runs'
Lfun.make_dir(outdir0)

if len(fn_list_WWTP) > 1:
    # prepare a directory for results if making a movie
    outdir = outdir0
    Lfun.make_dir(outdir / 'binary', clean=True)
    Lfun.make_dir(outdir / 'pcolormesh', clean=True)

###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

# # where to put output figures
# out_dir = Ldir['LOo'] / 'AL_custom_plots'
# Lfun.make_dir(out_dir)

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values

# # get model output
# fp_model1 = Ldir['roms_out'] / model1 / ('f'+date) / filetype
# fp_model2 = Ldir['roms_out'] / model2 / ('f'+date) / filetype
# ds_model1 = xr.open_dataset(fp_model1)
# ds_model2 = xr.open_dataset(fp_model2)

###################################################################
##                   Binary differences movie                    ##  
################################################################### 

for i,fn_WWTP in enumerate(fn_list_WWTP):

    # get model output
    fn_noWWTP = fn_list_noWWTP[i]
    ds_model1 = xr.open_dataset(fn_WWTP)
    ds_model2 = xr.open_dataset(fn_noWWTP)

    # Get data, and get rid of ocean_time dim (because this is at a single time)
    v1 = ds_model1[vn].squeeze()
    v2 = ds_model2['dye_02'].squeeze()

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
    # plt.suptitle('Locations where dye_01 and dye_02 differs between runs at any s-level\nWithin same run',
    #              fontsize=14,fontweight='bold')
    plt.suptitle('Locations where dye_01 and dye_02 differs between runs at any s-level\nBetween two different runs',
                 fontsize=14,fontweight='bold')
    plt.xlabel('Lon', fontsize=12)
    plt.ylabel('Lat', fontsize=12)
    pfun.dar(ax)

    ax.text(0.7, 0.95, 'Timestep {}'.format(i*timestep_interval), color='black', fontweight='bold', fontsize=12,
            transform=ax.transAxes)

    # # Salish Sea
    # ax.set_ylim([46.5,50])
    # ax.set_xlim([-126.5,-122])

    # West Point
    ax.set_ylim([47.4,48])
    ax.set_xlim([-122.9,-122.1])

    plt.tight_layout()
    # plt.show()

    # prepare a directory for results
    nouts = ('0000' + str(i))[-4:]
    outname = 'plot_' + nouts + '.png'
    outfile = outdir /'binary' / outname
    print('Plotting ' + str(fn_WWTP))
    sys.stdout.flush()
    plt.savefig(outfile)
    plt.close()

# make movie
if len(fn_list_WWTP) > 1:
    cmd_list = ['ffmpeg','-r','3','-i', str(outdir / 'binary')+'/plot_%04d.png', '-vcodec', 'libx264',
        '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir / 'binary')+'/movie.mp4']
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    if len(stdout) > 0:
        print('\n'+stdout.decode())
    if len(stderr) > 0:
        print('\n'+stderr.decode())

###################################################################
##                    Colormap differences                       ##  
################################################################### 

for i,fn_WWTP in enumerate(fn_list_WWTP):

    print(i)

    # get model output
    fn_noWWTP = fn_list_noWWTP[i]
    ds_model1 = xr.open_dataset(fn_WWTP)
    ds_model2 = xr.open_dataset(fn_noWWTP)

    # Get data, and get rid of ocean_time dim (because this is at a single time)
    v1 = ds_model1[vn].squeeze()
    v2 = ds_model2['dye_02'].squeeze()

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
        vmin = -0.00001
        vmax =  0.00001
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
    elif vn == 'dye_01':
        vn_name = vn
        vmin = -0.00001
        vmax =  0.00001
    else:
        print('vmin and vmax not provided for '+ vn)

    # Get model1 data
    surf_vn_model1 = ds_model1[vn][0,-1,:,:].values
    bott_vn_model1 = ds_model1[vn][0,0,:,:].values
    # Get model2 data
    surf_vn_model2 = ds_model2['dye_02'][0,-1,:,:].values
    bott_vn_model2 = ds_model2['dye_02'][0,0,:,:].values
    # Get difference
    surf_diff = (surf_vn_model1 - surf_vn_model2)
    bott_diff = (bott_vn_model1 - bott_vn_model2)

    # Initialize figure
    fig = plt.figure(figsize=(16,8)) # 15,11 for Puget sound and 18,8 for Salish Sea
    plt.tight_layout()

    subplotnums = [121,122]
    stexts = ['Surface','Bottom']
    values = [surf_diff,bott_diff]

    newcmap = cmocean.tools.crop_by_percent(cmocean.cm.balance_r, 20, which='both', N=None)

    # loop through all of the plots we need to make
    for j,stext in enumerate(stexts):

        # add water/land
        ax = fig.add_subplot(subplotnums[j])

        # plot values
        cs = ax.pcolormesh(lon, lat,values[j],vmin=vmin, vmax=vmax, cmap=newcmap)
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
        ax.set_title('Dye difference at ' + stext + '[kg/m3]', fontsize=16)
        fig.suptitle('dye_01 minus dye_02\nBetween two different runs' + date,
                    fontsize=14, fontweight='bold')
        
        if j == 1:
            ax.text(0.6, 0.1, 'Hour {}'.format(i), color='black', fontweight='bold', fontsize=12,
            transform=ax.transAxes)
        
    # prepare a directory for results
    nouts = ('0000' + str(i))[-4:]
    outname = 'plot_' + nouts + '.png'
    outfile = outdir /'pcolormesh' / outname
    print('Plotting ' + str(fn_WWTP))
    sys.stdout.flush()
    plt.savefig(outfile)
    plt.close()

# make movie
if len(fn_list_WWTP) > 1:
    cmd_list = ['ffmpeg','-r','3','-i', str(outdir / 'pcolormesh')+'/plot_%04d.png', '-vcodec', 'libx264',
        '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir / 'pcolormesh')+'/movie.mp4']
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    if len(stdout) > 0:
        print('\n'+stdout.decode())
    if len(stderr) > 0:
        print('\n'+stderr.decode())