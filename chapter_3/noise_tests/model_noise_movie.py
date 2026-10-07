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

##################################################################
# helper function from LO_tools

# format used for naming day folders
ds_fmt = '%Y.%m.%d'

def get_fn_list(list_type, Ldir, ds0, ds1, roms_out, his_num=2):
    """
    INPUT:
    A function for getting lists of history files.
    List items are Path objects
    
    NEW 2023.10.05: for list_type = 'hourly', if you pass his_num = 1
    it will start with ocean_his_0001.nc on the first day instead of the default which
    is to start with ocean_his_0025.nc on the day before.

    NEW 2025.06.20: for list_type = 'hourly0'
    which will start with ocean_his_0001.nc on the first day instead of the default which
    is to start with ocean_his_0025.nc on the day before.
    This is identical to passing his_num = 1, but may be more convenient, especially
    as we move to "continuation" start_type, which always writes an 0001 file.
    """
    dt0 = datetime.strptime(ds0, ds_fmt)
    dt1 = datetime.strptime(ds1, ds_fmt)
    dir0 = roms_out / Ldir['gtagex']
    if list_type == 'snapshot':
        # a single file name in a list
        his_string = ('0000' + str(his_num))[-4:]
        fn_list = [dir0 / ('f' + ds0) / ('ocean_his_' + his_string + '.nc')]
    elif list_type == 'hourly':
        # list of hourly files over a date range
        fn_list = Lfun.fn_list_utility(dt0,dt1,Ldir,his_num=his_num)
    elif list_type == 'hourly0':
        # list of hourly files over a date range, starting with 0001 of dt0.
        fn_list = Lfun.fn_list_utility(dt0,dt1,Ldir,his_num=1)
    elif list_type == 'daily':
        # list of history file 21 (Noon PST) over a date range
        fn_list = []
        date_list = Lfun.date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_his_0021.nc'
            fn_list.append(fn)
    elif list_type == 'lowpass':
        # list of lowpassed files (Noon PST) over a date range
        fn_list = []
        date_list = Lfun.date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'lowpassed.nc'
            fn_list.append(fn)
    elif list_type == 'average':
        # list of daily averaged files (Noon PST) over a date range
        fn_list = []
        date_list = Lfun.date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_avg_0001.nc'
            fn_list.append(fn)
    elif list_type == 'weekly':
        # like "daily" but at 7-day intervals
        fn_list = []
        date_list = Lfun.date_list_utility(dt0, dt1, daystep=7)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_his_0021.nc'
            fn_list.append(fn)
    elif list_type == 'allhours':
        # a list of all the history files in a directory
        # (this is the only list_type that actually finds files)
        in_dir = dir0 / ('f' + ds0)
        fn_list = [ff for ff in in_dir.glob('ocean_his*nc')]
        fn_list.sort()

    return fn_list

###################################################################
##                          User Inputs                          ##  
################################################################### 

vns = ['TIC','alkalinity']

d0= '2020.06.01'
d1 = '2020.10.31'
# for switching between adding alkalinity and not adding
d0end = '2020.06.30'
d1start = '2020.07.01'

list_type = 'average'


filetype = 'ocean_avg_0001.nc'
# filetype = 'ocean_his_0002.nc'

###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

model1 = 'Perturbation'
model2 = 'Baseline'

# where to put output figures
out_dir = Ldir['LOo'] / 'chapter_3' / 'figures' / 'noise_test'
Lfun.make_dir(out_dir)

# get his files
# ----------------------------------------------------------------

# gtagex of files to difference
Ldir_pert_1   = Lfun.Lstart(gridname='cas7', tag='t1dgeWB', ex_name='x11abd3monthscont')
Ldir_pert_2   = Lfun.Lstart(gridname='cas7', tag='t1dgeWB', ex_name='x11abd')
Ldir_base = Lfun.Lstart(gridname='cas7', tag='t1', ex_name='x11ab')

# get list of history files to plot (and skip ocean_his_0025 from previous day)
fn_list_pert_1   = Lfun.get_fn_list(list_type, Ldir_pert_1, d0, d0end, Ldir['roms_out'])
fn_list_pert_2   = Lfun.get_fn_list(list_type, Ldir_pert_2, d1start, d1, Ldir['roms_out'])
fn_list_pert = fn_list_pert_1 + fn_list_pert_2
fn_list_base = Lfun.get_fn_list(list_type, Ldir_base, d0, d1, Ldir['roms_out5'])

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values


###################################################################
##                      Binary differences                       ##  
################################################################### 

for vn in vns:

    for i,fn_pert in enumerate(fn_list_pert):

        # get model output
        fn_base = fn_list_base[i]
        ds_model1 = xr.open_dataset(fn_pert)
        ds_model2 = xr.open_dataset(fn_base)

        # Get data, and get rid of ocean_time dim (because this is at a single time)
        v1 = ds_model1[vn].squeeze()
        if vn == 'dye_01':
            v2 = v1 * 0 # there is no dye variable in the long hindcast
        else:
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

        # add alkalinity addition location
        inj_lon = -122.5674
        inj_lat = 48.1956
        ax.scatter(inj_lon,inj_lat,s=80, facecolors='none', edgecolors='deeppink')

        # format figure
        ax.text(0.96, 0.95, 'Differences', color='black', fontweight='bold', fontsize=12,
                transform=ax.transAxes, ha='right')
        ax.text(0.96, 0.92, 'No Differences', color='turquoise', fontweight='bold', fontsize=12,
                transform=ax.transAxes, ha='right')
        plt.suptitle('Locations where {} differs between runs at any s-level'.format(vn),fontsize=14,fontweight='bold')
        date = ds_model1["ocean_time"].dt.strftime("%Y.%m.%d").item()
        ax.set_title('{} and {} (daily avg on {})'.format(model1,model2,date))
        plt.xlabel('Lon', fontsize=12)
        plt.ylabel('Lat', fontsize=12)
        pfun.dar(ax)

        plt.tight_layout()

        # prepare a directory for results
        nouts = ('0000' + str(i))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = out_dir /('binary_' + vn)
        Lfun.make_dir(outfile)
        save_name = outfile / outname
        print('Plotting ' + str(fn_pert))
        sys.stdout.flush()
        plt.savefig(save_name)
        plt.close()

    # make movie
    if len(fn_list_pert) > 1:
        cmd_list = ['ffmpeg','-r','8','-i', str(out_dir / ('binary_' + vn) )+'/plot_%04d.png', '-vcodec', 'libx264',
            '-pix_fmt', 'yuv420p', '-crf', '25', str(out_dir / ('binary_' + vn))+'/movie.mp4']
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print('\n'+stdout.decode())
        if len(stderr) > 0:
            print('\n'+stderr.decode())

###################################################################
##                Vertical integral differences                  ##  
################################################################### 

for vn in vns:

    for i,fn_pert in enumerate(fn_list_pert):

        # get model output
        fn_base = fn_list_base[i]
        ds_model1 = xr.open_dataset(fn_pert)
        ds_model2 = xr.open_dataset(fn_base)

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
            vmin = -1e-6#-0.001
            vmax =  1e-6#0.001
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
        elif vn in ['alkalinity','dye_01','TIC']:
            vn_name = vn
            vmin = -0.1
            vmax =  0.1
            # vmin = -0.00001
            # vmax =  0.00001
        else:
            print('vmin and vmax not provided for '+ vn)

        # scale variable
        scale =  pinfo.fac_dict[vn_name]

        # Get model1 data
        vertint_vn_model1 = ds_model1[vn].sum(dim=vert_dim, skipna=False).values * scale
        vertint_vn_model2 = ds_model2[vn].sum(dim=vert_dim, skipna=False).values * scale

        # Get difference
        vertint_diff = vertint_vn_model1 - vertint_vn_model2

        # Initialize figure
        fig,ax = plt.subplots(1,1, figsize=(6,8))
        # plt.tight_layout()
        value = vertint_diff[0,:,:]
        # values = [surf_vn_model1,bott_vn_model1]

        newcmap = cmocean.tools.crop_by_percent(cmocean.cm.balance_r, 20, which='both', N=None)

        # plot values
        cs = ax.pcolormesh(lon, lat,value,vmin=vmin, vmax=vmax, cmap=newcmap)
        # cs = ax.pcolormesh(lon, lat,values[i],vmin=0, vmax=0.005)#, cmap=newcmap)
        cbar = fig.colorbar(cs)
        cbar.ax.tick_params(labelsize=12)
        cbar.outline.set_visible(False)
        # format figure
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.axis('off')
        ax.scatter(inj_lon,inj_lat,s=80, facecolors='none', edgecolors='deeppink')
        pfun.dar(ax)
        ax.set_title(vn + ' vertical integral difference ' + pinfo.units_dict[vn_name] + ' * m', fontsize=12)
        date = ds_model1["ocean_time"].dt.strftime("%Y.%m.%d").item()
        fig.suptitle('Perturbation minus Baseline ({})'.format(date),
                    fontsize=11, fontweight='bold')
        
        # prepare a directory for results
        nouts = ('0000' + str(i))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = out_dir /('vertint_' + vn)
        Lfun.make_dir(outfile)
        save_name = outfile / outname
        print('Plotting ' + str(fn_pert))
        sys.stdout.flush()
        plt.savefig(save_name)
        plt.close()

    # make movie
    if len(fn_list_pert) > 1:
        cmd_list = ['ffmpeg','-r','8','-i', str(out_dir / ('vertint_' + vn) )+'/plot_%04d.png', '-vcodec', 'libx264',
            '-pix_fmt', 'yuv420p', '-crf', '25', str(out_dir / ('vertint_' + vn))+'/movie.mp4']
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print('\n'+stdout.decode())
        if len(stderr) > 0:
            print('\n'+stderr.decode())

plt.close('all')
    