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

d0= '2020.05.31'
# d1 = '2020.05.31'
d1 = '2020.06.30'

list_type = 'average'

filetype = 'ocean_avg_0001.nc'
# filetype = 'ocean_his_0002.nc'

###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

# where to put output figures
out_dir = Ldir['LOo'] / 'chapter_3' / 'figures' / 'noise_test' / 'noise_diffs_og_twin'
Lfun.make_dir(out_dir)

# get his files
# ----------------------------------------------------------------

# twin variables model output
Ldir   = Lfun.Lstart(gridname='cas7', tag='t1twin', ex_name='x11aboaetwin')
fn_list_twinvar   = get_fn_list(list_type, Ldir, d0, d1, Ldir['roms_out'])

# og model output
Ldir_pert   = Lfun.Lstart(gridname='cas7', tag='t1alkNOISE', ex_name='x11ab')
Ldir_base = Lfun.Lstart(gridname='cas7', tag='t1', ex_name='x11ab')
fn_list_pert   = get_fn_list(list_type, Ldir_pert, d0, d1, Ldir['roms_out'])
fn_list_base = get_fn_list(list_type, Ldir_base, d0, d1, Ldir['roms_out5'])         # CHANGE ON APOGEE
# fn_list_base = get_fn_list(list_type, Ldir_base, d0, d1, Ldir['roms_out'])

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
plon, plat = pfun.get_plon_plat(lon,lat)
lons = lon[0,:]
lats = lat[:,0]

# get injection location
site_lon_test = -127.25
site_lat_test = 45.5
# get nearest lat/lon indices to random point selected on google maps
site_y = min(range(len(lats)), key=lambda i: abs(lats[i]-site_lat_test))
site_x = min(range(len(lons)), key=lambda i: abs(lons[i]-site_lon_test))
# get actual lat/lon
inj_lon = lons[site_x]
inj_lat = lats[site_y]

# format for gridded subplots
def generate_axes(fig):
    gridspec = fig.add_gridspec(nrows=12, ncols=12)
    ax = {}
    ax['his_og'] = fig.add_subplot(gridspec[0:4, 0:10])
    ax['his_tv'] = fig.add_subplot(gridspec[4:8, 0:10], sharex=ax['his_og'], sharey=ax['his_og'])
    ax['his_dif'] = fig.add_subplot(gridspec[8:12, 0:10], sharex=ax['his_og'], sharey=ax['his_og'])
    ax['map_og'] = fig.add_subplot(gridspec[0:4, 10:12])
    ax['map_tv'] = fig.add_subplot(gridspec[4:8, 10:12], sharex=ax['map_og'], sharey=ax['map_og'])
    ax['map_dif'] = fig.add_subplot(gridspec[8:12, 10:12], sharex=ax['map_og'], sharey=ax['map_og'])
    return ax

for vn in vns:

    for i,fn_twin in enumerate(fn_list_twinvar):

        # Initialize figure
        fig = plt.figure(figsize=(12,9))
        ax = generate_axes(fig)
        plt.suptitle(vn + ' vertical integral differences', fontsize=14, fontweight='bold')

###################################################################
##            Vertical integral differences (og model)           ##  
################################################################### 

        # get model output
        fn_base = fn_list_base[i]
        fn_pert = fn_list_pert[i]
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

        # set bounds
        vn_name = vn
        vmin = -0.1
        vmax = 0.1

        # Get model1 data
        vertint_vn_model1 = ds_model1[vn].sum(dim=vert_dim, skipna=False).values
        vertint_vn_model2 = ds_model2[vn].sum(dim=vert_dim, skipna=False).values

        # Get difference
        vertint_diff = vertint_vn_model1 - vertint_vn_model2
        vertint_diff_og = vertint_diff[0,:,:]
        # values = [surf_vn_model1,bott_vn_model1]

        newcmap = cmocean.tools.crop_by_percent(cmocean.cm.balance_r, 20, which='both', N=None)

        # plot values
        cs = ax['map_og'].pcolormesh(plon, plat,vertint_diff_og,vmin=vmin, vmax=vmax, cmap=newcmap)
        # cs = ax.pcolormesh(lon, lat,values[i],vmin=0, vmax=0.005)#, cmap=newcmap)
        # format figure
        ax['map_og'].set_yticklabels([])
        ax['map_og'].set_xticklabels([])
        ax['map_og'].scatter(inj_lon,inj_lat,s=80, facecolors='none', edgecolors='deeppink')
        pfun.dar(ax['map_og'])
        # ax['map_og'].set_title(vn + ' vertical integral difference (mmol/m2)', fontsize=12)
        date = ds_model1["ocean_time"].dt.strftime("%Y.%m.%d").item()
        # ax['map_og'].text(0.7,0.2,'{}'.format(date),transform=ax['map_og'].transAxes, fontsize=12, fontweight='bold', color='black')
        # fig.suptitle('Perturbation minus Baseline ({})'.format(date),
        #             fontsize=11, fontweight='bold')

###################################################################
##            Vertical integral differences (twin var)           ##  
################################################################### 

        if vn == 'TIC':
            vn_twin = 'dTIC'
        elif vn == 'alkalinity':
            vn_twin = 'dTA'

        # get model output
        ds_twinvar = xr.open_dataset(fn_twin)

        # Get data, and get rid of ocean_time dim (because this is at a single time)
        v = ds_twinvar[vn_twin].squeeze()

        # Get model data
        vertint_vn_model = ds_twinvar [vn_twin].sum(dim=vert_dim, skipna=False).values

        # Get difference
        vertint_diff = vertint_vn_model              
        vertint_diff_twinvar = vertint_diff[0,:,:]
        # values = [surf_vn_model1,bott_vn_model1]         

        # plot values
        cs = ax['map_tv'].pcolormesh(plon, plat, vertint_diff_twinvar,vmin=vmin, vmax=vmax, cmap=newcmap)
        # cs = ax.pcolormesh(lon, lat,values[i],vmin=0, vmax=0.005)#, cmap=newcmap)
        # format figure
        ax['map_tv'].set_yticklabels([])
        ax['map_tv'].set_xticklabels([])
        ax['map_tv'].scatter(inj_lon,inj_lat,s=80, facecolors='none', edgecolors='deeppink')
        pfun.dar(ax['map_tv'])
        # ax['map_tv'].set_title(vn + ' vertical integral difference (mmol/m2)', fontsize=12)
        date = ds_twinvar["ocean_time"].dt.strftime("%Y.%m.%d").item()
        # ax['map_tv'].text(0.7,0.2,'{}'.format(date),transform=ax['map_tv'].transAxes, fontsize=12, fontweight='bold', color='black')
        # fig.suptitle('Perturbation minus Baseline ({})'.format(date),
        #             fontsize=11, fontweight='bold')

###################################################################
##        Vertical integral differences (og minus twinvar)       ##  
################################################################### 


        # plot values
        cs = ax['map_dif'].pcolormesh(plon, plat, vertint_diff_og-vertint_diff_twinvar,vmin=vmin, vmax=vmax, cmap=newcmap)
        # format figure
        ax['map_dif'].set_yticklabels([])
        ax['map_dif'].set_xticklabels([])
        ax['map_dif'].scatter(inj_lon,inj_lat,s=80, facecolors='none', edgecolors='deeppink')
        pfun.dar(ax['map_dif'])
        # ax['map_tv'].set_title(vn + ' vertical integral difference (mmol/m2)', fontsize=12)
        date = ds_twinvar["ocean_time"].dt.strftime("%Y.%m.%d").item()
        # ax['map_tv'].text(0.7,0.2,'{}'.format(date),transform=ax['map_tv'].transAxes, fontsize=12, fontweight='bold', color='black')
        # fig.suptitle('Perturbation minus Baseline ({})'.format(date),
        #             fontsize=11, fontweight='bold')

        map_right = ax['map_dif'].get_position().x1
        plot_top = ax['map_og'].get_position().y1
        plot_bottom = ax['map_dif'].get_position().y0
        cax = fig.add_axes([map_right + 0.015, plot_bottom, 0.012, plot_top - plot_bottom])
        cbar = fig.colorbar(cs, cax=cax, orientation='vertical')
        cbar.ax.tick_params(labelsize=12)
        cbar.outline.set_visible(False)

###################################################################
##                       Histograms of diffs                     ##  
###################################################################

        if vn == 'TIC':
            his_min = -0.2
            his_max =  0.2
        elif vn == 'alkalinity':
            his_min = -20
            his_max =  20

        diff_og = np.asarray(vertint_diff_og)
        diff_tv = np.asarray(vertint_diff_twinvar)
        diff_dif = np.asarray(vertint_diff_og - vertint_diff_twinvar)
        hist_bins = np.linspace(his_min, his_max, 81)

        texts = ['Original Model','Twin Variables','Difference (Original - Twin)']

        hist_values = [
            (ax['his_og'], diff_og),
            (ax['his_tv'], diff_tv),
            (ax['his_dif'], diff_dif),
        ]
        for his_ax, values in hist_values:
            finite_values = values[np.isfinite(values)]
            his_ax.hist(
                finite_values[finite_values < 0],
                bins=hist_bins,
                color='crimson',
                alpha=0.5,
                edgecolor='none'
            )
            his_ax.hist(
                finite_values[finite_values >= 0],
                bins=hist_bins,
                color='royalblue',
                alpha=0.5,
                edgecolor='none'
            )
            his_ax.set_xlim(his_min, his_max)
            his_ax.set_yscale('log')
            his_ax.grid(True, color='silver', linestyle=':')
            his_ax.set_ylabel('Count', fontsize=12)
            his_ax.set_xlabel('Vertical Integral Difference (mmol/m2)', fontsize=12)
            his_ax.text(0.05, 0.95, texts[hist_values.index((his_ax, values))],
                        transform=his_ax.transAxes,
                        fontsize=14,
                        fontweight='bold',
                        color='black',
                        ha='left',
                        va='top')


###################################################################
##                  Make movie & format figure                   ##  
###################################################################
        
        # prepare a directory for results
        nouts = ('0000' + str(i))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = out_dir /('delta_' + vn)
        Lfun.make_dir(outfile)
        save_name = outfile / outname
        print('Plotting ' + str(fn_pert))
        sys.stdout.flush()
        plt.savefig(save_name)
        plt.close()
        # plt.show()

    # make movie
    if len(fn_list_twinvar) > 1:
        cmd_list = ['ffmpeg','-r','8','-i', str(out_dir / ('delta_' + vn) )+'/plot_%04d.png', '-vcodec', 'libx264',
            '-pix_fmt', 'yuv420p', '-crf', '25', str(out_dir / ('delta_' + vn))+'/movie.mp4']
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print('\n'+stdout.decode())
        if len(stderr) > 0:
            print('\n'+stderr.decode())

plt.close('all')
    