'''
Calculate hydrodynamic differencing plots for:
u_WWTP - u_noWWTP
'''

from matplotlib import markers
import numpy as np
import xarray as xr
import pickle
from datetime import datetime, timedelta
import pandas as pd
from cmocean import cm
import matplotlib.dates as mdates
import argparse
import math
import scipy.interpolate as interp
from scipy.optimize import curve_fit
from matplotlib.patches import Rectangle
from subprocess import Popen as Po
from subprocess import PIPE as Pi

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
import pinfo
from importlib import reload
reload(pfun)
reload(pinfo)

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu

pth = Path(__file__).absolute().parent.parent.parent / 'LO_user' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

import matplotlib as mpl
# mpl.use('QtAgg')
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.close('all')

# USER OPTIONS ----------------------------------------------------

d0= '2021.02.01'
d1 = '2021.03.01' #'2021.03.01'

list_type = 'hourly' #'snapshot', 'daily', 'hourly ', 'allhours'
his_num = 16

region = 'King County'

# ----------------------------------------------------------------

# gtagex of files to difference
Ldir_WWTP   = Lfun.Lstart(gridname='cas6', tag='trapsT00', ex_name='uu0mb')
Ldir_noWWTP = Lfun.Lstart(gridname='cas6', tag='trapsT01', ex_name='uu0mb')

# get list of history files to plot
fn_list_WWTP   = Lfun.get_fn_list(list_type, Ldir_WWTP,
    d0, d1, his_num)
fn_list_noWWTP = Lfun.get_fn_list(list_type, Ldir_noWWTP,
    d0, d1, his_num)

# Read in WWTP locations
wwtp_loc = pd.read_csv('../../LO_data/grids/cas6/wwtp_info.csv')

# Variables for plotting
vn_list = ['u', 'v', 'w', 'zeta']

# get plotting limits based on region
if region == 'King County':
    xmin = -122.55
    xmax = -122.3
    ymin = 47.58
    ymax = 47.73

# PLOTTING
outdir0 = Ldir['LOo'] / 'plots'
Lfun.make_dir(outdir0)

if len(fn_list_WWTP) > 1:
    # prepare a directory for results if making a movie
    outdir = outdir0 / (region + '_' + list_type+ '_withWWTPs_minus_noWWTP')
    Lfun.make_dir(outdir, clean=True)
#-------------------------------------------------------------------------------------

for i,fn_WWTP in enumerate(fn_list_WWTP):
    fs = 10
    # plt.tight_layout()
    fig = plt.figure(figsize=(16,8))
    fn_noWWTP = fn_list_noWWTP[i]
    ds_WWTP = xr.open_dataset(fn_WWTP)
    ds_noWWTP = xr.open_dataset(fn_noWWTP)

    # add map of Puget Sound
    ax = fig.add_subplot(1, 3, 1)
    pfun.add_coast(ax,color='gray')
    ax.set(xlim=(-123.5, -122), ylim=(47, 49))
    ax.tick_params('x',rotation=30)
    pfun.dar(ax)
    # draw region
    ax.add_patch(Rectangle((xmin, ymin), (xmax-xmin), (ymax-ymin),
        edgecolor = 'deeppink', facecolor = 'none', lw=2))
    pfun.add_info(ax, fn_WWTP, his_num=True)

    # Plot variables
    for j,vn in enumerate(vn_list):
        ii = j+1
        if 'lon_rho' in ds_WWTP[vn].coords:
            tag = 'rho'
        if 'lon_u' in ds_WWTP[vn].coords:
            tag = 'u'
        if 'lon_v' in ds_WWTP[vn].coords:
            tag = 'v'
        x = ds_WWTP['lon_'+tag].values
        y = ds_WWTP['lat_'+tag].values
        X = ds_WWTP['lon_rho'].values[0,:]
        Y = ds_WWTP['lat_rho'].values[:,0]
        px, py = pfun.get_plon_plat(x,y)
        if vn in ['u', 'v']:
            v = ds_WWTP[vn][0,-1,:,:].values - ds_noWWTP[vn][0,-1,:,:].values
            vmin = -2
            vmax = 2
            cmap=plt.get_cmap(cm.balance)
            vn = 'surface ' + str(vn) + ' (m/s)'
        elif vn in ['w']:
            v = ds_WWTP[vn][0,-1,:,:].values - ds_noWWTP[vn][0,-1,:,:].values
            vmin = -5e-4
            vmax = 5e-4
            cmap=plt.get_cmap(cm.balance)
            vn = 'surface w (m/s)'
        elif vn == 'zeta':
            v = ds_WWTP[vn][0,:,:].values - ds_noWWTP[vn][0,:,:].values
            mr = ds_WWTP.mask_rho.values
            v[mr==0] = np.nan
            vn = 'SSH anomaly (m)'
            vmin = -0.1
            vmax = 0.1
            cmap=plt.get_cmap(cm.balance)
        else:
            v = ds_WWTP[vn][0, -1,:,:].values
        if ii < 3:
            ax = fig.add_subplot(2, 3, ii+1)
        elif ii >=3:
            ax = fig.add_subplot(2, 3, ii+2)

        cs = ax.pcolormesh(px, py, v, cmap=cmap, vmin=vmin, vmax=vmax)
        pfun.add_coast(ax,color='gray')
        ax.axis(pfun.get_aa(ds_WWTP))
        ax.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        ax.set_xticks([])
        ax.set_yticks([])

        # plot location of wwtps
        lon_wwtps = [X[int(col)] for col in wwtp_loc['col_py']]
        lat_wwtps = [Y[int(row)] for row in wwtp_loc['row_py']]
        ax.scatter(lon_wwtps,lat_wwtps,s=20,marker='o',c='k',label='WWTPs')
            
        # add colorbar
        cbar = plt.colorbar(cs,ax=ax)
        cbar.ax.tick_params(labelsize=11, rotation=30)

        pfun.dar(ax)
        if ii == 1:
            ax.legend(loc='upper right', fontsize = 12)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_title(vn,fontsize=16)
        plt.suptitle('(with WWTPs) minus (no WWTPs)',fontsize=20)
    
# # FINISH
# if len(fn_list_WWTP) > 0:
#     plt.savefig(outdir0 / (list_type + '_'
#             + 'withWWTPs_minus_noWWTPs.png'))
# else:
#     pass

#-------------------------------------------------------

    if len(fn_list_WWTP) == 1:
        plt.savefig(outdir0 / (list_type + '_'
                + 'withWWTPs_minus_noWWTPs.png'))
    # # plot a single image to screen
    # fn = fn_list_WWTP[0]
    # in_dict['fn'] = fn
    # if Ldir['save_plot'] == True:
    #     in_dict['fn_out'] = outdir0 / (Ldir['list_type'] + '_'
    #         + Ldir['plot_type'] + '_' + Ldir['gtagex'] + '.png')
    # else:
    #     in_dict['fn_out'] = ''
    # whichplot(in_dict)
    
    elif len(fn_list_WWTP) > 1:
        # prepare a directory for results
        # outdir = outdir0 / (list_type + region + '_withWWTPs_minus_noWWTP')
        # Lfun.make_dir(outdir, clean=True)
        # plot to a folder of files
        # for fn in fn_list_WWTP:
        nouts = ('0000' + str(i))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir / outname
        print('Plotting ' + str(fn_WWTP))
        sys.stdout.flush()
        plt.savefig(outfile)
            # in_dict['fn'] = fn
            # in_dict['fn_out'] = outfile
            # whichplot(in_dict)
            # after the first plot we no longer change vlims
            # in_dict['auto_vlims'] = False
        # jj += 1
        plt.close()

# make movie
if len(fn_list_WWTP) > 1:
    cmd_list = ['ffmpeg','-r','8','-i', str(outdir)+'/plot_%04d.png', '-vcodec', 'libx264',
        '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir)+'/movie.mp4']
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    if len(stdout) > 0:
        print('\n'+stdout.decode())
    if len(stderr) > 0:
        print('\n'+stderr.decode())

