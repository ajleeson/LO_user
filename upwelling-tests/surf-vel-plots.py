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
from matplotlib import pyplot as plt
from subprocess import Popen as Po
from subprocess import PIPE as Pi


from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
# import pinfo
from importlib import reload
reload(pfun)
# reload(pinfo)

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

Ldir = Lfun.Lstart()

#---------------------------------------------------------


fn = 'results/roms_his_base.nc'
foldername = 'surface_base_upwelling'

# START
ds = xr.open_dataset(fn)
# print(list(ds.keys()))

# where to save files
outdir0 = Ldir['LOo'] / 'plots'
outdir = outdir0 / foldername
Lfun.make_dir(outdir, clean=True)

fs = 14
hgt = 10
# PLOT CODE
vn_list = ['u','v','w','omega','zeta','temp']

for t in range(len(ds.ocean_time)):
    print(t)
    pfun.start_plot(fs=fs, figsize=(18,8))
    fig = plt.figure()
    for i,vn in enumerate(vn_list):
        ii = i + 1
        ax = fig.add_subplot(1, len(vn_list), ii)
        if vn == 'w':
            x = ds['xi_rho'].values
            y = ds['eta_rho'].values
            v = ds[vn][t,-1,:,:].values
            cmap = cm.balance
            vmin = -1e-4
            vmax = 1e-4
            units = 'm/s'
        if vn == 'omega':
            x = ds['xi_rho'].values
            y = ds['eta_rho'].values
            v = ds[vn][t,-1,:,:].values
            cmap = cm.balance
            vmin = -1e-5
            vmax =  1e-5
            units = 'm3/s'
        elif vn == 'temp':
            x = ds['xi_rho'].values
            y = ds['eta_rho'].values
            v = ds[vn][t,-1,:,:].values
            cmap = cm.thermal
            vmin = 15
            vmax = 22
            units = 'C'
        elif vn == 'u':
            x = ds['xi_u'].values
            y = ds['eta_u'].values
            v = ds[vn][t,-1,:,:].values
            cmap = cm.balance
            vmin = -1
            vmax = 1
            units = 'm/s'
        elif vn == 'v':
            x = ds['xi_v'].values
            y = ds['eta_v'].values
            v = ds[vn][t,-1,:,:].values
            cmap = cm.balance
            vmin = -0.1
            vmax = 0.1
            units = 'm/s'
        elif vn == 'zeta':
            x = ds['xi_rho'].values
            y = ds['eta_rho'].values
            v = ds[vn][t,:,:].values
            cmap = cm.balance
            vmin = -0.5
            vmax = 0.5
            units = 'm'
        max = np.max(v)
        min = np.min(v)
        cs = ax.pcolormesh(x, y, v, cmap=cmap, vmin=vmin, vmax=vmax)
        cbar = plt.colorbar(cs,ax=ax, location='bottom')
        cbar.ax.tick_params(rotation=30)
        plt.locator_params(axis='x', nbins=3)
        pfun.dar(ax)
        ax.set_xlabel('E-W Distance (km)')

        ax.set_title('{} [{}] \n min = {:0.1e} \n max = {:0.1e}'.format(vn,units,min,max), fontsize=12)
        if ii == 1:
            ax.set_ylabel('N-S Distance (km)')

        else:
            ax.set_yticklabels([])
        # ii += 1

    # FINISH
    ds.close()
    # pfun.end_plot()
    nouts = ('0000' + str(t))[-4:]
    outname = 'plot_' + nouts + '.png'
    outfile = outdir / outname
    # print('Plotting ' + str(fn_WWTP))
    sys.stdout.flush()
    plt.savefig(outfile)
    plt.close()

# make movie
# if len(fn_list_WWTP) > 1:
cmd_list = ['ffmpeg','-r','6','-i', str(outdir)+'/plot_%04d.png', '-vcodec', 'libx264',
    '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir)+'/movie.mp4']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if len(stdout) > 0:
    print('\n'+stdout.decode())
if len(stderr) > 0:
    print('\n'+stderr.decode())

plt.close('all')
