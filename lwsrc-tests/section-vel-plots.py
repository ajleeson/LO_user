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

# USER OPTIONS

fn = 'his_V-less_LtracerSrcFT.nc' 
foldername = 'section_LtracerSrcFT_lwsrc_V-less' 
vn = 'omega' # options: temp, u, v, w, salt, omega

#---------------------------------------------------------

# START
ds = xr.open_dataset('../../LO_roms/lwsrc-test-results/results/'+fn)
in_dict = dict()
in_dict['fn'] = fn
# print(list(ds.keys()))
# print(ds)

# where to save files
outdir0 = Ldir['LOo'] / 'plots'
outdir = outdir0 / foldername
Lfun.make_dir(outdir, clean=True)

# print(ds.ocean_time)

fs = 14
hgt = 10

# PLOT CODE
for t in range(len(ds.ocean_time)):
    print(t)
    pfun.start_plot(fs=fs, figsize=(21,6))
    fig = plt.figure()

    h = ds['h'].values[:,21]

    if vn == 'w':
        x = ds['xi_rho'].values
        y = ds['eta_rho'].values
        horiz = x
        v_surf = ds[vn][t,-1,:,:].values
        v_sect = ds[vn][t,:,20,:]
        z_norm = ds['s_w'].values
        cmap = cm.balance
        vmin = -7e-4
        vmax =  7e-4
        units = 'm/s'
    if vn == 'omega':
        x = ds['xi_rho'].values
        y = ds['eta_rho'].values
        horiz = x
        v_surf = ds[vn][t,-1,:,:].values
        v_sect = ds[vn][t,:,20,:]
        z_norm = ds['s_w'].values
        cmap = cm.balance
        vmin = -1e-6
        vmax =  1e-6
        units = 'm3/s'
    elif vn == 'temp':
        x = ds['xi_rho'].values
        y = ds['eta_rho'].values
        horiz = x
        v_surf = ds[vn][t,-1,:,:].values
        v_sect = ds[vn][t,:,20,:]
        z_norm = ds['s_rho'].values
        cmap = cm.thermal
        vmin = 9
        vmax = 15
        units = 'C'
    elif vn == 'salt':
        x = ds['xi_rho'].values
        y = ds['eta_rho'].values
        horiz = x
        v_surf = ds[vn][t,-1,:,:].values
        v_sect = ds[vn][t,:,20,:]
        z_norm = ds['s_rho'].values
        cmap = cm.haline
        vmin = 34
        vmax = 35
        units = 'g/kg'
    elif vn == 'u':
        x = ds['xi_u'].values
        y = ds['eta_u'].values
        horiz = x
        v_surf = ds[vn][t,-1,:,:].values
        v_sect = ds[vn][t,:,20,:]
        z_norm = ds['s_rho'].values
        h = ds['h'].values[20,1::] # offset u grid
        cmap = cm.balance
        vmin = -2e-3
        vmax =  2e-3
        units = 'm/s'
    elif vn == 'v':
        x = ds['xi_v'].values
        y = ds['eta_v'].values
        horiz = y
        v_surf = ds[vn][t,-1,:,:].values
        v_sect = ds[vn][t,:,:,20]
        z_norm = ds['s_rho'].values
        h = ds['h'].values[1::,20] # offset v grid
        cmap = cm.balance
        vmin = -2e-3
        vmax =  2e-3
        units = 'm/s'

    # add surface profile
    ax = fig.add_subplot(1, 3, 1)
    cs = ax.pcolormesh(x*(5/39), y*(5/39), v_surf, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.locator_params(axis='x', nbins=3)
    pfun.dar(ax)
    ax.set_title('Surface {} [{}]'.format(vn,units))
    ax.set_xlabel('E-W Distance (km)')
    ax.set_ylabel('N-S Distance (km)')
    # add timestep
    ax.text(0.15,0.15,ds.ocean_time[t].values,fontsize=12)
    if vn == 'v':
        ax.axvline(np.max(x)/2*(5/39),0,5,linestyle=':', color='k')
        ax.scatter(np.max(x)/2*(5/39),0,marker='^',facecolor='k',edgecolors='k',s=120)
    else:
        ax.axhline(np.max(y)/2*(5/39),0,5,linestyle=':', color='k')
        ax.scatter(0,np.max(y)/2*(5/39),marker='>',facecolor='k',edgecolors='k',s=120)

    # add section profile
    ax = fig.add_subplot(1, 3, (2,3))
    z = [[sigma*depth for depth in h] for sigma in z_norm]
    max = np.max(v_sect)
    min = np.min(v_sect)
    # print(np.shape(y))
    # print(np.shape(z))
    # print(np.shape(v))
    cs = ax.pcolormesh(horiz*(5/39), z, v_sect, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title('Section {} [{}];  min = {:0.1e} , max = {:0.1e}'.format(vn,units,min,max))
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('Distance (km)')
    fig.colorbar(cs)
    ax.text(0.15,-19,fn,fontsize='12')

    # FINISH
    ds.close()
    # pfun.end_plot()
    nouts = ('0000' + str(t))[-4:]
    outname = 'plot_' + nouts + '.png'
    outfile = outdir / outname
    # print('Plotting ' + str(fn_WWTP))
    sys.stdout.flush()
    plt.savefig(outfile)
    # plt.show()
    plt.close()

# make movie
cmd_list = ['ffmpeg','-r','6','-i', str(outdir)+'/plot_%04d.png', '-vcodec', 'libx264',
    '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir)+'/movie.mp4']
proc = Po(cmd_list, stdout=Pi, stderr=Pi)
stdout, stderr = proc.communicate()
if len(stdout) > 0:
    print('\n'+stdout.decode())
if len(stderr) > 0:
    print('\n'+stderr.decode())