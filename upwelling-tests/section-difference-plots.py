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

fn = 'results/roms_his_botwwtp10.nc' #'results/roms_his_og.nc'
foldername = 'section_difference_botwwtp10'
vn = 'omega' # options: temp, u, v, w, omega

#---------------------------------------------------------

# get basecase file for differencing
ds0 = xr.open_dataset('results/roms_his_base.nc')

# START
ds = xr.open_dataset(fn)
in_dict = dict()
in_dict['fn'] = fn
# print(list(ds.keys()))
# print(ds)

# where to save files
outdir0 = Ldir['LOo'] / 'plots'
outdir = outdir0 / foldername
Lfun.make_dir(outdir, clean=True)

fs = 14
hgt = 10

# PLOT CODE
for t in range(len(ds.ocean_time)):
    print(t)
    pfun.start_plot(fs=fs, figsize=(15,8))
    fig = plt.figure()

    h = ds['h'].values[:,21]

    if vn == 'w':
        x = ds['xi_rho'].values
        y = ds['eta_rho'].values
        v = ds[vn][t,-1,:,:].values - ds0[vn][t,-1,:,:].values
        z_norm = ds['s_w'].values
        cmap = cm.balance
        vmin = -1e-5
        vmax =  1e-5
        units = 'm/s'
    elif vn == 'omega':
        x = ds['xi_rho'].values
        y = ds['eta_rho'].values
        v = ds[vn][t,-1,:,:].values - ds0[vn][t,-1,:,:].values
        z_norm = ds['s_w'].values
        cmap = cm.balance
        vmin = -3e-6
        vmax =  3e-6
        units = 'm3/s'
    elif vn == 'temp':
        x = ds['xi_rho'].values
        y = ds['eta_rho'].values
        v = ds[vn][t,-1,:,:].values - ds0[vn][t,-1,:,:].values
        z_norm = ds['s_rho'].values
        cmap = cm.balance
        vmin = -4e-3
        vmax =  4e-3
        units = 'C'
    elif vn == 'u':
        x = ds['xi_u'].values
        y = ds['eta_u'].values
        v = ds[vn][t,-1,:,:].values - ds0[vn][t,-1,:,:].values
        z_norm = ds['s_rho'].values
        cmap = cm.balance
        vmin = -3e-4
        vmax =  3e-4
        units = 'm/s'
    elif vn == 'v':
        x = ds['xi_v'].values
        y = ds['eta_v'].values
        v = ds[vn][t,-1,:,:].values - ds0[vn][t,-1,:,:].values
        z_norm = ds['s_rho'].values
        h = ds['h'].values[1::,21] # offset v grid
        cmap = cm.balance
        vmin = -3e-4
        vmax =  3e-4
        units = 'm/s'

    # add surface profile
    ax = fig.add_subplot(1, 3, 1)
    cs = ax.pcolormesh(x, y, v, cmap=cmap, vmin=vmin, vmax=vmax)
    plt.locator_params(axis='x', nbins=3)
    pfun.dar(ax)
    ax.set_title(r'Surface $\Delta$' + '{} ({})'.format(vn,units))
    ax.set_xlabel('E-W Distance (km)')
    ax.set_ylabel('N-S Distance (km)')
    ax.axvline(21,0,85,linestyle=':', color='k')
    ax.scatter(21,0,marker='^',facecolor='k',edgecolors='k',s=120)

    # add section profile
    ax = fig.add_subplot(1, 3, (2,3))
    z = [[sigma*depth for depth in h] for sigma in z_norm]
    v = ds[vn][t,:,:,21] - ds0[vn][t,:,:,21]
    max = np.max(v)
    min = np.min(v)
    # print(np.shape(y))
    # print(np.shape(z))
    # print(np.shape(v))
    cs = ax.pcolormesh(y, z, v, cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_title(r'Section $\Delta$' + '{} [{}]; min = {:0.1e} , max = {:0.1e}'.format(vn,units,min,max))
    ax.set_ylabel('Depth (m)')
    ax.set_xlabel('N-S Distance (km)')
    fig.colorbar(cs)


    # if ii == 1:
    #     pass
        # pfun.add_info(ax, in_dict['fn']) I commented this out so it is easier to see the point sources. Add back in later. --------------
        #pfun.add_windstress_flower(ax, ds)
        # pfun.add_bathy_contours(ax, ds, txt=False)

        # # plot wwtps if they exist
        # do_wwtp = False
        # wwtp_fn = Gr['wwtp_dir'] / 'wwtp_loc_info.csv'
        # # read wwtp lat lon info
        # if wwtp_fn.is_file():
        #     do_wwtp = True
        #     wwtp_df = pd.read_csv(wwtp_fn)
        #     # print(wwtp_df)
        # if do_wwtp:
        #     # plot wwtp locations on grid
        #     ax.scatter(wwtp_df['lon'],wwtp_df['lat'], color='black', label='wwtps')
        #     # print labels
        #     for i,wwtp in enumerate(wwtp_df['dname']):
        #         wwtp_lon = wwtp_df['lon'][i]
        #         wwtp_lat = wwtp_df['lat'][i]+0.05
        #         ax.text(wwtp_lon, wwtp_lat, wwtp, fontsize=14, horizontalalignment='center')

        # # plot point sources linked to the wwtp if the point sources have been created
        # do_ps = False
        # ps_fn = Ldir['data']/ 'grids'/ Gr['gridname'] / 'wwtp_info.csv'
        # # read point source location data
        # if ps_fn.is_file():
        #     do_ps = True
        #     ps_df = pd.read_csv(ps_fn)
        # if do_ps:
        #     # plot point source locations on grid
        #     lon = ds.lon_rho.values
        #     lat = ds.lat_rho.values
        #     X = lon[0,:]
        #     Y = lat[:,0]
        #     ps_lon = [X[int(ind)] for ind in ps_df['col_py']]
        #     ps_lat = [Y[int(ind)] for ind in ps_df['row_py']]
        #     ax.scatter(ps_lon,ps_lat, color='deeppink', marker='x', s=40, label='point sources')
        #     for i,ps in enumerate(ps_df['wname']):
        #         ax.plot([wwtp_df['lon'][i], ps_lon[i]],
        #         [wwtp_df['lat'][i], ps_lat[i]],
        #         color='deeppink', linewidth=1)
        #         ax.legend(loc='best',fontsize=14)

    # elif ii == 2:
    # ax.set_yticklabels([])
    # ii += 1

    plt.suptitle('[With WWTP] minus [No WWTP]')

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