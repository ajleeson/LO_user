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
vn_list = ['u','v','w','temp','zeta']

for t in range(len(ds.ocean_time)):
    print(t)
    pfun.start_plot(fs=fs, figsize=(16,8))
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
        cs = ax.pcolormesh(x, y, v, cmap=cmap, vmin=vmin, vmax=vmax)
        cbar = plt.colorbar(cs,ax=ax, location='bottom')
        cbar.ax.tick_params(rotation=30)
        plt.locator_params(axis='x', nbins=3)
        pfun.dar(ax)
        ax.set_xlabel('E-W Distance (km)')

        ax.set_title('Surface {} ({})'.format(vn,units))
        if ii == 1:
            ax.set_ylabel('N-S Distance (km)')
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
