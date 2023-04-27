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

fn = 'results/roms_his_botwwtp10.nc'
foldername = 'surface_height_botwwtp10'

#---------------------------------------------------------

# get basecase file for differencing
ds0 = xr.open_dataset('results/roms_his_base.nc')

# START
ds = xr.open_dataset(fn)
print(list(ds.keys()))
print(ds['x_psi'])

# where to save files
outdir0 = Ldir['LOo'] / 'plots'
outdir = outdir0 / foldername
Lfun.make_dir(outdir, clean=True)

fs = 11
# PLOT CODE
vn_list = ['zeta']

slopes = np.zeros((80,41))

time_s = np.linspace(0,432000,21)
# need to ignore values around the perimeter because
# rho grid extends one cell beyond domain in each direction
for row in range(1,81):
    for col in range(1,42):
            zeta = ds['zeta'][:,row,col]
            slope, intercept = np.polyfit(time_s, zeta, 1)
            slopes[row-1,col-1] = slope

# dt = 21600 # s
# for i in range(1,21):
#     for row in range(82):
#         for col in range(43):
#                 zeta0 = ds['zeta'][i-1,row,col]
#                 zeta1 = ds['zeta'][i,row,col]
#                 dzeta = zeta1-zeta0
#                 slope = dzeta/dt
#                 slopes[row,col] = slope
#     print(round(np.sum(slopes)*1000*1000,2))

print('SSH transport over entire domain = {} m3/s'.format(round(np.sum(slopes)*1000*1000,5)))
print('Source Flowrate = 10 m3/s')
# print('w transport in water column = {} (m3/s)'.format(round(np.sum(ds['w'][10,:,40,21].values)*1000*1000,2)))


# time_s = np.linspace(0,432000,21)
# zeta = ds['zeta'][:,1,1].values
# plt.figure(figsize = (5,5))
# plt.tight_layout
# plt.plot(time_s,zeta)
# plt.ylabel(r'$\zeta$ (m)')
# plt.xlabel('Time (s)')
# slope, intercept = np.polyfit(time_s, zeta, 1)
# plt.title('Slope = {:0.1e} m/s'.format(slope))
# plt.show()


# fs = 14
# hgt = 10
# for t in range(len(ds.ocean_time)):
#     print(t)
#     pfun.start_plot(fs=fs, figsize=(3,7))
#     fig = plt.figure()
#     for i,vn in enumerate(vn_list):
#         ii = i + 1
#         ax = fig.add_subplot(1, len(vn_list), ii)
#         if vn == 'w':
#             x = ds['xi_rho'].values
#             y = ds['eta_rho'].values
#             v = ds[vn][t,-1,:,:].values
#             v0 = ds0[vn][t,-1,:,:].values
#             vdiff = v - v0
#             cmap = cm.balance
#             vmin = -5e-6
#             vmax =  5e-6
#             units = 'm/s'
#             vn = r'$\Delta w$'
#         elif vn == 'omega':
#             x = ds['xi_rho'].values
#             y = ds['eta_rho'].values
#             v = ds[vn][t,-1,:,:].values
#             v0 = ds0[vn][t,-1,:,:].values
#             vdiff = v - v0
#             cmap = cm.balance
#             vmin = -1e-6
#             vmax =  1e-6
#             units = 'm3/s'
#             vn = r'$\Delta \omega$'
#         elif vn == 'temp':
#             x = ds['xi_rho'].values
#             y = ds['eta_rho'].values
#             v = ds[vn][t,-1,:,:].values
#             v0 = ds0[vn][t,-1,:,:].values
#             vdiff = v - v0
#             cmap = cm.balance
#             vmin = -5e-5
#             vmax =  5e-5
#             units = 'C'
#             vn = r'$\Delta T$'
#         elif vn == 'u':
#             x = ds['xi_u'].values
#             y = ds['eta_u'].values
#             v = ds[vn][t,-1,:,:].values
#             v0 = ds0[vn][t,-1,:,:].values
#             vdiff = v - v0
#             cmap = cm.balance
#             vmin = -5e-5
#             vmax =  5e-5
#             units = 'm/s'
#             vn = r'$\Delta u$'
#         elif vn == 'v':
#             x = ds['xi_v'].values
#             y = ds['eta_v'].values
#             v = ds[vn][t,-1,:,:].values
#             v0 = ds0[vn][t,-1,:,:].values
#             vdiff = v - v0
#             cmap = cm.balance
#             vmin = -5e-5
#             vmax =  5e-5
#             units = 'm/s'
#             vn = r'$\Delta v$'
#         elif vn == 'zeta':
#             x = ds['xi_rho'].values
#             y = ds['eta_rho'].values
#             v = ds[vn][t,:,:].values
#             v0 = ds0[vn][t,:,:].values
#             vdiff = v - v0
#             # subtract mean difference
#             vdiff = vdiff #- np.mean(vdiff)
#             cmap = cm.balance
#             vmin = -2e-3
#             vmax =  2e-3
#             units = 'm'
#             vn = r'$\Delta \zeta$'# - (\Delta \zeta)_{avg}$'
#         # calculate min and max values
#         max = np.max(vdiff)
#         min = np.min(vdiff)
#         cs = ax.pcolormesh(x, y, vdiff, cmap=cmap, vmin=vmin, vmax=vmax)
#         cbar = plt.colorbar(cs,ax=ax, location='bottom')
#         cbar.ax.tick_params(rotation=30, labelsize=12)
#         plt.locator_params(axis='x', nbins=3)
#         pfun.dar(ax)
#         ax.set_xlabel('E-W Distance (km)', fontsize=12)

#         ax.set_title('{} [{}] \n min = {:0.1e} \n max = {:0.1e}'.format(vn,units,min,max), fontsize=12)
#         if ii == 1:
#             ax.set_ylabel('N-S Distance (km)', fontsize=12)

#         else:
#             ax.set_yticklabels([])
#         # ii += 1

#         plt.suptitle('[With WWTP] minus [No WWTP] \n SSH',fontsize=11)
#         plt.subplots_adjust(top=0.8) 

#     # FINISH
#     ds.close()
#     # pfun.end_plot()
#     nouts = ('0000' + str(t))[-4:]
#     outname = 'plot_' + nouts + '.png'
#     outfile = outdir / outname
#     # print('Plotting ' + str(fn_WWTP))
#     sys.stdout.flush()
#     plt.savefig(outfile)
#     plt.close()

# # make movie
# # if len(fn_list_WWTP) > 1:
# cmd_list = ['ffmpeg','-r','6','-i', str(outdir)+'/plot_%04d.png', '-vcodec', 'libx264',
#     '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir)+'/movie.mp4']
# proc = Po(cmd_list, stdout=Pi, stderr=Pi)
# stdout, stderr = proc.communicate()
# if len(stdout) > 0:
#     print('\n'+stdout.decode())
# if len(stderr) > 0:
#     print('\n'+stderr.decode())
