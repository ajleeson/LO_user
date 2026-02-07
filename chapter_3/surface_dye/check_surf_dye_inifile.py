'''
Check history file created by add_dye_to_hisfile
to verify that I added surface dye correctly

The vertical integral should sum up to 1 kg/m3 * 5 m = 5 kg/m2

'''

import sys
import shutil
import argparse
from datetime import datetime, timedelta
from time import time
import xarray as xr
import numpy as np
import matplotlib.pylab as plt
import csv
import sys
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

###############################################################
# open ini file

fn = Ldir['LOo'] / 'forcing' / 'cas7' / 'f2020.01.02' / 'ocnG00d' / 'ocean_ini.nc'
ds = xr.open_dataset(fn)

###############################################################
# calculate vertical integral of dye
# we are aiming for 5 kg/m2 at every grid cell

# get dz
G, S, T = zrfun.get_basic_info(Ldir['roms_out'] / 'cas7_t1d_x11ad' / 'f2020.01.01' / 'ocean_his_0002.nc')
zr, zw = zrfun.get_z(G['h'],ds.zeta[0,:,:].to_numpy(),S)
dz = np.diff(zw, axis=0) # [m]

# zr_wrong, zw_wrong = zrfun.get_z(G['h'],0*G['h'],S)
# dz_wrong = np.diff(zw_wrong, axis=0) # [m]
# print(zw[:,600,353])
# print('\n')
# print(zw_wrong[:,600,353])

# get dye concentrations
dye_conc = ds['dye_01'][0,:,:,:].to_numpy() # [kg/m3]

# vertical integrals
dye_vert_int = (dye_conc * dz).sum(axis=0) # [kg/m2]

###############################################################
# plot

pfun.start_plot(fs=12, figsize=(11,7))
fig,axes = plt.subplots(1,2, sharex=True, sharey=True)
ax = axes.ravel()

# plot surface dye concentrations
x = ds['lon_rho'].values
y = ds['lat_rho'].values
px, py = pfun.get_plon_plat(x,y)
cs = ax[0].pcolormesh(px, py, ds['dye_01'][0,-1,:,:], cmap='rainbow')
# add colorbar
cbar = plt.colorbar(cs,ax=ax[0],location='left')
# format figure
pfun.dar(ax[0])
ax[0].set_title('Surface dye concentration\n' +  r'[kg m$^{-3}$]')

# plot vertical integral amount of dye
x = ds['lon_rho'].values
y = ds['lat_rho'].values
px, py = pfun.get_plon_plat(x,y)
cs = ax[1].pcolormesh(px, py, dye_vert_int, cmap='coolwarm',
                      vmin=4.9999, vmax=5.0001)
# add colorbar
cbar = plt.colorbar(cs,ax=ax[1],location='right')
# format figure
pfun.dar(ax[1])
ax[1].set_title('Vertically-integrated dye\n' +  r'[kg m$^{-2}$]')

plt.show()
