"""
Code to plot extract moor information from cas6_traps00_uu0mb,
the first LiveOcean run with TRAPS which blew up due to Oak Harbor Lagoon WWTP
"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import matplotlib as mpl
mpl.use('QtAgg')
import gsw

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun


# load mooring extraction data
grid_nc = '../../LO_output/extract/alpe2_itraps0_uu1k/moor/wwtp4/wwtp_2020.01.01_2020.01.01.nc'
ds = xr.open_dataset(grid_nc)

# get all velocities
u = ds['u'].transpose()
v = ds['v'].transpose()
w = ds['w'].transpose()

# get depth values
z_rho = ds['z_rho'].transpose() # depth of u and v-velocities
z_w   = ds['z_w'].transpose()   # depth of w-velocities
z_min = np.min(z_w.values)

# get other parameters
temp = ds['temp'].transpose() # potential temperature
salt = ds['salt'].transpose() # practical salinity
AKv = ds['AKv'].transpose()
AKs = ds['AKs'].transpose()

# Calculate density
# get pressure (dbar)
p = [gsw.p_from_z(depth,48.28559664274097) for depth in z_rho.values]
# calculate density
rho = np.ndarray((30,17))
for row in range(30):
    for col in range(17):
        # calculate absolute salinity from practical salinity
        salt_abs = gsw.conversions.SA_from_SP(salt[row][col], p[row][col], -122.60105555932371, 48.28559664274097)
        # calculate conservative temperature from potential temperature
        temp_cons = gsw.conversions.CT_from_pt(salt_abs, temp[row][col])
        # calculate density
        rho[row][col] = gsw.rho(salt_abs,temp_cons,p[row][col])

# create time variable
t = np.linspace(1,24,24)

# figure settings
fs = 14
ls = 12
ts = 16

plotting = True
if plotting:

    # Plot Velocities ---------------------------------------------------------------------------------
    fig, ax = plt.subplots(3,1,figsize = (10,8), sharex = True)
    # E-W
    cs = ax[0].pcolormesh(t, z_rho, u, vmin=-0.1, vmax=0.1, cmap=plt.get_cmap(cmocean.cm.balance))
    plt.colorbar(cs,ax=ax[0])
    ax[0].set_ylabel('Z (m)', fontsize = fs)
    ax[0].set_title('u (m/s)', fontsize = fs)
    ax[0].tick_params(axis='both', which='major', labelsize=ls)
    ax[0].set_ylim((z_min,6))
    # N-S
    cs = ax[1].pcolormesh(t, z_rho, v, vmin=-0.3, vmax=0.3, cmap=plt.get_cmap(cmocean.cm.balance))
    plt.colorbar(cs,ax=ax[1])
    ax[1].set_ylabel('Z (m)', fontsize = fs)
    ax[1].set_title('v (m/s)', fontsize = fs)
    ax[1].tick_params(axis='both', which='major', labelsize=ls)
    ax[1].set_ylim((z_min,6))
    # vertical
    cs = ax[2].pcolormesh(t, z_w, w, vmin=-0.0003, vmax=0.0003, cmap=plt.get_cmap(cmocean.cm.balance))
    plt.colorbar(cs,ax=ax[2])
    ax[2].set_ylabel('Z (m)', fontsize = fs)
    ax[2].set_title('w (m/s)', fontsize = fs)
    ax[2].set_xlabel('Hour', fontsize = fs)
    ax[2].tick_params(axis='both', which='major', labelsize=ls)
    ax[2].xaxis.set_ticks(np.arange(4, 25, 5))
    ax[2].set_ylim((z_min,6))
    # title
    fig.suptitle('Moor Extraction Velocities at WWTP4', fontsize = ts)
    plt.show()


# print(list(ds.keys()))
