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

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun


# load mooring extraction data
grid_nc = '../../LO_output/extract/cas6_TRAPS2_uu0mb/moor/oakharborlagoon_wwtp/wwtp_2021.02.19_2021.02.19.nc'
ds = xr.open_dataset(grid_nc)

# get all velocities
u = ds['u'].transpose()
v = ds['v'].transpose()
w = ds['w'].transpose()

# get depth values
z_rho = ds['z_rho'].transpose() # depth of u and v-velocities
z_w   = ds['z_w'].transpose()   # depth of w-velocities

# get other parameters
no3 = ds['NO3'].transpose() / 71.4 # (convert mmol/m3 to mg/L)
nh4 = ds['NH4'].transpose() / 71.4 # (convert mmol/m3 to mg/L)
temp = ds['temp'].transpose()
salt = ds['salt'].transpose()
oxygen = ds['oxygen'].transpose() / 31.26 # (convert mmol/m3 to mg/L)
alk = ds['alkalinity'].transpose()

# create time variable
t = np.linspace(1,20,20)
print(z_w)

# figure settings
fs = 14
ls = 12
ts = 16

# Plot Velocities ---------------------------------------------------------------------------------
fig, ax = plt.subplots(3,1,figsize = (10,8), sharex = True)
# E-W
cs = ax[0].pcolormesh(t, z_rho, u, vmin=-10, vmax=10, cmap=plt.get_cmap(cmocean.cm.balance))
plt.colorbar(cs,ax=ax[0])
ax[0].set_ylabel('Z (m)', fontsize = fs)
ax[0].set_title('u (m/s)', fontsize = fs)
ax[0].tick_params(axis='both', which='major', labelsize=ls)
# N-S
cs = ax[1].pcolormesh(t, z_rho, v, vmin=-10, vmax=10, cmap=plt.get_cmap(cmocean.cm.balance))
plt.colorbar(cs,ax=ax[1])
ax[1].set_ylabel('Z (m)', fontsize = fs)
ax[1].set_title('v (m/s)', fontsize = fs)
ax[1].tick_params(axis='both', which='major', labelsize=ls)
# vertical
cs = ax[2].pcolormesh(t, z_w, w, vmin=-0.1, vmax=0.1, cmap=plt.get_cmap(cmocean.cm.balance))
plt.colorbar(cs,ax=ax[2])
ax[2].set_ylabel('Z (m)', fontsize = fs)
ax[2].set_title('w (m/s)', fontsize = fs)
ax[2].set_xlabel('Hour', fontsize = fs)
ax[2].tick_params(axis='both', which='major', labelsize=ls)
ax[2].xaxis.set_ticks(np.arange(5, 21, 5))
# title
fig.suptitle('Moor Extraction Velocities at Oak Harbor Lagoon WWTP', fontsize = ts)
plt.show()

# Other variables -------------------------------------------------------------------------------
fig, ax = plt.subplots(3,2,figsize = (15,8), sharex = True, sharey = True)
fig.suptitle('Moor Extraction State Variables at Oak Harbor Lagoon WWTP', fontsize = ts)

# Plot Salinity ---------------------------------------------------------------------------------
cs = ax[0,0].pcolormesh(t, z_rho, salt, cmap=plt.get_cmap(cmocean.cm.haline))
plt.colorbar(cs,ax=ax[0,0])
ax[0,0].set_ylabel('Z (m)', fontsize = fs)
ax[0,0].set_title('Salinity (g/kg)', fontsize = fs)
ax[0,0].tick_params(axis='both', which='major', labelsize=ls)

# Plot Temperature ---------------------------------------------------------------------------------
cs = ax[0,1].pcolormesh(t, z_rho, temp, cmap=plt.get_cmap(cmocean.cm.thermal))
plt.colorbar(cs,ax=ax[0,1])
ax[0,1].set_title('Temperature (C)', fontsize = fs)
ax[0,1].tick_params(axis='both', which='major', labelsize=ls)

# Plot Nitrate ---------------------------------------------------------------------------------
cs = ax[1,0].pcolormesh(t, z_rho, no3, cmap=plt.get_cmap(cmocean.cm.deep))
plt.colorbar(cs,ax=ax[1,0])
ax[1,0].set_ylabel('Z (m)', fontsize = fs)
ax[1,0].set_title('Nitrate (mg/L)', fontsize = fs)
ax[1,0].tick_params(axis='both', which='major', labelsize=ls)

# Plot Oxygen ---------------------------------------------------------------------------------
cs = ax[1,1].pcolormesh(t, z_rho, oxygen, cmap=plt.get_cmap(cmocean.cm.deep))
plt.colorbar(cs,ax=ax[1,1])
ax[1,1].set_title('Oxygen (mg/L)', fontsize = fs)
ax[1,1].tick_params(axis='both', which='major', labelsize=ls)

# Plot Ammonium ---------------------------------------------------------------------------------
cs = ax[2,0].pcolormesh(t, z_rho, nh4, cmap=plt.get_cmap(cmocean.cm.deep))
plt.colorbar(cs,ax=ax[2,0])
ax[2,0].set_ylabel('Z (m)', fontsize = fs)
ax[2,0].set_title('Ammonium (mg/L)', fontsize = fs)
ax[2,0].set_xlabel('Hour', fontsize = fs)
ax[2,0].tick_params(axis='both', which='major', labelsize=ls)
ax[2,0].xaxis.set_ticks(np.arange(5, 21, 5))

# Plot Alkalinity ---------------------------------------------------------------------------------
cs = ax[2,1].pcolormesh(t, z_rho, alk, cmap=plt.get_cmap(cmocean.cm.deep))
plt.colorbar(cs,ax=ax[2,1])
ax[2,1].set_title('Alkalinity (mEq/m3)', fontsize = fs)
ax[2,1].set_xlabel('Hour', fontsize = fs)
ax[2,1].tick_params(axis='both', which='major', labelsize=ls)
ax[2,1].xaxis.set_ticks(np.arange(5, 21, 5))
plt.show()

# print what the sigma layers look like (so we can see whether they are too thin)
fig, ax = plt.subplots(1,1,figsize = (8,6))
print(z_w.shape)
for i in range(0,31):
    plt.plot(t,z_w[i],color='k')
plt.title('Sigma Layers at Oak Harbor Lagoon WWTP',fontsize=ts)
plt.xlabel('Hour',fontsize=fs)
plt.xlim(1,20)
plt.ylim(np.min(z_w),)
plt.ylabel('Depth of Sigma Layer Boundary (m)', fontsize=fs)
ax.tick_params(axis='both', which='major', labelsize=ls)
plt.xticks(np.arange(5, 21, 5))
plt.show()

print(list(ds.keys()))