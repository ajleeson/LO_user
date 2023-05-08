"""
Code to plot extract moor information for a 24 hour period.
User needs to specify some inputs at beginning.
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
import gsw

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun


#   USER DEFINED INPUTS 
# lat/lon for calculating density
lat = 48.28559664274097
lon = -122.60105555932371
site_name = 'Oak Harbor Lagoon WWTP'
date = '2021.02.19'
z_max = 6 # vertical limit

# Oak Harbor Lagoon prior LO mixing
grid_nc = '../../LO_output/extract/cas6_TRAPS2_uu1mb/moor/oakharborlagoon_wwtp/wwtp_2021.02.19_2021.02.19.nc'

# Option user defined colormap bounds
u_vmin = -2
u_vmax =  2
v_vmin = -2
v_vmax =  2
w_vmin = -0.1
w_vmax =  0.1
salt_vmin = 23
salt_vmax = 30
temp_vmin = 5.5
temp_vmax = 7.5
no3_vmin = 0.01
no3_vmax = 0.032
oxygen_vmin = 9.8
oxygen_vmax = 11
nh4_vmin = 0
nh4_vmax = 0.07
alk_vmin = 1800
alk_vmax = 2025
rho_vmin = 1018
rho_vmax = 1024
AKv_vmin = 0
AKv_vmax = 1
AKs_vmin = 0
AKs_vmax = 1


# ---------------------------------------------------------------------------

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
no3 = ds['NO3'].transpose() / 71.4 # (convert mmol/m3 to mg/L)
nh4 = ds['NH4'].transpose() / 71.4 # (convert mmol/m3 to mg/L)
temp = ds['temp'].transpose() # potential temperature
salt = ds['salt'].transpose() # practical salinity
oxygen = ds['oxygen'].transpose() / 31.26 # (convert mmol/m3 to mg/L)
alk = ds['alkalinity'].transpose()
AKv = ds['AKv'].transpose()
AKs = ds['AKs'].transpose()

# create time variable
t_end = u.values.shape[1]
t = np.linspace(1,t_end,t_end)

# Calculate density
# get pressure (dbar)
p = [gsw.p_from_z(depth,lat) for depth in z_rho.values]
# calculate density
rho = np.ndarray((30,t_end))
for row in range(30):
    for col in range(t_end):
        # calculate absolute salinity from practical salinity
        salt_abs = gsw.conversions.SA_from_SP(salt[row][col], p[row][col], lon, lat)
        # calculate conservative temperature from potential temperature
        temp_cons = gsw.conversions.CT_from_pt(salt_abs, temp[row][col])
        # calculate density
        rho[row][col] = gsw.rho(salt_abs,temp_cons,p[row][col])

# figure settings
fs = 14
ls = 12
ts = 16

plotting = True
if plotting:

    # Plot Velocities ---------------------------------------------------------------------------------
    fig, ax = plt.subplots(3,1,figsize = (10,8), sharex = True)
    # E-W
    cs = ax[0].pcolormesh(t, z_rho, u, vmin=u_vmin, vmax=u_vmax, cmap=plt.get_cmap(cmocean.cm.balance))
    plt.colorbar(cs,ax=ax[0])
    ax[0].set_ylabel('Z (m)', fontsize = fs)
    ax[0].set_title('u (m/s)', fontsize = fs)
    ax[0].tick_params(axis='both', which='major', labelsize=ls)
    ax[0].set_ylim((z_min,z_max))
    # N-S
    cs = ax[1].pcolormesh(t, z_rho, v, vmin=v_vmin, vmax=v_vmax, cmap=plt.get_cmap(cmocean.cm.balance))
    plt.colorbar(cs,ax=ax[1])
    ax[1].set_ylabel('Z (m)', fontsize = fs)
    ax[1].set_title('v (m/s)', fontsize = fs)
    ax[1].tick_params(axis='both', which='major', labelsize=ls)
    ax[1].set_ylim((z_min,z_max))
    # vertical
    cs = ax[2].pcolormesh(t, z_w, w, vmin=w_vmin, vmax=w_vmax, cmap=plt.get_cmap(cmocean.cm.balance))
    plt.colorbar(cs,ax=ax[2])
    ax[2].set_ylabel('Z (m)', fontsize = fs)
    ax[2].set_title('w (m/s)', fontsize = fs)
    ax[2].set_xlabel('Hour', fontsize = fs)
    ax[2].tick_params(axis='both', which='major', labelsize=ls)
    ax[2].xaxis.set_ticks(np.arange(4, 25, 4))
    ax[2].set_ylim((z_min,z_max))
    # title
    fig.suptitle('Moor Extraction Velocities at {} ({})'.format(site_name,date), fontsize = ts)
    plt.show()

    # Other variables -------------------------------------------------------------------------------
    fig, ax = plt.subplots(3,2,figsize = (15,8), sharex = True, sharey = True)
    fig.suptitle('Moor Extraction State Variables at {} ({})'.format(site_name,date), fontsize = ts)

    # Plot Salinity ---------------------------------------------------------------------------------
    cs = ax[0,0].pcolormesh(t, z_rho, salt, vmin=salt_vmin, vmax=salt_vmax, cmap=plt.get_cmap(cmocean.cm.haline))
    plt.colorbar(cs,ax=ax[0,0])
    ax[0,0].set_ylabel('Z (m)', fontsize = fs)
    ax[0,0].set_title('Salinity (g/kg)', fontsize = fs)
    ax[0,0].tick_params(axis='both', which='major', labelsize=ls)
    ax[0,0].set_ylim((z_min,z_max))

    # Plot Temperature ---------------------------------------------------------------------------------
    cs = ax[0,1].pcolormesh(t, z_rho, temp, vmin=temp_vmin, vmax=temp_vmax, cmap=plt.get_cmap(cmocean.cm.thermal))
    plt.colorbar(cs,ax=ax[0,1])
    ax[0,1].set_title('Temperature (C)', fontsize = fs)
    ax[0,1].tick_params(axis='both', which='major', labelsize=ls)
    ax[0,1].set_ylim((z_min,z_max))

    # Plot Nitrate ---------------------------------------------------------------------------------
    cs = ax[1,0].pcolormesh(t, z_rho, no3, vmin=no3_vmin, vmax=no3_vmax, cmap=plt.get_cmap(cmocean.cm.deep))
    plt.colorbar(cs,ax=ax[1,0])
    ax[1,0].set_ylabel('Z (m)', fontsize = fs)
    ax[1,0].set_title('Nitrate (mg/L)', fontsize = fs)
    ax[1,0].tick_params(axis='both', which='major', labelsize=ls)
    ax[1,0].set_ylim((z_min,z_max))

    # Plot Oxygen ---------------------------------------------------------------------------------
    cs = ax[1,1].pcolormesh(t, z_rho, oxygen, vmin=oxygen_vmin, vmax=oxygen_vmax, cmap=plt.get_cmap(cmocean.cm.deep))
    plt.colorbar(cs,ax=ax[1,1])
    ax[1,1].set_title('Oxygen (mg/L)', fontsize = fs)
    ax[1,1].tick_params(axis='both', which='major', labelsize=ls)
    ax[1,1].set_ylim((z_min,z_max))

    # Plot Ammonium ---------------------------------------------------------------------------------
    cs = ax[2,0].pcolormesh(t, z_rho, nh4, vmin=nh4_vmin, vmax=nh4_vmax, cmap=plt.get_cmap(cmocean.cm.deep))
    plt.colorbar(cs,ax=ax[2,0])
    ax[2,0].set_ylabel('Z (m)', fontsize = fs)
    ax[2,0].set_title('Ammonium (mg/L)', fontsize = fs)
    ax[2,0].set_xlabel('Hour', fontsize = fs)
    ax[2,0].tick_params(axis='both', which='major', labelsize=ls)
    ax[2,0].xaxis.set_ticks(np.arange(4, 25, 4))
    ax[2,0].set_ylim((z_min,z_max))

    # Plot Alkalinity ---------------------------------------------------------------------------------
    cs = ax[2,1].pcolormesh(t, z_rho, alk, vmin=alk_vmin, vmax=alk_vmax, cmap=plt.get_cmap(cmocean.cm.deep))
    plt.colorbar(cs,ax=ax[2,1])
    ax[2,1].set_title('Alkalinity (mEq/m3)', fontsize = fs)
    ax[2,1].set_xlabel('Hour', fontsize = fs)
    ax[2,1].tick_params(axis='both', which='major', labelsize=ls)
    ax[2,1].xaxis.set_ticks(np.arange(4, 25, 4))
    ax[2,1].set_ylim((z_min,z_max))
    plt.show()

    # print what the sigma layers look like (so we can see whether they are too thin)
    fig, ax = plt.subplots(1,1,figsize = (9,5))
    # print(z_w.shape)
    for i in range(0,31):
        plt.plot(t,z_w[i],color='k')
    plt.title('Sigma Layers at {} ({})'.format(site_name,date),fontsize=ts)
    plt.xlabel('Hour',fontsize=fs)
    plt.xlim(1,20)
    plt.ylim((z_min,z_max))
    plt.ylabel('Depth of Sigma Layer Boundary (m)', fontsize=fs)
    ax.tick_params(axis='both', which='major', labelsize=ls)
    plt.xticks(np.arange(4, 25, 4))
    plt.show()

    # Plot CFL Number ---------------------------------------------------------------------------------
    fig, ax = plt.subplots(3,1,figsize = (10,8), sharex = True)
    dt = 40 # timestep (seconds)
    # E-W
    # calculate courant number
    Cx = [(abs(U)*dt)/500 for U in u]
    cs = ax[0].pcolormesh(t, z_rho, Cx, vmin=0, vmax=1, cmap=plt.get_cmap(cmocean.cm.amp))
    plt.colorbar(cs,ax=ax[0])
    ax[0].set_ylabel('Z (m)', fontsize = fs)
    ax[0].set_title('C (x-direction)', fontsize = fs)
    ax[0].tick_params(axis='both', which='major', labelsize=ls)
    ax[0].set_ylim((z_min,z_max))
    # N-S
    # calculate courant number
    Cy = [(abs(U)*dt)/500 for U in v]
    cs = ax[1].pcolormesh(t, z_rho, Cy, vmin=0, vmax=1, cmap=plt.get_cmap(cmocean.cm.amp))
    plt.colorbar(cs,ax=ax[1])
    ax[1].set_ylabel('Z (m)', fontsize = fs)
    ax[1].set_title('C (y-direction)', fontsize = fs)
    ax[1].tick_params(axis='both', which='major', labelsize=ls)
    ax[1].set_ylim((z_min,z_max))
    # vertical
    # calculate courant number (using average vertical velocity between edges (to get w at rho point), time thickness of sigma layer)
    dz = np.diff(z_w.values, axis=0)
    w_avgs = (w.values[1::]+w.values[0:-1])/ 2
    Cz = (abs(w_avgs)*dt)/dz
    cs = ax[2].pcolormesh(t, z_rho, Cz, vmin=0, vmax=1, cmap=plt.get_cmap(cmocean.cm.amp))
    plt.colorbar(cs,ax=ax[2])
    ax[2].set_ylabel('Z (m)', fontsize = fs)
    ax[2].set_title('C (z-direction)', fontsize = fs)
    ax[2].set_xlabel('Hour', fontsize = fs)
    ax[2].tick_params(axis='both', which='major', labelsize=ls)
    ax[2].xaxis.set_ticks(np.arange(4, 25, 4))
    ax[2].set_ylim((z_min,z_max))
    # title
    fig.suptitle('CFL Number at {} ({})'.format(site_name,date), fontsize = ts)
    plt.show()

    # Plot Density and Mixing ---------------------------------------------------------------------------------
    fig, ax = plt.subplots(3,1,figsize = (10,8), sharex = True)
    # density
    cs = ax[0].pcolormesh(t, z_rho, rho, vmin = rho_vmin, vmax = rho_vmax, cmap=plt.get_cmap(cmocean.cm.dense))
    plt.colorbar(cs,ax=ax[0])
    ax[0].set_ylabel('Z (m)', fontsize = fs)
    ax[0].set_title(r'$\rho$ (kg $m^3$)', fontsize = fs)
    ax[0].tick_params(axis='both', which='major', labelsize=ls)
    ax[0].set_ylim((z_min,z_max))
    # vertical eddy diffusivity (temperature)
    cs = ax[1].pcolormesh(t, z_w, AKv, vmin = AKv_vmin, vmax = AKv_vmax, cmap=plt.get_cmap(cmocean.cm.matter))
    plt.colorbar(cs,ax=ax[1])
    ax[1].set_ylabel('Z (m)', fontsize = fs)
    ax[1].set_title(r'Vertical Eddy Viscosity $(m^2 \ s^{-1})$', fontsize = fs)
    ax[1].tick_params(axis='both', which='major', labelsize=ls)
    ax[1].set_ylim((z_min,z_max))
    # vertical eddy diffusivity (salt)
    cs = ax[2].pcolormesh(t, z_w, AKs, vmin = AKs_vmin, vmax = AKs_vmax, cmap=plt.get_cmap(cmocean.cm.matter))
    plt.colorbar(cs,ax=ax[2])
    ax[2].set_ylabel('Z (m)', fontsize = fs)
    ax[2].set_title(r'Salt Vertical Eddy Diffusivity $(m^2 \ s^{-1})$', fontsize = fs)
    ax[2].set_xlabel('Hour', fontsize = fs)
    ax[2].tick_params(axis='both', which='major', labelsize=ls)
    ax[2].xaxis.set_ticks(np.arange(4, 25, 4))
    ax[2].set_ylim((z_min,z_max))
    # title
    fig.suptitle('Density and Vertical Mixing at {} ({})'.format(site_name,date), fontsize = ts)
    plt.show()


# print(list(ds.keys()))
