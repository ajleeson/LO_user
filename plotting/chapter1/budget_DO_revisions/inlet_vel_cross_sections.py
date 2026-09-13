"""
Generate 2017 time-average velocity cross section plots
at the mouth of each inlet
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
import cmcrameri.cm as cmc
import matplotlib.pylab as plt
import gsw
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

# reload to make editing easier
from importlib import reload
reload(pinfo)

Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t1_x11b'

inlets = ['sinclair','quartermaster','lynchcove','penn','crescent','dyes',
          'case','holmes','elliot','carr','portsusan','commencement','dabob']

##########################################################
##                      Plotting                        ##
##########################################################

plt.close('all')

vmin = -0.12
vmax =  0.12

fig,axes = plt.subplots(5,3,figsize=(10,8))
ax = axes.ravel()

for i,inlet in enumerate(inlets):
    ds = xr.open_dataset('../../../../LO_output/extract/' + gtagex +
                     '/tef2/c21/extractions_2017.01.01_2017.12.31/' +
                     inlet + '.nc')

    # get time-averaged velocity and depth and zeta
    vel_tavg  = np.nanmean(ds['vel'].values, axis=0)   # (nz, nx)
    dz_tavg = np.nanmean(ds['DZ'].values,  axis=0)     # (nz, nx) vertical thickness
    zeta_tavg = np.nanmean(ds['zeta'].values, axis=0)  # (nx,) surface elevation
    lat_thick = ds['dd'].values                        # (nx,) horizontal thickness

    # get size of the inlet
    nz, nx = vel_tavg.shape

    # --- x edges (nx+1) from horizontal thickness ---
    x_edges = np.concatenate([[0.0], np.cumsum(lat_thick)])/1000   # (nx+1,) [km]

    # --- z edges per column (nz+1, nx) from vertical thickness ---
    # top boundary set by time-averaged zeta, multiply dz by -1 because depth goes down
    z_edges_col = np.vstack([zeta_tavg,np.cumsum(dz_tavg, axis=0) * -1])  # (nz+1, nx) [m]

    # --- convert column-based edges to corner grid (nz+1, nx+1) ---
    # by taking the average of all of the edges in each column,
    # except for the left and right edges which are set by the first and last column
    Z_edges = np.empty((nz+1, nx+1))
    Z_edges[:, 1:-1] = 0.5 * (z_edges_col[:, :-1] + z_edges_col[:, 1:])
    Z_edges[:, 0]    = z_edges_col[:, 0]
    Z_edges[:, -1]   = z_edges_col[:, -1]

    # X corner grid
    # repeat the x edges in the z dimensions
    # since the lateral thickness is unchanging
    X_edges = np.broadcast_to(x_edges, (nz+1, nx+1))

    # Plot
    cm = ax[i].pcolormesh(X_edges, Z_edges, vel_tavg, cmap=cmc.vik, vmin=vmin,vmax=vmax,shading='flat')
    # create colorbarlegend
    if i == 0:
        cbar_ax = fig.add_axes([0.4, 0.1, 0.5, 0.03])
        cbar = fig.colorbar(cm, cax=cbar_ax,orientation='horizontal')
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label('Outflow                         [m/s]                             Inflow',
                       fontsize=14, fontweight='bold')
        cbar.outline.set_visible(False)
    # format figure
    ax[i].set_title(inlet,fontsize=14,loc='left', fontweight='bold')
    ax[len(inlets)].set_visible(False)
    ax[len(inlets)+1].set_visible(False)
    if i >= 10:
        ax[i].set_xlabel('Distance along cross-section [km]')
    if i in [0,3,6,9,12]:
        ax[i].set_ylabel('z [m]')

    plt.tight_layout()
