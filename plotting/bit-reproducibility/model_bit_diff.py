"""
Plot difference between surface/bottom values of specified state variable.
Calculates difference between two different runs
(Written to compare long model1 to N-less run)

From ipython: run model_field_diff
Figures saved in LO_output/AL_custom_plots/[vn]_diff.png

"""

###################################################################
##                       import packages                         ##  
###################################################################      

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import pandas as pd
import cmocean
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()
Ldir = Lfun.Lstart()

plt.close('all')

###################################################################
##                          User Inputs                          ##  
################################################################### 

vns = ['NH4'] # u, v, w, NO3, NH4 oxygen
date = '2012.10.07'

filetype = 'ocean_avg_0001.nc'
# filetype = 'ocean_his_0002.nc'

###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

# model1 = 'cas7_newtraps00_debugx11ab'
# model2 = 'cas7_newtraps01_debugx11ab'
model1 = 'cas7_newtraps_x11ab'
model2 = 'cas7_newtrapsnoN_x11ab'

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values

# get model output
fp_model1 = Ldir['roms_out'] / model1 / ('f'+date) / filetype
fp_model2 = Ldir['roms_out'] / model2 / ('f'+date) / filetype
ds_model1 = xr.open_dataset(fp_model1)
ds_model2 = xr.open_dataset(fp_model2)


###################################################################
##                     Calculate differences                     ##  
################################################################### 

for vn in vns:

    # Get data, and get rid of ocean_time dim (because this is at a single time)
    v1 = ds_model1[vn].squeeze()
    v2 = ds_model2[vn].squeeze()

    # Identify vertical and horizontal dims
    if 's_rho' in v1.dims:
        vert_dim = 's_rho'
    elif 's_w' in v1.dims:
        vert_dim = 's_w'
    else:
        raise ValueError(f"No vertical dimension found for {vn}")

    if ('eta_rho' in v1.dims) and ('xi_rho' in v1.dims):
        h_dims = ('eta_rho', 'xi_rho')
        lon = ds_model1['lon_rho']
        lat = ds_model1['lat_rho']
    elif ('eta_u' in v1.dims) and ('xi_u' in v1.dims):
        h_dims = ('eta_u', 'xi_u')
        lon = ds_model1['lon_u']
        lat = ds_model1['lat_u']
    elif ('eta_v' in v1.dims) and ('xi_v' in v1.dims):
        h_dims = ('eta_v', 'xi_v')
        lon = ds_model1['lon_v']
        lat = ds_model1['lat_v']
    elif ('eta_psi' in v1.dims) and ('xi_psi' in v1.dims):
        h_dims = ('eta_psi', 'xi_psi')
        lon = ds_model1['lon_psi']
        lat = ds_model1['lat_psi']
    else:
        raise ValueError(f"Unknown grid type for variable '{vn}'.")
    
    print(list(ds_model1.keys()))

    # Compute strict difference (True where different, False where equal)
    diff_mask = (v1 != v2) | (v1.isnull() != v2.isnull())
    diff_2d = diff_mask.any(dim=vert_dim)

    # Mask out locations where both are NaN at all depths
    both_nan = v1.isnull() & v2.isnull()
    diff_2d = diff_2d.where(~both_nan.any(dim=vert_dim))

    # Convert to numeric: 0 = same, 1 = diff, NaN = both missing
    plot_data = diff_2d.astype(float)

    # Set up colormap: black = 0, lightblue = 1, white = NaN
    cmap = mcolors.ListedColormap(['lightblue', 'black'])
    bounds = [-0.5, 0.5, 1.5]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    ###################################################################
    ##                          Plotting                             ##  
    ################################################################### 

    # Initialize figure
    fig, ax = plt.subplots(1,1, figsize=(10, 8))

    # plot
    map = plt.pcolormesh(lon, lat, plot_data, cmap=cmap, norm=norm, shading='auto')
    cbar = plt.colorbar(map, ticks=[0, 1])
    cbar.ax.set_yticklabels(['No differences', 'Differences'],fontsize=12)

    # format figure
    plt.suptitle('Locations where {} differs between runs at any s-level'.format(vn),fontsize=14,fontweight='bold')
    ax.set_title('{} and {} ({})'.format(model1,model2, filetype))
    plt.xlabel('Lon', fontsize=12)
    plt.ylabel('Lat', fontsize=12)
    pfun.dar(ax)
    ax.set_ylim([46.5,50])
    ax.set_xlim([-126.5,-122])
    plt.tight_layout()
    plt.show()