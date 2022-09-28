"""
Plots a comparison of bottom DO for two identical grids with different oxygen value.

This is a custom function for a particular experiment, but code can be adapted for other use cases in the future.

From the terminal: python bottom_DO_comparison.py

"""

# import things
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
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

# ----------------------------------------------------------------------

# Provide information about models to compare
bc_gtagex = 'alpe2_npzd0nonutr_bpart'
c1_gtagex = 'alpe2_npzd0withnutr_bpart'
gtagexes = [bc_gtagex, c1_gtagex]

hr = '0025'
date = '2020.03.31'

# Variables to compare
vn_list = ['phytoplankton', 'oxygen']

# Loop through variable to compare
for vn in vn_list:
    if vn == 'phytoplankton':
        slev = -1
        stext = 'Surface'
    elif vn == 'oxygen':
        slev = 0
        stext = 'Bottom'

    # Initialize figure
    fs = 10
    pfun.start_plot(fs=fs, figsize=(6,9))
    fig = plt.figure()

    # loop through and plot both conditions
    for i,gtagex in enumerate(gtagexes):

        # get data
        fp = Ldir['roms_out'] / gtagex / ('f' + date) / ('ocean_his_'+ hr +'.nc')
        ds = xr.open_dataset(fp)

        # Plot individual map field
        ax = fig.add_subplot(3, 2, i+1)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
                    cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                    vlims_fac=pinfo.range_dict[vn], do_mask_edges=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        plt.suptitle('%s %s %s' % (stext, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.5*fs)
        ax.set_xlabel('Longitude')
        if i == 0:
            ax.set_title('No WWTP Nutrients')
            ax.set_ylabel('Latitude')
            # save dataset for later use
            bc_ds = ds
        elif i ==1:
            ax.set_title('With WWTP Nutrients')
            # fig.colorbar(cs)
            # add colorbar
            cb_ax = fig.add_axes([.91,.67,.03,.23])
            fig.colorbar(cs,orientation='vertical',cax=cb_ax)
            # save dataset for later use
            c1_ds = ds

    # plot the difference
    plt.subplots_adjust(hspace=0.5)
    ax = fig.add_subplot(3, 2, (3,6))
    G = zrfun.get_basic_info('/home/aleeson/LO_data/grids/alpe2/grid.nc', only_G=True)
    diff = c1_ds[vn] - bc_ds[vn]
    px, py = pfun.get_plon_plat(G['lon_rho'],G['lat_rho'])
    vmin = np.min(diff[0,slev,:,:])
    vmax = np.max(diff[0,slev,:,:])*0.1
    cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)
    cs = ax.pcolormesh(px,py,diff[0,slev,:,:], vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_title('(With Nutrients Condition) - (No Nutrients Condition)')
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')
    cb_ax = fig.add_axes([.91,.1,.03,.5])
    fig.colorbar(cs,orientation='vertical',cax=cb_ax)




    # Generate plot
    plt.show()