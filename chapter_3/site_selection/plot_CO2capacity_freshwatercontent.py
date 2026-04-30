"""
Get CO2 uptake capacity from model history file
and freshwater content
and other useful values for understanding the drivers of CO2 uptake capacity

.nc files are saved in LO_output/chapter_3/data
"""

# import things
import numpy as np
import xarray as xr
import csv
import matplotlib.pylab as plt
from pathlib import Path
from datetime import datetime
import gsw
import PyCO2SYS as pyco2
import cmcrameri.cm as cmc
from lo_tools import Lfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

plt.close('all')

##############################################################
##                       USER INPUTS                        ##
##############################################################

years = ['2019']#['2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025']

# which  model run to look at?
gtagex = 'cas7_t1_x11ab'

##############################################################
##                        GET DATA                          ##
##############################################################

for year in years:

    # get the dataset
    ds = xr.open_dataset(Ldir['LOo'] / 'chapter_3' / 'data' / (gtagex + '_' + year + '_freshwatercontent_CO2uptake.nc'))

    yrday = 87

    print(ds['ocean_time'][yrday].values)

    # get pcolormesh values
    CO2_capacity_ideal = ds['CO2_capacity_ideal'][yrday,:,:].values
    CO2_flux_actual = ds['CO2_flux_actual'][yrday,:,:].values
    surf_salt = ds['surf_salt'][yrday,:,:].values
    Fs_30 = ds['Fs_30'][yrday,:,:].values

    # make list of vars
    vars = [CO2_capacity_ideal,CO2_flux_actual, surf_salt,Fs_30]

    # list of vmins and vmax
    vmins = [-15,-50,25,0]
    vmaxs = [15,50,32,5]

    # get colormaps
    cmaps = [cmc.vik,cmc.vik,cmc.imola,cmc.batlowW_r]

    # titles
    titles = ['CO2 uptake capacity based on\n'+r'solubility and $\Delta$ pCO2 [mmol CO2 m$^{-3}$]',
            'CO2 flux over one day\n'+r'accounting for wind [mmol CO2 m$^{-2}$ day$^{-1}$]',
            'Surface practical salinity',
            'Freshwater content\n'+r'(s$_0$ = 30) [m]']

    # Get grid data
    G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
    grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
    lon = grid_ds.lon_rho.values
    lat = grid_ds.lat_rho.values
    lon_u = grid_ds.lon_u.values
    lat_u = grid_ds.lat_u.values
    lon_v = grid_ds.lon_v.values
    lat_v = grid_ds.lat_v.values
    px, py = pfun.get_plon_plat(lon,lat)

    # lon/lat limits (Study Domain)
    xmin = -126
    xmax = -122
    ymin = 45.5
    ymax = 50.5

    for var, vmin, vmax, cmap, title in zip(vars, vmins, vmaxs, cmaps, titles):

        # Initialize figure
        fig,ax = plt.subplots(1,1, figsize=(7,9))

        # plot values
        cs = ax.pcolormesh(px,py,var,vmin=vmin, vmax=vmax, cmap=cmap)

        # Add Puget Sound Inset
        # [x0, y0, width, height]
        axins = ax.inset_axes([0.71, 0.0, 0.45, 0.6])
        # plot values in inset
        axins.pcolormesh(px, py, var, vmin=vmin, vmax=vmax, cmap=cmap)
        # Puget Sound limits
        axins.set_xlim(-123.2, -122.1)
        axins.set_ylim(46.95, 48.4)
        axins.tick_params(left=False, bottom=False)
        # format
        axins.set_xticklabels([])
        axins.set_yticklabels([])
        for spine in axins.spines.values():
            spine.set_edgecolor('grey')
            spine.set_linewidth(2)

        # add colorbar
        cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
        cbar.ax.tick_params(labelsize=14, rotation=30)
        cbar.outline.set_visible(False)

        # format figure
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.tick_params(left=False, bottom=False)
        pfun.add_coast(ax, color='silver')
        pfun.add_coast(axins, color='silver')
        pfun.dar(ax)
        pfun.dar(axins)
        ax.set_title(title, fontsize=14,
                    loc='Left', fontweight='bold')

        # Generate plot
        plt.tight_layout
        plt.subplots_adjust(bottom=0.001, top=0.9)
        plt.show()