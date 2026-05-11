"""
Map monthly climatological means in parameter space
"""

# import things
import numpy as np
import xarray as xr
import pandas as pd
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

out_dir = Ldir['LOo'] / 'chapter_3' / 'data'


#############################################################
#                        PLOT DATA                         ##
#############################################################

# vns = ['CO2_capacity_ideal', 'CO2_flux_actual', 'surf_temp', 'surf_salt',
# 'surf_alk', 'surf_TIC', 'wind2', 'delta_pCO2', 'Fs_31', 'Fs_30', 'Fs_29']

stats = ['mean', 'min', 'max']

# initialize figure
fig,ax = plt.subplots(3,4, figsize=(10,8), sharex=True, sharey=True)
ax = ax.ravel()

locations=['Columbia River Plume',
        'Saratoga Passage',
        'Hood Canal',
        'Van Island Coast',
        'Quadra Island',
        'S. SoG (Fraser plume)',
        'S. of San Juan Islands']
natural_ys = [464,875,716,1047,1199,1075,898]
natural_xs = [360,578,504,186, 233, 468, 513]


# loop through months
for m,month in enumerate(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']):

    # open dataset for this month
    ds_month = xr.open_dataset(out_dir / f'monthly_climatology_freshwaterCO2_{month}.nc')
    ds_month_sml = xr.open_dataset(out_dir / f'monthly_climatology_SML_{month}.nc')

    # get freshwater content, wind speed, and CO2 capacity
    Fs = ds_month['Fs_30_mean'].values
    wind = np.sqrt(ds_month['wind2_mean'].values)
    CO2_capacity = ds_month['CO2_capacity_ideal_mean'].values
    sml = ds_month_sml['SML_thickness_mean'].values

    # loop through locations
    for l,(x,y,loc) in enumerate(zip(natural_xs,natural_ys,locations)):

        if loc in ['Quadra Island', 'S. of San Juan Islands']:
            edgecolor = 'deeppink'
        else:
            edgecolor = 'silver'
        edgecolor = 'none'

        # plot location in parameter space
        cs = ax[m].scatter(Fs[y,x], wind[y,x], s=80, edgecolor=edgecolor, zorder=5,
                      c=CO2_capacity[y,x], cmap=cmc.vik, vmin=-30, vmax=30)
        # cs = ax[m].scatter(Fs[y,x], sml[y,x], s=80, edgecolor=edgecolor, zorder=5,
        #               c=CO2_capacity[y,x], cmap=cmc.vik, vmin=-30, vmax=30)
        # cs = ax[m].scatter(sml[y,x], wind[y,x], s=80, edgecolor=edgecolor, zorder=5,
        #               c=CO2_capacity[y,x], cmap=cmc.vik, vmin=-30, vmax=30)

    # create colorbarlegend
    if m == 0 :
        cbar_ax = fig.add_axes([0.1, 0.13, 0.81, 0.03])
        cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=12, rotation=30)
        cbar.outline.set_visible(False)
        cbar.set_label(r'CO2 uptake capacity [mmol m$^{-3}$]', fontsize=12)

    # format figure
    ax[m].set_xlim([0,4])
    ax[m].set_ylim([0,8])
    ax[m].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    if month in ['Feb','May','Aug','Nov']:
        color = 'yellowgreen'
    else:
        color='black'
    ax[m].set_title(month,fontsize=12, fontweight='bold',loc='left', color=color)
    # plt.suptitle('(2015 - 2024) '+ month + ' - ' + title,
    #                 fontsize=14, fontweight='bold')

    # label axis
    if m in [0,4,8]:
        ax[m].set_ylabel('Wind speed [m s$^{-1}$]', fontsize=12)
    if m >= 8:
        ax[m].set_xlabel(r'Fs (s$_0$ = 30) [m]', fontsize=12)

    # Generate plot
    plt.tight_layout
    plt.subplots_adjust(bottom=0.25, top=0.92)
    plt.show()