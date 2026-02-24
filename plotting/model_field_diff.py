"""
Plot difference between surface/bottom values of specified state variable.
Calculates difference between two different runs
(Written to compare long hindcast to N-less run)

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
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
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

vns = ['salt','u', 'v','temp']#['DIN','u','v','oxygen'] # u, v, w, DIN
# date = '2013.04.04'
# date = '2012.10.31'
# date = '2012.10.07'
# date = '2014.01.09'
# date = '2018.06.12'
# date = '2017.02.12'
date = '2020.05.31'

###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

# gtagex_longhindcast = 'cas7_t0noN_x4b'
# gtagex_noN = 'cas7_t0noN_x4b_perf10'

# gtagex_longhindcast = 'cas7_t1_x11ab'#'cas7_t1_x11ab'
# gtagex_noN = 'cas7_t1_x11b'#'cas7_t1noDIN_x11ab'

gtagex_longhindcast = 'cas7_t1d_x11ad'#'cas7_t1_x11ab'
gtagex_noN = 'cas7_t1noDIN_x11ab'#'cas7_t1noDIN_x11ab'

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

# get WWTP locations
wwtp_ind_df = pd.read_csv(Ldir['data'] / 'grids/cas7/wwtp_info.csv', index_col='rname')
wwtp_col = [int(col) for col in wwtp_ind_df['col_py']]
wwtp_row = [int(row) for row in wwtp_ind_df['row_py']]
wwtp_lon = [lon[0,i] for i in wwtp_col]
wwtp_lat = [lat[i,0] for i in wwtp_row]

# layer, text, axis limits
slev_surf = -1
slev_bott = 0
stext_surf = 'Surface'
stext_bott = 'Bottom'

# lon/lat limits (Puget Sound)
xmin = -123.2
xmax = -122.1
ymin = 46.93
ymax = 48.45

# # lon/lat limits (Salish Sea)
# xmin = -125
# xmax = -121
# ymin = 46.9
# ymax = 50

# # lon/lat limits (Full grid)
# xmin = -130
# xmax = -121.5
# ymin = 42
# ymax = 52

# get model output
# fp_hindcast = Ldir['roms_out'] / gtagex_longhindcast / ('f'+date) / 'ocean_his_0025.nc'
# fp_noN = Ldir['roms_out'] / gtagex_noN / ('f'+date) / 'ocean_his_0025.nc'
fp_hindcast = Ldir['roms_out'] / gtagex_longhindcast / ('f'+date+'_192cores') / 'ocean_his_0002.nc'
fp_noN = Ldir['roms_out'] / gtagex_noN / ('f'+date) / 'ocean_his_0002.nc'
ds_hindcast = xr.open_dataset(fp_hindcast)
ds_noN = xr.open_dataset(fp_noN)

# ###################################################################
# ##                     Calculate differences                     ##  
# ################################################################### 

for vn in vns:

    # get coordinates for pcolormesh
    if vn == 'u':
        px, py = pfun.get_plon_plat(lon_u,lat_u)
    elif vn == 'v':
        px, py = pfun.get_plon_plat(lon_v,lat_v)
    else:
        px, py = pfun.get_plon_plat(lon,lat)

    # sum nitrate and ammonium for DIN
    if vn == 'DIN':
        vn_name = 'NO3'
        vn_no3 = 'NO3'
        vn_nh4 = 'NH4'
        # vmin = -5
        # vmax =  5
        vmin = -0.01 #-0.00001
        vmax =  0.01 #0.00001
    elif vn == 'u' or vn == 'v':
        vn_name = vn
        vmin = -0.00001#-0.01
        vmax =  0.00001#0.01
    elif vn == 'oxygen':
        vn_name = vn
        vmin = -0.01 #-0.001
        vmax =  0.01 #0.001
    elif vn == 'salt':
        vn_name = vn
        vmin = -0.00001
        vmax =  0.00001
    elif vn == 'temp':
        vn_name = vn
        vmin = -0.00001
        vmax =  0.00001
    else:
        print('vmin and vmax not provided for '+ vn)

    # scale variable
    scale =  pinfo.fac_dict[vn_name]

    # Get hindcast data
    if vn == 'DIN':
        surf_no3_hindcast = ds_hindcast[vn_no3][0,slev_surf,:,:].values
        surf_nh4_hindcast = ds_hindcast[vn_nh4][0,slev_surf,:,:].values
        bott_no3_hindcast = ds_hindcast[vn_no3][0,slev_bott,:,:].values
        bott_nh4_hindcast = ds_hindcast[vn_nh4][0,slev_bott,:,:].values
    else:
        surf_vn_hindcast = ds_hindcast[vn][0,slev_surf,:,:].values
        bott_vn_hindcast = ds_hindcast[vn][0,slev_bott,:,:].values

    # Get noN data
    if vn == 'DIN':
        surf_no3_noN = ds_noN[vn_no3][0,slev_surf,:,:].values
        surf_nh4_noN = ds_noN[vn_nh4][0,slev_surf,:,:].values
        bott_no3_noN = ds_noN[vn_no3][0,slev_bott,:,:].values
        bott_nh4_noN = ds_noN[vn_nh4][0,slev_bott,:,:].values
    else:
        surf_vn_noN = ds_noN[vn][0,slev_surf,:,:].values
        bott_vn_noN = ds_noN[vn][0,slev_bott,:,:].values

    # Calculate DIN
    if vn == 'DIN':
        surf_DIN_hindcast = surf_no3_hindcast + surf_nh4_hindcast
        bott_DIN_hindcast = bott_no3_hindcast + bott_nh4_hindcast
        surf_DIN_noN = surf_no3_noN + surf_nh4_noN
        bott_DIN_noN = bott_no3_noN + bott_nh4_noN

    # Get difference
    if vn == 'DIN':
        surf_diff = (surf_DIN_hindcast - surf_DIN_noN) * scale
        bott_diff = (bott_DIN_hindcast - bott_DIN_noN) * scale
    else:
        surf_diff = (surf_vn_hindcast - surf_vn_noN) * scale
        bott_diff = (bott_vn_hindcast - bott_vn_noN) * scale

    # ###################################################################
    # ##                  Plotting and saving figure                   ##  
    # ################################################################### 

    # Initialize figure
    fig = plt.figure(figsize=(12,9)) # 15,11 for Puget sound and 18,8 for Salish Sea
    plt.tight_layout()

    subplotnums = [121,122]
    stexts = [stext_surf,stext_bott]
    values = [surf_diff,bott_diff]

    newcmap = cmocean.cm.balance_r

    # loop through all of the plots we need to make
    for i,stext in enumerate(stexts):

        # add water/land
        ax = fig.add_subplot(subplotnums[i])

        # plot values
        cs = ax.pcolormesh(px,py,values[i],vmin=vmin, vmax=vmax, cmap=newcmap)
        cbar = fig.colorbar(cs)
        cbar.ax.tick_params(labelsize=14)
        cbar.outline.set_visible(False)
        # format figure
        # ax.set_xlim([xmin,xmax])
        # ax.set_ylim([ymin,ymax])
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.axis('off')
        # pfun.add_coast(ax, color='k')
        pfun.dar(ax)
        ax.set_title(vn + ' difference at ' + stext + pinfo.units_dict[vn_name], fontsize=16)
        fig.suptitle('{} minus {}\n'.format(gtagex_longhindcast,gtagex_noN) + date + ' ocean_avg_0001',
                    fontsize=18, fontweight='bold')

        # # add 10 km bar
        # lat0 = 47
        # lon0 = -122.4
        # lat1 = lat0
        # lon1 = -122.27
        # distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
        # x_dist_km = round(distances_m[0]/1000)
        # # ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=6)
        # # ax.text((lon0+lon1)/2,lat0+0.02,'{} km'.format(x_dist_km),color='k',
        # #         horizontalalignment='center', fontsize=15)
        
        # # add WWTP locations
        # ax.scatter(wwtp_lon,wwtp_lat,s=30,alpha=0.5,
        #         facecolors='none',edgecolors='deeppink')

    # Generate plot
    plt.tight_layout
    plt.subplots_adjust(wspace=0.05)
    plt.savefig(out_dir / (date+'_'+vn+'_diff.png'))