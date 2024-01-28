"""
Plots a comparison of bottom DO for two identical grids with different oxygen value.

This is a custom function for a particular experiment,
but code can be adapted for other use cases in the future.

From the terminal: python money_plot.py

"""

# import things
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

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# ----------------------------------------------------------------------

# Provide information about models to compare
# basecase (N-less run)
bc_gtagex = 'cas6_traps3_x2b'
# test condition (long hindcast)
c1_gtagex = 'cas6_traps2_x2b'
gtagexes = [bc_gtagex, c1_gtagex]

#'Puget Sound','Whidbey Basin','North King County','Lynch Cove'
region = 'Puget Sound'

# Variables to compare
vn_list = ['oxygen']#['NO3','phytoplankton','oxygen']

# Show WWTP locations?
WWTP_loc = True

##########################################################
# helper funcions

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = '../../LO_data/trapsD00/LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_LO = repeatrivs_df.loc[repeatrivs_df['SSM_rname'] == rname, 'LO_rname'].values[0]

    return rname_LO

def LO2SSM_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = Ldir['data'] / 'trapsD00' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_SSM = repeatrivs_df.loc[repeatrivs_df['LO_rname'] == rname, 'SSM_rname'].values[0]

    return rname_SSM

# riv fun function to get biology
def get_bio_vec(vn, rn, yd_ind):
    ndt = len(yd_ind)
    yd = yd_ind.values
    ovec = np.ones(ndt)
    if vn == 'NO3':
        if rn == 'fraser':
            vv = 2 + (13/2) + (13/2)*np.cos(2*np.pi*((yd-30)/366))
        elif rn == 'columbia':
            vv = 5 + (35/2) + (35/2)*np.cos(2*np.pi*((yd)/366))
        else:
            vv = 5 * ovec
    elif vn == 'Oxyg':
        vv = 350 * ovec
    elif vn in ['TAlk', 'TIC']:
        if rn in ['columbia', 'deschutes', 'duwamish']:
            vv = 1000 * ovec
        else:
            vv = 300 * ovec
    else:
        vv = 0 * ovec # all others filled with zeros
    return vv

##########################################################
if WWTP_loc == True:
    # set up the time index for the record
    Ldir = Lfun.Lstart()
    dsf = Ldir['ds_fmt']
    dt0 = datetime.strptime('2020.01.01',dsf)
    dt1 = datetime.strptime('2020.12.31',dsf)
    days = (dt0, dt1)
        
    # pandas Index objects
    dt_ind = pd.date_range(start=dt0, end=dt1)
    yd_ind = pd.Index(dt_ind.dayofyear)

    # Get LiveOcean grid info --------------------------------------------------

    # get the grid data
    ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc')
    z = -ds.h.values
    mask_rho = np.transpose(ds.mask_rho.values)
    lon = ds.lon_rho.values
    lat = ds.lat_rho.values
    X = lon[0,:] # grid cell X values
    Y = lat[:,0] # grid cell Y values
    plon, plat = pfun.get_plon_plat(lon,lat)
    # make a version of z with nans where masked
    zm = z.copy()
    zm[np.transpose(mask_rho) == 0] = np.nan
    zm[np.transpose(mask_rho) != 0] = -1

    # Point Sources -------------------------------------------------------------
    # Prepare data for spatial summary plots

    # get flow, nitrate, and ammonium values
    fp_wwtps = '../../LO_output/pre/trapsP00/point_sources/lo_base/Data_historical/'
    flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow.p')    # m3/s
    no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3.p')      # mmol/m3
    nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4.p')      # mmol/m3

    # calculate total DIN concentration in mg/L
    dindf_wwtps = (no3df_wwtps + nh4df_wwtps)/71.4    # mg/L

    # calculate daily loading timeseries in kg/d
    dailyloaddf_wwtps = 86.4*dindf_wwtps*flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s

    # calculate average daily load over the year (kg/d)
    avgload_wwtps = dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

    # add row and col index for plotting on LiveOcean grid
    griddf0_wwtps = pd.read_csv('../../LO_data/grids/cas6/wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    avgload_wwtps = avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    avgload_wwtps = avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    # get point source lat and lon
    lon_wwtps = [X[int(col)] for col in avgload_wwtps['col_py']]
    lat_wwtps = [Y[int(row)] for row in avgload_wwtps['row_py']]
    
    # define marker sizes (minimum size is 10 so dots don't get too small)
    sizes_wwtps = [max(0.3*load,30) for load in avgload_wwtps['avg-daily-load(kg/d)']]

##########################################################

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# get plotting limits based on region
if region == 'North King County':
    xmin = -122.6
    xmax = -122.3
    ymin = 47.55
    ymax = 47.75
elif region == 'Lynch Cove':
    xmin = -123.2
    xmax = -122.8
    ymin = 47.3
    ymax = 47.5
elif region == 'Whidbey Basin':
    xmin = -122.8
    xmax = -122.1
    ymin = 48
    ymax = 48.4
elif region == 'Puget Sound':
    xmin = -123.2
    xmax = -122.1
    ymin = 46.93
    ymax = 48.45

# Get grid data
G = zrfun.get_basic_info('/home/aleeson/LO_data/grids/cas6/grid.nc', only_G=True)
grid_ds = xr.open_dataset('/home/aleeson/LO_data/grids/cas6/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values

# Loop through variable to compare
for vn in vn_list:
    if vn == 'NO3':
        slev = -1
        stext = 'Surface'
        vmin = 0
        vmax = 44
        cmap = cmocean.cm.matter
    elif vn == 'NH4':
        slev = -1
        stext = 'Surface'
        vmin = 0
        vmax = 5
        cmap = cmocean.cm.matter
    elif vn == 'phytoplankton':
        slev = -1
        stext = 'Surface'
        vmin = 0
        vmax = 10
        cmap = cmocean.cm.algae
    elif vn == 'oxygen':
        slev = 0
        stext = 'Bottom'
        vmin = 0
        vmax = 10
        cmap = cmocean.cm.oxy

    # scale variable
    scale =  pinfo.fac_dict[vn]

    # Initialize figure
    fs = 10
    pfun.start_plot(fs=fs, figsize=(36,27))
    fig = plt.figure()
    # gs = fig.add_gridspec(nrows=4, ncols=3, left=0.05, right=0.48, wspace=0.05)
    gs = fig.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.95, wspace=0.05, hspace=0.05)

    # loop through and plot both conditions
    for i,gtagex in enumerate(gtagexes):

        # get data
        fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / 'prelimDO_2017.08.01_2017.08.31.nc'
        ds = xr.open_dataset(fp)

        # Get coordinates for pcolormesh
        px, py = pfun.get_plon_plat(ds.coords['lon_rho'].values,ds.coords['lat_rho'].values)

        # Plot basecase map field
        if i == 0:
            ax = fig.add_subplot(1,2,1)
            v = ds[vn][:,slev,:,:].values * scale
            v = np.nanmean(v,axis=0)
            cs = ax.pcolormesh(px,py,v, vmin=vmin, vmax=vmax, cmap=cmap)
            # add colorbar
            cbar = fig.colorbar(cs, location='left')
            cbar.ax.tick_params(labelsize=32)#,length=10, width=2)
            cbar.outline.set_visible(False)
            # format figure
            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.axis('off')
            # pfun.add_coast(ax)
            pfun.dar(ax)
            ax.set_title('(a) N-less run', fontsize=36)
            # save dataset for later use
            bc_ds = ds
            bc = v

            # add 10 km bar
            lat0 = 46.94
            lon0 = -123.05
            lat1 = lat0
            lon1 = -122.91825
            distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
            x_dist_km = round(distances_m[0]/1000)
            ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=8)
            ax.text(lon0,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=28)

            # add puget sound map
            inset_map = plt.imread('puget_sound.png')
            imagebox = OffsetImage(inset_map)#, zoom = 0.15)
            ab = AnnotationBbox(imagebox, (-122.3, 47.05), frameon = False)
            ax.add_artist(ab)

        elif i == 1:
            # save test condition
            v = ds[vn][:,slev,:,:].values * scale
            v = np.nanmean(v,axis=0)
            c1_ds = ds
            c1 = v

    # plot the pcolormesh difference 
    plt.subplots_adjust(wspace=0.01)
    ax = fig.add_subplot(1,2,2)
    diff = (c1 - bc)
    vmin = np.nanmin(diff)
    vmax = np.nanmax(diff)
    # make sure the colorbar is always centered about zero
    cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)
    cs = ax.pcolormesh(px,py,diff, vmin=vmin, vmax=vmax, cmap=cmap)
    cbar = fig.colorbar(cs, location='right')
    cbar.ax.tick_params(labelsize=32)
    cbar.outline.set_visible(False)
    # format everything else
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_title('(b) Anomaly (Hindcast minus N-less run)', fontsize=36)
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.axis('off')
    # pfun.add_coast(ax)
    pfun.dar(ax)

    # add wwtp locations
    if WWTP_loc == True:
        ax.scatter(lon_wwtps,lat_wwtps,color='none', edgecolors='k', linewidth=3, s=sizes_wwtps, label='WWTPs')
        leg_szs = [100, 1000, 10000]
        szs = [0.3*(leg_sz) for leg_sz in leg_szs]
        l0 = plt.scatter([],[], s=szs[0], color='none', edgecolors='k', linewidth=3)
        l1 = plt.scatter([],[], s=szs[1], color='none', edgecolors='k', linewidth=3)
        l2 = plt.scatter([],[], s=szs[2], color='none', edgecolors='k', linewidth=3)
        labels = ['< 100', '1,000', '10,000']
        legend = ax.legend([l0, l1, l2], labels, fontsize = 24, markerfirst=False,
            title='WWTP N loading \n'+r' (kg N d$^{-1}$)',loc='lower right', labelspacing=1, borderpad=0.8)
        plt.setp(legend.get_title(),fontsize=28)

                           
    # Add colormap title
    plt.suptitle('August 2017 Average Bottom DO [mg/L]',
                 fontsize=44, fontweight='bold')#, x=0.64, y=0.94)

    # Generate plot
    plt.tight_layout
    plt.savefig(out_dir / (stext+'_'+vn+'_month_avg_diff.png'))
    # plt.show()