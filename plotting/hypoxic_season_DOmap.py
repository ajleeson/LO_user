"""
Plots a comparison of bottom DO for two identical grids with different oxygen value.

This is a custom function for a particular experiment,
but code can be adapted for other use cases in the future.

From the terminal: python month_avg_comparison.py

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

##############################################################
##                       USER INPUTS                        ##
##############################################################

year = '2013'

# Variables to compare
vn_list = ['oxygen'] #['oxygen','NO3','NH4','phytoplankton','SdetritusN','LdetritusN']

# Provide information about models to compare
# basecase (natural-- N-less run)
natural_gtagex = 'cas7_t0noN_x4b'
# long hindcast (anthropogenic-- existing conditions)
anthro_gtagex = 'cas7_t0_x4b'
gtagexes = [natural_gtagex, anthro_gtagex]

#'Puget Sound'
region = 'Puget Sound'

# Show WWTP locations?
WWTP_loc = True

###########################################################
# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots' / (anthro_gtagex + '_MINUS_' + natural_gtagex) / 'hypoxic_season'
Lfun.make_dir(out_dir)

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
    ds = xr.open_dataset('../../LO_data/grids/cas7/grid.nc')
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

cbar_pad = 0.02

# get plotting limits based on region
if region == 'Puget Sound':
    # # box extracion limits: [-123.29, -122.1, 46.95, 48.93]
    # xmin = -123.29
    # xmax = -122.1 + 0.1 # to make room for legend key
    # ymin = 46.95 - 0.1 # to make room for legend key
    # ymax = 48.93
    # old box extracion limits
    xmin = -123.2
    xmax = -122.1 + 0.1
    ymin = 46.93 - 0.1
    ymax = 48.45

# # Get grid data
# G = zrfun.get_basic_info('/home/aleeson/LO_data/grids/cas6/grid.nc', only_G=True)
# grid_ds = xr.open_dataset('/home/aleeson/LO_data/grids/cas6/grid.nc')
# lon = grid_ds.lon_rho.values
# lat = grid_ds.lat_rho.values
# lon_u = grid_ds.lon_u.values
# lat_u = grid_ds.lat_u.values
# lon_v = grid_ds.lon_v.values
# lat_v = grid_ds.lat_v.values

# # get month
# if month == 'July':
#     day0 = '07.01'
#     day1 = '07.31'
# elif month == 'August':
#     day0 = '08.01'
#     day1 = '08.31'
# elif month == 'September':
#     day0 = '09.01'
#     day1 = '09.30'

# Loop through variable to compare
for vn in vn_list:

    # make plots for both surface and bottom water
    for stext in ['surface','bottom']:

        if stext == 'surface':
            slev = -1
        elif stext == 'bottom':
            slev = 0

        # set axes range for different state variables
        if vn == 'NO3':
            vmin = 0
            vmax = 40
            cmap = cmocean.cm.matter
        elif vn == 'NH4':
            vmin = 0
            vmax = 6
            cmap = cmocean.cm.matter
        elif vn == 'phytoplankton':
            vmin = 0
            vmax = 30
            cmap = cmocean.cm.algae
        elif vn == 'oxygen':
            vmin = 0
            vmax = 10
            cmap = plt.cm.get_cmap('rainbow_r', 10) #cmocean.cm.thermal #cmocean.cm.thermal#cmocean.cm.oxy
        elif vn == 'SdetritusN':
            vmin = 0
            vmax = 5
            cmap = cmocean.cm.matter
        elif vn == 'LdetritusN':
            vmin = 0
            vmax = 0.1
            cmap = cmocean.cm.matter

        # scale variable & get units
        scale =  pinfo.fac_dict[vn]
        units = pinfo.units_dict[vn]

        # Initialize figure
        fs = 10
        pfun.start_plot(fs=fs, figsize=(31,25))
        fig = plt.figure()
        gs = fig.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.95, wspace=0.05, hspace=0.05)

        # loop through and plot both conditions
        for i,gtagex in enumerate(gtagexes):

            # get data
            fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_'+year+'.01.01_'+year+'.12.31.nc')
            ds = xr.open_dataset(fp)

            # crop to just hypoxic season (may 1 - nov 30)
            ds = ds.sel(ocean_time=slice(np.datetime64(year+'-05-01'),np.datetime64(year+'-12-01')))

            # Get coordinates for pcolormesh
            # get lat and lon
            # if vn == 'u':
            #     lons = ds.coords['lon_u'].values
            #     lats = ds.coords['lat_u'].values
            # elif vn == 'v':
            #     lons = ds.coords['lon_v'].values
            #     lats = ds.coords['lat_v'].values
            # else:
            lons = ds.coords['lon_rho'].values
            lats = ds.coords['lat_rho'].values
            px, py = pfun.get_plon_plat(lons,lats)

            # Plot natural map field
            if i == 0:
                ax = fig.add_subplot(1,2,1)
                v = ds[vn][:,slev,:,:].values * scale
                v = np.nanmean(v,axis=0)
                cs = ax.pcolormesh(px,py,v, vmin=vmin, vmax=vmax, cmap=cmap)
                # add colorbar
                cbar = fig.colorbar(cs, location='left', pad = cbar_pad)
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
                ax.set_title('(a) Natural', fontsize=38)
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
                ax.text(lon0-0.01,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=28)

                # add puget sound map
                inset_map = plt.imread('puget_sound.png')
                imagebox = OffsetImage(inset_map, zoom = 1.05)
                ab = AnnotationBbox(imagebox, (-122.3, 46.98), frameon = False)
                ax.add_artist(ab)

                # add wwtp locations
                if WWTP_loc == True:
                    ax.scatter(lon_wwtps,lat_wwtps,color='none', edgecolors='k', linewidth=3, s=sizes_wwtps, label='WWTPs')

            elif i == 1:
                # save test condition
                v = ds[vn][:,slev,:,:].values * scale
                v = np.nanmean(v,axis=0)
                c1_ds = ds
                c1 = v

        # plot the pcolormesh difference 
        plt.subplots_adjust(wspace=0.05)
        ax = fig.add_subplot(1,2,2)
        diff = (c1 - bc)
        mindiff = np.nanmin(diff)
        maxdiff = np.nanmax(diff)
        # make sure colorbar axis contains zero
        if mindiff > 0 and maxdiff > 0:
            mindiff = maxdiff*-1.01
        if mindiff < 0 and maxdiff < 0:
            maxdiff = mindiff*-1.01
        # don't let colorbar axis scale get too large
        if maxdiff > vmax:
            maxdiff = vmax
        # make sure the colorbar is always centered about zero
        cmap = cmocean.tools.crop(cmocean.cm.balance_r, mindiff, maxdiff, 0)
        cs = ax.pcolormesh(px,py,diff, vmin=mindiff, vmax=maxdiff, cmap=cmap)
        cbar = fig.colorbar(cs, location='right', pad = cbar_pad)
        cbar.ax.tick_params(labelsize=32)
        cbar.outline.set_visible(False)
        # format everything else
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.set_title('(b) Anthropogenic - Natural', fontsize=38)
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
                title='WWTP loading \n'+r' (kg N d$^{-1}$)',loc='lower right', labelspacing=1, borderpad=0.8)
            plt.setp(legend.get_title(),fontsize=28)

                            
        # Add colormap title
        plt.suptitle('Hypoxic season ' + year + ' average ' + stext + ' ' + vn + ' ' + units,
                    fontsize=44, fontweight='bold', y=0.95)

        # Generate plot
        plt.tight_layout
        plt.savefig(out_dir / (vn+'_'+stext+'_avg_diff_may2nov.png'))
        # plt.show()