"""
Compare average bottom DO between multiple years
(Set up to run for 6 years)

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
import csv
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator
import matplotlib.patches as patches
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
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

##############################################################
##                       USER INPUTS                        ##
##############################################################

# DO threshold for calculating days of bottom DO < threshold,
# and bottom area with DO < threshold
DO_thresh = 2 # [mg/L]

# Show WWTP locations?
WWTP_loc = False

remove_straits = True

vn = 'oxygen'

years =  ['2014','2015','2016','2017','2018','2019']

# which  model run to look at?
gtagex = 'cas7_t0_x4b' # long hindcast (anthropogenic)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures'
Lfun.make_dir(out_dir)

region = 'Puget Sound'

plt.close('all')

##############################################################
##                    HELPER FUNCTIONS                      ##
##############################################################

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = '../../../LO_data/trapsD00/LiveOcean_SSM_rivers.xlsx'
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
    ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
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

    # get flow, nitrate, and ammonium values
    fp_wwtps = '../../../LO_output/pre/trapsP00/point_sources/lo_base/Data_historical/'
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
    griddf0_wwtps = pd.read_csv('../../../LO_data/grids/cas6/wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    avgload_wwtps = avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    avgload_wwtps = avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    # get point source lat and lon
    lon_wwtps = [X[int(col)] for col in avgload_wwtps['col_py']]
    lat_wwtps = [Y[int(row)] for row in avgload_wwtps['row_py']]
    
    # define marker sizes (minimum size is 10 so dots don't get too small)
    sizes_wwtps = [max(0.3*load,30) for load in avgload_wwtps['avg-daily-load(kg/d)']]

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# open datasets
if remove_straits:
    straits = 'noStraits'
else:
    straits = 'withStraits'

# get the grid data
ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
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

# initialize empty dictionaries
DO_bot_dict = {} # dictionary with DO_bot values
hyp_thick_dict = {}

for year in years:
    # add ds to dictionary
    ds = xr.open_dataset(Ldir['LOo'] / 'pugetsound_DO' / 'data' / (year + '_DO_info_' + straits + '.nc'))
    DO_bot = ds['DO_bot'].values
    hyp_thick = ds['hyp_thick'].values
    # if not a leap year, add a nan on feb 29 (julian day 60 - 1 because indexing from 0)
    if np.mod(int(year),4) != 0: 
        DO_bot = np.insert(DO_bot,59,'nan',axis=0)
        hyp_thick = np.insert(hyp_thick,59,'nan',axis=0)
    DO_bot_dict[year] = DO_bot
    hyp_thick_dict[year] = hyp_thick

# calculate average of all of the arrays
v_avg = sum(DO_bot_dict.values())/len(DO_bot_dict)
# add average to dictionary
DO_bot_dict['avg'] = v_avg

# initialize new dictionary for number of days with bottom DO < threshold
DO_days = {}
DO_bot_threshold_dict = {} # and dictionary with 1 if DO < thresh, else nan
# get days with DO < threshold (for average year, and difference in days for individual years)
for year in ['avg'] + years:
    # Get boolean array. True if DO < threshold mg/L, otherwise, nan
    DO_bot_threshold = np.where(DO_bot_dict[year] < DO_thresh, 1, np.nan)
    DO_bot_threshold_dict[year] = DO_bot_threshold
    # Sum over time to compress into just lat/lon dimension, with values indicating days with bottom DO < 2mg/L
    DO_days_threshold = np.nansum(DO_bot_threshold, axis=0)
    # save difference from average for all years
    if year != 'avg':
        DO_days_threshold = DO_days_threshold - DO_days['avg']
    # convert zeros to nans in dataset
    DO_days_threshold = np.where(DO_days_threshold == 0, np.nan, DO_days_threshold)
    DO_days[year] = DO_days_threshold

# initialize dictionary for hypoxic area
hyp_area = {}
for year in years:
    # calculate bottom hypoxic area
    hyp_area_arr = 0.5 * 0.5 * DO_bot_threshold_dict[year] # area of grid cell = 0.5 km by 0.5 km times number of gridcells with hypoxia
    # get timeseries of bottom hypoxic area (by summing over spatial dimensions)
    hyp_area_timeseries = np.nansum(hyp_area_arr,axis=1)
    hyp_area_timeseries = np.nansum(hyp_area_timeseries,axis=1)
    hyp_area[year] = hyp_area_timeseries

# get grid cell area
fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_2013.01.01_2013.12.31.nc')
ds_2013 = xr.open_dataset(fp)
DX = (ds_2013.pm.values)**-1
DY = (ds_2013.pn.values)**-1
DA = DX*DY*(1/1000)*(1/1000) # get area, but convert from m^2 to km^2

# initialize dictionary for hypoxic volume [km3]
hyp_vol = {}
for year in years:
    # get hypoxic thickness
    hyp_thick = hyp_thick_dict[year]/1000 # [km]
    hyp_vol_timeseries = np.sum(hyp_thick * DA, axis=(1, 2)) # km^3
    hyp_vol[year] = hyp_vol_timeseries

##############################################################
##                         MAKE MAP                         ##
##############################################################

cbar_pad = 0.02

letters= ['(a)','(b)','(c)',
          '(d)','(e)','(f)',
          '(g)']

# get plotting limits based on region
if region == 'Puget Sound':
    # box extracion limits: [-123.29, -122.1, 46.95, 48.93]
    xmin = -123.29
    xmax = -122.1
    ymin = 46.95
    ymax = 48.93

# Initialize figure
fs = 10
pfun.start_plot(fs=fs, figsize=(24,24))
fig = plt.figure()
# ax = axes.ravel()

gs = fig.add_gridspec(2,6)
ax0 = fig.add_subplot(gs[:, 0:3])

# get lat and lon
fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_2013.01.01_2013.12.31.nc')
ds = xr.open_dataset(fp)
lons = ds.coords['lon_rho'].values
lats = ds.coords['lat_rho'].values
px, py = pfun.get_plon_plat(lons,lats)


# loop through and plot both conditions
for i,year in enumerate(['avg'] + years):

    # get days
    DO_bot = DO_days[year]
                
    if i == 0:
        ax = ax0
        ax.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-1.25, vmax=0, cmap=plt.get_cmap('Greys'))
        cs = ax.pcolormesh(px,py,DO_bot, vmin=0, vmax=np.nanmax(DO_bot), cmap='rainbow')
        cbar = fig.colorbar(cs, location='left')
        cbar.ax.tick_params(labelsize=28)#,length=10, width=2)
        cbar.outline.set_visible(False)
        ax.set_title(letters[0] + ' six-year avg.', fontsize=28, loc='left')

        # add wwtp locations
        if WWTP_loc == True:
            ax.scatter(lon_wwtps,lat_wwtps,color='none', edgecolors='k', linewidth=3, s=sizes_wwtps, label='WWTPs')
            leg_szs = [100, 1000, 10000]
            szs = [0.3*(leg_sz) for leg_sz in leg_szs]
            l0 = plt.scatter([],[], s=szs[0], color='none', edgecolors='k', linewidth=3)
            l1 = plt.scatter([],[], s=szs[1], color='none', edgecolors='k', linewidth=3)
            l2 = plt.scatter([],[], s=szs[2], color='none', edgecolors='k', linewidth=3)
            labels = ['< 100', '1,000', '10,000']
            legend = ax.legend([l0, l1, l2], labels, fontsize = 18, markerfirst=False,
                title='WWTP loading \n'+r' (kg N d$^{-1}$)',loc='lower right', labelspacing=1, borderpad=0.8)
            plt.setp(legend.get_title(),fontsize=20)

        # add 10 km bar
        lat0 = 47.0
        lon0 = -123.05
        lat1 = lat0
        lon1 = -122.91825
        distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
        x_dist_km = round(distances_m[0]/1000)
        ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=8)
        ax.text(lon0-0.04,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=24)

    else: 
        if i in [1,2,3]:
            ax = fig.add_subplot(gs[0,i+2])
        elif i in [4,5,6]:
            ax = fig.add_subplot(gs[1,i-1])
        # Create map
        ax.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-1.25, vmax=0, cmap=plt.get_cmap('Greys'))
        cs_all = ax.pcolormesh(px,py,DO_bot, vmin=-60, vmax=60, cmap=cmocean.cm.balance)
        if i == 5:#6:
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
            cbar = fig.colorbar(cs_all, cax=cbar_ax)
        #     cbar = fig.colorbar(cs, location='right', anchor=(2,0.5), #anchor=(1.8,0.5),
        #                 ax=[ax[0],ax[1],ax[2],ax[3],ax[4],ax[5],ax[6]])
            cbar.ax.tick_params(labelsize=28)#,length=10, width=2)
            cbar.outline.set_visible(False)
        ax.set_title(letters[i] + ' ' + year + ' - avg.', fontsize=24, loc = 'left')

    # format figure
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.axis('off')
    pfun.dar(ax)
    # pfun.add_coast(ax)
                                
# Add colormap title
plt.suptitle('Days with bottom DO < {} [mg/L]'.format(str(DO_thresh)),
            fontsize=36, y=0.95)

# Generate plot
plt.tight_layout
plt.subplots_adjust(left=0.05, right=0.90, top=0.85, wspace=0.02)
plt.savefig(out_dir / ('days_DO_less_than_{}.png'.format(str(DO_thresh))))

##############################################################
##                HYPOXIC AREA TIMESERIES                   ##
##############################################################

# initialize figure
plt.close('all)')
pfun.start_plot(figsize=(12,5))
fig,ax = plt.subplots(1,1)


# create time vector
startdate = '2020.01.01'
enddate = '2020.12.31'
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

# colors = ['red','darkorange','gold','green','blue','purple','deeppink']
colors = ['#F9627D','#62B6CB','#A8C256','#96031A','#957FEF','#476ad1','darkorange']

# # plot timeseries
# for i,year in enumerate(years):
#     # plot hypoxic area timeseries
#     plt.plot(dates_local,hyp_area[year],color=colors[i],
#              linewidth=3,alpha=0.5,label=year)

# # format figure
# ax.grid(visible=True, color='w')
# # format background color
# ax.set_facecolor('#EEEEEE')
# for border in ['top','right','bottom','left']:
#     ax.spines[border].set_visible(False)
# ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
# plt.legend(loc='best')
# plt.title('Bottom Area with DO < {} mg/L '.format(DO_thresh) + r'[km$^2$]')
# ax.set_xlim([dates_local[0],dates_local[-1]])

# plt.savefig(out_dir / 'area_with_DO_lt{}'.format(DO_thresh))

##############################################################
##                HYPOXIC VOLUME TIMESERIES                 ##
##############################################################

# initialize figure
plt.close('all)')
pfun.start_plot(figsize=(12,4))
f, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 4]})

# format figure
ax0.set_xlim([xmin,xmax])
ax0.set_ylim([ymin,ymax])
# ax0.set_yticklabels([])
# ax0.set_xticklabels([])
# ax0.axis('off')
ax0.set_ylabel('Latitude', fontsize=12)
ax0.set_xlabel('Longitude', fontsize=12)
ax0.tick_params(axis='both', labelsize=12)
ax0.pcolormesh(plon, plat, zm, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
pfun.dar(ax0)
# pfun.add_coast(ax0, color='gray')
# Create a Rectangle patch to omit Straits
if remove_straits:
    # get lat and lon
    fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_2013.01.01_2013.12.31.nc')
    ds = xr.open_dataset(fp)
    lons = ds.coords['lon_rho'].values
    lats = ds.coords['lat_rho'].values
    lon = lons[0,:]
    lat = lats[:,0]
    # Straits
    lonmax = -122.76
    lonmin = xmin
    latmax = ymax
    latmin = 48.14
    # convert lat/lon to eta/xi
    diff = np.absolute(lon-lonmin)
    ximin = diff.argmin()
    diff = np.absolute(lon-lonmax)
    ximax = diff.argmin()
    diff = np.absolute(lat-latmin)
    etamin = diff.argmin()
    diff = np.absolute(lat-latmax)
    etamax = diff.argmin()
    rect = patches.Rectangle((lon[ximin], lat[etamin]), lon[ximax]-lon[ximin], lat[etamax]-lat[etamin],
                            edgecolor='none', facecolor='white', alpha=0.9)
    # Add the patch to the Axes
    ax0.add_patch(rect)
    ax0.text(-123.2,48.25,'Straits\nomitted', rotation=90, fontsize=12)
    ax0.set_title('(a) Region', fontsize = 14, loc='left')

# Puget Sound volume
if remove_straits:
    PS_vol = 195.2716230839466 # [km^3]

# create time vector
startdate = '2020.01.01'
enddate = '2020.12.31'
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

# plot timeseries
for i,year in enumerate(years):
    # plot hypoxic area timeseries
    # ax1.plot(dates_local,hyp_vol[year],color=colors[i],
    #          linewidth=3,alpha=0.5,label=year)
    if year == '2017':
        ax1.plot(dates_local,hyp_vol[year],color='white',
                linewidth=4,zorder=4)
        ax1.plot(dates_local,hyp_vol[year],color='black',
                linewidth=3.5,label=year,zorder=4)
    else:
        ax1.plot(dates_local,hyp_vol[year],color='white',
                linewidth=2.5)
        ax1.plot(dates_local,hyp_vol[year],color=colors[i],
                linewidth=2,label=year)

# get median hypoxic volume
med_vol = np.nanmedian(list(hyp_vol.values()), axis=0)
ax1.plot(dates_local,med_vol,color='k',
        linestyle='--',linewidth=2,label='median')

# format figure
# ax1.grid(visible=True, axis='x', color='w')
ax1.grid(visible=True, axis='both', color='silver', linestyle='--')
# format background color
# ax1.set_facecolor('#EEEEEE')
# for border in ['top','right','bottom','left']:
#     ax1.spines[border].set_visible(False)
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax1.tick_params(axis='both', labelsize=12)
ax1.set_ylabel(r'Hypoxic volume [km$^3$]', fontsize=12)
plt.legend(loc='upper left', fontsize=12)
plt.title('(b) Puget Sound hypoxic volume (DO < 2 mg/L)', fontsize = 14, loc='left')
ax1.set_xlim([dates_local[0],dates_local[-1]])
ax1.set_ylim([0,13])

# create hypoxic volume y-axis
# convert hypoxic volume to percent hypoxic volume
percent = lambda hyp_vol: hyp_vol/PS_vol*100
# get left axis limits
ymin, ymax = ax1.get_ylim()
# match ticks
ax2 = ax1.twinx()
ax2.set_ylim((percent(ymin),percent(ymax)))
ax2.plot([],[])
for border in ['top','right','bottom','left']:
    ax2.spines[border].set_visible(False)
ax2.set_ylabel(r'Percent of regional volume [%]', fontsize=12)

plt.savefig(out_dir / ('volume_with_DO_lt2_'+straits+'.png'))