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
pth = Path(__file__).absolute().parent.parent.parent.parent.parent / 'LO' / 'pgrid'
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

years =  ['2015','2016','2017','2018','2019','2020']

# which  model run to look at?
gtagex = 'cas7_t1_x11ab' # long hindcast (anthropogenic)

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
    repeatrivs_fn = '../../../../LO_data/trapsD00/LiveOcean_SSM_rivers.xlsx'
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
    ds = xr.open_dataset('../../../../LO_data/grids/cas7/grid.nc')
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
    fp_wwtps = '../../../../LO_output/pre/trapsP00/point_sources/lo_base/Data_historical/'
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
    griddf0_wwtps = pd.read_csv('../../../../LO_data/grids/cas6/wwtp_info.csv')
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

# get the grid data
ds = xr.open_dataset('../../../../LO_data/grids/cas7/grid.nc')
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

# read in masks
basin_mask_ds = xr.open_dataset('../../../../LO_output/chapter_2/data/basin_masks_from_pugetsoundDObox.nc')
mask_rho = basin_mask_ds.mask_rho.values
mask_ps = basin_mask_ds.mask_pugetsound.values

# initialize empty dictionaries
DO_bot_dict = {} # dictionary with DO_bot values
hyp_thick_dict = {}

for year in years:
    # add ds to dictionary
    # ds = xr.open_dataset(Ldir['LOo'] / 'pugetsound_DO' / 'data' / (year + '_DO_info_' + straits + '.nc'))
    ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / (gtagex + '_pugetsoundDO_' + year + '_DO_info.nc'))
    DO_bot = ds['DO_bot'].values
    hyp_thick = ds['hyp_thick'].values
    # if not a leap year, add a nan on feb 29 (julian day 60 - 1 because indexing from 0)
    if np.mod(int(year),4) != 0: 
        DO_bot = np.insert(DO_bot,59,'nan',axis=0)
        hyp_thick = np.insert(hyp_thick,59,'nan',axis=0)
    DO_bot_dict[year] = DO_bot  * mask_ps
    hyp_thick_dict[year] = hyp_thick * mask_ps

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

# get grid cell area
DX = (basin_mask_ds.pm.values)**-1
DY = (basin_mask_ds.pn.values)**-1
DA = DX*DY*(1/1000)*(1/1000) # get area, but convert from m^2 to km^2

# initialize dictionary for hypoxic volume [km3]
hyp_vol = {}
for year in years:
    # get hypoxic thickness
    hyp_thick = hyp_thick_dict[year]/1000 # [km]
    hyp_vol_timeseries = np.sum(hyp_thick * DA, axis=(1, 2)) # km^3
    hyp_vol[year] = hyp_vol_timeseries

# get plotting limits based on region
if region == 'Puget Sound':
    # box extracion limits: [-123.29, -122.1, 46.95, 48.93]
    xmin = -123.29
    xmax = -122.1
    ymin = 46.95
    # ymax = 48.93
    ymax=48.5

colors = ['#62B6CB','#A8C256','#96031A','#957FEF','#F9627D','#476ad1','darkorange']

##############################################################
##                HYPOXIC VOLUME TIMESERIES                 ##
##############################################################

# initialize figure
plt.close('all)')
f, (ax0, ax1) = plt.subplots(1,2,figsize = (12,5.5),gridspec_kw={'width_ratios': [1, 3]})

# format figure
ax0.set_xlim([xmin,xmax])
ax0.set_ylim([ymin,ymax])
# ax0.set_yticklabels([])
# ax0.set_xticklabels([])
# ax0.axis('off')
ax0.set_ylabel('Latitude', fontsize=12)
ax0.set_xlabel('Longitude', fontsize=12)
ax0.tick_params(axis='both', labelsize=12)
# ax0.pcolormesh(plon, plat, zm, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
lon = basin_mask_ds.lon_rho.values
lat = basin_mask_ds.lat_rho.values
plon, plat = pfun.get_plon_plat(lon,lat)
ax0.pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho),
                vmin=0, vmax=10, cmap='Blues')
ax0.pcolormesh(plon, plat, np.where(mask_ps == 0, np.nan, mask_ps),
                vmin=0, vmax=2, cmap='Blues')
pfun.dar(ax0)

ax0.set_title('(a)', fontsize = 14, loc='left', fontweight='bold')

# Puget Sound volume
# Puget Sound volume
PS_vol = np.nansum(basin_mask_ds['h'].values/1000 * DA * mask_ps) # [km^3]
print('Puget Sound volume: {} km3'.format(round(PS_vol,1)))

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
# med_vol = np.nanmedian(list(hyp_vol.values()), axis=0)
# ax1.plot(dates_local,med_vol,color='k',
#         linestyle='--',linewidth=2,label='median')

mean_vol = np.nanmean(list(hyp_vol.values()), axis=0)
ax1.plot(dates_local,mean_vol,color='k',
        linestyle='--',linewidth=2,label='mean')

# format figure
ax1.grid(visible=True, axis='both', color='silver', linestyle='--')
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax1.tick_params(axis='both', labelsize=12)
ax1.set_ylabel(r'Hypoxic volume [km$^3$]', fontsize=12)
plt.legend(loc='upper left', fontsize=12)
plt.title('(b)', fontsize = 14, loc='left', fontweight='bold')
ax1.set_xlim([dates_local[0],dates_local[-1]])
ax1.set_ylim([0,6])

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

# save to csv file
hyp_vol_df = pd.DataFrame.from_dict(hyp_vol)
hyp_vol_df.insert(0, 'yearday', list(range(1,367)))
hyp_vol_df.to_csv('../../../../terminal_inlet_DO_rev2/daily_hypoxic_volume_km3.csv', index=False)
print(hyp_vol_df)