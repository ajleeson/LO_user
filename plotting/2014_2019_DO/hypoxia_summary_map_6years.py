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

remove_straits = False

vn = 'oxygen'

years =  ['2014','2015','2016','2017','2018','2019']

# which  model run to look at?
gtagex = 'cas7_t0_x4b' # long hindcast (anthropogenic)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures'
Lfun.make_dir(out_dir)

region = 'Puget Sound'

# start date
start = '08-01'
end = '09-30'

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
##               PLOT DAYS WITH BOTTOM HYPOXIA              ##
##############################################################

cbar_pad = 0.02

# get plotting limits based on region
if region == 'Puget Sound':
    # box extracion limits: [-123.29, -122.1, 46.95, 48.93]
    xmin = -123.29
    xmax = -122.1
    ymin = 46.95
    ymax = 48.93

# Initialize figure
fs = 10
pfun.start_plot(fs=16, figsize=(7,9))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)


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

# Create map
plt.close('all')
fig = plt.figure(figsize=(5,9))
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-1.5, vmax=0, cmap=plt.get_cmap('Greys'))

# plot hypoxic days

# get lat and lon
fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_2013.01.01_2013.12.31.nc')
ds = xr.open_dataset(fp)
lons = ds.coords['lon_rho'].values
lats = ds.coords['lat_rho'].values
px, py = pfun.get_plon_plat(lons,lats)

# get days
DO_days = DO_days['avg']
            
cs = ax.pcolormesh(px,py,DO_days, vmin=0, vmax=np.nanmax(DO_days), cmap='rainbow')
cbar = fig.colorbar(cs, location='left', pad=0.05)
cbar.ax.tick_params(labelsize=14)#,length=10, width=2)
cbar.outline.set_visible(False)
# ax.set_title('six-year avg.', fontsize=14, loc='left')

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
lat0 = 47#46.94
lon0 = -123.05 + 0.7
lat1 = lat0
lon1 = -122.91825 + 0.7
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=2)
ax.text(lon0-0.04,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=12)

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
            fontsize=16, fontweight='bold', y=0.95)

# Generate plot
plt.tight_layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.05, wspace=0.02)
plt.show()

##############################################################
##             PLOT MEAN BOTTOM DO CONCENTRATION            ##
##############################################################

cbar_pad = 0.02

# get plotting limits based on region
if region == 'Puget Sound':
    # box extracion limits: [-123.29, -122.1, 46.95, 48.93]
    xmin = -123.29
    xmax = -122.1
    ymin = 46.95
    ymax = 48.93

# Initialize figure
fs = 10
pfun.start_plot(fs=16, figsize=(5,9))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

# initialize empty dictionary
ds_dict = {}
# add ds to dictionary
for year in years:
    ds = xr.open_dataset(Ldir['LOo'] / 'pugetsound_DO' / 'data' / (year + '_DO_info_' + straits + '.nc'))
    ds_dict[year] = ds

# store all values in new dictionary (the averages over hypoxic season)
val_dict = {}
# plot average year
for year in years:
    ds = ds_dict[year]
    # crop to just hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))
    v = ds['DO_bot'].values
    # take average over season
    v = np.nanmean(v,axis=0)
    val_dict[year] = v

cmap = plt.cm.get_cmap('rainbow_r', 10)
vmin = 0
vmax = 10
# calculate average of all of the arrays
v_avg = sum(val_dict.values())/len(val_dict)
cs = ax.pcolormesh(px,py,v_avg, vmin=vmin, vmax=vmax, cmap=cmap)
cbar = fig.colorbar(cs, location='left',pad=0.05)
cbar.ax.tick_params(labelsize=15)#,length=10, width=2)
cbar.outline.set_visible(False)
            

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
lat0 = 47#46.94
lon0 = -123.05 + 0.7
lat1 = lat0
lon1 = -122.91825 + 0.7
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=2)
ax.text(lon0-0.04,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=12)

# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.axis('off')
pfun.dar(ax)
# pfun.add_coast(ax)
                                
# Add colormap title
plt.suptitle(start+' to '+end+' mean bottom DO [mg/L]',
            fontsize=16, fontweight='bold', y=0.95)

# Generate plot
plt.tight_layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.05, wspace=0.02)
plt.show()
