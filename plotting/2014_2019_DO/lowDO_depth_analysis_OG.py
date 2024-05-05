"""
Analyzes 2013 run for relationship between watercolumn depth and hypoxia

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
import matplotlib as mpl
from matplotlib.dates import MonthLocator
from matplotlib import cm
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import pandas as pd
import cmocean
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
from matplotlib.colors import Normalize 
from scipy.interpolate import interpn
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

DO_thresh = 2 # mg/L DO threshold
remove_straits = True

year = '2013'

# Variables to compare
vn = 'oxygen'

# Provide information about models to compare
# natural
natural_gtagex = 'cas7_t0noN_x4b'
# anthropogenic (long hindcast)
anthropogenic_gtagex = 'cas7_t0_x4b'

# Show WWTP locations?
WWTP_loc = True

###########################################################
# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots' / 'DO_quantification'
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

##############################################################
##                      PROCESS DATA                        ##
##############################################################

print('Processing started...')

# get data
fp_natural = Ldir['LOo'] / 'extract' / natural_gtagex / 'box' / 'pugetsoundDO_2013.01.01_2013.12.31.nc'
ds_natural = xr.open_dataset(fp_natural)
fp_anthropogenic = Ldir['LOo'] / 'extract' / anthropogenic_gtagex / 'box' / 'pugetsoundDO_2013.01.01_2013.12.31.nc'
ds_anthropogenic = xr.open_dataset(fp_anthropogenic)

if remove_straits:
    # set values in strait of juan de fuca and strait of georgia are nan (so they don't interfere with analysis)
    lat_threshold = 48.14
    lon_threshold = -122.76
    # Create a mask for latitudes and longitudes that meet the condition
    mask_natural = (ds_natural['lat_rho'] > lat_threshold) & (ds_natural['lon_rho'] < lon_threshold)
    mask_anthropogenic = (ds_anthropogenic['lat_rho'] > lat_threshold) & (ds_anthropogenic['lon_rho'] < lon_threshold)
    # Expand mask dimensions to match 'oxygen' dimensions
    expanded_mask_natural = mask_natural.expand_dims(ocean_time=len(ds_natural['oxygen']), s_rho=len(ds_natural['s_rho']))
    expanded_mask_anthropogenic = mask_anthropogenic.expand_dims(ocean_time=len(ds_anthropogenic['oxygen']), s_rho=len(ds_anthropogenic['s_rho']))
    # Apply the mask to the 'oxygen' variable
    ds_natural['oxygen'] = xr.where(expanded_mask_natural, np.nan, ds_natural['oxygen'])
    ds_anthropogenic['oxygen'] = xr.where(expanded_mask_anthropogenic, np.nan, ds_anthropogenic['oxygen'])

# get s-rho of the lowest DO (array with dimensions of (ocean_time: 365, eta_rho: 441, xi_rho: 177))
srho_min_natural = ds_natural['oxygen'].idxmin(dim='s_rho', skipna=True).values
srho_min_anthropogenic = ds_anthropogenic['oxygen'].idxmin(dim='s_rho', skipna=True).values

# get depths
depths = ds_natural['h'].values
# reshape
depths_reshape = depths.reshape((1,441,177))

# convert srho to depths
# natural
depth_min_natural = depths_reshape * srho_min_natural
# anthropogenic
depth_min_anthropogenic = depths_reshape * srho_min_anthropogenic

# # add new dimension with s-levels
s_level = np.linspace(0,29,30)
# natural
ds_natural['oxygen'] = ds_natural['oxygen'].assign_coords({'s_level': ('s_rho',s_level)})
ds_natural['oxygen'] = ds_natural['oxygen'].swap_dims({'s_rho': 's_level'})
# anthropogenic
ds_anthropogenic['oxygen'] = ds_anthropogenic['oxygen'].assign_coords({'s_level': ('s_rho',s_level)})
ds_anthropogenic['oxygen'] = ds_anthropogenic['oxygen'].swap_dims({'s_rho': 's_level'})

# get s-level of the lowest DO (array with dimensions of (ocean_time: 365, eta_rho: 441, xi_rho: 177))
slev_min_natural = ds_natural['oxygen'].idxmin(dim='s_level', skipna=True).values
slev_min_anthropogenic = ds_anthropogenic['oxygen'].idxmin(dim='s_level', skipna=True).values

# get corresponding DO minima concentration (mg/L)
DO_min_natural = pinfo.fac_dict['oxygen'] * ds_natural['oxygen'].min(dim='s_level', skipna=True).values
DO_min_anthropogenic = pinfo.fac_dict['oxygen'] * ds_anthropogenic['oxygen'].min(dim='s_level', skipna=True).values

# print(s_min_natural.shape)
# print(DO_min_natural.shape)


#############################################################
##      S-LEVEL VS. NATURAL DO (MIN OF WATER COLUMN)       ## (2D histogram)
#############################################################


# compress spatial and time dimensions
# natural
slev_min_all_natural = np.reshape(slev_min_natural,-1)
DO_min_all_natural = np.reshape(DO_min_natural,-1)
# anthropogenic
slev_min_all_anthropogenic = np.reshape(slev_min_anthropogenic,-1)
DO_min_all_anthropogenic = np.reshape(DO_min_anthropogenic,-1)

# rename variables so its easier to manipulate
# natural
x_natural = DO_min_all_natural
y_natural = slev_min_all_natural
# anthropogenic
x_anthropogenic = DO_min_all_anthropogenic
y_anthropogenic = slev_min_all_anthropogenic
# get rid of nans in dataset
# natural
bad_indices = np.isnan(x_natural) | np.isnan(y_natural)
good_indices = ~bad_indices
good_x_natural = x_natural[good_indices]
good_y_natural = y_natural[good_indices]
# anthropogenic
bad_indices = np.isnan(x_anthropogenic) | np.isnan(y_anthropogenic)
good_indices = ~bad_indices
good_x_anthropogenic = x_anthropogenic[good_indices]
good_y_anthropogenic = y_anthropogenic[good_indices]

# calculate difference between two histograms
# natural
plt.figure(1)
h, xedges, yedges, image0 = plt.hist2d(good_x_natural,good_y_natural, bins=(100, 30), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
# anthropogenic
plt.figure(2)
h1, xedges1, yedges1, image1 = plt.hist2d(good_x_anthropogenic,good_y_anthropogenic, bins=(xedges, yedges), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)

# initialize figure
plt.close('all')
pfun.start_plot(figsize=(17,8))
fig,ax = plt.subplots(1,2, sharey = True)

# plot 2d histogram to get colored scatter by point density (natural)
cs = ax[0].hist2d(good_x_natural, good_y_natural, bins = [100,30], norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
cbar = fig.colorbar(cs[3])
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)

# plot difference between colormaps
# anthropogenic minus natural
vmin = np.min(h1-h)
vmax = np.max(h1-h)
cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)
cs = ax[1].pcolormesh(xedges, yedges, (h1-h).T, cmap=cmap)
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)


# format natural figure
ax[0].set_xlabel('DO minima concentration [mg/L]')
ax[0].set_title('(a) Natural')
ax[0].grid(visible=True, color='w')
# add labels
ax[0].set_ylabel('S-Level')
# format background color
ax[0].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[0].spines[border].set_visible(False)

# format difference figure
ax[1].set_xlabel('DO minima concentration [mg/L]')
ax[1].set_title('(b) Anthropogenic - Natural')
ax[1].grid(visible=True, color='w')
# add labels
ax[0].set_ylabel('S-Level')
# format background color
ax[1].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[1].spines[border].set_visible(False)


plt.tight_layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.80, wspace=0.02)
plt.suptitle(r'2D Histogram: S-Level vs. DO Minima in Water Column' + '\n for all cells and days in Puget Sound')
plt.savefig(out_dir / 'slevel_vs_minDO')
plt.close('all')

#############################################################
##        DEPTH VS. NATURAL DO (MIN OF WATER COLUMN)       ## (2D histogram)
#############################################################


# compress spatial and time dimensions
# natural
depth_min_all_natural = np.reshape(depth_min_natural,-1)
# anthropogenic
depth_min_all_anthropogenic = np.reshape(depth_min_anthropogenic,-1)

# rename variables so its easier to manipulate
# natural
x_natural = DO_min_all_natural
y_natural = depth_min_all_natural
# anthropogenic
x_anthropogenic = DO_min_all_anthropogenic
y_anthropogenic = depth_min_all_anthropogenic
# get rid of nans in dataset
# natural
bad_indices = np.isnan(x_natural) | np.isnan(y_natural)
good_indices = ~bad_indices
good_x_natural = x_natural[good_indices]
good_y_natural = y_natural[good_indices]
# anthropogenic
bad_indices = np.isnan(x_anthropogenic) | np.isnan(y_anthropogenic)
good_indices = ~bad_indices
good_x_anthropogenic = x_anthropogenic[good_indices]
good_y_anthropogenic = y_anthropogenic[good_indices]

# calculate difference between two histograms
# natural
plt.figure(1)
h, xedges, yedges, image0 = plt.hist2d(good_x_natural,good_y_natural, bins=(100, 100), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
# anthropogenic
plt.figure(2)
h1, xedges1, yedges1, image1 = plt.hist2d(good_x_anthropogenic,good_y_anthropogenic, bins=(xedges, yedges), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)

# initialize figure
plt.close('all')
pfun.start_plot(figsize=(17,8))
fig,ax = plt.subplots(1,2, sharey = True)

# plot 2d histogram to get colored scatter by point density (natural)
cs = ax[0].hist2d(good_x_natural, good_y_natural, bins = [100,100], norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
cbar = fig.colorbar(cs[3])
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)

# plot difference between colormaps
# anthropogenic minus natural
vmin = np.min(h1-h)
vmax = np.max(h1-h)
cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)
cs = ax[1].pcolormesh(xedges, yedges, (h1-h).T, cmap=cmap)
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)

# format natural figure
ax[0].set_xlabel('DO minima concentration [mg/L]')
ax[0].set_title('(a) Natural')
ax[0].grid(visible=True, color='w')
# add labels
ax[0].set_ylabel('Depth [m]')
# format background color
ax[0].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[0].spines[border].set_visible(False)

# format difference figure
ax[1].set_xlabel('DO minima concentration [mg/L]')
ax[1].set_title('(b) Anthropogenic - Natural')
ax[1].grid(visible=True, color='w')
# format background color
ax[1].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[1].spines[border].set_visible(False)


plt.tight_layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.80, wspace=0.02)
plt.suptitle(r'2D Histogram: Depth vs. DO Minima in Water Column' + '\n for all cells and days in Puget Sound')
plt.savefig(out_dir / 'depth_vs_minDO')
plt.close('all')

#############################################################
##                   S-LEVEL VS. YEARDAY                   ## (2D histogram)
#############################################################

# create time vector
startdate = '2013.01.01'
enddate = '2013.12.31'
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]
# duplicate time dimension
day = np.linspace(1,365,365)
dates_timeonly = np.repeat(day[:, np.newaxis], 78057, axis=1)

# compress dimensions, but keep time dimension (lost spatial resolution)
# natural
slev_min_timeonly_natural = np.reshape(slev_min_natural,(365,78057))
# anthropogenic
slev_min_timeonly_anthropogenic = np.reshape(slev_min_anthropogenic,(365,78057))


# rename variables so its easier to manipulate
# natural
x_natural = dates_timeonly
y_natural = slev_min_timeonly_natural
# anthropogenic
x_anthropogenic = dates_timeonly
y_anthropogenic = slev_min_timeonly_anthropogenic
# get rid of nans in dataset
# natural
bad_indices = np.isnan(x_natural) | np.isnan(y_natural)
good_indices = ~bad_indices
good_x_natural = x_natural[good_indices]
good_y_natural = y_natural[good_indices]
# anthropogenic
bad_indices = np.isnan(x_anthropogenic) | np.isnan(y_anthropogenic)
good_indices = ~bad_indices
good_x_anthropogenic = x_anthropogenic[good_indices]
good_y_anthropogenic = y_anthropogenic[good_indices]

# calculate difference between two histograms
# natural
plt.figure(1)
h, xedges, yedges, image0 = plt.hist2d(good_x_natural,good_y_natural, bins=(365, 30), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
# anthropogenic
plt.figure(2)
h1, xedges1, yedges1, image1 = plt.hist2d(good_x_anthropogenic,good_y_anthropogenic, bins=(xedges, yedges), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)

# initialize figure
plt.close('all')
pfun.start_plot(figsize=(17,8))
fig,ax = plt.subplots(1,2, sharey = True)

# plot 2d histogram to get colored scatter by point density (natural)
cs = ax[0].hist2d(good_x_natural, good_y_natural, bins = [365,30], norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
cbar = fig.colorbar(cs[3])
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)

# plot difference between colormaps
# anthropogenic minus natural
vmin = np.min(h1-h)
vmax = np.max(h1-h)
cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)
cs = ax[1].pcolormesh(xedges, yedges, (h1-h).T, cmap=cmap)
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)


# format natural figure
ax[0].set_xlabel('Year day')
ax[0].set_title('(a) Natural')
ax[0].grid(visible=True, color='w')
# add labels
ax[0].set_ylabel('S-Level')
# format background color
ax[0].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[0].spines[border].set_visible(False)

# format difference figure
ax[1].set_xlabel('Year day')
ax[1].set_title('(b) Anthropogenic - Natural')
ax[1].grid(visible=True, color='w')
# add labels
ax[0].set_ylabel('S-Level')
# format background color
ax[1].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[1].spines[border].set_visible(False)


plt.tight_layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.80, wspace=0.02)
plt.suptitle(r'2D Histogram: S-Level of DO minima vs. day of year' + '\n for all cells in Puget Sound')
plt.savefig(out_dir / 'slevel_vs_time')
plt.close('all')

#############################################################
##                    DEPTH VS. YEARDAY                    ## (2D histogram)
#############################################################

# create time vector
startdate = '2013.01.01'
enddate = '2013.12.31'
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]
# duplicate time dimension
day = np.linspace(1,365,365)
dates_timeonly = np.repeat(day[:, np.newaxis], 78057, axis=1)

# compress dimensions, but keep time dimension (lost spatial resolution)
# natural
depth_min_timeonly_natural = np.reshape(depth_min_natural,(365,78057))
# anthropogenic
depth_min_timeonly_anthropogenic = np.reshape(depth_min_anthropogenic,(365,78057))


# rename variables so its easier to manipulate
# natural
x_natural = dates_timeonly
y_natural = depth_min_timeonly_natural
# anthropogenic
x_anthropogenic = dates_timeonly
y_anthropogenic = depth_min_timeonly_anthropogenic
# get rid of nans in dataset
# natural
bad_indices = np.isnan(x_natural) | np.isnan(y_natural)
good_indices = ~bad_indices
good_x_natural = x_natural[good_indices]
good_y_natural = y_natural[good_indices]
# anthropogenic
bad_indices = np.isnan(x_anthropogenic) | np.isnan(y_anthropogenic)
good_indices = ~bad_indices
good_x_anthropogenic = x_anthropogenic[good_indices]
good_y_anthropogenic = y_anthropogenic[good_indices]

# calculate difference between two histograms
# natural
plt.figure(1)
h, xedges, yedges, image0 = plt.hist2d(good_x_natural,good_y_natural, bins=(365, 100), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
# anthropogenic
plt.figure(2)
h1, xedges1, yedges1, image1 = plt.hist2d(good_x_anthropogenic,good_y_anthropogenic, bins=(xedges, yedges), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)

# initialize figure
plt.close('all')
pfun.start_plot(figsize=(17,8))
fig,ax = plt.subplots(1,2, sharey = True)

# plot 2d histogram to get colored scatter by point density (natural)
cs = ax[0].hist2d(good_x_natural, good_y_natural, bins = [365,100], norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
cbar = fig.colorbar(cs[3])
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)

# plot difference between colormaps
# anthropogenic minus natural
vmin = np.min(h1-h)
vmax = np.max(h1-h)
cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)
cs = ax[1].pcolormesh(xedges, yedges, (h1-h).T, cmap=cmap)
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)


# format natural figure
ax[0].set_xlabel('Year day')
ax[0].set_title('(a) Natural')
ax[0].grid(visible=True, color='w')
# add labels
ax[0].set_ylabel('Depth [m]')
# format background color
ax[0].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[0].spines[border].set_visible(False)

# format difference figure
ax[1].set_xlabel('Year day')
ax[1].set_title('(b) Anthropogenic - Natural')
ax[1].grid(visible=True, color='w')
# format background color
ax[1].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[1].spines[border].set_visible(False)


plt.tight_layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.80, wspace=0.02)
plt.suptitle(r'2D Histogram: Depth of DO minima vs. day of year' + '\n for all cells in Puget Sound')
plt.savefig(out_dir / 'depth_vs_time')
plt.close('all')

#############################################################
##                 DO MINIMA VS. YEARDAY                   ## (2D histogram)
#############################################################


# compress dimensions, but keep time dimension (lost spatial resolution)
DO_min_all_natural = np.reshape(DO_min_natural,(365,78057))
# anthropogenic
DO_min_all_anthropogenic = np.reshape(DO_min_anthropogenic,(365,78057))

# rename variables so its easier to manipulate
# natural
x_natural = dates_timeonly
y_natural = DO_min_all_natural
# anthropogenic
x_anthropogenic = dates_timeonly
y_anthropogenic = DO_min_all_anthropogenic
# get rid of nans in dataset
# natural
bad_indices = np.isnan(x_natural) | np.isnan(y_natural)
good_indices = ~bad_indices
good_x_natural = x_natural[good_indices]
good_y_natural = y_natural[good_indices]
# anthropogenic
bad_indices = np.isnan(x_anthropogenic) | np.isnan(y_anthropogenic)
good_indices = ~bad_indices
good_x_anthropogenic = x_anthropogenic[good_indices]
good_y_anthropogenic = y_anthropogenic[good_indices]

# calculate difference between two histograms
# natural
plt.figure(1)
h, xedges, yedges, image0 = plt.hist2d(good_x_natural,good_y_natural, bins=(365, 100), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
# anthropogenic
plt.figure(2)
h1, xedges1, yedges1, image1 = plt.hist2d(good_x_anthropogenic,good_y_anthropogenic, bins=(xedges, yedges), norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)

# initialize figure
plt.close('all')
pfun.start_plot(figsize=(17,8))
fig,ax = plt.subplots(1,2, sharey = True)

# plot 2d histogram to get colored scatter by point density (natural)
cs = ax[0].hist2d(good_x_natural, good_y_natural, bins = [365,100], norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
cbar = fig.colorbar(cs[3])
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)

# plot difference between colormaps
# anthropogenic minus natural
vmin = np.min(h1-h)
vmax = np.max(h1-h)
cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)
cs = ax[1].pcolormesh(xedges, yedges, (h1-h).T, cmap=cmap)
cbar = fig.colorbar(cs)
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)


# format natural figure
ax[0].set_xlabel('Year day')
ax[0].set_title('(a) Natural')
ax[0].grid(visible=True, color='w')
# add labels
ax[0].set_ylabel('DO minima concentration [mg/L]')
# format background color
ax[0].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[0].spines[border].set_visible(False)

# format difference figure
ax[1].set_xlabel('Year day')
ax[1].set_title('(b) Anthropogenic - Natural')
ax[1].grid(visible=True, color='w')
# format background color
ax[1].set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax[1].spines[border].set_visible(False)


plt.tight_layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.80, wspace=0.02)
plt.suptitle(r'2D Histogram: Concentration of DO minima vs. day of year' + '\n for all cells in Puget Sound')
plt.savefig(out_dir / 'minDO_vs_time')
plt.close('all')

print('Done.')
