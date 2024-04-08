"""
Quantifies DO differences between anthropogenic and natural run.

This is a custom function for a particular experiment,
but code can be adapted for other use cases in the future.

From the terminal: DO_diff_quantification.py

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

DO_thresh = 6 # mg/L DO threshold

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

# get bottom DO  (in mg/L)
DO_bot_natural = pinfo.fac_dict[vn] * ds_natural[vn][:,0,:,:].values # s_rho = 0 for bottom
DO_bot_anthropogenic = pinfo.fac_dict[vn] * ds_anthropogenic[vn][:,0,:,:].values # s_rho = 0 for bottom

# Get boolean array. True if DO < 2 mg/L, otherwise, nan
DO_bot_lt2_natural = np.where(DO_bot_natural < DO_thresh, 1, np.nan)
DO_bot_lt2_anthropogenic = np.where(DO_bot_anthropogenic < DO_thresh, 1, np.nan)

# Sum over time to compress into just lat/lon dimension, with values indicating days with bottom DO < 2mg/L
DO_days_lt2_natural = np.nansum(DO_bot_lt2_natural, axis=0)
DO_days_lt2_anthropogenic = np.nansum(DO_bot_lt2_anthropogenic, axis=0)
DO_days_lt2_diff = DO_days_lt2_anthropogenic - DO_days_lt2_natural
# get convert zeros to nans in dataset
DO_days_lt2_natural = np.where(DO_days_lt2_natural==0, np.nan, DO_days_lt2_natural)
DO_days_lt2_anthropogenic = np.where(DO_days_lt2_anthropogenic==0, np.nan, DO_days_lt2_anthropogenic)
DO_days_lt2_diff = np.where(DO_days_lt2_diff==0, np.nan, DO_days_lt2_diff)

# calculate bottom hypoxic area
hyp_area_natural = 0.5 * 0.5 * DO_bot_lt2_natural # area of grid cell = 05 km by 0.5 km
hyp_area_anthropogenic = 0.5 * 0.5 * DO_bot_lt2_anthropogenic # area of grid cell = 05 km by 0.5 km

# get timeseries of bottom hypoxia area
hyp_area_timeseries_natural = np.nansum(hyp_area_natural,axis=1)
hyp_area_timeseries_natural = np.nansum(hyp_area_timeseries_natural,axis=1)
hyp_area_timeseries_anthropogenic = np.nansum(hyp_area_anthropogenic,axis=1)
hyp_area_timeseries_anthropogenic = np.nansum(hyp_area_timeseries_anthropogenic,axis=1)

print('Data processing done...')

##############################################################
##                     PLOTTING DO MAPS                     ##
##############################################################

# get plotting limits of box extraction
xmin = -123.29
xmax = -122.1
ymin = 46.95
ymax = 48.93

# get lat/lon
lons = ds_natural.coords['lon_rho'].values
lats = ds_natural.coords['lat_rho'].values
px, py = pfun.get_plon_plat(lons,lats)

# initialize figure
plt.close('all')
pfun.start_plot(figsize=(32,27))
fig = plt.figure()
gs = fig.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.95, wspace=0.05, hspace=0.05)

# plot natural run
ax = fig.add_subplot(1,2,1)
cs = ax.pcolormesh(px,py,DO_days_lt2_natural, cmap='rainbow')
# add colorbar
cbar = fig.colorbar(cs, location='left')
cbar.ax.tick_params(labelsize=32)#,length=10, width=2)
cbar.outline.set_visible(False)
# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
# ax.axis('off')
pfun.add_coast(ax)
pfun.dar(ax)
ax.set_title('(a) Natural', fontsize=38)

vmax = int(np.nanmax(DO_days_lt2_diff))

# plot anthropogenic minus natural run
ax = fig.add_subplot(1,2,2)
cs = ax.pcolormesh(px,py,DO_days_lt2_diff, vmin=-1, vmax=vmax,
                   cmap=cmocean.tools.crop(cmocean.cm.balance_r, -1, vmax, 0))
# add colorbar
cbar = fig.colorbar(cs, location='right')
cbar.ax.tick_params(labelsize=32)#,length=10, width=2)
cbar.outline.set_visible(False)
# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
# ax.axis('off')
pfun.add_coast(ax)
pfun.dar(ax)
ax.set_title('(b) Anthropogenic - Natural', fontsize=38)

# Add title
plt.suptitle('Days with Bottom DO < {} mg/L'.format(DO_thresh),
            fontsize=44, fontweight='bold', y=0.95)

# save figure
plt.savefig(out_dir / 'map_days_do_lt_{}'.format(DO_thresh))

print('DO maps done...')

##############################################################
##           PLOTTING HYPOXIC AREA TIMESERIES               ##
##############################################################

# initialize figure
plt.close('all)')
pfun.start_plot(figsize=(12,5))
fig,ax = plt.subplots(1,1)


# create time vector
startdate = '2013.01.01'
enddate = '2013.12.31'
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

# plot hypoxic area timeseries
plt.plot(dates_local,hyp_area_timeseries_natural,color='darkturquoise',linewidth=3,label='Natural')
plt.plot(dates_local,hyp_area_timeseries_anthropogenic,color='k',linewidth=1,linestyle='--',label='Anthropogenic')
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
plt.legend(loc='best')
plt.title('Bottom Area with DO < {} mg/L '.format(DO_thresh) + r'[km$^2$]')
ax.set_xlim([dates_local[0],dates_local[-1]])

plt.savefig(out_dir / 'timeseries_area_with_DO_lt{}'.format(DO_thresh))

print('Hypoxic time series done...')


# #############################################################
# ##              CHANGE IN DO VS. NATURAL DO                ## (with histograms)
# #############################################################

# # initialize figure
# plt.close('all)')
# pfun.start_plot(figsize=(10,8))
# fig = plt.figure()

# # set up main scatter and histogram margin
# gs = GridSpec(4,4)
# ax_main = fig.add_subplot(gs[1:4,0:3])
# ax_hist_right = fig.add_subplot(gs[1:4,3])
# ax_hist_top = fig.add_subplot(gs[0,0:3])

# # get points with DO < threshold, and get difference between model runs
# DO_bot_lt6_natural = np.where(DO_bot_natural < DO_thresh, DO_bot_natural, np.nan)
# DO_bot_lt6_anthropogenic = np.where(DO_bot_anthropogenic < DO_thresh, DO_bot_anthropogenic, np.nan)
# # compress spatial and time dimensions
# all_bot_DO_natural = np.reshape(DO_bot_lt6_natural,(28490805))
# all_bot_DO_anthropogenic = np.reshape(DO_bot_lt6_anthropogenic,(28490805))
# all_bot_DO_diff = all_bot_DO_anthropogenic - all_bot_DO_natural

# # colors
# colors = np.where(all_bot_DO_diff<0,'firebrick','royalblue')
# # add points
# ax_main.scatter(all_bot_DO_natural,all_bot_DO_diff,color=colors,
#          marker='o',alpha=0.1, edgecolor='none')
# # add y = -1 line
# ax_main.fill_between([0,0.8],[0,-0.8], -0.8, color='white')
# # histogram
# ax_hist_right.hist(all_bot_DO_diff,orientation='horizontal',bins=120, # 120 bins means each bid width is 0.01 mg/L
#              color='silver', log=True) 
# ax_hist_top.hist(all_bot_DO_natural,bins=120, # 120 bins
#              color='silver', log=True) 

# # add text
# ax_main.text(0.12, 0.1, 'slope = -1', horizontalalignment='left',
#      verticalalignment='center', transform=ax_main.transAxes, fontsize=12)
# ax_main.text(0.35, 0.25, 'Anthropogenic DO Lower', horizontalalignment='left',
#         color='firebrick',verticalalignment='center', fontweight='bold', fontsize=12)
# ax_main.text(0.35, 0.3, 'Anthropogenic DO Higher', horizontalalignment='left',
#         color='royalblue',verticalalignment='center', fontweight='bold', fontsize=12)

# # format figure
# plt.setp(ax_hist_right.get_yticklabels(), visible=False)
# plt.setp(ax_hist_top.get_xticklabels(), visible=False)
# ax_main.set_xlabel('Natural DO [mg/L]')
# # format grid and labels
# ax_main.set_xlim([0.01,DO_thresh])
# ax_main.set_ylim([-0.8,0.4])
# ax_hist_right.set_ylim([-0.8,0.4])
# ax_hist_top.set_xlim([0.01,DO_thresh])
# ax_main.grid(visible=True, color='w')
# ax_hist_right.grid(visible=True, color='w')
# ax_hist_top.grid(visible=True, color='w')
# # put grid behind points and histogram
# ax_main.set_axisbelow(True)
# ax_hist_right.set_axisbelow(True)
# ax_hist_top.set_axisbelow(True)
# # add labels
# ax_main.set_ylabel('Anthropogenic - Natural [mg/L]')
# ax_hist_right.set_xlabel('Count')
# ax_hist_top.set_ylabel('Count')
# plt.suptitle(r'$\Delta$ DO vs. Natural DO (Bottom)' + '\n for all cells/days in Puget Sound where DO < {} mg/L'.format(DO_thresh))
# # format figure color
# ax_main.set_facecolor('#EEEEEE')
# for border in ['top','right','bottom','left']:
#     ax_main.spines[border].set_visible(False)
# ax_hist_right.set_facecolor('#EEEEEE')
# for border in ['top','right','bottom','left']:
#     ax_hist_right.spines[border].set_visible(False)
# ax_hist_top.set_facecolor('#EEEEEE')
# for border in ['top','right','bottom','left']:
#     ax_hist_top.spines[border].set_visible(False)


# plt.savefig(out_dir / 'Delta_DO_vs_natural_DO_lt_{}'.format(DO_thresh))

# print('Delta DO scatter done...')

# print('Done.')

#############################################################
##              CHANGE IN DO VS. NATURAL DO                ## (density colormap)
#############################################################

# initialize figure
plt.close('all)')
pfun.start_plot(figsize=(10,8))
fig,ax = plt.subplots(1,1)

# get points with DO < threshold, and get difference between model runs
DO_bot_lt6_natural = np.where(DO_bot_natural < DO_thresh, DO_bot_natural, np.nan)
DO_bot_lt6_anthropogenic = np.where(DO_bot_anthropogenic < DO_thresh, DO_bot_anthropogenic, np.nan)
# compress spatial and time dimensions
all_bot_DO_natural = np.reshape(DO_bot_lt6_natural,(28490805))
all_bot_DO_anthropogenic = np.reshape(DO_bot_lt6_anthropogenic,(28490805))
all_bot_DO_diff = all_bot_DO_anthropogenic - all_bot_DO_natural

# rename variables so its easier to manipulate
x = all_bot_DO_natural
y = all_bot_DO_diff
# get rid of nans in dataset
bad_indices = np.isnan(x) | np.isnan(y)
good_indices = ~bad_indices
good_x = x[good_indices]
good_y = y[good_indices]
# plot 2d histogram to get colored scatter by point density
cs = ax.hist2d(good_x, good_y, bins = [100,100], norm=mpl.colors.LogNorm(), cmap=cmocean.cm.thermal)
cbar = fig.colorbar(cs[3])
cbar.ax.set_ylabel('Count')
cbar.outline.set_visible(False)

# add y = -1 line
ax.fill_between([0,0.8],[0,-0.8], -0.8, color='white')

# add text
ax.text(0.12, 0.1, 'slope = -1', horizontalalignment='left',
     verticalalignment='center', transform=ax.transAxes, fontsize=12)

# format figure
ax.set_xlabel('Natural DO [mg/L]')
# format grid and labels
ax.set_xlim([0,DO_thresh])
ax.set_ylim([-0.8,0.4])
ax.grid(visible=True, color='w')
# put grid behind points and histogram
# ax.set_axisbelow(True)
# add labels
ax.set_ylabel('Anthropogenic - Natural [mg/L]')
plt.suptitle(r'$\Delta$ DO vs. Natural DO (Bottom)' + '\n for all cells/days in Puget Sound where DO < {} mg/L'.format(DO_thresh))
# format figure color
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)


plt.savefig(out_dir / 'Delta_DO_vs_natural_DO_lt_{}'.format(DO_thresh))

print('Delta DO scatter done...')

print('Done.')