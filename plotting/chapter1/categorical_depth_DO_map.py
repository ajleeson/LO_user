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
from matplotlib.colors import LinearSegmentedColormap
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
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
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
remove_straits = False

year = '2014'

# Variables to compare
vn = 'oxygen'

# model gtagex
gtagex = 'cas7_t0noN_x4b'

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures'
Lfun.make_dir(out_dir)

##############################################################
##                      PROCESS DATA                        ##
##############################################################

print('Processing started...')

# get data
fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_'+year+'.01.01_'+year+'.12.31.nc')
ds_box = xr.open_dataset(fp)

if remove_straits:
    # set values in strait of juan de fuca and strait of georgia are nan (so they don't interfere with analysis)
    lat_threshold = 48.14
    lon_threshold = -122.76
    # Create a mask for latitudes and longitudes that meet the condition
    mask_natural = (ds_box['lat_rho'] > lat_threshold) & (ds_box['lon_rho'] < lon_threshold)
    # Expand mask dimensions to match 'oxygen' dimensions
    expanded_mask_natural = mask_natural.expand_dims(ocean_time=len(ds_box['oxygen']), s_rho=len(ds_box['s_rho']))
    # Apply the mask to the 'oxygen' variable
    ds_box['oxygen'] = xr.where(expanded_mask_natural, np.nan, ds_box['oxygen'])

# get bottom DO  (in mg/L)
DO_bot_natural = pinfo.fac_dict[vn] * ds_box[vn][:,0,:,:].values # s_rho = 0 for bottom

# get depth
depths = -1 * ds_box['h'].values
# duplicate depths for all times
depths = np.repeat(depths[np.newaxis, :, :], 365, axis=0)

# Get boolean array. True if DO < threshold mg/L, otherwise, nan
DO_bot_lt2 = np.where(DO_bot_natural < DO_thresh, 1, np.nan)

# Sum over time to compress into just lat/lon dimension, with values indicating days with bottom DO < 2mg/L
DO_days_lt2 = np.nansum(DO_bot_lt2, axis=0)
# get convert zeros to nans in dataset
DO_days_lt2 = np.where(DO_days_lt2==0, np.nan, DO_days_lt2)

print('Data processing done...')

# get plotting limits of box extraction
xmin = -123.29
xmax = -122.1
ymin = 46.95
ymax = 48.93

# get lat/lon
lons = ds_box.coords['lon_rho'].values
lats = ds_box.coords['lat_rho'].values
px, py = pfun.get_plon_plat(lons,lats)

#############################################################
##                HYPOXIA AND DEPTH  MAPS                  ## 
#############################################################

deep_cutoff = 125
shallow_cutoff = 10

# Get boolean array. True if DO < 2 mg/L, otherwise, nan
DO_bot_hyp_natural = np.where(DO_bot_natural < 2, 1, np.nan)
# compress time dimension, so we get true values if gridcell sees hypoxia at least once in the year
DO_bot_hyp_natural = np.nansum(DO_bot_hyp_natural,axis=0)
# convert zeros to nans, and make any positive value = 1
DO_bot_hyp_natural[DO_bot_hyp_natural==0] = np.nan
DO_bot_hyp_natural[DO_bot_hyp_natural > 0] = 1
# DO_bot_hyp_natural and DO_bot_hyp_anthropogenic are now arrays in space that are nan if there is no hypoxia,
# and 1 if the bottom grid cell experiences hypoxia at least one time during the year

# Now lets get boolean array for depth. True if -10 m  > depth > -50 m, otherwsie, nan
depths = -1 * ds_box['h'].values
shallow_depths = (depths > -deep_cutoff) & (depths < -shallow_cutoff)
# set true and false to one and nan
shallow_depths = np.where(shallow_depths == False, np.nan, 1)

# Now lets get boolean array for depth. True if depth > -10 m, otherwsie, nan
shallowest_depths = (depths >= -shallow_cutoff) & (depths < -4) # need to remove the -4 depth mask
# set true and false to one and nan
shallowest_depths = np.where(shallowest_depths == False, np.nan, 1)

# land mask
land_depths = (depths >= -4)
# set true and false to one and nan
land_depths = np.where(land_depths == False, np.nan, 1)

#--------------------------------------------------------------------------------

# (6) Top
# Get locations where depths are surface-shallow, and the bottom becomes hypoxia
surface_hypoxic_natural = (DO_bot_hyp_natural == 1) & (shallowest_depths == 1)
# set true and false to one and nan
surface_hypoxic_natural = np.where(surface_hypoxic_natural == False, np.nan, 6)

# (5)
# Get locations where depths are surface-shallow, and bottomd does NOT get hypoxic
surface_natural = (DO_bot_hyp_natural != 1) & (shallowest_depths == 1)
# set true and false to one and nan
surface_natural = np.where(surface_natural == False, np.nan, 5)

# (4)
# Get locations where depths are shallow (10 - 50 m deep), and the bottom becomes hypoxia
lt_50_hypoxic_natural = (DO_bot_hyp_natural == 1) & (shallow_depths == 1)
# set true and false to one and nan
lt_50_hypoxic_natural = np.where(lt_50_hypoxic_natural == False, np.nan, 4)

# (3)
# Get locations where depths are shallow (10 - 50 m deep), and bottomd does NOT get hypoxic
lt_50_natural = (DO_bot_hyp_natural != 1) & (shallow_depths == 1)
# set true and false to one and nan
lt_50_natural = np.where(lt_50_natural == False, np.nan, 3)

# (2)
# Get locations where depths are deep, and bottom does get hypoxic
deep_hypoxic_natural = (DO_bot_hyp_natural == 1) & (shallow_depths != 1) & (shallowest_depths != 1)
# set true and false to one and nan
deep_hypoxic_natural = np.where(deep_hypoxic_natural == False, np.nan, 2)

# (1) Bottom
# Get locations where depths are deep, and bottom does not get hypoxic
deep_natural = (DO_bot_hyp_natural != 1) & (shallow_depths != 1) & (shallowest_depths != 1) & (land_depths != 1)
# set true and false to one and nan
deep_natural = np.where(deep_natural == False, np.nan, 1)

# sum up all of the categorical values, and set zero to nan
depth_hypoxia_category_natural = np.nansum(np.dstack((lt_50_hypoxic_natural,lt_50_natural,surface_hypoxic_natural,
                                                      deep_hypoxic_natural,surface_natural,deep_natural)),2)
depth_hypoxia_category_natural[depth_hypoxia_category_natural==0] = np.nan

colors = ['cornflowerblue',
          'crimson',
          'powderblue',
          'plum',
          'orange',
          'black']
cmap = LinearSegmentedColormap.from_list('categorycolor', colors, N=6)


# get plotting limits of box extraction
xmin = -123.29
xmax = -122.1
ymin = 46.95
ymax = 48.93

# get lat/lon
lons = ds_box.coords['lon_rho'].values
lats = ds_box.coords['lat_rho'].values
px, py = pfun.get_plon_plat(lons,lats)

# initialize figure
plt.close('all')
pfun.start_plot(figsize=(12,18))
fig = plt.figure()
# gs = fig.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.95, wspace=0.05, hspace=0.05)

# plot natural run
ax = fig.add_subplot(1,1,1)
cs = ax.pcolormesh(px,py,depth_hypoxia_category_natural, cmap=cmap)
cbar = fig.colorbar(cs, location='left')
cbar.ax.tick_params(labelsize=20)#,length=10, width=2)
cbar.outline.set_visible(False)
# label colorbar
yticks = np.linspace(*cbar.ax.get_ylim(), cmap.N+1)[:-1]
yticks += (yticks[1] - yticks[0]) / 2
cbar.set_ticks(yticks, labels=['z > {} m\nNOT hypoxic'.format(str(deep_cutoff)), # bottom label
                               'z > {} m\nhypoxic'.format(str(deep_cutoff)),
                               '{} m < z < {} m\nNOT hypoxic'.format(str(shallow_cutoff),str(deep_cutoff)),
                               '{} m < z < {} m\nhypoxic'.format(str(shallow_cutoff), str(deep_cutoff)),
                               'z < {} m\nNOT hypoxic'.format(str(shallow_cutoff)),
                               'z < {} m\nhypoxic'.format(str(shallow_cutoff))]) # top label
# color hypoxic labels red
cbar.ax.get_yticklabels()[1].set_color('crimson')
cbar.ax.get_yticklabels()[3].set_color('crimson')
cbar.ax.get_yticklabels()[5].set_color('crimson')
# remove tick lines
cbar.ax.tick_params(length=0)


# format figure
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.axis('off')
# pfun.add_coast(ax)
pfun.dar(ax)

# Add title
plt.suptitle(year + ' depth and hypoxia (< 2 mg/L)',
            fontweight='bold', fontsize=24, y=0.92, x=0.65)
ax.set_title('[hypoxic if DO < 2 mg/L for at least one day of year]',
             fontsize = 20)

# save figure
plt.subplots_adjust(left=0.21, right=0.98, bottom = 0.08) # wspace=0.02, 
plt.savefig(out_dir / 'depth_hypoxia_categorical_map')

print('hypoxia depth map done...')