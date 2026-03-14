"""
Make boxplot of bottom DO concentration in 
Southern Hood Canal

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches
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

remove_straits = True

years =  ['2014','2015','2016','2017','2018','2019']

# which  model run to look at?
gtagex = 'cas7_t0_x4b' # long hindcast (anthropogenic)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures'
Lfun.make_dir(out_dir)

# start date
start = '08-01'
end = '09-30'

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# open dataset for every year, and add to dictionary, with year as key

# open datasets
if remove_straits:
    straits = 'noStraits'
else:
    straits = 'withStraits'
# initialize empty dictionary
ds_dict = {}
# add ds to dictionary
for year in years:
    ds = xr.open_dataset(Ldir['LOo'] / 'pugetsound_DO' / 'data' / (year + '_DO_info_' + straits + '.nc'))
    # crop to just hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))
    ds_dict[year] = ds

#############################################################
##              CREATE BOXPLOT OF BOTTOM DO                ##
#############################################################

# get lat and lon
fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_2013.01.01_2013.12.31.nc')
ds = xr.open_dataset(fp)
lons = ds.coords['lon_rho'].values
lats = ds.coords['lat_rho'].values
lon = lons[0,:]
lat = lats[:,0]

# sourthern hood canal
lonmax = -122.846790
lonmin = -123.171810
latmax = 47.630114 #47.5
latmin = 47.339921

# convert lat/lon to eta/xi
diff = np.absolute(lon-lonmin)
ximin = diff.argmin()
diff = np.absolute(lon-lonmax)
ximax = diff.argmin()
diff = np.absolute(lat-latmin)
etamin = diff.argmin()
diff = np.absolute(lat-latmax)
etamax = diff.argmin()

# create dictionary of bottom DO concentration (reshpaed to 1D array)
# initialize empty dataframe
DO_bot_dict = {}
# add ds to dictionary
for year in years:
    # get min DO values
    DO_bot = ds_dict[year].DO_bot
    DO_bot = DO_bot.values[:,etamin:etamax,ximin:ximax]
    # compress spatial and time dimensions
    DO_bot = np.reshape(DO_bot,-1)
    # remove nans
    DO_bot = DO_bot[~np.isnan(DO_bot)]
    DO_bot_dict[year] = DO_bot
    print(year)
    print(np.nanmin(DO_bot))

# initialize figure
plt.close('all')
pfun.start_plot(figsize=(12,7))
f, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3]})

# plot map
xmin = -123.29
xmax = -122.1
ymin = 46.95
ymax = 48.93
# format figure
ax0.set_xlim([xmin,xmax])
ax0.set_ylim([ymin,ymax])
ax0.set_yticklabels([])
ax0.set_xticklabels([])
ax0.axis('off')
pfun.dar(ax0)
pfun.add_coast(ax0)
# Create a Rectangle patch
rect = patches.Rectangle((lon[ximin], lat[etamin]), lon[ximax]-lon[ximin], lat[etamax]-lat[etamin],
                         edgecolor='none', facecolor='mediumorchid', alpha=0.5)
# Add the patch to the Axes
ax0.add_patch(rect)


# Plot boxplot
dot = dict(markerfacecolor='mediumorchid', marker='.',alpha=0.2,markeredgecolor='none')
bp = ax1.boxplot(DO_bot_dict.values(), patch_artist = True, flierprops=dot)
plt.setp(bp['medians'], color='mediumorchid')
for patch in bp['boxes']:
    patch.set(facecolor='mediumorchid', alpha=0.3)

# format figure  
ax1.set_xticklabels(DO_bot_dict.keys())
ax1.set_ylabel('DO [mg/L]')
ax1.set_xlabel('Year')
ax1.grid(visible=True, color='w')
# format background color
ax1.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax1.spines[border].set_visible(False)
ax1.set_title('Distribution of bottom DO in Southern Hood Canal\nbetween ' + start + ' to ' + end)

plt.savefig(out_dir / ('Boxplot_HC_botDO.png'))
plt.close('all')
