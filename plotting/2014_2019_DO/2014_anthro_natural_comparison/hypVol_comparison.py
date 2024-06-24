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
pth = Path(__file__).absolute().parent.parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

##############################################################
##                       USER INPUTS                        ##
##############################################################

remove_straits = True

vn = 'oxygen'

years =  ['2014noN','2014']

# which  model run to look at?
gtagexes = ['cas7_t0noN_x4b','cas7_t0_x4b'] # long hindcast (anthropogenic)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures'
Lfun.make_dir(out_dir)

region = 'Puget Sound'

plt.close('all')

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
    if np.mod(int(years[1]),4) != 0: 
        DO_bot = np.insert(DO_bot,59,'nan',axis=0)
        hyp_thick = np.insert(hyp_thick,59,'nan',axis=0)
    DO_bot_dict[year] = DO_bot
    hyp_thick_dict[year] = hyp_thick

# initialize dictionary for hypoxic volume [km3]
hyp_vol = {}
for year in years:
    # get hypoxic thickness
    hyp_thick = hyp_thick_dict[year] # [m]
    # get timeseries of hypoxic volume (by summing over spatial dimensions)
    hyp_thick_timeseries = np.nansum(hyp_thick,axis=1)
    hyp_thick_timeseries = np.nansum(hyp_thick_timeseries,axis=1)
    # calculate hypoxic volume
    # area of grid cell = 0.5 km by 0.5 km times thickness of hypoxic layer (converted from m to km)
    hyp_vol_timeseries = 0.5 * 0.5 * (hyp_thick_timeseries/1000) 
    hyp_vol[year] = hyp_vol_timeseries
    print(np.nansum(hyp_vol_timeseries))

# get plotting limits based on region
if region == 'Puget Sound':
    # box extracion limits: [-123.29, -122.1, 46.95, 48.93]
    xmin = -123.29
    xmax = -122.1 + 0.1 # to make room for legend key
    ymin = 46.95 - 0.1 # to make room for legend key
    ymax = 48.93

##############################################################
##                HYPOXIC VOLUME TIMESERIES                 ##
##############################################################

# initialize figure
plt.close('all)')
pfun.start_plot(figsize=(16,8))
f, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 3]})

# format figure
ax0.set_xlim([xmin,xmax])
ax0.set_ylim([ymin,ymax])
ax0.set_yticklabels([])
ax0.set_xticklabels([])
# ax0.axis('off')
pfun.dar(ax0)
pfun.add_coast(ax0)
# Create a Rectangle patch to omit Straits
if remove_straits:
    # get lat and lon
    fp = Ldir['LOo'] / 'extract' / gtagexes[1] / 'box' / ('pugetsoundDO_2013.01.01_2013.12.31.nc')
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
                            edgecolor='none', facecolor='gray', alpha=0.5)
    # Add the patch to the Axes
    ax0.add_patch(rect)
    ax0.text(-123.2,48.25,'Straits\nomitted')

# create time vector
startdate = '2020.01.01'
enddate = '2020.12.31'
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

colors = ['darkturquoise','black']
linewidths = [5,2]
alphas = [0.5,1]
linestyles = ['-',':']
labels=['Natural','Anthropogenic']

# plot timeseries
for i,year in enumerate(years):
    # plot hypoxic area timeseries
    ax1.plot(dates_local,hyp_vol[year],color=colors[i],linestyle=linestyles[i],
             linewidth=linewidths[i],alpha=alphas[i],label=labels[i])

# format figure
ax1.grid(visible=True, color='w')
# format background color
ax1.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax1.spines[border].set_visible(False)
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
plt.legend(loc='best')
plt.title('2014 Puget Sound hypoxic volume (DO < 2 mg/L) ' + r'[km$^3$]')
ax1.set_xlim([dates_local[0],dates_local[-1]])

plt.savefig(out_dir / 'anthro_natural_hypVol')