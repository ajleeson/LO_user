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

remove_straits = True

vn = 'oxygen'

years =  ['2014']

# which  model run to look at?
gtagexes = ['cas7_t0_x4b'] # long hindcast (anthropogenic)

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
    if np.mod(int(years[0]),4) != 0: 
        DO_bot = np.insert(DO_bot,59,'nan',axis=0)
        hyp_thick = np.insert(hyp_thick,59,'nan',axis=0)
    DO_bot_dict[year] = DO_bot
    hyp_thick_dict[year] = hyp_thick

# get grid cell area
fp = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
PSbox_ds = xr.open_dataset(fp)
DX = (PSbox_ds.pm.values)**-1
DY = (PSbox_ds.pn.values)**-1
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
    xmax = -122.1 + 0.1 # to make room for legend key
    ymin = 46.95 - 0.1 # to make room for legend key
    ymax = 48.93

# get grid data
grid_ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
z = -grid_ds.h.values
mask_rho = np.transpose(grid_ds.mask_rho.values)
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
plon, plat = pfun.get_plon_plat(lon,lat)
# make a version of z with nans where masked
# this gives us a binary map of land and water cells
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
zm[np.transpose(mask_rho) != 0] = -1

##############################################################
##             Plot hypoxic volume time series              ##
##############################################################

# Puget Sound bounds
xmin = -123.29
xmax = -122.1
ymin = 46.95
ymax = 48.93

# initialize figure
fig, (ax0, ax1) = plt.subplots(1,2,figsize = (12,5.5),gridspec_kw={'width_ratios': [1, 4]})

# format figure
ax0.set_xlim([xmin,xmax])
ax0.set_ylim([ymin,ymax])
ax0.set_ylabel('Latitude', fontsize=12)
ax0.set_xlabel('Longitude', fontsize=12)
ax0.tick_params(axis='both', labelsize=12)
ax0.pcolormesh(plon, plat, zm, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
pfun.dar(ax0)
# Create a Rectangle patch to omit Straits
# get lat and lon
# fp = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
# PSbox_ds = xr.open_dataset(fp)
lons = PSbox_ds.coords['lon_rho'].values
lats = PSbox_ds.coords['lat_rho'].values
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
ax0.set_title('(a)', fontsize = 14, loc='left', fontweight='bold')

# Puget Sound volume with straits omitted
PS_vol = 195.2716230839466 # [km^3]

# create time vector
startdate = '2020.01.01'
enddate = '2020.12.31'
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

colors = ['#62B6CB','#A8C256','#96031A']
linewidths = [4,4,4]
alphas = [0.9,0.9,0.9]
linestyles = ['-','-','-']
labels=gtagexes

# plot timeseries
for i,year in enumerate(years):
    # plot hypoxic area timeseries
    ax1.plot(dates_local,hyp_vol[year],color=colors[i],linestyle=linestyles[i],
             linewidth=linewidths[i],alpha=alphas[i],label=labels[i])

# format figure
# ax1.grid(visible=True, axis='x', color='w')
ax1.grid(visible=True, axis='both', color='silver', linestyle='--')
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax1.tick_params(axis='both', labelsize=12)
ax1.set_ylabel(r'Hypoxic volume [km$^3$]', fontsize=12)
plt.legend(loc='upper left', fontsize=12)
plt.title('(b)', fontsize = 14, loc='left', fontweight='bold')
ax1.set_xlim([dates_local[0],dates_local[-1]])
ax1.set_ylim([0,10])

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
ax2.set_ylabel(r'Percent of regional volume [%]', fontsize=14)