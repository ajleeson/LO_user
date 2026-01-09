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

year =  '2014'

# which  model run to look at?
gtagexes = ['cas7_t0_x4b','cas7_t1_x11ab','cas7_t1noDIN_x11ab'] 

# where to put output figures
out_dir = Ldir['LOo'] / 'chapter_2' / 'figures'
Lfun.make_dir(out_dir)

plt.close('all')

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# open datasets
if remove_straits:
    straits = 'noStraits'
else:
    straits = 'withStraits'

# region = 'pugetsoundDO'
regions = ['pugetsoundDO','HC_up','HC_low','SS_and_HC_low']

# initialize empty dictionaries
DO_bot_dict = {} # dictionary with DO_bot values
hyp_thick_dict = {}

for region in regions:
    for gtagex in gtagexes:
        # add ds to dictionary
        if region == 'pugetsoundDO':
            # only include strait infor when looking at the whole of puget sound
            ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / (gtagex + '_' + region + '_' + year + '_DO_info_' + straits + '.nc'))
        else:
            if gtagex == 'cas7_t0_x4b':
                # don't compare old version of LiveOcean in the sub-basins
                continue
            else:
                ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / (gtagex + '_' + region + '_' + year + '_DO_info.nc'))
        DO_bot = ds['DO_bot'].values
        hyp_thick = ds['hyp_thick'].values
        # if not a leap year, add a nan on feb 29 (julian day 60 - 1 because indexing from 0)
        if np.mod(int(year),4) != 0: 
            DO_bot = np.insert(DO_bot,59,'nan',axis=0)
            hyp_thick = np.insert(hyp_thick,59,'nan',axis=0)
        DO_bot_dict[gtagex+region] = DO_bot
        hyp_thick_dict[gtagex+region] = hyp_thick

# initialize empty dictionary to get grid cell areas
DA = {}
for region in regions:
    fp = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'box' / (region + '_2014.01.01_2014.12.31.nc')
    box_ds = xr.open_dataset(fp)
    if region == 'pugetsoundDO':
        PSbox_ds = box_ds # save Puget Sound bounds for later
    DX = (box_ds.pm.values)**-1
    DY = (box_ds.pn.values)**-1
    DA[region] = DX*DY*(1/1000)*(1/1000) # get area, but convert from m^2 to km^2

# initialize dictionary for hypoxic volume [km3]
hyp_vol = {}
for region in regions:
    for gtagex in gtagexes:
        # don't use old version of LiveOcean for the sub-basins
        if gtagex == 'cas7_t0_x4b' and region != 'pugetsoundDO':
            continue
        else:
            # get hypoxic thickness
            hyp_thick = hyp_thick_dict[gtagex+region]/1000 # [km]
            # subtract the bottom row of grid cells for the HC_up region
            # (so we don't double count cells when we add them to HC_low)
            if region == 'HC_up':
                hyp_thick[:, 0, :] = 0 # dimensions of (t,y,x), so y=0 is the lowest latitude of HC_up
            hyp_vol_timeseries = np.sum(hyp_thick * DA[region], axis=(1, 2)) # km^3
            hyp_vol[gtagex+region] = hyp_vol_timeseries

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
##   Puget Sound hypoxic volume & comparison to old model   ##
##############################################################

# Puget Sound bounds
xmin = -123.29
xmax = -122.1
ymin = 46.95
ymax = 48.93

# initialize figure
fig, (ax0, ax1) = plt.subplots(1,2,figsize = (12,5.5),gridspec_kw={'width_ratios': [1, 2]})

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

colors = ['gray','deeppink','black']
linewidths = [4,3,2]
alphas = [0.5,0.5,1]
linestyles = ['-','-','--']
labels=gtagexes

# plot timeseries
for i,gtagex in enumerate(gtagexes):
    # plot hypoxic area timeseries
    ax1.plot(dates_local,hyp_vol[gtagex+'pugetsoundDO'],color=colors[i],linestyle=linestyles[i],
             linewidth=linewidths[i],alpha=alphas[i],label=labels[i])

# format figure
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

plt.show()

##############################################################
##                  Sub-basin hypoxic volume                ##
##############################################################

gtagexes = ['cas7_t1_x11ab','cas7_t1noDIN_x11ab'] 

# Puget Sound bounds
xmin = -123.29
xmax = -122.1
ymin = 46.95
ymax = 48.93

# initialize figure
fig, (ax0, ax1) = plt.subplots(1,2,figsize = (12,5.5),gridspec_kw={'width_ratios': [1, 2]})

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

# upper hood canal
hcminlat = 47.5
hcmaxlat = 47.86
hcminlon = -123.24
hcmaxlon = -122.72
rect = patches.Rectangle((hcminlon, hcminlat), hcmaxlon-hcminlon, hcmaxlat-hcminlat,
                        edgecolor='black', facecolor='turquoise', alpha=0.3)
# Add the patch to the Axes
ax0.add_patch(rect)

# lower hood canal
hclowminlat = 47.305
hclowmaxlat = hcminlat
hclowminlon = hcminlon
hclowmaxlon = -122.88
rect = patches.Rectangle((hclowminlon, hclowminlat), hclowmaxlon-hclowminlon, hclowmaxlat-hclowminlat,
                        edgecolor='black', facecolor='turquoise', alpha=0.3)
# Add the patch to the Axes
ax0.add_patch(rect)

# south sound
ssminlat = 46.96
ssmaxlat = hclowmaxlat
ssminlon = hclowminlon
ssmaxlon = -122.57
rect = patches.Rectangle((ssminlon, ssminlat), ssmaxlon-ssminlon, ssmaxlat-ssminlat,
                        edgecolor='black', facecolor='hotpink', alpha=0.3)
# Add the patch to the Axes
ax0.add_patch(rect)

# create time vector
startdate = '2020.01.01'
enddate = '2020.12.31'
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

colors = ['black','turquoise','hotpink']
linewidths = [3,2]
alphas = [0.5,1]
linestyles = ['-','--']

# plot timeseries
basins = ['Puget Sound','Hood Canal','South Sound']
for j,basin in enumerate(basins):
    for i,gtagex in enumerate(gtagexes):
        # plot hypoxic area timeseries
        if basin == 'Puget Sound':
            hypoxic_volume = hyp_vol[gtagex+'pugetsoundDO']
        elif basin == 'Hood Canal':
            hypoxic_volume = hyp_vol[gtagex+'HC_up'] + hyp_vol[gtagex+'HC_low']
        elif basin == 'South Sound':
            hypoxic_volume = hyp_vol[gtagex+'SS_and_HC_low'] - hyp_vol[gtagex+'HC_low']

        if basin == 'Puget Sound':
            if i == 0:
                loading = '; Loading'
            else:
                loading = '; No-Loading'
            ax1.plot(dates_local,hypoxic_volume,color=colors[j],linestyle=linestyles[i],
                    linewidth=linewidths[i],alpha=alphas[i],label=basin+loading)
        else:
            if i == 0:
                ax1.plot(dates_local,hypoxic_volume,color=colors[j],linestyle=linestyles[i],
                    linewidth=linewidths[i],alpha=alphas[i],label=basin)
            else:
                ax1.plot(dates_local,hypoxic_volume,color=colors[j],linestyle=linestyles[i],
                    linewidth=linewidths[i],alpha=alphas[i])

# format figure
ax1.grid(visible=True, axis='both', color='silver', linestyle='--')
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
ax1.tick_params(axis='both', labelsize=12)
ax1.set_ylabel(r'Hypoxic volume [km$^3$]', fontsize=12)
plt.legend(loc='upper left', fontsize=12)
plt.title('(b)', fontsize = 14, loc='left', fontweight='bold')
ax1.set_xlim([dates_local[0],dates_local[-1]])
ax1.set_ylim([0,1.25])

#  # convert hypoxic volume to percent hypoxic volume
# percent = lambda hyp_vol: hyp_vol/PS_vol*100
# # get left axis limits
# ymin, ymax = ax1.get_ylim()
# # match ticks
# ax2 = ax1.twinx()
# ax2.set_ylim((percent(ymin),percent(ymax)))
# ax2.plot([],[])
# for border in ['top','right','bottom','left']:
#     ax2.spines[border].set_visible(False)
# ax2.set_ylabel(r'Percent of regional volume [%]', fontsize=14)

plt.show()