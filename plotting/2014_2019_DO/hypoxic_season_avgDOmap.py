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
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

##############################################################
##                       USER INPUTS                        ##
##############################################################

# Show WWTP locations?
WWTP_loc = True

remove_straits = False

vn = 'oxygen'

years =  ['2014','2015','2016','2017','2018','2019'] # ['2014','2015','2016','2017','2018','2019']

# which  model run to look at?
gtagex = 'cas7_t0_x4b' # long hindcast (anthropogenic)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures'
Lfun.make_dir(out_dir)

region = 'Puget Sound'

# start date
start = '08-01'
end = '09-30'

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
    ds_dict[year] = ds

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
    xmax = -122.1 + 0.1 # to make room for legend key
    ymin = 46.95 - 0.1 # to make room for legend key
    ymax = 48.93

# set axes range for different state variables
if vn == 'NO3':
    vmin = 0
    vmax = 40
    cmap = cmocean.cm.matter
elif vn == 'NH4':
    vmin = 0
    vmax = 6
    cmap = cmocean.cm.matter
elif vn == 'phytoplankton':
    vmin = 0
    vmax = 30
    cmap = cmocean.cm.algae
elif vn == 'oxygen':
    vmin = 0
    vmax = 10
    cmap = plt.cm.get_cmap('rainbow_r', 10)
    stext = 'bottom'
    slev = '0'
    var = 'DO_bot'
elif vn == 'SdetritusN':
    vmin = 0
    vmax = 5
    cmap = cmocean.cm.matter
elif vn == 'LdetritusN':
    vmin = 0
    vmax = 0.1
    cmap = cmocean.cm.matter

# scale variable & get units
# scale =  pinfo.fac_dict[vn]
units = pinfo.units_dict[vn]

# Initialize figure
fs = 10
pfun.start_plot(fs=fs, figsize=(51,21))
fig,axes = plt.subplots(1,len(years)+1)
ax = axes.ravel()

# get lat and lon
fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_2013.01.01_2013.12.31.nc')
ds = xr.open_dataset(fp)
lons = ds.coords['lon_rho'].values
lats = ds.coords['lat_rho'].values
px, py = pfun.get_plon_plat(lons,lats)

# store all values in new dictionary (the averages over hypoxic season)
val_dict = {}
# plot average year
for year in years:
    ds = ds_dict[year]
    # crop to just hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))
    v = ds[var].values
    # take average over season
    v = np.nanmean(v,axis=0)
    val_dict[year] = v

# calculate average of all of the arrays
v_avg = sum(val_dict.values())/len(val_dict)
cs = ax[0].pcolormesh(px,py,v_avg, vmin=vmin, vmax=vmax, cmap=cmap)
cbar = fig.colorbar(cs, location='left', anchor=(-1,0.5),
                     ax=[ax[0],ax[1],ax[2],ax[3],ax[4],ax[5],ax[6]])
cbar.ax.tick_params(labelsize=32)#,length=10, width=2)
cbar.outline.set_visible(False)
ax[0].set_title(letters[0] + ' Average', fontsize=38, loc='left')

 # add wwtp locations
if WWTP_loc == True:
    ax[0].scatter(lon_wwtps,lat_wwtps,color='none', edgecolors='k', linewidth=3, s=sizes_wwtps, label='WWTPs')
    leg_szs = [100, 1000, 10000]
    szs = [0.3*(leg_sz) for leg_sz in leg_szs]
    l0 = plt.scatter([],[], s=szs[0], color='none', edgecolors='k', linewidth=3)
    l1 = plt.scatter([],[], s=szs[1], color='none', edgecolors='k', linewidth=3)
    l2 = plt.scatter([],[], s=szs[2], color='none', edgecolors='k', linewidth=3)
    labels = ['< 100', '1,000', '10,000']
    legend = ax[0].legend([l0, l1, l2], labels, fontsize = 18, markerfirst=False,
        title='WWTP loading \n'+r' (kg N d$^{-1}$)',loc='lower right', labelspacing=1, borderpad=0.8)
    plt.setp(legend.get_title(),fontsize=20)

# add 10 km bar
lat0 = 46.94
lon0 = -123.05
lat1 = lat0
lon1 = -122.91825
distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
x_dist_km = round(distances_m[0]/1000)
ax[0].plot([lon0,lon1],[lat0,lat1],color='k',linewidth=8)
ax[0].text(lon0-0.04,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=24)

# format figure
ax[0].set_xlim([xmin,xmax])
ax[0].set_ylim([ymin,ymax])
ax[0].set_yticklabels([])
ax[0].set_xticklabels([])
ax[0].axis('off')
pfun.dar(ax[0])


# loop through and plot both conditions
for j,year in enumerate(years):

    i = j+1
                
    v = val_dict[year]
    cs = ax[i].pcolormesh(px,py,v-v_avg, vmin=-1.5, vmax=1.5, cmap=cmocean.cm.balance_r)

    if i == 6:
        cbar = fig.colorbar(cs, location='right', anchor=(1.8,0.5),
                     ax=[ax[0],ax[1],ax[2],ax[3],ax[4],ax[5],ax[6]])
        cbar.ax.tick_params(labelsize=32)#,length=10, width=2)
        cbar.outline.set_visible(False)

    # format figure
    ax[i].set_xlim([xmin,xmax])
    ax[i].set_ylim([ymin,ymax])
    ax[i].set_yticklabels([])
    ax[i].set_xticklabels([])
    ax[i].axis('off')
    pfun.dar(ax[i])
    ax[i].set_title(letters[i] + ' ' + year + ' - Avg', fontsize=38, loc = 'left')
                                
# Add colormap title
plt.suptitle(start+' to '+end+' average ' + stext + ' ' + vn + ' ' + units,
            fontsize=44, fontweight='bold', y=0.95)

# Generate plot
plt.tight_layout
plt.subplots_adjust(left=0.05, right=0.95, top=0.85, wspace=0.02)
plt.savefig(out_dir / (vn+'_'+stext+'_avg_'+ start + 'THRU'+end+'.png'))