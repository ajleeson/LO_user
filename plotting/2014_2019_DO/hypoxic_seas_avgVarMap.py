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
import csv
import gsw
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

remove_straits = False

years =  ['2014','2015','2016','2017','2018','2019'] 

# which  model run to look at?
gtagex = 'cas7_t0_x4b' # long hindcast (anthropogenic)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures'
Lfun.make_dir(out_dir)

region = 'Puget Sound'

# start date
start = '08-01'
end = '10-31'

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
ds_dict_DO = {}
ds_dict_vars = {}
# add ds to dictionary
for year in years:
    ds_DO = xr.open_dataset(Ldir['LOo'] / 'pugetsound_DO' / 'data' / (year + '_DO_info_' + straits + '.nc'))
    ds_vars = xr.open_dataset(Ldir['LOo'] / 'pugetsound_DO' / 'data' / (year + '_vars_' + straits + '.nc'))
    # crop to just hypoxic season
    ds_DO = ds_DO.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))
    ds_vars = ds_vars.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))
    # save dataset
    ds_dict_DO[year] = ds_DO
    ds_dict_vars[year] = ds_vars

##############################################################
##                         MAKE MAP                         ##
##############################################################

year = '2014'

vns = ['oxygen','stratification','DIN','phytoplankton','zooplankton','SdetritusN','LdetritusN']

labels = ['Bottom DO\n[mg L' + r'$^{-1}$' + ']',
      r'$\Delta\rho$' + '\n[kg m' + r'$^{-3}$' + ']',
      'Surface DIN (NO3+NH4)\n[mmol m'+ r'$^{-3}$' + ']',
      r'$\int_{H}^{\zeta}$' +' phytoplankton\n[kg Chl]',
      r'$\int_{H}^{\zeta}$' +' zooplankton\n[kg N]',
      r'$\int_{H}^{\zeta}$' +' S detritus\n[kg N]',
      r'$\int_{H}^{\zeta}$' +' L detritus\n[kg N]',
      ]

# get lat and lon
fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_2013.01.01_2013.12.31.nc')
ds = xr.open_dataset(fp)
lons = ds.coords['lon_rho'].values
lats = ds.coords['lat_rho'].values
px, py = pfun.get_plon_plat(lons,lats)
# get cell thickness
# get S for the whole grid
Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
reader = csv.DictReader(open(Sfp))
S_dict = {}
for row in reader:
    S_dict[row['ITEMS']] = row['VALUES']
S = zrfun.get_S(S_dict)
h = ds['h'].values # height of water column
z_rho, z_w = zrfun.get_z(h, 0*h, S) 
z_rho_surf = z_rho[-1,:,:]
z_rho_bot = z_rho[0,:,:]

# get plotting limits based on region
xmin = -123.29
xmax = -122.1
ymin = 46.95
ymax = 48.93

# Initialize figure
plt.close('all')
fs = 20
pfun.start_plot(fs=fs, figsize=(35,15))
fig,axes = plt.subplots(1,7)
ax = axes.ravel()

# plot each variable
for i, vn in enumerate(vns):

    if vn == 'oxygen':
        # get bottom DO (already in units of mg/L)
        bot_DO = ds_dict_DO[year]['DO_bot'].values
        # time average
        var = bot_DO.mean(axis=0)
        # get colorbar options
        vmin = 0
        vmax = 10
        cmap = plt.cm.get_cmap('rainbow_r', 10)
    elif vn == 'phytoplankton': 
        # get cell thickness * phytoplankton (mmol/m2)
        thick_phyto = ds_dict_vars[year]['intphyto'].values
        # multiply by horizontal area
        int_phyto = thick_phyto * 500 * 500
        # convert to kg from mmol
        int_phyto = int_phyto * pinfo.fac_dict['phytoplankton'] * (1/1000) * (1/1000) # 2.5 mg/mmol * 1/1000 g/mg * 1/1000 kg/g
        # time average
        var = int_phyto.mean(axis=0)
        # replace zeros with nans
        var[var == 0] = 'nan'
        vmin = 0
        vmax = 100
        cmap = cmocean.cm.algae
    elif vn == 'stratification':  
        # get surface and bottom temp and salinity
        surfS = ds_dict_vars[year]['surfS'].values 
        botS = ds_dict_vars[year]['botS'].values
        surfT = ds_dict_vars[year]['surfT'].values 
        botT = ds_dict_vars[year]['botT'].values 
        # Calculate pressure
        p_surf = gsw.p_from_z(z_rho_surf,lats)
        p_bot = gsw.p_from_z(z_rho_bot,lats)
        # calculate absolute salinity from practical salinity
        salt_abs_surf = gsw.conversions.SA_from_SP(surfS, z_rho_surf, lons, lats)
        salt_abs_bot = gsw.conversions.SA_from_SP(botS, z_rho_bot, lons, lats)
        # calculate conservative temperature from potential temperature
        temp_cons_surf = gsw.conversions.CT_from_pt(salt_abs_surf, surfT)
        temp_cons_bot = gsw.conversions.CT_from_pt(salt_abs_bot, botT)
        # calculate density
        rho_surf = gsw.rho(salt_abs_surf,temp_cons_surf,p_surf)
        rho_bot = gsw.rho(salt_abs_bot,temp_cons_bot,p_bot)
        # calculate density difference
        delta_rho = rho_bot - rho_surf
        # take average over season
        var = delta_rho.mean(axis=0)
        vmin = 0
        vmax = 12
        cmap = cmocean.tools.crop_by_percent(cmocean.cm.dense, 15, which='min')
    elif vn == 'DIN':
        # get surface NO3 and NH4
        surfNO3 = ds_dict_vars[year]['surfNO3'].values
        surfNH4 = ds_dict_vars[year]['surfNH4'].values
        # sum to get DIN
        surfDIN = surfNO3 + surfNH4
        # time average
        var = surfDIN.mean(axis=0)
        # get colorbar options
        vmin = 0
        vmax = np.nanmax(var)
        cmap = cmocean.cm.matter
    elif vn == 'zooplankton': 
        # get cell thickness * zooplankton (mmol/m2)
        thick_zoop = ds_dict_vars[year]['intzoop'].values
        # multiply by horizontal area
        int_zoop = thick_zoop * 500 * 500
        # convert to kg N from mmol
        int_zoop = int_zoop * 14 * (1/1000) * (1/1000) # 14 g/mol * 1/1000 g/mg * 1/1000 kg/g
        # time average
        var = int_zoop.mean(axis=0)
        # replace zeros with nans
        var[var == 0] = 'nan'
        vmin = 0
        vmax = 50
        cmap = cmocean.tools.crop_by_percent(cmocean.cm.amp, 10, which='min')
    elif vn == 'SdetritusN': 
        # get cell thickness * setritus (mmol/m2)
        thick_sdet = ds_dict_vars[year]['intSdetN'].values
        # multiply by horizontal area
        int_sdet = thick_sdet * 500 * 500
        # convert to kg N from mmol
        int_sdet = int_sdet * 14 * (1/1000) * (1/1000) # 14 g/mol * 1/1000 g/mg * 1/1000 kg/g
        # time average
        var = int_sdet.mean(axis=0)
        # replace zeros with nans
        var[var == 0] = 'nan'
        vmin = 0
        vmax = 500
        cmap = cmocean.tools.crop_by_percent(cmocean.cm.turbid, 10, which='min')
    elif vn == 'LdetritusN': 
        # get cell thickness * setritus (mmol/m2)
        thick_ldet = ds_dict_vars[year]['intLdetN'].values
        # multiply by horizontal area
        int_ldet = thick_ldet * 500 * 500
        # convert to kg N from mmol
        int_ldet = int_ldet * 14 * (1/1000) * (1/1000) # 14 g/mol * 1/1000 g/mg * 1/1000 kg/g
        # time average
        var = int_ldet.mean(axis=0)
        # replace zeros with nans
        var[var == 0] = 'nan'
        vmin = 0
        vmax = 30
        cmap = cmocean.tools.crop_by_percent(cmocean.cm.turbid, 10, which='min')
    else:
        continue

    # plot values
    cs = ax[i].pcolormesh(px,py,var, vmin=vmin, vmax=vmax, cmap=cmap)
    cbar = fig.colorbar(cs, location='bottom', pad = 0.01)
    cbar.ax.tick_params(rotation=30)
    cbar.outline.set_visible(False)

    # add title
    ax[i].set_title(labels[i])

    # add 10 km bar
    lat0 = 46.99
    lon0 = -122.4
    lat1 = lat0
    lon1 = lon0 + 0.13175
    distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
    x_dist_km = round(distances_m[0]/1000)
    ax[0].plot([lon0,lon1],[lat0,lat1],color='k',linewidth=8)
    ax[0].text(lon0-0.04,lat0+0.02,'{} km'.format(x_dist_km),color='k',fontsize=18)

    # format figure
    ax[i].set_xlim([xmin,xmax])
    ax[i].set_ylim([ymin,ymax])
    ax[i].set_yticklabels([])
    ax[i].set_xticklabels([])
    # ax[i].axis('off')
    pfun.dar(ax[i])

plt.suptitle(year + '\n' + start + ' to ' + end + ' average values', fontsize = 26)


# Generate plot
plt.tight_layout
plt.subplots_adjust(left=0.01, right=0.99, top=0.85, bottom=0.05, wspace=0.15)
plt.savefig(out_dir / (year+'_allvars_avg_'+ start + 'THRU'+end+'.png'))