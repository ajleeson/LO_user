"""
Compare LiveOcean (cas7_t0_x4b) values to KC ctd observations in Whidbey Basin
with particular focus in Port Susan
"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
from matplotlib.patches import Rectangle
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import gsw
import pinfo
import pickle
from importlib import reload
reload(pinfo)

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
# jobname = 'twentyoneinlets'
# startdate = '2017.01.01'
# enddate = '2017.12.31'
year = '2023' # for making a date label

# jobname = 'twentyoneinlets'
# # find job lists from the extract moor
# job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')
# # Get stations:
# sta_dict = job_lists.get_sta_dict(jobname)

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# # find job lists from the extract moor
# job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# # Get mooring stations:
# sta_dict = job_lists.get_sta_dict(jobname)

# get observations
in_dir = Ldir['parent'] / 'LO_output' / 'obsmod' / 'kc_whidbey'
in_fn = in_dir / ('multi_ctd_' + year + '.p')
df_dict = pickle.load(open(in_fn, 'rb'))
# only look at ecology stations
source = 'kc_whidbey'
df_obs = df_dict['obs'].loc[df_dict['obs'].source==source,:]
df_model = df_dict[gtagex].loc[df_dict[gtagex].source==source,:]

# get observational and model information for Port Susan
df_ob_PS = df_obs.loc[df_obs.lat>=48.1,:].loc[df_obs.lon>=-122.45,:]
df_mo_PS = df_model.loc[df_model.lat>=48.1,:].loc[df_model.lon>=-122.45,:]

# get observational and model information for all other stations
df_ob_filtered = pd.merge(df_obs, df_ob_PS, on=['cid', 'time', 'lat', 'lon', 'z', 'cruise', 'name', 'CT', 'SA', 'DO (uM)', 'Chl (mg m-3)', 'source'],
            how='outer', indicator=True).query("_merge != 'both'").drop('_merge', axis=1).reset_index(drop=True)
df_mo_filtered = pd.merge(df_model, df_mo_PS, on=['cid', 'time', 'lat', 'lon', 'z', 'cruise', 'name', 'CT', 'SA', 'DO (uM)', 'Chl (mg m-3)', 'source'],
            how='outer', indicator=True).query("_merge != 'both'").drop('_merge', axis=1).reset_index(drop=True)

##########################################################
##                 Get station locations                ##
##########################################################

plt.close('all')

# Port Susan stations
cids = df_ob_PS['cid'].unique()
ctd_lat_PS = [df_ob_PS.loc[df_ob_PS['cid'] == x, 'lat'].iloc[0] for x in cids]
ctd_lon_PS = [df_ob_PS.loc[df_ob_PS['cid'] == x, 'lon'].iloc[0] for x in cids]

# all other stations
cids = df_ob_filtered['cid'].unique()
ctd_lat_other = [df_ob_filtered.loc[df_ob_filtered['cid'] == x, 'lat'].iloc[0] for x in cids]
ctd_lon_other = [df_ob_filtered.loc[df_ob_filtered['cid'] == x, 'lon'].iloc[0] for x in cids]

##########################################################
##          Plotting terminal inlets                    ##
##########################################################

letter = ['(a)','(b)','(c)','(d)','(e)']

# initialize figure
fig = plt.figure(figsize=(14, 7))
ax_map = plt.subplot(1,2,1)
ax_ct = plt.subplot(2,4,3)
ax_sa = plt.subplot(2,4,4)
ax_chl = plt.subplot(2,4,7)
ax_do = plt.subplot(2,4,8)

##########################################################
##           Terminal inlet location map                ##
##########################################################

# create map of Puget Sound and inlets
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
# zm[np.transpose(mask_rho) != 0] = -1

# Create map
lon_low = -122.75
lon_high =  -122.2
lat_low = 47.9
lat_high = 48.3
cs = ax_map.pcolormesh(plon, plat, zm, linewidth=0.5, alpha=0.5, vmin=-200, vmax=0, cmap=plt.get_cmap(cmocean.cm.thermal))
cbar = fig.colorbar(cs)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_ylabel(r'Depth [m]', rotation=90, fontsize=12)
cbar.outline.set_visible(False)

# format
ax_map.set_xlim(lon_low,lon_high) # Puget Sound
ax_map.set_ylim(lat_low, lat_high) # Puget Sound
pfun.dar(ax_map)

# add cast locations
ax_map.plot(ctd_lon_other,ctd_lat_other,linestyle='none',marker='o',markersize=9,
    color='navy',markeredgecolor='white',markeredgewidth=0.8)

# highlight port susan casts
ax_map.plot(ctd_lon_PS,ctd_lat_PS,linestyle='none',marker='o',markersize=9,
    color='deeppink',markeredgecolor='white',markeredgewidth=0.8)

ax_map.set_title('(a) Inlet obs comparison locations',fontsize=12,loc='left')

##########################################################
##     Terminal inlets Property-property plots          ##
##########################################################

axes = [ax_ct,ax_sa,ax_chl,ax_do]
vars = [r'Cons Temp [$^\circ$C]',r'Abs Salinity [g kg$^{-1}$]',
        r'Chl [mg m$^{-3}$]', r'DO [mg L$^{-1}$]']
vns = ['CT','SA','Chl (mg m-3)','DO (uM)']

# format grid
for i,ax in enumerate(axes):
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax.tick_params(axis='both', labelsize=10)
    ax.set_title(letter[i+1] + ' ' + vars[i],loc='left',fontsize=12)

# plot property-property plots
for i,vn in enumerate(vns):
    # initialize arrays to save values
    deep_obs_PS = np.array([])
    deep_mod_PS = np.array([])
    deep_obs_other = np.array([])
    deep_mod_other = np.array([])
    if vn == 'DO (uM)':
        obsvals_PS = df_ob_PS[vn].values * 32/1000
        modvals_PS = df_mo_PS[vn].values * 32/1000
        obsvals_other = df_ob_filtered[vn].values * 32/1000
        modvals_other = df_mo_filtered[vn].values * 32/1000
    else:
        obsvals_PS = df_ob_PS[vn].values 
        modvals_PS = df_mo_PS[vn].values
        obsvals_other = df_ob_filtered[vn].values 
        modvals_other = df_mo_filtered[vn].values
    # get max limits
    if vn == 'CT':
        maxval = 25
    if vn == 'SA':
        maxval = 40
    if vn == 'Chl (mg m-3)':
        maxval = 35
    if vn == 'DO (uM)':
        maxval = 15
    # plot data and y = x line
    axes[i].scatter(obsvals_other,modvals_other, color='navy', alpha=0.3, s=6, zorder=10, edgecolor='none')
    axes[i].scatter(obsvals_PS,modvals_PS, color='deeppink', alpha=0.3, s=6, zorder=10, edgecolor='none')
    axes[i].plot([0,maxval*1.01], [0,maxval*1.01], color='grey', linestyle='-', zorder=5)
    axes[i].set_ylim([0,maxval])
    axes[i].set_xlim([0,maxval])
    axes[i].set_aspect('equal', adjustable='box')
    axes[i].xaxis.set_major_locator(plt.MaxNLocator(5))
    axes[i].yaxis.set_major_locator(plt.MaxNLocator(5))
    axes[i].set_ylabel('Modeled',fontsize=12)
    axes[i].set_xlabel('Observed',fontsize=12)
    # add data to arrays
    deep_obs_other = np.concatenate((deep_obs_other,obsvals_other))
    deep_mod_other = np.concatenate((deep_mod_other,modvals_other))
    deep_obs_PS = np.concatenate((deep_obs_PS,obsvals_PS))
    deep_mod_PS = np.concatenate((deep_mod_PS,modvals_PS))
    # # add bias and rmse label
    # calculate bias and rmse
    bias_other = np.nanmean(deep_mod_other-deep_obs_other)
    rmse_other = np.sqrt(np.nanmean((deep_mod_other-deep_obs_other)**2))
    bias_PS = np.nanmean(deep_mod_PS-deep_obs_PS)
    rmse_PS = np.sqrt(np.nanmean((deep_mod_PS-deep_obs_PS)**2))
    t01 = axes[i].text(0.17, 0.95, 'Bias', fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        transform=axes[i].transAxes, fontsize=10, color = 'black')
    t02 = axes[i].text(0.37, 0.95, 'RMSE', fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        transform=axes[i].transAxes, fontsize=10, color = 'black')
    t1 = axes[i].text(0.17, 0.87, '{:.2f}'.format(round(bias_other,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'navy')
    t2 = axes[i].text(0.37, 0.87, '{:.2f}'.format(round(rmse_other,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'navy')
    t3 = axes[i].text(0.17, 0.79, '{:.2f}'.format(round(bias_PS,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'deeppink')
    t4 = axes[i].text(0.37, 0.79, '{:.2f}'.format(round(rmse_other,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'deeppink')
    t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    t3.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    t4.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

plt.suptitle('{}'.format(year), fontsize=16)
plt.subplots_adjust(wspace=-0.2)      
plt.tight_layout()
