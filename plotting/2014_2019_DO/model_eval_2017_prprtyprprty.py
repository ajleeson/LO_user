"""
Generate depth vs. time property plots using mooring extraction data. 
Used to look at 21 inlets in Puget Sound, and compare to Ecology monitoring stations, if available
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
jobname = 'twentyoneinlets'
startdate = '2017.01.01'
enddate = '2017.12.31'
year = '2017' # for making a date label

# vn_list = ['temp','salt','phytoplankton','oxygen']
# rows = len(vn_list)

# # figure settings
# fs = 12 # figure font size
# ls = 11 # label size
# ts = 14 # title size

# # initialize arrays for all obs and model values
# # to calculate summary bias and rmse
# allobsvals_SA = np.array([])
# allmodvals_SA = np.array([])
# allobsvals_CT = np.array([])
# allmodvals_CT = np.array([])
# allobsvals_Chl = np.array([])
# allmodvals_Chl = np.array([])
# allobsvals_DO = np.array([])
# allmodvals_DO = np.array([])

jobname = 'twentyoneinlets'
# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')
# Get stations:
sta_dict = job_lists.get_sta_dict(jobname)

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)

# get observations
in_dir = Ldir['parent'] / 'LO_output' / 'obsmod'
in_fn = in_dir / ('multi_ctd_' + year + '.p')
df_dict = pickle.load(open(in_fn, 'rb'))
# only look at ecology stations
source = 'ecology_nc'
df_obs = df_dict['obs'].loc[df_dict['obs'].source==source,:]
df_model = df_dict[gtagex].loc[df_dict[gtagex].source==source,:]

##########################################################
##                 Get station locations                ##
##########################################################

plt.close('all')

# dictionary of ecology stations in terminal inlets
terminl_stns = {
    'elliot': 'ELB015',
    'lynchcove': 'HCB007',
    'lynchcove2': 'HCB004',
    'commencement': 'CMB003',
    'hammersley': 'OAK004',
    # 'totten': 'TOT002',
    'budd': 'BUD005',
    'carr': 'CRR001',
    'sinclair': 'SIN001',
}

# dictionary of ecology station colors
stn_color = {
    'sinclair': '#A8C256',#black',
    'elliot': '#62B6CB',#'red',
    'lynchcove': '#F9627D',#'blue',
    'lynchcove2': '#F9627D',
    'commencement': 'dimgrey',#'purple',
    'hammersley': '#957FEF',#'green',
    'budd': '#476ad1',#'chocolate'
    'carr': '#96031A'
}

# dictionary of shallow inlets
shallow_stns = {
    'hammersley': 'OAK004',
    'budd': 'BUD005',
}

# dictionary of deep inlets
deep_stns = {
    'elliot': 'ELB015',
    'lynchcove': 'HCB007',
    'lynchcove2': 'HCB004',
    'commencement': 'CMB003',
    'carr': 'CRR001',
    'sinclair': 'SIN001',
}

# # get all station names
# all_stns = set(df_obs['name'].values)
# # subtract out ecology stations (not inlet stations remain)
# notinl_stns = [x for x in all_stns if x not in list(terminl_stns.values())]

# get Puget Sound stations that aren't in terminal inlets
PS_stns = ['ADM002','PTH005','ADM001',
           'ADM003','PSB003','EAP001']


# # print station depths
# for station in terminl_stns:
#     print('========='+station+'=========')
#     fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
#     ds_moor = xr.open_dataset(fn)
#     h = ds_moor.h.values
#     print(h)

##########################################################
##          Plotting terminal inlets                    ##
##########################################################

letter = ['(a)','(b)','(c)','(d)','(e)']

# initialize figure
plt.figure(figsize=(12, 8))
ax_map = plt.subplot(1,3,1)
ax_ct = plt.subplot(2,3,2)
ax_sa = plt.subplot(2,3,3)
ax_chl = plt.subplot(2,3,5)
ax_do = plt.subplot(2,3,6)

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
zm[np.transpose(mask_rho) != 0] = -1

# Create map
lon_low = -123.42
lon_high =  -122
lat_low = 46.93
lat_high = 48.46
ax_map.add_patch(Rectangle((lon_low, lat_low), lon_high-lon_low,lat_high-lat_low, facecolor='#EEEEEE'))
ax_map.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

# format
# ax_map.axes.xaxis.set_visible(False)
# ax_map.axes.yaxis.set_visible(False)
ax_map.set_xlim(-123.29, -122.1) # Puget Sound
ax_map.set_ylim(lat_low, lat_high) # Puget Sound
pfun.dar(ax_map)
# # remove border
# for border in ['top','right','bottom','left']:
#     ax_map.spines[border].set_visible(False)

# # add cast locations
# for sta in terminl_stns:
#     sta_lon = sta_dict[sta][0]
#     sta_lat = sta_dict[sta][1]
#     color = stn_color[sta]
#     ecol = '\n('+ terminl_stns[sta] +')'
#     ax_map.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=12,
#     color=color,markeredgecolor='white',markeredgewidth=0.8)
#     # add inlet label
#     # move text to make things more visible
#     if sta in ['sinclair']:
#         lon_off = -0.03
#         lat_off = 0
#         ha = 'right'
#     if sta in ['elliot']:
#         lon_off = 0.05
#         lat_off = 0
#         ha = 'left'
#     if sta in ['commencement']:
#         lat_off = -0.06
#         lon_off = -0.06
#     if sta in ['budd']:
#         lon_off = -0.05
#         lat_off = -0.05
#         ha = 'left'
#     if sta in ['totten']:
#         lon_off = -0.07
#     if sta in ['lynchcove']:
#         lat_off = 0.05
#         lon_off = -0.1
#     if sta in ['lynchcove2']:
#         lat_off = -0.048
#         lon_off = -0.06
#     if sta in ['hammersley']:
#         lat_off = -0.05
#         lon_off = -0.2
#         ha = 'left'
#     if sta in ['carr']:
#         lat_off = 0.05
#     ax_map.text(sta_lon+lon_off,sta_lat+lat_off,sta+ecol,va='center',
#             ha=ha,color=color,fontweight='bold',fontsize=9)
    
# add cast locations
for sta in deep_stns:
    sta_lon = sta_dict[sta][0]
    sta_lat = sta_dict[sta][1]
    color = 'navy'
    ecol = '\n('+ deep_stns[sta] +')'
    ax_map.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=12,
    color=color,markeredgecolor='white',markeredgewidth=0.8)
    # add inlet label
    # move text to make things more visible
    if sta in ['sinclair']:
        lon_off = -0.03
        lat_off = 0
        ha = 'right'
    if sta in ['elliot']:
        lon_off = 0.05
        lat_off = 0
        ha = 'left'
    if sta in ['commencement']:
        lat_off = -0.06
        lon_off = -0.06
    if sta in ['lynchcove']:
        lat_off = 0.05
        lon_off = -0.1
    if sta in ['lynchcove2']:
        lat_off = -0.048
        lon_off = -0.06
    if sta in ['carr']:
        lat_off = 0.05
    ax_map.text(sta_lon+lon_off,sta_lat+lat_off,sta+ecol,va='center',
            ha=ha,color=color,fontweight='bold',fontsize=9)
for sta in shallow_stns:
    sta_lon = sta_dict[sta][0]
    sta_lat = sta_dict[sta][1]
    color = 'deeppink'
    ecol = '\n('+ shallow_stns[sta] +')'
    ax_map.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=12,
    color=color,markeredgecolor='white',markeredgewidth=0.8)
    # add inlet label
    # move text to make things more visible
    if sta in ['budd']:
        lon_off = -0.05
        lat_off = -0.05
        ha = 'left'
    if sta in ['hammersley']:
        lat_off = -0.05
        lon_off = -0.2
        ha = 'left'
    ax_map.text(sta_lon+lon_off,sta_lat+lat_off,sta+ecol,va='center',
            ha=ha,color=color,fontweight='bold',fontsize=9)

# types of stations
t1 = ax_map.text(0.03, 0.98, 'Shallow (mean inlet depth < 10 m)',
            verticalalignment='top', horizontalalignment='left',fontweight='bold',
            transform=ax_map.transAxes, fontsize=10, color = 'deeppink')
t2 = ax_map.text(0.03, 0.94, 'Deep (mean inlet depth > 10 m)',
            verticalalignment='top', horizontalalignment='left',fontweight='bold',
            transform=ax_map.transAxes, fontsize=10, color = 'navy')
t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

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
    # ax.set_facecolor('#EEEEEE')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    # for border in ['top','right','bottom','left']:
    #     ax.spines[border].set_visible(False)
    ax.tick_params(axis='both', labelsize=10)
    ax.set_title(letter[i+1] + ' ' + vars[i],loc='left',fontsize=12)

# # plot property-property plots
# for i,vn in enumerate(vns):
#     for stn,station in enumerate(terminl_stns):
#         # get observational and model information
#         df_ob_stn = df_obs.loc[df_obs.name==terminl_stns[station],:]
#         df_mo_stn = df_model.loc[df_model.name==terminl_stns[station],:]
#         # # get depth and time of observations
#         # z = df_ob_stn['z']
#         # time = [pd.Timestamp(x) for x in df_ob_stn['time']]
#         if vn == 'DO (uM)':
#             obsvals = df_ob_stn[vn].values * 32/1000
#             modvals = df_mo_stn[vn].values * 32/1000
#         else:
#             obsvals = df_ob_stn[vn].values 
#             modvals = df_mo_stn[vn].values
#         # calculate bias and rmse
#         bias = np.nanmean(modvals-obsvals)
#         rmse = np.sqrt(np.nanmean((modvals-obsvals)**2))
#         # get max limits
#         if vn == 'CT':
#             maxval = 25
#         if vn == 'SA':
#             maxval = 40
#         if vn == 'Chl (mg m-3)':
#             maxval = 120
#         if vn == 'DO (uM)':
#             maxval = 17
#         # plot data and y = x line
#         axes[i].scatter(obsvals,modvals, color=stn_color[station], alpha=0.3, s=8, zorder=10, edgecolor='none')
#         axes[i].plot([0,maxval*1.01], [0,maxval*1.01], color='grey', linestyle='-', zorder=5)
#         axes[i].set_ylim([0,maxval])
#         axes[i].set_xlim([0,maxval])
#         axes[i].set_aspect('equal', adjustable='box')
#         axes[i].xaxis.set_major_locator(plt.MaxNLocator(5))
#         axes[i].yaxis.set_major_locator(plt.MaxNLocator(5))
#         # add bias and rmse label
#         t01 = axes[i].text(0.17, 0.95, 'Bias', fontweight='bold',
#             verticalalignment='top', horizontalalignment='right',
#             transform=axes[i].transAxes, fontsize=10, color = 'black')
#         t02 = axes[i].text(0.35, 0.95, 'RMSE', fontweight='bold',
#             verticalalignment='top', horizontalalignment='right',
#             transform=axes[i].transAxes, fontsize=10, color = 'black')
#         t1 = axes[i].text(0.17, 0.89-0.06*stn, '{:.2f}'.format(round(bias,2)),
#             verticalalignment='top', horizontalalignment='right',fontweight='bold',
#             transform=axes[i].transAxes, fontsize=10, color = stn_color[station])
#         t2 = axes[i].text(0.35, 0.89-0.06*stn, '{:.2f}'.format(round(rmse,2)),
#             verticalalignment='top', horizontalalignment='right',fontweight='bold',
#             transform=axes[i].transAxes, fontsize=10, color = stn_color[station])
#         t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
#         t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

# plot property-property plots for deep terminal inlets
for i,vn in enumerate(vns):
    # initialize arrays to save values
    deep_obs = np.array([])
    deep_mod = np.array([])
    for stn,station in enumerate(deep_stns):
        # get observational and model information
        df_ob_stn = df_obs.loc[df_obs.name==deep_stns[station],:]
        df_mo_stn = df_model.loc[df_model.name==deep_stns[station],:]
        if vn == 'DO (uM)':
            obsvals = df_ob_stn[vn].values * 32/1000
            modvals = df_mo_stn[vn].values * 32/1000
        else:
            obsvals = df_ob_stn[vn].values 
            modvals = df_mo_stn[vn].values
        # get max limits
        if vn == 'CT':
            maxval = 25
        if vn == 'SA':
            maxval = 40
        if vn == 'Chl (mg m-3)':
            maxval = 120
        if vn == 'DO (uM)':
            maxval = 17
        # plot data and y = x line
        axes[i].scatter(obsvals,modvals, color='navy', alpha=0.3, s=6, zorder=10, edgecolor='none')
        axes[i].plot([0,maxval*1.01], [0,maxval*1.01], color='grey', linestyle='-', zorder=5)
        axes[i].set_ylim([0,maxval])
        axes[i].set_xlim([0,maxval])
        axes[i].set_aspect('equal', adjustable='box')
        axes[i].xaxis.set_major_locator(plt.MaxNLocator(5))
        axes[i].yaxis.set_major_locator(plt.MaxNLocator(5))
        # add data to arrays
        deep_obs = np.concatenate((deep_obs,obsvals))
        deep_mod = np.concatenate((deep_mod,modvals))
    # # add bias and rmse label
    # calculate bias and rmse
    bias = np.nanmean(deep_mod-deep_obs)
    rmse = np.sqrt(np.nanmean((deep_mod-deep_obs)**2))
    t01 = axes[i].text(0.17, 0.95, 'Bias', fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        transform=axes[i].transAxes, fontsize=10, color = 'black')
    t02 = axes[i].text(0.35, 0.95, 'RMSE', fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        transform=axes[i].transAxes, fontsize=10, color = 'black')
    t1 = axes[i].text(0.17, 0.89, '{:.2f}'.format(round(bias,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'navy')
    t2 = axes[i].text(0.35, 0.89, '{:.2f}'.format(round(rmse,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'navy')
    t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

# plot property-property plots for shallow terminal inlets
for i,vn in enumerate(vns):
    # initialize arrays to save values
    shallow_obs = np.array([])
    shallow_mod = np.array([])
    for stn,station in enumerate(shallow_stns):
        # get observational and model information
        df_ob_stn = df_obs.loc[df_obs.name==shallow_stns[station],:]
        df_mo_stn = df_model.loc[df_model.name==shallow_stns[station],:]
        if vn == 'DO (uM)':
            obsvals = df_ob_stn[vn].values * 32/1000
            modvals = df_mo_stn[vn].values * 32/1000
        else:
            obsvals = df_ob_stn[vn].values 
            modvals = df_mo_stn[vn].values
        # get max limits
        if vn == 'CT':
            maxval = 25
        if vn == 'SA':
            maxval = 40
        if vn == 'Chl (mg m-3)':
            maxval = 120
        if vn == 'DO (uM)':
            maxval = 17
        # plot data and y = x line
        axes[i].scatter(obsvals,modvals, color='deeppink', alpha=0.5, s=6, zorder=10, edgecolor='none')
        axes[i].plot([0,maxval*1.01], [0,maxval*1.01], color='grey', linestyle='-', zorder=5)
        axes[i].set_ylim([0,maxval])
        axes[i].set_xlim([0,maxval])
        axes[i].set_aspect('equal', adjustable='box')
        axes[i].xaxis.set_major_locator(plt.MaxNLocator(5))
        axes[i].yaxis.set_major_locator(plt.MaxNLocator(5))
        # add data to arrays
        shallow_obs = np.concatenate((shallow_obs,obsvals))
        shallow_mod = np.concatenate((shallow_mod,modvals))
    # # add bias and rmse label
    # calculate bias and rmse
    bias = np.nanmean(shallow_mod-shallow_obs)
    rmse = np.sqrt(np.nanmean((shallow_mod-shallow_obs)**2))
    t1 = axes[i].text(0.17, 0.83, '{:.2f}'.format(round(bias,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'deeppink')
    t2 = axes[i].text(0.35, 0.83, '{:.2f}'.format(round(rmse,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'deeppink')
    t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

    axes[i].set_ylabel('Modeled',fontsize=12)
    axes[i].set_xlabel('Observed',fontsize=12)

plt.suptitle('Model vs. Observations property-property plots (terminal inlets 2017)', fontsize=15, x=0.65)
plt.subplots_adjust(wspace=-0.2)      
plt.tight_layout()


#####################################################################################################################################

##########################################################
##       Plotting all Ecology CTD locations             ##
##########################################################

letter = ['(a)','(b)','(c)','(d)','(e)']

# initialize figure
plt.figure(figsize=(12, 8))
ax_map = plt.subplot(1,3,1)
ax_ct = plt.subplot(2,3,2)
ax_sa = plt.subplot(2,3,3)
ax_chl = plt.subplot(2,3,5)
ax_do = plt.subplot(2,3,6)

##########################################################
##              Puget Sound location map                ##
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
zm[np.transpose(mask_rho) != 0] = -1

# Create map
lon_low = -123.42
lon_high =  -122
lat_low = 46.93
lat_high = 48.46
ax_map.add_patch(Rectangle((lon_low, lat_low), lon_high-lon_low,lat_high-lat_low, facecolor='white'))#, facecolor='#EEEEEE'))
ax_map.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

# format
# ax_map.axes.xaxis.set_visible(False)
# ax_map.axes.yaxis.set_visible(False)
ax_map.tick_params(axis='x', labelrotation=30)
ax_map.set_ylabel('Latitude',fontsize=12)
ax_map.set_xlabel('Longitude',fontsize=12)
ax_map.set_xlim(-123.29, -122.1) # Puget Sound
ax_map.set_ylim(lat_low, lat_high) # Puget Sound
pfun.dar(ax_map)
# # remove border
# for border in ['top','right','bottom','left']:
#     ax_map.spines[border].set_visible(False)
    

# add terminal inlets cast locations
for sta in terminl_stns:
    sta_lon = sta_dict[sta][0]
    sta_lat = sta_dict[sta][1]
    color = stn_color[sta]
    ecol = '\n('+ terminl_stns[sta] +')'
    ax_map.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=8,
    color='deeppink',markeredgecolor='white',markeredgewidth=0.8)
    # add inlet label
    lon_off = 0.12
    lat_off = 0
    ha = 'center'
    ax_map.text(sta_lon+lon_off,sta_lat+lat_off,terminl_stns[sta],va='center',
                ha=ha,color='deeppink',fontsize=9, fontweight='bold')
    
# add the remaining cast locations
for sta in PS_stns:
    # sta_lon = sta_dict[sta][0]
    # sta_lat = sta_dict[sta][1]
    sta_lat = df_obs.loc[df_obs['name'] == sta, 'lat'].values[0]
    sta_lon = df_obs.loc[df_obs['name'] == sta, 'lon'].values[0]
    ax_map.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=8,
    color='navy',markeredgecolor='white',markeredgewidth=0.8)
    # add inlet label
    lon_off = 0.12
    lat_off = 0
    ha = 'center'
    ax_map.text(sta_lon+lon_off,sta_lat+lat_off,sta,va='center',
                ha=ha,color='navy',fontsize=9, fontweight='bold')

# types of stations
t1 = ax_map.text(0.03, 0.98, 'Terminal Inlet stations',
            verticalalignment='top', horizontalalignment='left',fontweight='bold',
            transform=ax_map.transAxes, fontsize=10, color = 'deeppink')
t2 = ax_map.text(0.03, 0.94, 'Main Basin stations',
            verticalalignment='top', horizontalalignment='left',fontweight='bold',
            transform=ax_map.transAxes, fontsize=10, color = 'navy')
t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

ax_map.set_title('(a) Inlet obs comparison locations',fontsize=12,loc='left')

##########################################################
##         Puget Sound Property-property plots          ##
##########################################################

axes = [ax_ct,ax_sa,ax_chl,ax_do]
vars = [r'Cons Temp [$^\circ$C]',r'Abs Salinity [g kg$^{-1}$]',
        r'Chl [mg m$^{-3}$]', r'DO [mg L$^{-1}$]']
vns = ['CT','SA','Chl (mg m-3)','DO (uM)']

# format grid
for i,ax in enumerate(axes):
    # ax.set_facecolor('#EEEEEE')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    # ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    # for border in ['top','right','bottom','left']:
    #     ax.spines[border].set_visible(False)
    ax.tick_params(axis='both', labelsize=10)
    ax.set_title(letter[i+1] + ' ' + vars[i],loc='left',fontsize=12)

# plot property-property plots for terminal inlets
for i,vn in enumerate(vns):
    # initialize arrays to save values
    terminl_obs = np.array([])
    terminl_mod = np.array([])
    for stn,station in enumerate(terminl_stns):
        # get observational and model information
        df_ob_stn = df_obs.loc[df_obs.name==terminl_stns[station],:]
        df_mo_stn = df_model.loc[df_model.name==terminl_stns[station],:]
        if vn == 'DO (uM)':
            obsvals = df_ob_stn[vn].values * 32/1000
            modvals = df_mo_stn[vn].values * 32/1000
        else:
            obsvals = df_ob_stn[vn].values 
            modvals = df_mo_stn[vn].values
        # get max limits
        if vn == 'CT':
            maxval = 25
        if vn == 'SA':
            maxval = 40
        if vn == 'Chl (mg m-3)':
            maxval = 120
        if vn == 'DO (uM)':
            maxval = 17
        # plot data and y = x line
        axes[i].scatter(obsvals,modvals, color='deeppink', alpha=0.3, s=6, zorder=10, edgecolor='none')
        axes[i].plot([0,maxval*1.01], [0,maxval*1.01], color='grey', linestyle='-', zorder=5)
        axes[i].set_ylim([0,maxval])
        axes[i].set_xlim([0,maxval])
        axes[i].set_aspect('equal', adjustable='box')
        axes[i].xaxis.set_major_locator(plt.MaxNLocator(5))
        axes[i].yaxis.set_major_locator(plt.MaxNLocator(5))
        # add data to arrays
        terminl_obs = np.concatenate((terminl_obs,obsvals))
        terminl_mod = np.concatenate((terminl_mod,modvals))
    # # add bias and rmse label
    # calculate bias and rmse
    bias = np.nanmean(terminl_mod-terminl_obs)
    rmse = np.sqrt(np.nanmean((terminl_mod-terminl_obs)**2))
    nse = 1 - (
        np.nansum( (terminl_mod-terminl_obs)**2 )/
        np.nansum( (terminl_obs- np.nanmean(terminl_obs) )**2 )
        )
    # nse = 1 - (
    #     np.nansum( np.abs(terminl_mod-terminl_obs) )/
    #     np.nansum( np.abs(terminl_obs- np.nanmean(terminl_obs) ) )
    #     )
    t01 = axes[i].text(0.17, 0.95, 'Bias', fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        transform=axes[i].transAxes, fontsize=10, color = 'black')
    t02 = axes[i].text(0.37, 0.95, 'RMSE', fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        transform=axes[i].transAxes, fontsize=10, color = 'black')
    t03 = axes[i].text(0.57, 0.95, 'NSE', fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        transform=axes[i].transAxes, fontsize=10, color = 'black')
    t1 = axes[i].text(0.17, 0.89, '{:.2f}'.format(round(bias,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'deeppink')
    t2 = axes[i].text(0.37, 0.89, '{:.2f}'.format(round(rmse,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'deeppink')
    t3 = axes[i].text(0.57, 0.89, '{:.2f}'.format(round(nse,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'deeppink')
    t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    t3.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

    axes[i].set_ylabel('Modeled',fontsize=12)
    axes[i].set_xlabel('Observed',fontsize=12)


    # initialize arrays to save values
    PSstns_obs = np.array([])
    PSstns_mod = np.array([])
    for stn,station in enumerate(PS_stns):
        # get observational and model information
        df_ob_stn = df_obs.loc[df_obs.name==station,:]
        df_mo_stn = df_model.loc[df_model.name==station,:]
        # # get depth and time of observations
        # z = df_ob_stn['z']
        # time = [pd.Timestamp(x) for x in df_ob_stn['time']]
        if vn == 'DO (uM)':
            obsvals = df_ob_stn[vn].values * 32/1000
            modvals = df_mo_stn[vn].values * 32/1000
        else:
            obsvals = df_ob_stn[vn].values 
            modvals = df_mo_stn[vn].values
        # calculate bias and rmse
        bias = np.nanmean(modvals-obsvals)
        rmse = np.sqrt(np.nanmean((modvals-obsvals)**2))
        # get max limits
        if vn == 'CT':
            maxval = 25
        if vn == 'SA':
            maxval = 40
        if vn == 'Chl (mg m-3)':
            maxval = 120
        if vn == 'DO (uM)':
            maxval = 17
        # plot data and y = x line
        axes[i].scatter(obsvals,modvals, color='navy', alpha=0.1, s=6, zorder=10, edgecolor='none')
        axes[i].plot([0,maxval*1.01], [0,maxval*1.01], color='grey', linestyle='-', zorder=5)
        axes[i].set_ylim([0,maxval])
        axes[i].set_xlim([0,maxval])
        axes[i].set_aspect('equal', adjustable='box')
        axes[i].xaxis.set_major_locator(plt.MaxNLocator(5))
        axes[i].yaxis.set_major_locator(plt.MaxNLocator(5))
    # add data to arrays
        PSstns_obs = np.concatenate((PSstns_obs,obsvals))
        PSstns_mod = np.concatenate((PSstns_mod,modvals))
    # # add bias and rmse label
    # calculate bias and rmse and nash-sutcliffe model efficiency index
    bias = np.nanmean(PSstns_mod-PSstns_obs)
    rmse = np.sqrt(np.nanmean((PSstns_mod-PSstns_obs)**2))

    nse = 1 - (
        np.nansum((PSstns_mod-PSstns_obs)**2)/
        np.nansum((PSstns_obs-np.nanmean(PSstns_obs))**2)
        )
    
    if vn == 'SA':
        print((PSstns_mod-PSstns_obs)[210:215])
        # print(np.nansum((PSstns_obs-np.nanmean(PSstns_obs))**2))
    
    # nse = 1 - (
    #     np.nansum( np.abs(PSstns_mod-PSstns_obs) )/
    #     np.nansum( np.abs(PSstns_obs-np.nanmean(PSstns_obs)) )
    #     )
    t1 = axes[i].text(0.17, 0.83, '{:.2f}'.format(round(bias,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'navy')
    t2 = axes[i].text(0.37, 0.83, '{:.2f}'.format(round(rmse,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'navy')
    t3 = axes[i].text(0.57, 0.83, '{:.2f}'.format(round(nse,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=10, color = 'navy')
    t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    t3.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))


plt.suptitle('Model vs. Observations property-property plots (Puget Sound 2017)', fontsize=15, x=0.65)
plt.subplots_adjust(wspace=-0.2)      
plt.tight_layout()


############################################################################################################################
## Selected time series

# create time vector
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

# dictionary of selected stations in terminal inlets for time series
selected_stns = {
    'elliot': 'ELB015',
    'lynchcove': 'HCB007',
    # 'hammersley': 'OAK004',
}

# letter = ['(a)','(b)','(c)',
#           '(d)','(e)','(f)',
#           '(g)','(h)','(i)',
#           '(j)','(k)','(l)']

letter = ['(a)','(b)','(c)',
          '(d)','(e)','(f)',
          '(g)','(h)']

# initialize figure
fig, ax = plt.subplots(4,2,figsize = (10,9.5),sharey='row',sharex='col')

# add data
for i,vn in enumerate(vns):
    for stn, station in enumerate(selected_stns):
        # add a title
        if i == 0:
            ax[i,stn].set_title(station + ' ('+selected_stns[station]+')',fontsize=12,fontweight='bold')
        # get axis
        axis = ax[i,stn]

        # calculate lat/lon for station
        lon = sta_dict[station][0]
        lat = sta_dict[station][1]

        # get observational information
        df_ob_stn = df_obs.loc[df_obs.name==selected_stns[station],:]
        df_mo_stn = df_model.loc[df_model.name==selected_stns[station],:]
        # get depth and time of observations
        z = df_ob_stn['z']
        time = [pd.Timestamp(x) for x in df_ob_stn['time']]

        # get station depth
        fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds_moor = xr.open_dataset(fn)
        h = ds_moor.h.values
        # get depth values
        z_w   = ds_moor['z_w'].transpose()   # depth of w-velocities
        z_min = np.min(z_w.values)
        
        # get all variables and convert to correct units
        if vn == 'DO (uM)': # already in uM
            ds_moor['DO (mg/L)'] = ds_moor['oxygen'] * 32/1000
            val = ds_moor['DO (mg/L)'].transpose()
        if vn == 'Chl (mg m-3)': # need to multiply by 2.5
            ds_moor['Chl (mg m-3)'] = ds_moor['phytoplankton'] * 2.5
            val = ds_moor['Chl (mg m-3)'].transpose()
        if vn == 'CT': # convert to conservative temprature
            ds_moor['SA'] = gsw.conversions.SA_from_SP(ds_moor['salt'], ds_moor['z_rho'], lon, lat)
            ds_moor['CT'] = gsw.conversions.CT_from_pt(ds_moor['SA'], ds_moor['temp'])
            val = ds_moor['CT'].transpose()
        if vn == 'SA': # convert to absolute salinity
            ds_moor['SA'] = gsw.conversions.SA_from_SP(ds_moor['salt'], ds_moor['z_rho'], lon, lat)
            val = ds_moor['SA'].transpose()
            # print(val.values[0,:])

        # if deeper than 10 m, split into top 5 m and bottom 5 m layer
        d = 10
        if h > d:

            # model
            surf_z_ds = ds_moor.where((ds_moor.z_w >= -(d/2)))
            bott_z_ds = ds_moor.where((ds_moor.z_w <= (-1*h) + (d/2)))
            surf_mod_ds = ds_moor.where((ds_moor.z_rho >= -(d/2)))
            bott_mod_ds = ds_moor.where((ds_moor.z_rho <= (-1*h) + (d/2)))
            # get depth of each vertical layer
            z_thick_surf = np.diff(surf_z_ds['z_w'],axis=1) # 365,30
            z_thick_bott = np.diff(bott_z_ds['z_w'],axis=1) # 365,30
            # get model data 
            if vn == 'DO (uM)':
                surf_mod = surf_mod_ds['DO (mg/L)'].values
                bott_mod = bott_mod_ds['DO (mg/L)'].values
            else:
                surf_mod = surf_mod_ds[vn].values
                bott_mod = bott_mod_ds[vn].values
            # check if thickness array is all nan
            # (this occurs if the first z_w is already greater than the threshold, so we don't have two z_w values to diff)
            # in which case, average value is just the single value
            if np.isnan(z_thick_surf).all():
                surf_mod_avg = np.nansum(surf_mod,axis=1)
            # take weighted average given depth
            # first, multiply by thickness of each layer and sum in z, then divide by layer thickness
            else:
                surf_mod_avg = np.nansum(surf_mod * z_thick_surf, axis=1)/np.nansum(z_thick_surf,axis=1)
            # repeat check for surface
            if np.isnan(z_thick_bott).all():
                bott_mod_avg = np.nansum(bott_mod,axis=1)
            else:
                bott_mod_avg = np.nansum(bott_mod * z_thick_bott, axis=1)/np.nansum(z_thick_bott,axis=1)
            # plot model output
            axis.plot(dates_local, surf_mod_avg, color='deeppink', linewidth=2, alpha=0.5, zorder=5)
            axis.plot(dates_local, bott_mod_avg, color='navy', linewidth=2, alpha=0.5, zorder=5,label='model')

            # observation
            surf_obs_df = df_ob_stn.where(z >= -(d/2))
            bott_obs_df = df_ob_stn.where(z <= (-1*h) + (d/2))
            # get average value
            if vn == 'DO (uM)':
                surf_obs_avg = surf_obs_df.groupby('time')[vn].mean() * 32/1000
                bott_obs_avg = bott_obs_df.groupby('time')[vn].mean() * 32/1000
            else:
                surf_obs_avg = surf_obs_df.groupby('time')[vn].mean()
                bott_obs_avg = bott_obs_df.groupby('time')[vn].mean()
            # plot observations
            unique_time = [pd.Timestamp(x) for x in df_ob_stn['time'].unique()] # one point per timestampe
            axis.scatter(unique_time, surf_obs_avg, color='deeppink', s=20,zorder=10)
            axis.scatter(unique_time, bott_obs_avg, color='navy', s=20,zorder=10, label='obs')


            # label
            if i == 0 and stn == 0:
                axis.text(0.75, 0.9, 'surface {} m'.format(str(round(d/2))),
                    verticalalignment='top', horizontalalignment='left',
                    transform=axis.transAxes, fontsize=12, color = 'deeppink',fontweight='bold')
                axis.text(0.75, 0.8, 'bottom {} m'.format(str(round(d/2))),
                    verticalalignment='top', horizontalalignment='left',
                    transform=axis.transAxes, fontsize=12, color = 'navy',fontweight='bold')
                axis.legend(loc='lower right',fontsize=11, frameon=False, handletextpad=0.1,
                            handlelength=1, markerfirst = False)
                

        # otherwise just one layer
        else:

            # get depth of each vertical layer
            z_thick = np.diff(z_w,axis=0).transpose() # 365,30

            # model
            mod_ds = ds_moor
            # get model data 
            if vn == 'DO (uM)':
                mod = mod_ds['DO (mg/L)']
            else:
                mod = mod_ds[vn]
            # take weighted average given depth
            # first, multiply by thickness of each layer and sum in z, then divide by water column depth
            mod_avg = np.nansum(mod * z_thick, axis=1)/np.nansum(z_thick,axis=1)
            # plot model output
            axis.plot(dates_local, mod_avg, color='mediumorchid', linewidth=2, alpha=0.5, zorder=5, label='model')

            # observation
            obs_df = df_ob_stn
            # get average value
            if vn == 'DO (uM)':
                obs_avg = obs_df.groupby('time')[vn].mean() * 32/1000
            else:
                obs_avg = obs_df.groupby('time')[vn].mean()
            # plot observations
            unique_time = [pd.Timestamp(x) for x in df_ob_stn['time'].unique()]
            # axis.scatter(time, obs_df[var], color='k',s=3)
            axis.scatter(unique_time, obs_avg, color='mediumorchid', s=20, zorder=10, label='obs')

        # set y-axis title and limits
        if i == 0:
            axis.set_ylim([0,25])
        if i == 1:
            axis.set_ylim([0,36])
        if i == 2:
            axis.set_ylim([0,25])
        if i == 3:
            axis.set_ylim([0,14])
        if stn == 0:
            axis.set_ylabel(vars[i],fontsize=12)


# format grids and add titles
ax = ax.ravel()
for i,axis in enumerate(ax):
    axis.set_xlim([dates_local[0],dates_local[-1]])
    # axis.set_facecolor('#EEEEEE')
    # axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    axis.grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    axis.tick_params(axis='x', labelrotation=30)
    # for border in ['top','right','bottom','left']:
    #     axis.spines[border].set_visible(False)
    axis.text(0.03, 0.95, letter[i],
            verticalalignment='top', horizontalalignment='left',
            transform=axis.transAxes, fontsize=12, color = 'k')

plt.suptitle('Model vs. Observations time series (2017)', fontsize=15)
plt.subplots_adjust(wspace=-0.2)      
plt.tight_layout()