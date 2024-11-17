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

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures' / 'twentyone' / '2017_inlet_obs_comparison'
Lfun.make_dir(out_dir)

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
    'sinclair': 'SIN001',
    'elliot': 'ELB015',
    'lynchcove': 'HCB007',
    'lynchcove2': 'HCB004',
    'commencement': 'CMB003',
    'hammersley': 'OAK004',
    # 'totten': 'TOT002',
    'budd': 'BUD005',
    'carr': 'CRR001'
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

# # get all station names
# all_stns = set(df_obs['name'].values)
# # subtract out ecology stations (not inlet stations remain)
# notinl_stns = [x for x in all_stns if x not in list(terminl_stns.values())]

# get Puget Sound stations that aren't in terminal inlets
PS_stns = ['ADM002','PTH005','ADM001',
           'ADM003','PSB003','EAP001']

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
# define colormaps for coloring the inlets
cmaps = ['Greys','Reds','Blues','Purples','Greens','autumn']
# # get terminal inlet locations
# for stn,station in enumerate(ecol_stn): # stations: 
#     # get segment information
#     seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
#     seg_df = pd.read_pickle(seg_name)
#     ji_list = seg_df[station+'_p']['ji_list']
#     jj = [x[0] for x in ji_list]
#     ii = [x[1] for x in ji_list]
#     # set everything to nan that is not in the inlet
#     # first make everything a nan
#     inlet_loc = np.full(zm.shape, np.nan) 
#     # set values of 1 for everything that is in the inlet
#     inlet_loc[jj,ii] = 70
#     # add inlet locations
#     ax_map.pcolormesh(plon, plat, inlet_loc, linewidth=0.5, vmin=0, vmax=100, cmap=plt.get_cmap(cmaps[stn]))

# format
ax_map.axes.xaxis.set_visible(False)
ax_map.axes.yaxis.set_visible(False)
ax_map.set_xlim(-123.29, -122.1) # Puget Sound
ax_map.set_ylim(lat_low, lat_high) # Puget Sound
pfun.dar(ax_map)
# remove border
for border in ['top','right','bottom','left']:
    ax_map.spines[border].set_visible(False)

# add cast locations
for sta in terminl_stns:
    sta_lon = sta_dict[sta][0]
    sta_lat = sta_dict[sta][1]
    color = stn_color[sta]
    ecol = '\n('+ terminl_stns[sta] +')'
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
    if sta in ['budd']:
        lon_off = -0.05
        lat_off = -0.05
        ha = 'left'
    if sta in ['totten']:
        lon_off = -0.07
    if sta in ['lynchcove']:
        lat_off = 0.05
        lon_off = -0.1
    if sta in ['lynchcove2']:
        lat_off = -0.048
        lon_off = -0.06
    if sta in ['hammersley']:
        lat_off = -0.05
        lon_off = -0.2
        ha = 'left'
    if sta in ['carr']:
        lat_off = 0.05
    ax_map.text(sta_lon+lon_off,sta_lat+lat_off,sta+ecol,va='center',
            ha=ha,color=color,fontweight='bold',fontsize=9)
    
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
    ax.set_facecolor('#EEEEEE')
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    ax.tick_params(axis='both', labelsize=10)
    ax.set_title(letter[i+1] + ' ' + vars[i],loc='left',fontsize=12)

# plot property-property plots
for i,vn in enumerate(vns):
    for stn,station in enumerate(terminl_stns):
        if station == 'carr':
            continue
        else:
            # get observational and model information
            df_ob_stn = df_obs.loc[df_obs.name==terminl_stns[station],:]
            df_mo_stn = df_model.loc[df_model.name==terminl_stns[station],:]
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
            axes[i].scatter(obsvals,modvals, color=stn_color[station], alpha=0.6, s=8, zorder=10)
            axes[i].plot([0,maxval*1.01], [0,maxval*1.01], color='grey', linestyle='-', zorder=5)
            axes[i].set_ylim([0,maxval])
            axes[i].set_xlim([0,maxval])
            axes[i].set_aspect('equal', adjustable='box')
            axes[i].xaxis.set_major_locator(plt.MaxNLocator(5))
            axes[i].yaxis.set_major_locator(plt.MaxNLocator(5))
            # add bias and rmse label
            t01 = axes[i].text(0.17, 0.95, 'Bias', fontweight='bold',
                verticalalignment='top', horizontalalignment='right',
                transform=axes[i].transAxes, fontsize=10, color = 'black')
            t02 = axes[i].text(0.35, 0.95, 'RMSE', fontweight='bold',
                verticalalignment='top', horizontalalignment='right',
                transform=axes[i].transAxes, fontsize=10, color = 'black')
            t1 = axes[i].text(0.17, 0.89-0.06*stn, '{:.2f}'.format(round(bias,2)),
                verticalalignment='top', horizontalalignment='right',fontweight='bold',
                transform=axes[i].transAxes, fontsize=10, color = stn_color[station])
            t2 = axes[i].text(0.35, 0.89-0.06*stn, '{:.2f}'.format(round(rmse,2)),
                verticalalignment='top', horizontalalignment='right',fontweight='bold',
                transform=axes[i].transAxes, fontsize=10, color = stn_color[station])
            t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
            t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))


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
ax_map.add_patch(Rectangle((lon_low, lat_low), lon_high-lon_low,lat_high-lat_low, facecolor='#EEEEEE'))
ax_map.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

# # get terminal inlet locations
# for stn,station in enumerate(sta_dict): # stations: 
#     if station != 'lynchcove2':
#         # get segment information
#         seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
#         seg_df = pd.read_pickle(seg_name)
#         ji_list = seg_df[station+'_p']['ji_list']
#         jj = [x[0] for x in ji_list]
#         ii = [x[1] for x in ji_list]
#         # set everything to nan that is not in the inlet
#         # first make everything a nan
#         inlet_loc = np.full(zm.shape, np.nan) 
#         # set values of 1 for everything that is in the inlet
#         inlet_loc[jj,ii] = 20
#         # add inlet locations
#         # plt.pcolormesh(plon, plat, inlet_loc, linewidth=0.5, vmin=0, vmax=65, cmap=plt.get_cmap('coolwarm'))
#         ax_map.pcolormesh(plon, plat, inlet_loc, linewidth=0.5, vmin=0, vmax=35, cmap=plt.get_cmap(cmocean.cm.ice))

# format
ax_map.axes.xaxis.set_visible(False)
ax_map.axes.yaxis.set_visible(False)
ax_map.set_xlim(-123.29, -122.1) # Puget Sound
ax_map.set_ylim(lat_low, lat_high) # Puget Sound
pfun.dar(ax_map)
# remove border
for border in ['top','right','bottom','left']:
    ax_map.spines[border].set_visible(False)
    

# add terminal inlets cast locations
for sta in terminl_stns:
    sta_lon = sta_dict[sta][0]
    sta_lat = sta_dict[sta][1]
    color = stn_color[sta]
    ecol = '\n('+ terminl_stns[sta] +')'
    ax_map.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=8,
    color='k',markeredgecolor='white',markeredgewidth=0.8)
    
# add the remaining cast locations
for sta in PS_stns:
    # sta_lon = sta_dict[sta][0]
    # sta_lat = sta_dict[sta][1]
    sta_lat = df_obs.loc[df_obs['name'] == sta, 'lat'].values[0]
    sta_lon = df_obs.loc[df_obs['name'] == sta, 'lon'].values[0]
    # print('\n==============='+sta+'================')
    # print('lat: {}'.format(sta_lat))
    # print('lon: {}'.format(sta_lon))
    ax_map.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=8,
    color='orchid',markeredgecolor='white',markeredgewidth=0.8)
    # add inlet label
    lon_off = 0.14
    lat_off = 0
    ha = 'center'
    ax_map.text(sta_lon+lon_off,sta_lat+lat_off,sta,va='center',
                ha=ha,color='mediumorchid',fontsize=9, fontweight='bold')

# types of stations
t1 = ax_map.text(0.03, 0.98, 'Terminal Inlet stations',
            verticalalignment='top', horizontalalignment='left',fontweight='bold',
            transform=ax_map.transAxes, fontsize=10, color = 'k')
t2 = ax_map.text(0.03, 0.94, 'Main Basin stations',
            verticalalignment='top', horizontalalignment='left',fontweight='bold',
            transform=ax_map.transAxes, fontsize=10, color = 'orchid')
t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

ax_map.set_title('(a) Inlet obs comparison locations',fontsize=12,loc='left')

# ##########################################################
# ##         Puget Sound Property-property plots          ##
# ##########################################################

# axes = [ax_ct,ax_sa,ax_chl,ax_do]
# vars = [r'Cons Temp [$^\circ$C]',r'Abs Salinity [g kg$^{-1}$]',
#         r'Chl [mg m$^{-3}$]', r'DO [mg L$^{-1}$]']
# vns = ['CT','SA','Chl (mg m-3)','DO (uM)']

# # format grid
# for i,ax in enumerate(axes):
#     ax.set_facecolor('#EEEEEE')
#     ax.tick_params(axis='x', labelrotation=30)
#     ax.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#     for border in ['top','right','bottom','left']:
#         ax.spines[border].set_visible(False)
#     ax.tick_params(axis='both', labelsize=10)
#     ax.set_title(letter[i+1] + ' ' + vars[i],loc='left',fontsize=12)

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
#         axes[i].scatter(obsvals,modvals, color=stn_color[station], alpha=0.3, s=6, zorder=10)
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


plt.suptitle('Model vs. Observations property-property plots (Puget Sound 2017)', fontsize=15, x=0.65)
plt.subplots_adjust(wspace=-0.2)      
plt.tight_layout()
