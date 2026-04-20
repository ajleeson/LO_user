"""
Bottom minus top density differences
"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
import numpy as np
import xarray as xr
from scipy.stats import pearsonr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import gsw
import pinfo
import pickle
from matplotlib.colors import LinearSegmentedColormap

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun


Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t1_x11ab'
jobname = 'twentyoneinlets'
startdate = '2017.01.01'
enddate = '2017.12.31'
year = '2017' # for making a date label

# figure settings
fs = 12 # figure font size
ls = 11 # label size
ts = 12 # title size

# # hypoxic season
# start = '08-01'
# end = '09-30'

# annual
start = '01-01'
end = '12-31'

inlet_label = False

# basin color list

color_list = ['limegreen',  # whidbey
            'hotpink',    # hood canal
            'deepskyblue',  # main basin
            'blueviolet',       # south sound
            'black']        # admiralty inlet

##########################################################
##                  Helper Functions                    ##
##########################################################

def get_surf_average(vn,d,ds):
    '''
    calculate average value of vn in surface d meters
    '''
    z_ds = ds.where((ds.z_w >= -d))
    var_ds = ds.where((ds.z_rho >= -d))
    # get thickness of each vertical layer
    z_thick = np.diff(z_ds['z_w'],axis=1) # 365,30
    # get variable data 
    var = var_ds[vn].values
    # check if thickness array is all nan
    # (this occurs if the first z_w is already greater than the threshold, so we don't have two z_w values to diff)
    # in which case, average value is just the single value
    if np.isnan(z_thick).all():
        val = np.nansum(var,axis=1)
    # take weighted average given depth
    # first, multiply by thickness of each layer and sum in z, then divide by layer thickness
    else:
        val = np.nansum(var * z_thick, axis=1)/np.nansum(z_thick,axis=1)
    return val

def get_bott_average(vn,d,h,ds):
    '''
    calculate average value of vn in bottom d meters
    '''
    z_ds = ds.where((ds.z_w <= (-1*h) + d))
    var_ds = ds.where((ds.z_rho <= (-1*h) + d))
    # get thickness of each vertical layer
    z_thick = np.diff(z_ds['z_w'],axis=1) # 365,30
    # get variable data 
    var = var_ds[vn].values
    # check if thickness array is all nan
    # (this occurs if the first z_w is already greater than the threshold, so we don't have two z_w values to diff)
    # in which case, average value is just the single value
    if np.isnan(z_thick).all():
        val = np.nansum(var,axis=1)
    # take weighted average given depth
    # first, multiply by thickness of each layer and sum in z, then divide by layer thickness
    else:
        val = np.nansum(var * z_thick, axis=1)/np.nansum(z_thick,axis=1)
    return val

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
del sta_dict['lynchcove2']
# remove shallow inlets (< 10 m deep)
del sta_dict['hammersley']
del sta_dict['henderson']
del sta_dict['oak']
del sta_dict['totten']
del sta_dict['similk']
del sta_dict['budd']
del sta_dict['eld']
del sta_dict['killsut']

# inlets in different basins
whidbey = ['similk','oak','crescent','penn','portsusan','holmes']
hoodcanal = ['dabob','lynchcove']
mainbasin = ['dyes','sinclair','elliot','quartermaster','commencement']
southsound = ['case','carr','hammersley','totten','eld','budd','henderson']
admiraltysill = ['killsut']
    
# Calculate average bottom DO and assign colors to different basins
bott_DO = np.zeros(len(sta_dict))
colors = np.zeros(len(sta_dict))
for i,station in enumerate(sta_dict):
    # download .nc files
    fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))
    # get average bottom DO
    bott_DO_all = ds['oxygen'].values[:,0] * pinfo.fac_dict['oxygen'] # mg/L
    # average over hypoxic season
    bott_DO_avg = np.nanmean(bott_DO_all)
    bott_DO[i] = bott_DO_avg
    # print(station)
    # print(bott_DO_avg)
    # assign colors to different basins
    if station in whidbey:
        colors[i] = 1
    elif station in hoodcanal:
        colors[i] = 2
    elif station in mainbasin:
        colors[i] = 3
    elif station in southsound:
        colors[i] = 4
    else:
        colors[i] = 5

cmap = LinearSegmentedColormap.from_list('categorycolor', color_list, N=5)

##########################################################
##             Average bottom DO vs. depth              ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
depth = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):

    # download .nc files
    fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    # get depth at this point
    depth[i] =  ds['h'].values
    
# create scatter plot
# ax.plot(depth,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(depth,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# print('Mean Depth: {}'.format(np.mean(depth)))

# calculate correlation coefficient (Pearson)
r,p = pearsonr(depth, bott_DO)
r_depth = r
ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=ts, fontweight='bold')
ax.text(0.1, 0.79, r'$p =$' + str(round(p,4)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=ts, fontweight='bold')

# add basin label
ax.text(0.95, 0.95, 'Whidbey Basin',color=color_list[0],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.95, 0.91, 'Hood Canal',color=color_list[1],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.95, 0.87, 'Main Basin',color=color_list[2],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.95, 0.83, 'South Sound',color=color_list[3],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.95, 0.79, 'Admiralty Inlet',color=color_list[4],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')

# format labels
ax.set_title('Average bottom DO vs. Inlet depth\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel('Depth [m]',fontsize=ls)
ax.set_ylabel('Avg. bottom DO [mg/L]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

# add inlet label
if inlet_label:
    for i, txt in enumerate(sta_dict):
        ax.annotate(txt, (depth[i],bott_DO[i]), fontsize=10)

plt.show()

##########################################################
##         Average bottom DO vs. delta rho              ##
##########################################################


pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
strat = np.zeros(len(sta_dict))

print('Delta rho -------------------')
for i,station in enumerate(sta_dict):
    
    # calculate lat/lon for station
    lon = sta_dict[station][0]
    lat = sta_dict[station][1]

    # download .nc files
    fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    # get stratification at this point
    # calculate absolute salinity from practical salinity
    salt_abs_top_all = gsw.conversions.SA_from_SP(ds['salt'].values[:,-1], ds['z_rho'].values[:,-1], lon, lat)
    salt_abs_bott_all = gsw.conversions.SA_from_SP(ds['salt'].values[:,0], ds['z_rho'].values[:,0], lon, lat)
    # calculate conservative temperature from potential temperature
    cons_temp_top_all = gsw.conversions.CT_from_pt(salt_abs_top_all, ds['temp'].values[:,-1])
    cons_temp_bott_all = gsw.conversions.CT_from_pt(salt_abs_bott_all, ds['temp'].values[:,0])
    # calculate potential density
    rho_top_all = gsw.density.sigma0(salt_abs_top_all,cons_temp_top_all)
    rho_bott_all = gsw.density.sigma0(salt_abs_bott_all,cons_temp_bott_all)
    # calculate density difference
    rho_diff_all = rho_bott_all - rho_top_all
    # get average and divide by depth
    rho_diff_avg = np.nanmean(rho_diff_all)
    strat[i] = rho_diff_avg

    if station in ['dabob','lynchcove','sinclair','carr','commencement']:
        print('    {} = {}'.format(station,rho_diff_avg))
    
# create scatter plot
# ax.plot(strat,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(strat,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(strat, bott_DO)
r_drho = r
ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=ts, fontweight='bold')
ax.text(0.1, 0.79, r'$p =$' + str(round(p,4)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=ts, fontweight='bold')

# add basin label
ax.text(0.95, 0.95, 'Whidbey Basin',color=color_list[0],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.95, 0.91, 'Hood Canal',color=color_list[1],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.95, 0.87, 'Main Basin',color=color_list[2],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.95, 0.83, 'South Sound',color=color_list[3],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
ax.text(0.95, 0.79, 'Admiralty Inlet',color=color_list[4],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts-3, fontweight='bold')

# format labels
ax.set_title(r'Average bottom DO vs. $\Delta\rho$' + '\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'$\rho_{bottom} - \rho_{surface}$ [kg m$^{-3}$]',fontsize=ls)
ax.set_ylabel('Avg. bottom DO [mg/L]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

# add inlet label
if inlet_label:
    for i, txt in enumerate(sta_dict):
        ax.annotate(txt, (strat[i],bott_DO[i]), fontsize=10)

plt.show()
