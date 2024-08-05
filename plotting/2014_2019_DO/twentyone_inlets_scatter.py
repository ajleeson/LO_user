"""
Generate depth vs. time property plots using mooring extraction data. 
Used to look at 21 inlets in Puget Sound, and compare to Ecology monitoring stations, if available
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

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'extract' / 'tef2'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import tef_fun

Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
jobname = 'twentyoneinlets'
startdate = '2014.01.01'
enddate = '2014.12.31'
year = '2014' # for making a date label

# figure settings
fs = 12 # figure font size
ls = 11 # label size
ts = 12 # title size

# hypoxic season
start = '08-01'
end = '09-30'

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

# # remove totten and hammersley, to see if correlation still exists
# del sta_dict['totten']
# del sta_dict['hammersley']

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures' / 'twentyone' / 'scatter_R2'
Lfun.make_dir(out_dir)

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
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))
    # get average bottom DO
    bott_DO_all = ds['oxygen'].values[:,0] * pinfo.fac_dict['oxygen'] # mg/L
    # average over hypoxic season
    bott_DO_avg = np.nanmean(bott_DO_all)
    bott_DO[i] = bott_DO_avg
    print(station)
    print(bott_DO_avg)
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
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    # get depth at this point
    depth[i] =  ds['h'].values
    
# create scatter plot
# ax.plot(depth,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(depth,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

print('Mean Depth: {}'.format(np.mean(depth)))

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

plt.savefig(out_dir / 'botDO_vs_depth.png')

##########################################################
##         Average bottom DO vs. delta rho              ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
strat = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):
    
    # calculate lat/lon for station
    lon = sta_dict[station][0]
    lat = sta_dict[station][1]

    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
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

plt.savefig(out_dir / 'botDO_vs_strat.png')

##########################################################
##       Average bottom DO vs. tidal currents           ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
tcurr = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):

    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '_hourly/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)

    # get average tidal current
    ubar = ds['ubar'].values
    vbar = ds['vbar'].values
    # calculate U_rms
    #U_rms = sqrt( Godin( ubar^2 + vbar^2 ) ) # this equation got Parker's OK
    U_rms = np.sqrt( zfun.lowpass(ubar**2 + vbar**2, f='godin')) 
    tcurr[i] = np.nanmean(U_rms) # average over time
    
# create scatter plot
# ax.plot(tcurr,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(tcurr,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(tcurr, bott_DO)
r_Urms = r
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
ax.set_title('Average bottom DO vs. Tidal currents\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel('Annual average tidal currents [m/s]',fontsize=ls)
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
        ax.annotate(txt, (tcurr[i],bott_DO[i]), fontsize=10)

plt.savefig(out_dir / 'botDO_vs_tidalcurr.png')


##########################################################
##      Average bottom DO vs. bottom large detritus     ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
botLdetN = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):

    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    d = 5 # bottom 5 m
    botLdetN[i] = np.nanmean(get_bott_average('LdetritusN',d,ds.h,ds))
    
# create scatter plot
# ax.plot(botLdetN,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(botLdetN,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(botLdetN, bott_DO)
rLdetN = r
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
ax.set_title('Average bottom DO vs. bottom 5 m large detritus\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'Avg. bottom 5 m large detritus [mmol/m$^3$]',fontsize=ls)
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
        ax.annotate(txt, (botLdetN[i],bott_DO[i]), fontsize=10)

plt.savefig(out_dir / 'botDO_vs_botLdet.png')

##########################################################
##      Average bottom DO vs. surface small detritus    ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
surfSdetN = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):

    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    d = 80 # surface 80 m
    surfSdetN[i] = np.nanmean(get_surf_average('SdetritusN',d,ds))
    
# create scatter plot
# ax.plot(surfSdetN,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(surfSdetN,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(surfSdetN, bott_DO)
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
ax.set_title('Average bottom DO vs. surface 80 m small detritus\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'Avg. surface 80 m small detritus [mmol/m$^3$]',fontsize=ls)
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
        ax.annotate(txt, (surfSdetN[i],bott_DO[i]), fontsize=10)

plt.savefig(out_dir / 'botDO_vs_surfSdet.png')

##########################################################
##      Average bottom DO vs. surface phytoplankton     ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
surfphyto = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):

    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    d = 20 # surface 20 m
    surfphyto[i] = np.nanmean(get_surf_average('phytoplankton',d,ds))
    
# create scatter plot
# ax.plot(surfphyto,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(surfphyto,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(surfphyto, bott_DO)
rphyto = r
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
ax.set_title('Average bottom DO vs. surface 20 m phytoplankton\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'Avg. surface 20 m phytoplankton [mmol/m$^3$]',fontsize=ls)
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
        ax.annotate(txt, (surfphyto[i],bott_DO[i]), fontsize=10)

plt.savefig(out_dir / 'botDO_vs_surfphyto.png')

##########################################################
##      Average bottom DO vs. surface zooplankton       ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
surfzoop = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):

    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    d = 20 # surface 20 m
    surfzoop[i] = np.nanmean(get_surf_average('zooplankton',d,ds))
    
# create scatter plot
# ax.plot(surfzoop,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(surfzoop,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(surfzoop, bott_DO)
rzoop = r
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
ax.set_title('Average bottom DO vs. surface 20 m zooplankton\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'Avg. surface 20 m zooplankton [mmol/m$^3$]',fontsize=ls)
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
        ax.annotate(txt, (surfzoop[i],bott_DO[i]), fontsize=10)

plt.savefig(out_dir / 'botDO_vs_surfzoop.png')

##########################################################
##          Average bottom DO vs. surface NO3           ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
surfNO3 = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):

    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    d = 20 # surface 20 m
    surfNO3[i] = np.nanmean(get_surf_average('NO3',d,ds))
    
# create scatter plot
# ax.plot(surfNO3,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(surfNO3,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(surfNO3, bott_DO)
rNO3 = r
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
ax.set_title('Average bottom DO vs. surface 20 m NO3\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'Avg. surface 20 m NO3 [mmol/m$^3$]',fontsize=ls)
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
        ax.annotate(txt, (surfNO3[i],bott_DO[i]), fontsize=10)

plt.savefig(out_dir / 'botDO_vs_surfNO3.png')

##########################################################
##          Average bottom DO vs. surface NH4           ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
surfNH4 = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):

    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    d = 20 # surface 20 m
    surfNH4[i] = np.nanmean(get_surf_average('NH4',d,ds))
    
# create scatter plot
# ax.plot(surfNH4,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(surfNH4,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(surfNH4, bott_DO)
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
ax.set_title('Average bottom DO vs. surface 20 m NH4\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'Avg. surface 20 m NH4 [mmol/m$^3$]',fontsize=ls)
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
        ax.annotate(txt, (surfNH4[i],bott_DO[i]), fontsize=10)

plt.savefig(out_dir / 'botDO_vs_surfNH4.png')

##########################################################
##           Average bottom DO vs. V/Qprism            ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# open pickle file to get volume
seg_pickle = '../../../LO_output/extract/tef2/vol_df_cas7_c21.p'
with open(seg_pickle, 'rb') as f:
    segments = pickle.load(f)

# initialize arrays for plotting
Tflush = np.zeros(len(sta_dict))

# get inlet volume
for i,station in enumerate(sta_dict):
    # get inlet info # m3
    volume = segments.loc[(station + '_p')]['volume m3']
    # get tidal prism
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(Ldir['LOo'] / 'extract' /'cas7_t0_x4b' / 'tef2' / 'bulk_2014.01.01_2014.12.31', station)
    Qprism_alltime = tef_df['qprism'] # m3/s, with one value per day
    # annual average Qprism
    avg_Qprism = np.nanmean(Qprism_alltime)
    # calculate Volume/Qprism
    Tflush[i] = volume/avg_Qprism * (1/60) * (1/60) * (1/24) # seconds, then convert to days (1/60 min/sec) * (1/60 hr/min) * (1/24 day/hr)

# create scatter plot
# ax.plot(Tflush,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(Tflush,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(Tflush, bott_DO)
r_Tflush = r
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
ax.set_title('Average bottom DO vs. flushing time\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'Annual avg $T_{flush} = V/Q_{prism}$ [days]',fontsize=ls)
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
        ax.annotate(txt, (Tflush[i],bott_DO[i]), fontsize=10)

plt.savefig(out_dir / 'botDO_vs_Tflush.png')

##########################################################
##      Average bottom DO vs. yearday of peak bloom     ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# open pickle file to get volume
seg_pickle = '../../../LO_output/extract/tef2/vol_df_cas7_c21.p'
with open(seg_pickle, 'rb') as f:
    segments = pickle.load(f)

# initialize arrays for plotting
pbloom_day = np.zeros(len(sta_dict))

# get inlet volume
for i,station in enumerate(sta_dict):
    # get data
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    bott_ox_raw = ds['oxygen'].values[:,0] # 0 = bottom layer
    val = bott_ox_raw *  pinfo.fac_dict['oxygen'] # convert to mg/L
    # get index of max value-- this is the yearday minus one
    max_i = np.argmax(val)
    # yearday peak bloom
    pbloom_day[i] = max_i + 1


# create scatter plot
# ax.plot(Tflush,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
ax.scatter(pbloom_day,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(pbloom_day, bott_DO)
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
ax.set_title('Average bottom DO vs. yearday of peak bloom\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel('Yearday of peak bloom',fontsize=ls)
ax.set_ylabel('Avg. bottom DO [mg/L]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)


plt.savefig(out_dir / 'botDO_vs_pbloomday.png')

##########################################################
##            Depth vs. Tidal currents                  ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5.5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

cmap_oxy = plt.cm.get_cmap('rainbow_r', 10)
cs = ax.scatter(tcurr,depth,c=bott_DO,edgecolor='k',s=75,cmap=cmap_oxy,zorder=5)
cbar = fig.colorbar(cs)
cbar.ax.tick_params(labelsize=ls)
cbar.ax.set_ylabel('Aug/Sep Bottom DO [mg/L]', rotation=90, fontsize=ls)
cbar.outline.set_visible(False)

# format labels
ax.set_title('Depth vs. Tidal Currents',
             fontsize=ts)
ax.set_xlabel('Annual avg. ' + r'$U_{rms}$ [m s$^{-1}$]', fontsize=ls)
ax.set_ylabel('Depth [m]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig(out_dir / 'depth_v_Urms.png')


##########################################################
##          Stratification vs. Tidal currents           ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5.5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

cmap_oxy = plt.cm.get_cmap('rainbow_r', 10)
cs = ax.scatter(tcurr,strat,c=bott_DO,edgecolor='k',s=75,cmap=cmap_oxy,zorder=5)
cbar = fig.colorbar(cs)
cbar.ax.tick_params(labelsize=ls)
cbar.ax.set_ylabel('Aug/Sep Bottom DO [mg/L]', rotation=90, fontsize=ls)
cbar.outline.set_visible(False)

# format labels
ax.set_title(r'$\Delta\rho$ vs. Tidal Currents',
             fontsize=ts)
ax.set_xlabel('Annual avg. ' + r'$U_{rms}$ [m s$^{-1}$]', fontsize=ls)
ax.set_ylabel(r'Aug/Sep $\rho_{bottom} - \rho_{surface}$ [kg m$^{-3}$]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig(out_dir / 'deltarho_v_Urms.png')

##########################################################
##               Depth vs. delta rho                    ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5.5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

cs = ax.scatter(strat,depth,c=bott_DO,edgecolor='k',s=75,cmap=cmap_oxy,zorder=5)
cbar = fig.colorbar(cs)
cbar.ax.tick_params(labelsize=ls)
cbar.ax.set_ylabel('Aug/Sep Bottom DO [mg/L]', rotation=90, fontsize=ls)
cbar.outline.set_visible(False)

# format labels
ax.set_title(r'Depth vs. $\Delta\rho$',
             fontsize=ts)
ax.set_xlabel(r'Aug/Sep $\rho_{bottom} - \rho_{surface}$ [kg m$^{-3}$]', fontsize=ls)
ax.set_ylabel('Depth [m]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig(out_dir / 'depth_v_deltarho.png')

##########################################################
##            Tflushing   vs. Tidal currents            ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5.5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

cmap_oxy = plt.cm.get_cmap('rainbow_r', 10)
cs = ax.scatter(tcurr,Tflush,c=bott_DO,edgecolor='k',s=75,cmap=cmap_oxy,zorder=5)
cbar = fig.colorbar(cs)
cbar.ax.tick_params(labelsize=ls)
cbar.ax.set_ylabel('Aug/Sep Bottom DO [mg/L]', rotation=90, fontsize=ls)
cbar.outline.set_visible(False)

# format labels
ax.set_title(r'$T_{flush}$ vs. Tidal Currents',
             fontsize=ts)
ax.set_xlabel('Annual avg. ' + r'$U_{rms}$ [m s$^{-1}$]', fontsize=ls)
ax.set_ylabel(r'Annual avg $T_{flush} = V/Q_{prism}$ [days]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)


plt.savefig(out_dir / 'Tflush_v_Urms.png')

##########################################################
##          stratification vs. tidal flushing           ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5.5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

cmap_oxy = plt.cm.get_cmap('rainbow_r')
cs = ax.scatter(Tflush,strat,c=bott_DO,edgecolor='k',s=75,cmap=cmap_oxy,zorder=5)
cbar = fig.colorbar(cs)
cbar.ax.tick_params(labelsize=ls)
cbar.ax.set_ylabel('Aug/Sep Bottom DO [mg/L]', rotation=90, fontsize=ls)
cbar.outline.set_visible(False)


# calculate correlation coefficient (Pearson)
r,p = pearsonr(Tflush,strat)
rLdetN = r
ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=ts, fontweight='bold')
ax.text(0.1, 0.79, r'$p =$' + str(round(p,4)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=ax.transAxes, fontsize=ts, fontweight='bold')

# format labels
ax.set_title(r'$\Delta\rho$ vs. $T_{flush}$',
             fontsize=ts)
ax.set_xlabel(r'Annual avg $T_{flush} = V/Q_{prism}$ [days]',fontsize=ls)
ax.set_ylabel(r'Aug/Sep $\rho_{bottom} - \rho_{surface}$ [kg m$^{-3}$]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)


plt.savefig(out_dir / 'deltarho_v_Tflush.png')

##########################################################
##                 Scatter matrix                       ##
##########################################################

plt.close('all')
fig, axes = plt.subplots(1,1,figsize = (8,8))

# put all variables in pandas dataframe
d = {'Depth [m]': depth,
    r'$\Delta\rho$ [kg m$^{-3}$]': strat,
    r'$U_{rms}$ [m s$^{-1}$]': tcurr,
    r'$T_{flush}$ [days]': Tflush}
df = pd.DataFrame(data=d)

pd.plotting.scatter_matrix(df[['Depth [m]', r'$\Delta\rho$ [kg m$^{-3}$]',
                            r'$U_{rms}$ [m s$^{-1}$]', r'$T_{flush}$ [days]']],
                           figsize=(8,8), range_padding=0.2,
                           ax = axes,
                           hist_kwds={'ec':'gray',
                                      'facecolor':'silver',
                                      'alpha':0.5,
                                     'bins':8},
                           c=bott_DO, marker='.',s=180, #edgecolors='k'
                           alpha=1,cmap=cmap_oxy)

plt.suptitle('Physical Characteristic Scatter Matrix')

# add colorbar
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.88, 0.12, 0.03, 0.74])
cbar = fig.colorbar(cs, cax=cbar_ax)
cbar.outline.set_visible(False)
cbar.ax.set_ylabel('Aug/Sep Bottom DO [mg/L]', rotation=90, fontsize=ls)

plt.savefig(out_dir / 'scatter_matrix.png')

##########################################################
##               Plot all in subplots                   ##
##########################################################

# histrogram of bottom DO
fig, ax = plt.subplots(1,1,figsize = (6,5))

# create histogram
bins = [0,1,2,3,4,5,6,7,8,9,10]
ax.hist(bott_DO, bins=bins, zorder=5, facecolor='mediumorchid',
        alpha=0.6, edgecolor='darkviolet')
ax.set_xlim([0,10])
ax.set_xlabel(r'Avg. Aug/Sep $DO_{bot}$ [mg/L]')
ax.set_ylim([0,5])
ax.set_ylabel('Count')
ax.set_title(r'Distribution of $DO_{bot}$ in 21 inlets' + '\n(averaged over Aug/Sep)')

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)
plt.savefig(out_dir / 'bottDO_hist.png')

# --------------------------------------------------------

plt.close('all')

fig, axes = plt.subplots(1,4,figsize = (10,3.4),sharey=True)
axes[0].set_ylabel('Avg. bottom DO [mg/L]',fontsize=ls)
# add basin label
axes[3].text(0.95, 0.9, 'Whidbey Basin',color=color_list[0],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
axes[3].text(0.95, 0.81, 'Hood Canal',color=color_list[1],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
axes[3].text(0.95, 0.72, 'Main Basin',color=color_list[2],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
axes[3].text(0.95, 0.63, 'South Sound',color=color_list[3],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
axes[3].text(0.95, 0.54, 'Admiralty Inlet',color=color_list[4],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')

# add data
# Depth
axes[0].scatter(depth,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
axes[0].set_title(r'(a) vs. depth', loc='left', fontsize=ts)
axes[0].set_xlabel('Depth [m]', fontsize=ls)
axes[0].text(0.1, 0.84, r'$r =$' + str(round(r_depth,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=axes[0].transAxes, fontsize=ts)
# delta rho
axes[1].scatter(strat,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
axes[1].set_title(r'(b) vs. Aug/Sep $\Delta\rho$', loc='left', fontsize=ts)
axes[1].set_xlabel(r'$\rho_{bottom} - \rho_{surface}$ [kg m$^{3}$]', fontsize=ls)
axes[1].text(0.1, 0.84, r'$r =$' + str(round(r_drho,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=axes[1].transAxes, fontsize=ts)
# RMS Tidal currents
axes[2].scatter(tcurr,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
axes[2].set_title(r'(c) vs. annual avg. $U_{rms}$', loc='left', fontsize=ts)
axes[2].set_xlabel(r'$U_{rms}$ [m s$^{-1}$]', fontsize=ls)
axes[2].text(0.1, 0.84, r'$r =$' + str(round(r_Urms,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=axes[2].transAxes, fontsize=ts)
# flushing time
axes[3].scatter(Tflush,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
axes[3].set_title(r'(d) vs. annual avg. $T_{flush}$', loc='left', fontsize=ts)
axes[3].set_xlabel(r'$T_{flush} = V/Q_{prism}$ [days]', fontsize=ls)
axes[3].text(0.1, 0.84, r'$r =$' + str(round(r_Tflush,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=axes[3].transAxes, fontsize=ts)

for ax in axes:
    # format grid
    ax.tick_params(axis='both', which='major', labelsize=ls)
    ax.grid(True,color='white',linewidth=1)
    # format colors
    ax.set_facecolor('#EEEEEE')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)

plt.suptitle(r'Aug/Sep DO$_{Bot}$ vs. physical characteristics',fontsize=ts+1)

plt.subplots_adjust(wspace=0.05,top=0.82, bottom=0.16,right=0.99,left=0.06)
plt.savefig(out_dir / 'physics_summary.png')

# --------------------------------------------------------

plt.close('all')

fig, axes = plt.subplots(1,4,figsize = (10,3.4),sharey=True)
axes[0].set_ylabel('Avg. bottom DO [mg/L]',fontsize=ls)
# add basin label
axes[3].text(0.95, 0.9, 'Whidbey Basin',color=color_list[0],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
axes[3].text(0.95, 0.81, 'Hood Canal',color=color_list[1],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
axes[3].text(0.95, 0.72, 'Main Basin',color=color_list[2],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
axes[3].text(0.95, 0.63, 'South Sound',color=color_list[3],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
axes[3].text(0.95, 0.54, 'Admiralty Inlet',color=color_list[4],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[3].transAxes, fontsize=ls, fontweight='bold')

# add data
# NO3
axes[0].scatter(surfNO3,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
axes[0].set_title('(a) vs. NO3 (surf 20 m)', loc='left', fontsize=ts)
axes[0].set_xlabel(r'NO3 [mmol m$^{-3}$]', fontsize=ls)
axes[0].text(0.1, 0.84, r'$r =$' + str(round(rNO3,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=axes[0].transAxes, fontsize=ts)
# phytoplankton
axes[1].scatter(surfphyto,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
axes[1].set_title('(b) vs. phyto (surf 20 m)', loc='left', fontsize=ts)
axes[1].set_xlabel(r'phytoplankton [mmol m$^{-3}$]', fontsize=ls)
axes[1].text(0.1, 0.84, r'$r =$' + str(round(rphyto,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=axes[1].transAxes, fontsize=ts)
# zooplankton
axes[2].scatter(surfzoop,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
axes[2].set_title('(c) vs. zoo (surf 20 m)', loc='left', fontsize=ts)
axes[2].set_xlabel(r'zooplankton [mmol m$^{-3}$]', fontsize=ls)
axes[2].text(0.1, 0.84, r'$r =$' + str(round(rzoop,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=axes[2].transAxes, fontsize=ts)
# LdetN
axes[3].scatter(botLdetN,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
axes[3].set_title('(d) vs. L det (bott 5 m)', loc='left', fontsize=ts)
axes[3].set_xlabel(r'large detritus [mmol m$^{-3}$]', fontsize=ls)
axes[3].text(0.1, 0.84, r'$r =$' + str(round(rLdetN,2)) ,color='black',
                        verticalalignment='bottom', horizontalalignment='left',
                        transform=axes[3].transAxes, fontsize=ts)

for ax in axes:
    # format grid
    ax.tick_params(axis='both', which='major', labelsize=ls)
    ax.grid(True,color='white',linewidth=1)
    # format colors
    ax.set_facecolor('#EEEEEE')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)

plt.suptitle(r'Aug/Sep DO$_{Bot}$ vs. Aug/Sep biogeochemistry',fontsize=ts+1)

plt.subplots_adjust(wspace=0.05,top=0.82, bottom=0.16,right=0.99,left=0.06)
plt.savefig(out_dir / 'biogeochem_summary.png')

plt.close('all')