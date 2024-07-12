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

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

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
ts = 14 # title size

# hypoxic season
start = '08-01'
end = '09-30'

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
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures' / 'twentyone' / 'scatter_R2'
Lfun.make_dir(out_dir)

# Calculate average bottom DO
bott_DO = np.zeros(len(sta_dict))
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
ax.plot(depth,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(depth, bott_DO)
ax.text(0.85, 0.85, r'$r =$' + str(round(r,2)) ,color='navy',
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts)
ax.text(0.85, 0.78, r'$p =$' + str(round(p,3)) ,color='navy',
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts)

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
    # get pressure
    press_top = [gsw.p_from_z(z,lat) for z in ds['z_rho'].values[:,-1]]
    press_bott = [gsw.p_from_z(z,lat) for z in ds['z_rho'].values[:,0]]
    # calculate absolute salinity from practical salinity
    salt_abs_top_all = gsw.conversions.SA_from_SP(ds['salt'].values[:,-1], ds['z_rho'].values[:,-1], lon, lat)
    salt_abs_bott_all = gsw.conversions.SA_from_SP(ds['salt'].values[:,0], ds['z_rho'].values[:,0], lon, lat)
    # calculate conservative temperature from potential temperature
    cons_temp_top_all = gsw.conversions.CT_from_pt(salt_abs_top_all, ds['temp'].values[:,-1])
    cons_temp_bott_all = gsw.conversions.CT_from_pt(salt_abs_bott_all, ds['temp'].values[:,0])
    # calculate density
    rho_top_all = gsw.rho(salt_abs_top_all,cons_temp_top_all,press_top)
    rho_bott_all = gsw.rho(salt_abs_bott_all,cons_temp_bott_all,press_bott)
    # calculate density difference
    rho_diff_all = rho_bott_all - rho_top_all
    # get average and divide by depth
    rho_diff_avg = np.nanmean(rho_diff_all)
    strat[i] = rho_diff_avg
    
# create scatter plot
ax.plot(strat,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(strat, bott_DO)
ax.text(0.85, 0.85, r'$r =$' + str(round(r,2)) ,color='navy',
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts)
ax.text(0.85, 0.78, r'$p =$' + str(round(p,4)) ,color='navy',
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts)

# format labels
ax.set_title('Average bottom DO vs. Stratification\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'Stratification $\Delta\rho$ [kg m$^{-3}$]',fontsize=ls)
ax.set_ylabel('Avg. bottom DO [mg/L]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig(out_dir / 'botDO_vs_strat.png')

##########################################################
##         Average bottom DO vs. delta rho/h            ##
##########################################################

plt.close('all')

pfun.start_plot(figsize=(5,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# initialize arrays for plotting
strat_over_h = np.zeros(len(sta_dict))

for i,station in enumerate(sta_dict):

    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    # get depth at this point
    h =  ds['h'].values

    strat_over_h[i] = strat[i]/h
    
# create scatter plot
ax.plot(strat_over_h,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(strat_over_h, bott_DO)
ax.text(0.85, 0.85, r'$r =$' + str(round(r,2)) ,color='navy',
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts)
ax.text(0.85, 0.78, r'$p =$' + str(round(p,2)) ,color='navy',
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts)

# format labels
ax.set_title('Average bottom DO vs. Stratification/Depth\n('+year+'-'+start+' to '+year+'-'+end+')',
             fontsize=ts)
ax.set_xlabel(r'$\Delta\rho/H$ [kg m$^{-4}$]',fontsize=ls)
ax.set_ylabel('Avg. bottom DO [mg/L]',fontsize=ls)

# format grid
ax.tick_params(axis='both', which='major', labelsize=ls)
ax.grid(True,color='white',linewidth=1)
# format colors
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

plt.savefig(out_dir / 'botDO_vs_strat_over_h.png')


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
    U = np.sqrt(ubar**2 + vbar**2) # add in quadrature (RMS)
    tcurr[i] = np.nanmean(U)
    
# create scatter plot
ax.plot(tcurr,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)

# calculate correlation coefficient (Pearson)
r,p = pearsonr(tcurr, bott_DO)
ax.text(0.85, 0.85, r'$r =$' + str(round(r,2)) ,color='navy',
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts)
ax.text(0.85, 0.78, r'$p =$' + str(round(p,4)) ,color='navy',
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=ax.transAxes, fontsize=ts)

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

plt.savefig(out_dir / 'botDO_vs_tidalcurr.png')
