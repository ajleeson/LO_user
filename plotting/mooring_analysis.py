#%%
"""
Plots a timeseries of different variables at mooring stations to compare different model runs.

Inputs: list of gtagex, label names, moor job name (As defined in job_lists), start date, end date

The labels names should be the distinguishing factor between the gtagexes.

I need to define these inputs at line 38 in this document. 

Then, from the terminal, I can simply run: python mooring_timeseries.py

"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Ldir = Lfun.Lstart()

#%% DEFINE INPUTS ----------------------------------------------------------------------------

gtagexes = ['cas7_t0noN_x4b', 'cas7_t0_x4b']
labels = ['N-less', 'Hindcast']
jobname = 'noWWTPNtest'
startdate = '2013.01.01'
enddate = '2013.12.31'
dayafterend = '2014.01.01' # for making date arrays. The next day after we have no more data
year = '2013' # for making a date label

#%%  PLOT THE MOORING LOCATIONS -----------------------------------------

plt.close('all')

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots' / 'mooring_timeseries'
Lfun.make_dir(out_dir)

# parse gtagex
gtagexample = gtagexes[0]
gridname, tag, ex_name = gtagexample.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)

# load grid data
grid_nc = '../../LO_data/grids/' + gridname + '/grid.nc'
ds = xr.open_dataset(grid_nc)
z = -ds.h.values
mask_rho = ds.mask_rho.values

lon = ds.lon_rho.values
lat = ds.lat_rho.values

plon, plat = pfun.get_plon_plat(lon,lat)
# pad = 0.05*(plat[-1,0]-plat[0,0])
# ax_lims = (plon[0,0]-pad, plon[0,-1]+pad, plat[0,0]-pad, plat[-1,0]+pad)

# make a version of z with nans where masked
zm = z.copy()
zm[mask_rho == 0] = np.nan

# # Initialize plots
# plt.close('all')
# pfun.start_plot(figsize=(4,4))

# lon/lat limits (Puget Sound)
xmin = -123.2
xmax = -122.1
ymin = 46.93
ymax = 48.45
# plot bathymetry
fig = plt.figure(figsize=(8,9))
ax = fig.add_subplot(111)
cs = ax.pcolormesh(plon, plat, zm, vmin=-250, vmax=0, cmap=plt.get_cmap(cmocean.cm.deep_r))
cbar = fig.colorbar(cs, ax=ax)
cbar.ax.set_title('Depth (m)', fontsize = 12)
cbar.outline.set_visible(False)
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
ax.set_title('Mooring Locations')
ax.set_ylabel('Latitude', fontsize = 16)
ax.tick_params(axis='x', labelrotation=30)
ax.set_xlabel('Longitude', fontsize = 16)
pfun.dar(ax)

# plot mooring locations
moor_lats = [x[1] for x in sta_dict.values()]
moor_lons= [x[0] for x in sta_dict.values()]
plt.scatter(moor_lons,moor_lats, s=150, c='deeppink', marker='*')#, edgecolors='k')

# label mooring locations
for i,station in enumerate(sta_dict.keys()):
    ax.text(moor_lons[i],moor_lats[i]+0.05,station, fontsize=12,
             c='deeppink', horizontalalignment='center', fontweight = 'bold')


#%% Chlorophyll and DO --------------------------------------------------------------------

# create time vector
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

# Loop through all of the mooring stations
for i,station in enumerate(sta_dict.keys()):

    # initalize a new figure
    fig = plt.figure(figsize=(10,8))#,facecolor='none')
    plt.tight_layout()
    fig.suptitle(station)

    # pick colors
    phyto_color = 'darkcyan'
    oxy_color = 'orchid'

    # create all of the axis
    axN_0 = plt.subplot(311)
    axDO_0 = axN_0.twinx()
    axN_1 = plt.subplot(312)
    axDO_1 = axN_1.twinx()
    axN_diff = plt.subplot(313)
    axDO_diff = axN_diff.twinx()

    axN = [axN_0,axN_1,axN_diff]
    axDO = [axDO_0,axDO_1,axDO_diff]

    letters = ['(a) ', '(b) ']

    # loop through all of the model scenarios
    for j,gtagex in enumerate(gtagexes):

        # download .nc files
        fn = '../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds = xr.open_dataset(fn)

        # save datasets for later, and add titles
        if gtagex == 'cas7_t0noN_x4b':
            ds_noN = ds
            title = 'No-loading'
        elif gtagex == 'cas7_t0_x4b':
            ds_hindcast = ds
            title = 'With-loading'

        # add subtitles for each panel
        axN[j].text(0.02, 0.8, letters[j]+title, fontweight='bold',
            verticalalignment='bottom', horizontalalignment='left',
            transform=axN[j].transAxes, fontsize=12)
        axN[2].text(0.02, 0.8, '(c) With-loading minus No-loading', fontweight='bold',
            verticalalignment='bottom', horizontalalignment='left',
            transform=axN[2].transAxes, fontsize=12)

        # depth layer
        surf_lev = -1
        bott_lev = 0

        # plot surface phytoplankton concentration
        scale =  pinfo.fac_dict['phytoplankton']
        units = pinfo.units_dict['phytoplankton']
        phyto_val = ds.phytoplankton.values[:,surf_lev] * scale
        axN[j].plot(dates_local, phyto_val,color=phyto_color,linewidth=2)
        axN[j].tick_params(axis='y', labelcolor=phyto_color, labelsize=12, color=phyto_color)
        axN[j].set_ylim([0,35])
        axN[j].set_xlim([dates_local[0],dates_local[-1]])
        axN[j].set_xticklabels([])
        # add legend
        axN[0].text(0.7, 0.8, 'Surface Chl '+ units, color=phyto_color,
            verticalalignment='bottom', horizontalalignment='left',
            transform=axN[0].transAxes, fontsize=12)
        # plot difference
        if j == 1:
            # draw horizontal line at 0
            axN[j+1].axhline(y = 0, color = 'w', linewidth=2)
            # get values from both runs
            phyto_val_noN = ds_noN.phytoplankton.values[:,surf_lev] * scale
            phyto_val_hindcast = ds_hindcast.phytoplankton.values[:,surf_lev] * scale
            # calculate hindcast minus no N
            diff = phyto_val_hindcast - phyto_val_noN
            # plot difference
            axN[j+1].plot(dates_local, diff, color = phyto_color,linewidth=2)
            axN[j+1].tick_params(axis='y', labelcolor=phyto_color, labelsize=12, color=phyto_color)
            # pick limits
            lim = 1.2*(max(abs(min(diff)),max(diff)))
            axN[j+1].set_ylim([-1*lim,lim])
            axN[j+1].set_xlim([dates_local[0],dates_local[-1]])


        # plot bottom DO concentration
        scale =  pinfo.fac_dict['oxygen']
        units = pinfo.units_dict['oxygen']
        oxy_val = ds.oxygen.values[:,bott_lev] * scale
        axDO[j].plot(dates_local, oxy_val, color = oxy_color,linewidth=2)
        axDO[j].tick_params(axis='y', labelcolor=oxy_color, labelsize=12, color = oxy_color)
        axDO[j].set_ylim([0,12])
        axDO[j].set_xlim([dates_local[0],dates_local[-1]])
        axDO[j].set_xticklabels([])
        # add legend
        axDO[0].text(0.7, 0.68, 'Bottom DO '+ units, color=oxy_color,
            verticalalignment='bottom', horizontalalignment='left',
            transform=axN[0].transAxes, fontsize=12)
        # plot difference
        if j == 1:
            # get values from both runs
            oxy_val_noN = ds_noN.oxygen.values[:,bott_lev] * scale
            oxy_val_hindcast = ds_hindcast.oxygen.values[:,bott_lev] * scale
            # calculate hindcast minus no N
            diff = oxy_val_hindcast - oxy_val_noN
            # plot difference
            axDO[j+1].plot(dates_local, diff, color = oxy_color,linewidth=2)
            axDO[j+1].tick_params(axis='y', labelcolor=oxy_color, labelsize=12, color = oxy_color)
            # pick limits
            lim = 1.2*(max(abs(min(diff)),max(diff)))
            axDO[j+1].set_ylim([-1*lim,lim])
            axDO[j+1].set_xlim([dates_local[0],dates_local[-1]])

        # format date label
        axDO[2].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        axDO[2].set_xlabel('2013')#, color='#EEEEEE')
        axDO[2].tick_params(axis='x', labelrotation=30, labelsize=12) # colors='#EEEEEE'
        axN[2].tick_params(axis='x', labelrotation=30, labelsize=12) # colors='#EEEEEE'

        # format subplot colors
        for ax in axN + axDO:
            ax.grid(True,color='w',linewidth=2,axis='x')
            ax.set_facecolor('#EEEEEE')
            for border in ['top','right','bottom','left']:
                ax.spines[border].set_visible(False)

    plt.savefig(out_dir / (station+'_chl_do.png'))#, facecolor=ax.get_facecolor(), transparent=True)


#%% Total water column N --------------------------------------------------------------------
                
# Loop through all of the mooring stations
for i,station in enumerate(sta_dict.keys()):

    # get multiplier
    scale =  pinfo.fac_dict['NO3'] # mmol/m3

    # initalize a new figure
    fig = plt.figure(figsize=(10,6),facecolor='none')
    plt.tight_layout()
    fig.suptitle(station)

    # loop through all of the model scenarios
    for gtagex in gtagexes:
        # download .nc files
        fn = '../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds = xr.open_dataset(fn)
        # save datasets for later, and add titles
        if gtagex == 'cas7_t0noN_x4b':
            ds_noN = ds
        elif gtagex == 'cas7_t0_x4b':
            ds_hindcast = ds

    # get NO3 and NH4 arrays (in mmol/m3)
    NO3_noN = ds_noN['NO3'].values * scale # dimensions of time, s_rho
    NH4_noN = ds_noN['NH4'].values * scale # dimensions of time, s_rho
    DIN_noN = NO3_noN + NH4_noN # dimensions of time, s_rho
    NO3_hindcast = ds_hindcast['NO3'].values * scale # dimensions of time, s_rho
    NH4_hindcast = ds_hindcast['NH4'].values * scale # dimensions of time, s_rho
    DIN_hindcast = NO3_hindcast + NH4_hindcast # dimensions of time, s_rho

    # get depth of sigma layers
    zw_noN = ds_noN['z_w'] # dimensions of time, s_w
    zw_hindcast = ds_hindcast['z_w'] # dimensions of time, s_w

    # get total depth at mooring location
    depth_noN = zw_noN[:,-1] - zw_noN[:,0] # dimensions of time
    depth_hindcast = zw_hindcast[:,-1] - zw_hindcast[:,0] # dimensions of time

    # calculate thickness of sigma layers
    thick_noN = np.diff(zw_noN.values) # dimensions of time, s_rho
    thick_hindcast = np.diff(zw_hindcast.values) # dimensions of time, s_rho

    # multiply cell N concentration by cell thickness and area (500x500m) (mmol)
    molNO3_noN = NO3_noN * thick_noN * 500 * 500
    molNO3_hindcast = NO3_hindcast * thick_hindcast * 500 * 500
    molNH4_noN = NH4_noN * thick_noN * 500 * 500
    molNH4_hindcast = NH4_hindcast * thick_hindcast * 500 * 500
    molDIN_noN = DIN_noN * thick_noN * 500 * 500
    molDIN_hindcast = DIN_hindcast * thick_hindcast * 500 * 500
    
    # sum throughout water column to get total N in water column (kg where 1 mmol = 7.1378E-8 kg)
    TmolNO3_noN = np.sum(molNO3_noN,axis=1) * 7.1378E-8
    TmolNO3_hindcast = np.sum(molNO3_hindcast,axis=1) * 7.1378E-8
    TmolNH4_noN = np.sum(molNH4_noN,axis=1) * 7.1378E-8
    TmolNH4_hindcast = np.sum(molNH4_hindcast,axis=1) * 7.1378E-8
    TmolDIN_noN = np.sum(molDIN_noN,axis=1) * 7.1378E-8
    TmolDIN_hindcast = np.sum(molDIN_hindcast,axis=1) * 7.1378E-8
    
    # plot total water column N of the two runs
    ax = plt.subplot(211)
    # NO3
    # ax.plot(dates_local, TmolNO3_noN, color='deeppink',linewidth=0.6,label='NO3 N-less')
    # ax.plot(dates_local, TmolNO3_hindcast, color='deeppink',linewidth=2,alpha=0.6,label='NO3 with-N')
    # # NH4
    # ax.plot(dates_local, TmolNH4_noN, color='mediumorchid',linewidth=0.6,label='NH4 N-less')
    # ax.plot(dates_local, TmolNH4_hindcast, color='mediumorchid',linewidth=2,alpha=0.6,label='NH4 with-N')
    # DIN
    ax.plot(dates_local, TmolDIN_noN, color='darkcyan',linewidth=0.6,label='No-loading')
    ax.plot(dates_local, TmolDIN_hindcast, color='darkcyan',linewidth=2.5,alpha=0.4,label='With-loading')
    # format subplot
    ax.set_xlim([dates_local[0],dates_local[-1]])
    ax.grid(True,color='w',linewidth=2)
    ax.set_facecolor('#EEEEEE')
    ax.set_xticklabels([])
    ax.tick_params(axis='y', labelcolor='#EEEEEE', labelsize=12, color='#EEEEEE')
    ax.set_ylim([0,1.2*max(TmolDIN_hindcast)])
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    ax.legend(loc='lower left')
    ax.text(0.02, 0.8, '(a) Depth-integrated DIN [kg]', fontweight='bold',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=12)
    
    # plot the difference between the two runs
    ax = plt.subplot(212)
    # # NO3
    # ax.plot(dates_local, TmolNO3_hindcast-TmolNO3_noN, color='deeppink',linewidth=2,label='NO3')
    # # NH4
    # ax.plot(dates_local, TmolNH4_hindcast-TmolNH4_noN, color='mediumorchid',linewidth=2,label='NH4')
    # DIN
    ax.plot(dates_local, TmolDIN_hindcast-TmolDIN_noN, color='darkcyan',linewidth=2,label='DIN')
    # format subplot
    ax.set_xlim([dates_local[0],dates_local[-1]])
    ax.grid(True,color='w',linewidth=2)
    ax.tick_params(axis='y', labelcolor='#EEEEEE', labelsize=12, color='#EEEEEE')
    ax.set_facecolor('#EEEEEE')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    # ax.legend(loc='lower left')
    ax.text(0.02, 0.8, '(b) With-loading minus No-loading', fontweight='bold',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes, fontsize=12)

    plt.suptitle(station)
    # format date label
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.set_xlabel('2013', color='#EEEEEE',fontsize=12)
    ax.tick_params(axis='x', labelrotation=30, colors='#EEEEEE', labelsize=12)
                
    plt.savefig(out_dir / (station+'_depth_integrated_N.png'))


# #%% GENERATE PLOTS ------------------------------------------------------------------------------
# plt.show()