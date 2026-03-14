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

# initialize figure
plt.close('all')
# pfun.start_plot(figsize=(5,5))
# fig = plt.figure()
# ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)
fig, axes = plt.subplots(1,2,figsize = (8,4.5),sharey=True)
axes[0].set_ylabel('Bottom DO [mg/L]',fontsize=ls)


# Loop through all of the basins
    
# Calculate average bottom DO and assign colors to different basins
bott_DO = np.zeros(len(sta_dict))
for i,station in enumerate(sta_dict):
    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    # calculate lat/lon for station
    lon = sta_dict[station][0]
    lat = sta_dict[station][1]

    # get average bottom DO during hypoxic season
    bott_DO_all = ds['oxygen'].values[:,0] * pinfo.fac_dict['oxygen'] # mg/L

    # get stratification at this station
    # calculate absolute salinity from practical salinity
    salt_abs_top_all = gsw.conversions.SA_from_SP(ds['salt'].values[:,-1], ds['z_rho'].values[:,-1], lon, lat)
    salt_abs_bott_all = gsw.conversions.SA_from_SP(ds['salt'].values[:,0], ds['z_rho'].values[:,0], lon, lat)
    # calculate conservative temperature from potential temperature
    cons_temp_top_all = gsw.conversions.CT_from_pt(salt_abs_top_all, ds['temp'].values[:,-1])
    cons_temp_bott_all = gsw.conversions.CT_from_pt(salt_abs_bott_all, ds['temp'].values[:,0])
    # calculate potential density
    rho_top_all = gsw.density.sigma0(salt_abs_top_all,cons_temp_top_all)
    rho_bott_all = gsw.density.sigma0(salt_abs_bott_all,cons_temp_bott_all)
    # calculate density difference on every day of hypoxic season
    rho_diff_all = rho_bott_all - rho_top_all

    # download hourly .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '_hourly/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # crop to hypoxic season
    ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))
    # get barotropic velocities
    ubar = ds['ubar'].values
    vbar = ds['vbar'].values
    #U_rms = sqrt( Godin( ubar^2 + vbar^2 ) ) # this equation got Parker's OK
    U_rms = np.sqrt( zfun.lowpass(ubar**2 + vbar**2, f='godin')) 
    # average every day (24 hours)
    U_rms_daily = np.average(U_rms[0:-1].reshape(-1, 24), axis=1)
    
    # assign colors to different basins
    if station in whidbey:
        color = color_list[0]
    elif station in hoodcanal:
        color = color_list[1]
    elif station in mainbasin:
        color = color_list[2]
    elif station in southsound:
        color = color_list[3]
    else:
        color = color_list[4]

    # plot scribble
    # stratification
    axes[0].plot(rho_diff_all,bott_DO_all,alpha=0.5,linewidth=2,color=color)
    axes[0].set_xlabel(r'$\rho_{bottom} - \rho_{surface}$ [kg m$^{-3}$]',fontsize=ls)
    axes[0].set_title(r'(a) vs. daily $\Delta\rho$', loc='left', fontsize=ts)
    # tidal currents
    axes[1].plot(U_rms_daily,bott_DO_all,alpha=0.5,linewidth=2,color=color)
    axes[1].set_xlabel(r'$U_{rms}$ [m s$^{-1}$]',fontsize=ls)
    axes[1].set_title(r'(b) vs. daily avg. $U_{rms}$', loc='left', fontsize=ts)

# format grid
for ax in axes:
    ax.tick_params(axis='both', which='major', labelsize=ls)
    ax.grid(True,color='white',linewidth=1)
    # format colors
    ax.set_facecolor('#EEEEEE')
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)

# add basin label
axes[0].text(0.95, 0.9, 'Whidbey Basin',color=color_list[0],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[0].transAxes, fontsize=ls, fontweight='bold')
axes[0].text(0.95, 0.85, 'Hood Canal',color=color_list[1],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[0].transAxes, fontsize=ls, fontweight='bold')
axes[1].text(0.95, 0.80, 'Main Basin',color=color_list[2],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[0].transAxes, fontsize=ls, fontweight='bold')
axes[0].text(0.95, 0.75, 'South Sound',color=color_list[3],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[0].transAxes, fontsize=ls, fontweight='bold')
axes[0].text(0.95, 0.70, 'Admiralty Inlet',color=color_list[4],
                        verticalalignment='bottom', horizontalalignment='right',
                        transform=axes[0].transAxes, fontsize=ls, fontweight='bold')

plt.suptitle(r'Daily Aug/Sep DO$_{Bot}$ vs. physical mechanisms',fontsize=ts+1)

plt.subplots_adjust(wspace=0.1,top=0.82, bottom=0.16,right=0.9,left=0.1)
plt.savefig(out_dir / 'scribble_plot.png')


# ##########################################################
# ##       Average bottom DO vs. tidal currents           ##
# ##########################################################

# plt.close('all')

# pfun.start_plot(figsize=(5,5))
# fig = plt.figure()
# ax = fig.add_subplot(1,1,1)
# plt.subplots_adjust(wspace=0, hspace=0.1)

# # initialize arrays for plotting
# tcurr = np.zeros(len(sta_dict))

# for i,station in enumerate(sta_dict):

#     # download .nc files
#     fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '_hourly/' + station + '_' + startdate + '_' + enddate + '.nc'
#     ds = xr.open_dataset(fn)

#     # get average tidal current
#     ubar = ds['ubar'].values
#     vbar = ds['vbar'].values
#     # calculate U_rms
#     #U_rms = sqrt( Godin( ubar^2 + vbar^2 ) ) # this equation got Parker's OK
#     U_rms = np.sqrt( zfun.lowpass(ubar**2 + vbar**2, f='godin')) 
#     tcurr[i] = np.nanmean(U_rms) # average over time
    
# # create scatter plot
# # ax.plot(tcurr,bott_DO,linestyle='none',marker='o',color='navy',alpha=0.5,markersize=10)
# ax.scatter(tcurr,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)

# # calculate correlation coefficient (Pearson)
# r,p = pearsonr(tcurr, bott_DO)
# r_Urms = r
# ax.text(0.1, 0.85, r'$r =$' + str(round(r,2)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=ax.transAxes, fontsize=ts, fontweight='bold')
# ax.text(0.1, 0.79, r'$p =$' + str(round(p,4)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=ax.transAxes, fontsize=ts, fontweight='bold')

# # add basin label
# ax.text(0.95, 0.95, 'Whidbey Basin',color=color_list[0],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
# ax.text(0.95, 0.91, 'Hood Canal',color=color_list[1],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
# ax.text(0.95, 0.87, 'Main Basin',color=color_list[2],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
# ax.text(0.95, 0.83, 'South Sound',color=color_list[3],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=ax.transAxes, fontsize=ts-3, fontweight='bold')
# ax.text(0.95, 0.79, 'Admiralty Inlet',color=color_list[4],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=ax.transAxes, fontsize=ts-3, fontweight='bold')

# # format labels
# ax.set_title('Average bottom DO vs. Tidal currents\n('+year+'-'+start+' to '+year+'-'+end+')',
#              fontsize=ts)
# ax.set_xlabel('Annual average tidal currents [m/s]',fontsize=ls)
# ax.set_ylabel('Avg. bottom DO [mg/L]',fontsize=ls)

# # format grid
# ax.tick_params(axis='both', which='major', labelsize=ls)
# ax.grid(True,color='white',linewidth=1)
# # format colors
# ax.set_facecolor('#EEEEEE')
# for border in ['top','right','bottom','left']:
#     ax.spines[border].set_visible(False)

# # add inlet label
# if inlet_label:
#     for i, txt in enumerate(sta_dict):
#         ax.annotate(txt, (tcurr[i],bott_DO[i]), fontsize=10)

# plt.savefig(out_dir / 'botDO_vs_tidalcurr.png')


# ##########################################################
# ##               Plot all in subplots                   ##
# ##########################################################

# # histrogram of bottom DO
# fig, ax = plt.subplots(1,1,figsize = (6,5))

# # create histogram
# bins = [0,1,2,3,4,5,6,7,8,9,10]
# ax.hist(bott_DO, bins=bins, zorder=5, facecolor='mediumorchid',
#         alpha=0.6, edgecolor='darkviolet')
# ax.set_xlim([0,10])
# ax.set_xlabel(r'Avg. Aug/Sep $DO_{bot}$ [mg/L]')
# ax.set_ylim([0,5])
# ax.set_ylabel('Count')
# ax.set_title(r'Distribution of $DO_{bot}$ in 21 inlets' + '\n(averaged over Aug/Sep)')

# # format grid
# ax.tick_params(axis='both', which='major', labelsize=ls)
# ax.grid(True,color='white',linewidth=1)
# # format colors
# ax.set_facecolor('#EEEEEE')
# for border in ['top','right','bottom','left']:
#     ax.spines[border].set_visible(False)
# plt.savefig(out_dir / 'bottDO_hist.png')

# # --------------------------------------------------------

# plt.close('all')

# fig, axes = plt.subplots(1,4,figsize = (10,3.4),sharey=True)
# axes[0].set_ylabel('Avg. bottom DO [mg/L]',fontsize=ls)
# # add basin label
# axes[3].text(0.95, 0.9, 'Whidbey Basin',color=color_list[0],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
# axes[3].text(0.95, 0.81, 'Hood Canal',color=color_list[1],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
# axes[3].text(0.95, 0.72, 'Main Basin',color=color_list[2],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
# axes[3].text(0.95, 0.63, 'South Sound',color=color_list[3],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=axes[3].transAxes, fontsize=ls, fontweight='bold')
# axes[3].text(0.95, 0.54, 'Admiralty Inlet',color=color_list[4],
#                         verticalalignment='bottom', horizontalalignment='right',
#                         transform=axes[3].transAxes, fontsize=ls, fontweight='bold')

# # add data
# # Depth
# axes[0].scatter(depth,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
# axes[0].set_title(r'(a) vs. depth', loc='left', fontsize=ts)
# axes[0].set_xlabel('Depth [m]', fontsize=ls)
# axes[0].text(0.1, 0.84, r'$r =$' + str(round(r_depth,2)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=axes[0].transAxes, fontsize=ts)
# # delta rho
# axes[1].scatter(strat,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
# axes[1].set_title(r'(b) vs. Aug/Sep $\Delta\rho$', loc='left', fontsize=ts)
# axes[1].set_xlabel(r'$\rho_{bottom} - \rho_{surface}$ [kg m$^{3}$]', fontsize=ls)
# axes[1].text(0.1, 0.84, r'$r =$' + str(round(r_drho,2)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=axes[1].transAxes, fontsize=ts)
# # RMS Tidal currents
# axes[2].scatter(tcurr,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
# axes[2].set_title(r'(c) vs. annual avg. $U_{rms}$', loc='left', fontsize=ts)
# axes[2].set_xlabel(r'$U_{rms}$ [m s$^{-1}$]', fontsize=ls)
# axes[2].text(0.1, 0.84, r'$r =$' + str(round(r_Urms,2)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=axes[2].transAxes, fontsize=ts)
# # flushing time
# axes[3].scatter(Tflush,bott_DO,c=colors,alpha=0.5,s=100,cmap=cmap,zorder=5)
# axes[3].set_title(r'(d) vs. annual avg. $T_{flush}$', loc='left', fontsize=ts)
# axes[3].set_xlabel(r'$T_{flush} = V/Q_{prism}$ [days]', fontsize=ls)
# axes[3].text(0.1, 0.84, r'$r =$' + str(round(r_Tflush,2)) ,color='black',
#                         verticalalignment='bottom', horizontalalignment='left',
#                         transform=axes[3].transAxes, fontsize=ts)

# for ax in axes:
#     # format grid
#     ax.tick_params(axis='both', which='major', labelsize=ls)
#     ax.grid(True,color='white',linewidth=1)
#     # format colors
#     ax.set_facecolor('#EEEEEE')
#     for border in ['top','right','bottom','left']:
#         ax.spines[border].set_visible(False)

# plt.suptitle(r'Aug/Sep DO$_{Bot}$ vs. physical characteristics',fontsize=ts+1)

# plt.subplots_adjust(wspace=0.05,top=0.82, bottom=0.16,right=0.99,left=0.06)
# plt.savefig(out_dir / 'physics_summary.png')