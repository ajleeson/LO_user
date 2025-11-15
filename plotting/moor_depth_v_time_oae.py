"""
Generate depth vs. time property plots using mooring extraction data. 
Used to compare the with-loading run to the no-loading run.
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
import gsw
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex_nomodule = 'cas7_t1jxoae_x11bjx'
gtagex_module   = 'cas7_t1jxoae_x11ecb'
jobname = 'oae_mod_test'
startdate = '2013.07.01'
enddate = '2013.07.01'
enddate_plus1 = '2013.07.02'
year = '2013' # for making a date label

# gtagexes = gtagexes[0:1]

vn_list = ['alkalinity', 'TIC']
rows = len(vn_list)

# figure settings
fs = 14 # figure font size
ls = 12 # label size
ts = 16 # title size

# look at full year or only spring bloom?
spring_bloom_only  = False

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# parse gtagex
gtagexample = gtagex_nomodule
gridname, tag, ex_name = gtagexample.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)

###########################################################
# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots' / 'oae_test'
Lfun.make_dir(out_dir)


##########################################################
##                      Plotting                        ##
##########################################################

plt.close('all')

z_max = 5 # upper vertical limit (to see top of sea surface)

station = list(sta_dict.keys())[0]

# calculate lat/lon for station
lon = sta_dict[station][0]
lat = sta_dict[station][1]

# Initialize Figure
fig, ax = plt.subplots(rows,3,figsize = (10,rows*3.5), sharex = True, sharey = True)
fig.suptitle(station, fontsize = ts)
    
# loop through different state variables
for j,vn in enumerate(vn_list):

    # get scale and units
    scale =  pinfo.fac_dict[vn]
    units = pinfo.units_dict[vn]
    vlims = pinfo.vlims_dict[vn]
    if vn == 'alkalinity':
        vmin = 0
        vmax = 600
    if vn == 'TIC':
        vmin = 0
        vmax = 20#2000
    cmap = 'inferno_r'#'rainbow_r'

    # create time vector
    dates = pd.date_range(start= startdate, end= enddate_plus1, freq= 'h')
    dates_local = [pfun.get_dt_local(x) for x in dates]

    # no oae module run ----------------------------------------------------------------
    # add title
    ax[0,0].set_title('No OAE module\n(minus ambient)',fontsize=ts)
    # download .nc files
    fn = '../../LO_output/extract/' + gtagex_nomodule + '/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds_noModule = xr.open_dataset(fn)
    # get dataset
    val = ds_noModule[vn].transpose() * scale
    # subtract ambient alkalinity concentrations
    fn = '../../LO_output/extract/' + 'cas7_t1jxoae_x11b' + '/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds_ambient = xr.open_dataset(fn)
    ambient = ds_ambient[vn].transpose() * scale
    val = val - ambient
    # get depth values
    z_rho = ds_noModule['z_rho'].transpose() # depth of u and v-velocities
    z_w   = ds_noModule['z_w'].transpose()   # depth of w-velocities
    z_min = np.min(z_w.values)

    # Plot pcolormesh
    cs = ax[j,0].pcolormesh(dates_local, z_rho, val, vmin=vmin, vmax=vmax, cmap=cmap)#cmocean.cm.haline))
    cbar = fig.colorbar(cs, ax=ax[j,0], aspect=10)
    cbar.outline.set_visible(False)
    ax[j,0].set_ylabel('z (m)', fontsize = fs)
    ax[j,0].text(0.05, 0.05, vn + units, fontweight='bold',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax[j,0].transAxes, fontsize=ls, color = 'k')
    ax[j,0].tick_params(axis='both', which='major', labelsize=ls)
    ax[j,0].set_ylim((z_min,z_max))
    ax[j,0].grid(True,color='k',linewidth=1,linestyle=':',axis='x')

    # with oae module run ----------------------------------------------------------------
    # add title
    ax[0,1].set_title('With OAE module\n(minus ambient)',fontsize=ts)
    # download .nc files
    fn = '../../LO_output/extract/' + gtagex_module + '/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds_Module = xr.open_dataset(fn)
    # get dataset
    val = ds_Module[vn].transpose() * scale
    val = val - ambient
    # get depth values
    z_rho = ds_Module['z_rho'].transpose() # depth of u and v-velocities
    z_w   = ds_Module['z_w'].transpose()   # depth of w-velocities
    z_min = np.min(z_w.values)

    # Plot pcolormesh
    cs = ax[j,1].pcolormesh(dates_local, z_rho, val, vmin=vmin, vmax=vmax, cmap=cmap)#cmocean.cm.haline))
    cbar = fig.colorbar(cs, ax=ax[j,1], aspect=10)
    cbar.outline.set_visible(False)
    ax[j,1].set_ylabel('z (m)', fontsize = fs)
    ax[j,1].text(0.05, 0.05, vn + units, fontweight='bold',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax[j,1].transAxes, fontsize=ls, color = 'k')
    ax[j,1].tick_params(axis='both', which='major', labelsize=ls)
    ax[j,1].set_ylim((z_min,z_max))
    ax[j,1].grid(True,color='k',linewidth=1,linestyle=':',axis='x')


    # difference between runs ----------------------------------------------------------------
    # add title
    ax[0,2].set_title('With OAE module minus No Module',fontsize=ts)

    # caculatate difference between runs
    ds_Module = ds_Module.assign(t0_minus_t0noN=(ds_Module[vn]-ds_noModule[vn]))
    val = ds_Module['t0_minus_t0noN'].transpose() * scale
    # val = ds_withModule[vn].transpose() * scale
    
    # get min and max for plotting
    # (we use the average min/max value in time multiplied by a scale because using
    # the straight min/max values obscures small differences-- since min/max are large)
    factor = 4#3
    vmin = np.nanmin(val.values,axis=0)
    vmin = factor*np.nanmean(vmin)
    vmax = np.nanmax(val.values,axis=0)
    vmax = factor*np.nanmean(vmax)
    # make sure colorbar axis contains zero
    if vmin > 0 and vmax > 0:
        vmin = vmax*-1.01
    if vmin < 0 and vmax < 0:
        vmax = vmin*-1.01
    if vmin == 0 and vmax == 0:
        vmin = -0.11
        vmax = 0.1
    vmin = -10
    vmax = 10
    # cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)
    cmap = cmocean.tools.crop_by_percent(cmocean.cm.balance_r, 20, which='both', N=None)

    # Plot pcolormesh
    cs = ax[j,2].pcolormesh(dates_local, z_rho, val, vmin=vmin, vmax=vmax, cmap=cmap)#cmocean.cm.haline))
    cbar = fig.colorbar(cs, ax=ax[j,2], aspect=10)
    cbar.outline.set_visible(False)
    ax[j,2].set_ylabel('z (m)', fontsize = fs)
    ax[j,2].text(0.05, 0.05, vn + units, fontweight='bold',
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax[j,2].transAxes, fontsize=ls, color = 'k')
    ax[j,2].tick_params(axis='both', which='major', labelsize=ls)
    ax[j,2].set_ylim((z_min,z_max))
    ax[j,2].grid(True,color='k',linewidth=1,linestyle=':',axis='x')

plt.tight_layout
# plt.subplots_adjust(left=0.05, right=0.95, top=0.90, wspace=0.02)
plt.show()


##########################################################
##                   Depth Profile                      ##
##########################################################


# Initialize Figure
fig, ax = plt.subplots(1,2,figsize = (8,6))
fig.suptitle(station, fontsize = ts)

hour = 12

ax[0].set_title('Alkalinity profile [meq/m3]')
ax[0].plot(ds_noModule['alkalinity'][hour,:],z_rho[:,hour],
           color='hotpink', linewidth=3, linestyle='-',label='No Module')
ax[0].plot(ds_Module['alkalinity'][hour,:],z_rho[:,hour],
           color='black', linewidth=3, linestyle='--', label='With Module',)

ax[0].legend(loc='best')

ax[1].set_title('TIC profile [mmol/m3]')
ax[1].plot(ds_noModule['TIC'][hour,:],z_rho[:,hour],
           color='hotpink', linewidth=3, linestyle='-',label='No Module')
ax[1].plot(ds_Module['TIC'][hour,:],z_rho[:,hour],
           color='black', linewidth=3, linestyle='--', label='With Module',)

plt.show()