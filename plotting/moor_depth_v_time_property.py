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

gtagexes = ['cas7_t0noN_x4b', 'cas7_t0_x4b']
jobname = 'noWWTPNtest'
startdate = '2013.01.01'
enddate = '2013.12.31'
year = '2013' # for making a date label

# gtagexes = gtagexes[0:1]

vn_list = ['rho', 'NO3', 'NH4', 'phytoplankton','zooplankton', 'SdetritusN', 'LdetritusN', 'oxygen']
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
gtagexample = gtagexes[0]
gridname, tag, ex_name = gtagexample.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)

###########################################################
# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots' / 'depth_v_time_property'
Lfun.make_dir(out_dir)


##########################################################
##                      Plotting                        ##
##########################################################

plt.close('all')

z_max = 6 # upper vertical limit (to see top of sea surface)

# Loop through all of the mooring stations
for i,station in enumerate(sta_dict.keys()):

    # calculate lat/lon for station
    lon = sta_dict[station][0]
    lat = sta_dict[station][1]

    # Initialize Figure
    fig, ax = plt.subplots(rows,2,figsize = (18,rows*2), sharex = True, sharey = True)
    fig.suptitle(station, fontsize = ts)
    
    # loop through different state variables
    for j,vn in enumerate(vn_list):

        # loop through two model runs
        for gtagex in gtagexes:

            # no loading run
            if gtagex == 'cas7_t0noN_x4b':
                title = 'Natural'
                # download .nc files
                fn = '../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
                ds = xr.open_dataset(fn)
                # get depth values
                z_rho = ds['z_rho'].transpose() # depth of u and v-velocities
                z_w   = ds['z_w'].transpose()   # depth of w-velocities
                z_min = np.min(z_w.values)
                # column number
                col = 0
                # calculate density
                if vn == 'rho':
                    # Calculate density
                    ds = ds.assign(p=gsw.p_from_z(ds['z_rho'],lat))
                    # calculate absolute salinity from practical salinity
                    ds = ds.assign(salt_abs=gsw.conversions.SA_from_SP(ds['salt'], ds['z_rho'], lon, lat))
                    # calculate conservative temperature from potential temperature
                    ds = ds.assign(temp_cons=gsw.conversions.CT_from_pt(ds['salt_abs'], ds['temp']))
                    # calculate density
                    ds = ds.assign(rho=gsw.rho(ds['salt_abs'],ds['temp_cons'],ds['p']))
                    # set scale and units
                    scale = 1
                    units = ' $(kg\ m^{-3})$'
                    vmin = 1015
                    vmax = 1025
                    cmap = 'gist_stern'
                else:
                    # get scale and units
                    scale =  pinfo.fac_dict[vn]
                    units = pinfo.units_dict[vn]
                    vlims = pinfo.vlims_dict[vn]
                    vmin = vlims[0]
                    vmax = vlims[1]
                    cmap = pinfo.cmap_dict[vn]
                    # cmap = 'rainbow_r'
                # get dataset
                val = ds[vn].transpose() * scale
                # autoscale nutrient colorbar
                if vn == 'NO3' or vn == 'NH4':
                    vmin = 0.9*np.nanmin(val)
                    vmax = 1.1*np.nanmax(val)

            elif gtagex == 'cas7_t0_x4b':
                title = 'Anthropogenic minus Natural'
                # download .nc files
                fn = '../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
                ds_withLoading = xr.open_dataset(fn)
                # get depth values
                z_rho = ds_withLoading['z_rho'].transpose() # depth of u and v-velocities
                z_w   = ds_withLoading['z_w'].transpose()   # depth of w-velocities
                z_min = np.min(z_w.values)
                # column number
                col = 1
                # calculate density
                if vn == 'rho':
                    # Calculate density
                    ds_withLoading = ds_withLoading.assign(p=gsw.p_from_z(ds_withLoading['z_rho'],lat))
                    # calculate absolute salinity from practical salinity
                    ds_withLoading = ds_withLoading.assign(salt_abs=gsw.conversions.SA_from_SP(ds_withLoading['salt'], ds_withLoading['z_rho'], lon, lat))
                    # calculate conservative temperature from potential temperature
                    ds_withLoading = ds_withLoading.assign(temp_cons=gsw.conversions.CT_from_pt(ds_withLoading['salt_abs'], ds_withLoading['temp']))
                    # calculate density
                    ds_withLoading = ds_withLoading.assign(rho=gsw.rho(ds_withLoading['salt_abs'],ds_withLoading['temp_cons'],ds_withLoading['p']))
                    # set scale and units
                    scale = 1
                    units = ' $(kg\ m^{-3})$'
                else:
                    # get scale and units
                    scale =  pinfo.fac_dict[vn]
                    units = pinfo.units_dict[vn]
                    vlims = pinfo.vlims_dict[vn]
                # caculatate difference between runs
                ds_withLoading = ds_withLoading.assign(t0_minus_t0noN=(ds_withLoading[vn]-ds[vn]))
                val = ds_withLoading['t0_minus_t0noN'].transpose() * scale
                
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
                cmap = cmocean.tools.crop(cmocean.cm.balance_r, vmin, vmax, 0)

            # need white text to see some of the labels on natural run (first column)
            if (vn == 'NH4' or vn == 'zooplankton' or vn == 'SdetritusN'
                or vn == 'LdetritusN' or vn == 'phytoplankton') and col == 0:
                font_color = 'white'
            else:
                font_color = 'black'

            # create time vector
            dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
            dates_local = [pfun.get_dt_local(x) for x in dates]

            # Plot pcolormesh
            cs = ax[j,col].pcolormesh(dates_local, z_rho, val, vmin=vmin, vmax=vmax, cmap=cmap)#cmocean.cm.haline))
            cbar = fig.colorbar(cs, ax=ax[j,col], aspect=10)
            cbar.outline.set_visible(False)
            if col == 0:
                ax[j,col].set_ylabel('z (m)', fontsize = fs)
            ax[j,col].text(0.05, 0.05, vn + units, fontweight='bold',
                    verticalalignment='bottom', horizontalalignment='left',
                    transform=ax[j,col].transAxes, fontsize=ls, color = font_color)
            ax[j,col].tick_params(axis='both', which='major', labelsize=ls)
            ax[j,col].set_ylim((z_min,z_max))
            ax[j,col].grid(True,color='k',linewidth=1,linestyle=':',axis='x')

            # add bottom axis
            if j == rows-1:
                ax[j,col].set_xlabel('2013', fontsize = fs)
                ax[j,col].tick_params(axis='both', which='major', labelsize=ls)
                ax[j,col].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
                ax[j,col].tick_params(axis='x', labelrotation=30, labelsize=ls)

            # add title
            ax[0,col].set_title(title,fontsize=ts)

            # add note about colorbar
            ax[0,1].text(0.95, 0.7, 'Anthropogenic higher',color='royalblue',
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax[0,1].transAxes, fontsize=fs)
            ax[0,1].text(0.95, 0.5, 'Natural higher',color='crimson',
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax[0,1].transAxes, fontsize=fs)

            # Look at onset of spring bloom (march 1 - may 1)
            if spring_bloom_only:
                ax[j,col].set_xlim(dates_local[59],dates_local[121])
                ax[rows-1,col].xaxis.set_major_formatter(mdates.DateFormatter("%b-%d"))

    plt.tight_layout
    plt.subplots_adjust(left=0.05, right=0.95, top=0.90, wspace=0.02)
    plt.savefig(out_dir / (station + '(' + str(round(lon,2)) + ',' + str(round(lat,2)) + ').png'))
