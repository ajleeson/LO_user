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

vn_list = ['rho','phytoplankton','oxygen']#['rho', 'NO3', 'NH4', 'phytoplankton','zooplankton', 'SdetritusN', 'LdetritusN', 'oxygen']
rows = len(vn_list)

# figure settings
fs = 12 # figure font size
ls = 11 # label size
ts = 14 # title size

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
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures' / 'twentyone' / 'depthVtime_property' / 'obsraw'
Lfun.make_dir(out_dir)

# get observations
in_dir = Ldir['parent'] / 'LO_output' / 'obsmod'
in_fn = in_dir / ('multi_ctd_' + year + '.p')
df_dict = pickle.load(open(in_fn, 'rb'))
# only look at ecology stations
source = 'ecology_nc'
df_obs = df_dict['obs'].loc[df_dict['obs'].source==source,:]

##########################################################
##                 Plot station locations               ##
##########################################################

plt.close('all')

# dictionary of ecology stations
ecol_stn = {
    'sinclair': 'SIN001',
    'elliot': 'ELB015',
    'lynchcove': 'HCB007',
    'commencement': 'CMB003',
    'hammersley': 'OAK004',
    'totten': 'TOT002',
    'budd': 'BUD005'
}

# Puget Sound only
lat_low = 46.95
lat_high = 48.5
lon_low = -123.5
lon_high = -122

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

pfun.start_plot(figsize=(6,9))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# Plot map
ax.add_patch(Rectangle((lon_low, lat_low), lon_high-lon_low,lat_high-lat_low, facecolor='#EEEEEE'))
plt.pcolormesh(plon, plat, zm, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
# pfun.add_coast(ax, color='gray')
ax.axis([lon_low,lon_high,lat_low,lat_high])
pfun.dar(ax)
ax.set_xlabel('')
ax.set_ylabel('')
plt.xticks(rotation=30)
for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)

# Plot cast locations
for sta in sta_dict:
    sta_lon = sta_dict[sta][0]
    sta_lat = sta_dict[sta][1]
    lon_off = 0
    lat_off = 0.04
    ha = 'center'
    ecol = ''
    # format depending on whether there is an ecology station or not
    if sta in ecol_stn:
         color = 'mediumorchid'
         ecol = '\n('+ ecol_stn[sta] +')'
    else:
         color = 'k'
    ax.plot(sta_lon,sta_lat,linestyle='none',marker='o',color=color,markeredgecolor = 'k',markeredgewidth=0.5)
    # move text to make things more visible
    if sta_lon > -122.61 or sta in ['henderson','budd']:
        lon_off = 0.03
        lat_off = 0
        ha = 'left'
    if sta_lon < -123 or sta in ['oak','sinclair','case','penn','dabob']:
        lon_off = -0.02
        lat_off = 0
        ha = 'right'
    if sta in ['budd','totten']:
        lat_off = -0.02
    ax.text(sta_lon+lon_off,sta_lat+lat_off,sta+ecol,va='center',ha=ha,color=color,fontsize = 11)

plt.title('Inlet Locations',fontsize = 16)
fig.tight_layout

# add labels
ax.add_patch(Rectangle((-123.47, 48.3), 0.65,0.15, facecolor='white',alpha=0.5))
ax.text(-123.45,48.4,'Ecology monitoring\nstation (ID)',va='center',ha='left',fontweight='bold',color='mediumorchid',fontsize = 12)
ax.text(-123.45,48.33,'No observations',va='center',ha='left',fontweight='bold',color='k',fontsize = 12)

plt.savefig(out_dir / ('inlet_extraction_locs.png'))

##########################################################
##                      Plotting                        ##
##########################################################

z_max = 5 # upper vertical limit (to see top of sea surface)

letter = ['(a)','(b)','(c)','(d)',
          '(e)','(f)','(g)','(h)',
          '(i)','(j)','(k)','(l)']

# Loop through all of the mooring stations
for i,station in enumerate(sta_dict): # enumerate(['budd','penn']):
    plt.close('all')

    print(station)

    # calculate lat/lon for station
    lon = sta_dict[station][0]
    lat = sta_dict[station][1]

    # check if station has ecology data
    if station in ecol_stn:
        cols = 2
    else:
        cols = 1

    # Initialize Figure
    scale = 3 #1.2
    fig, ax = plt.subplots(rows,cols,figsize = (9*cols,rows*scale), sharex = True, sharey = True)
    fig.suptitle(station + ' ' + year, fontsize = 18)
    
    # loop through different state variables
    for j,vn in enumerate(vn_list):

        # MODEL MODEL MODEL MODEL MODEL MODEL MODEL MODEL  ----------------------------------------------------------------
        # download .nc files
        fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds = xr.open_dataset(fn)
        # get depth values
        z_rho = ds['z_rho'].transpose() # depth of u and v-velocities
        z_rho = z_rho.values
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
            cmap = cmocean.cm.dense
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
        # vmin = np.nanmin(val)
        # vmax = np.nanmax(val)

        # autoscale nutrient colorbar
        # if vn == 'NO3' or vn == 'NH4':
        #     vmin = 0.9*np.nanmin(val)
        #     vmax = 1.1*np.nanmax(val)

        # index through axes differently if there are obs vs not
        if station in ecol_stn:
            axis = ax[j,col]
        else:
            axis = ax[j]

        # need white text to see some of the labels on natural run (first column)
        if vn in ['rho','NO3'] and col == 0:
            font_color = 'white'
        else:
            font_color = 'black'

        # create time vector
        dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
        dates_local = [pfun.get_dt_local(x) for x in dates]

        # plot
        cs = axis.pcolormesh(dates_local, z_rho, val, vmin=vmin, vmax=vmax, cmap=cmap)
        # place colorbar
        cbar = fig.colorbar(cs, ax=axis, aspect=10)
        cbar.outline.set_visible(False)
        if col == 0:
            axis.set_ylabel('z (m)', fontsize = fs)
        axis.text(0.02, 0.05, letter[j] + ' ' + vn + units, fontweight='bold',
                verticalalignment='bottom', horizontalalignment='left',
                transform=axis.transAxes, fontsize=ls, color = font_color)
        axis.tick_params(axis='both', which='major', labelsize=ls)
        axis.set_ylim((z_min,z_max))
        axis.grid(True,color='k',linewidth=1,linestyle=':',axis='x')

        # title
        if j == 0 and station in ecol_stn:
            axis.set_title('Model')

        # format colors
        axis.set_facecolor('#EEEEEE')
        for border in ['top','right','bottom','left']:
            axis.spines[border].set_visible(False)

        # add bottom axis
        if j == rows-1:
            axis.set_xlabel(year, fontsize = fs)
            axis.tick_params(axis='both', which='major', labelsize=ls)
            axis.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
            axis.tick_params(axis='x', labelrotation=30, labelsize=ls)

        
        # OBSERVATIONS OBSERVATIONS OBSERVATIONS OBSERVATIONS ----------------------------------------------------------------
        if station in ecol_stn:
            col = 1
            axis = ax[j,col]

            if j == 0:
                axis.set_title('Observations')

            # get current station
            df_ob_stn = df_obs.loc[df_obs.name==ecol_stn[station],:]

            # only plot if observations contain data
            if vn == 'oxygen':
                # get data
                var = df_ob_stn['DO (uM)']*pinfo.fac_dict['oxygen'] # convert to mg/L
                z = df_ob_stn['z']
                time = [pd.Timestamp(x) for x in df_ob_stn['time']]
                axis.scatter(time,z, c=var, vmin=vmin, vmax=vmax, cmap=cmap)
            if vn == 'phytoplankton':
                # get data
                var = df_ob_stn['Chl (mg m-3)']
                z = df_ob_stn['z']
                time = [pd.Timestamp(x) for x in df_ob_stn['time']]
                axis.scatter(time,z, c=var, vmin=vmin, vmax=vmax, cmap=cmap)
            if vn == 'rho':
                z = df_ob_stn['z']
                time = [pd.Timestamp(x) for x in df_ob_stn['time']]
                # Calculate density
                p = gsw.p_from_z(z,lat)
                # calculate density
                rho = gsw.rho(df_ob_stn['SA'],df_ob_stn['CT'],p)
                axis.scatter(time,z, c=rho, vmin=vmin, vmax=vmax, cmap=cmap)


            # format grid
            axis.tick_params(axis='both', which='major', labelsize=ls)
            axis.set_ylim((z_min,z_max))
            axis.grid(True,color='k',linewidth=1,linestyle=':',axis='x')

            # format colors
            axis.set_facecolor('#EEEEEE')
            for border in ['top','right','bottom','left']:
                axis.spines[border].set_visible(False)

            # add bottom axis
            if j == rows-1:
                axis.set_xlabel(year, fontsize = fs)
                axis.tick_params(axis='both', which='major', labelsize=ls)
                axis.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
                axis.tick_params(axis='x', labelrotation=30, labelsize=ls)

            # place colorbar
            cbar = fig.colorbar(cs, ax=axis, aspect=10)
            cbar.outline.set_visible(False)

    plt.tight_layout
    plt.subplots_adjust(hspace = 0.2, wspace=0.02, top=0.93)
    plt.savefig(out_dir / (station + '(' + str(round(lon,2)) + ',' + str(round(lat,2)) + ').png'))
