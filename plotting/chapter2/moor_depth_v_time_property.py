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
import cmcrameri.cm as cmc
import matplotlib.pylab as plt
import gsw
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

# reload to make editing easier
from importlib import reload
reload(pinfo)

Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagexes = ['cas7_t1noDIN_x11ab', 'cas7_t1_x11ab']
jobname = 'wwtp_tests'
# startdate = '2012.10.07'
# enddate = '2013.10.07'
startdate = '2017.01.01'
enddate = '2017.12.31'
# year = '2017' # for making a date label

WWTP_loc = True

# gtagexes = gtagexes[0:1]

# vn_list = ['rho', 'NO3', 'NH4', 'phytoplankton','zooplankton', 'SdetritusN', 'LdetritusN', 'oxygen']
vn_list = ['NO3', 'NH4', 'phytoplankton','zooplankton', 'SdetritusN', 'LdetritusN', 'oxygen']
vn_list_pecs = ['phytoplankton','SdetritusN','oxygen','NO3']
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

# ###########################################################
# # where to put output figures
# out_dir = Ldir['LOo'] / 'AL_custom_plots' / 'depth_v_time_property'
# Lfun.make_dir(out_dir)


##########################################################
##                      Plotting                        ##
##########################################################

plt.close('all')

def generate_axes(fig):
    gridspec = fig.add_gridspec(nrows=4, ncols=12, wspace=0.8, hspace=0.1)
    axes = {}
    axes['map'] = fig.add_subplot(gridspec[0:4, 0:2])
    axes['mainbasin']     = fig.add_subplot(gridspec[0:1, 2:7])
    axes['lynchcove']     = fig.add_subplot(gridspec[1:2, 2:7], sharex=axes['mainbasin'])
    axes['penncove']      = fig.add_subplot(gridspec[2:3, 2:7], sharex=axes['mainbasin'])
    axes['caseinlet']     = fig.add_subplot(gridspec[3:4, 2:7], sharex=axes['mainbasin'])
    axes['mainbasindiff'] = fig.add_subplot(gridspec[0:1, 7:12],sharex=axes['mainbasin'],sharey=axes['mainbasin'])
    axes['lynchcovediff'] = fig.add_subplot(gridspec[1:2, 7:12],sharex=axes['mainbasin'],sharey=axes['lynchcove'])
    axes['penncovediff']  = fig.add_subplot(gridspec[2:3, 7:12],sharex=axes['mainbasin'],sharey=axes['penncove'])
    axes['caseinletdiff'] = fig.add_subplot(gridspec[3:4, 7:12],sharex=axes['mainbasin'],sharey=axes['caseinlet'])
    return axes

# read in masks
basin_mask_ds = xr.open_dataset('../../../LO_output/chapter_2/data/basin_masks_from_pugetsoundDObox.nc')
mask_rho = basin_mask_ds.mask_rho.values
lon = basin_mask_ds['lon_rho'].values
lat = basin_mask_ds['lat_rho'].values
plon, plat = pfun.get_plon_plat(lon,lat)

# Loop through all of the mooring stations
for i,vn in enumerate(vn_list):

    # Initialize figure
    fig = plt.figure(figsize=(14,8))
    ax = generate_axes(fig)
    
    # loop through different state variables
    for j,station in enumerate(list(sta_dict.keys())):

        # calculate lat/lon for station
        lon = sta_dict[station][0]
        lat = sta_dict[station][1]

        # plot map with station locations
        ax['map'].pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho), vmin=0, vmax=1.2, cmap='bone' )
        ax['map'].scatter(lon,lat,marker='o',color='navy',s=10,zorder=5)
        ax['map'].tick_params(axis='both', labelbottom=False, labelleft=False)
        pfun.dar(ax['map'])

        # loop through two model runs
        for gtagex in gtagexes:

            # no loading run
            if gtagex == 'cas7_t1noDIN_x11ab':
                title = 'No-loading'
                # download .nc files
                fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
                ds = xr.open_dataset(fn)
                # get depth values
                z_rho = ds['z_rho'].transpose() # depth of u and v-velocities
                z_w   = ds['z_w'].transpose()   # depth of w-velocities
                z_min = np.min(z_w.values)
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

            elif gtagex == 'cas7_t1_x11ab':
                title = r'Loading $-$ No-loading'
                # download .nc files
                fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
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
                # get min and max for plotting
                vlims = pinfo.vlims_dict[vn]
                if vn == 'NH4':
                    scale = 1
                elif vn == 'NO3':
                    scale = 0.1
                else:
                    scale = 0.01
                vmin = vlims[1] * scale * -1
                vmax = vlims[1] * scale
                cmap = cmocean.cm.balance_r

            # create time vector
            dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
            dates_local = [pfun.get_dt_local(x) for x in dates]

            # Plot pcolormesh
            if gtagex == 'cas7_t1noDIN_x11ab':
                cs_noload = ax[station].pcolormesh(dates_local, z_rho, val, vmin=vmin, vmax=vmax, cmap=cmap)
            else:
                cs_diff = ax[station+'diff'].pcolormesh(dates_local, z_rho, val, vmin=vmin, vmax=vmax, cmap=cmap)

            # formatting
            ax[station].text(0.03, 0.08, station, fontweight='bold', transform=ax[station].transAxes, fontsize=ls, color = 'k')
            if station == 'caseinlet':
                ax[station].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
                ax[station+'diff'].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
            else:
                ax[station].tick_params(axis='x', which='both', labelbottom=False)
                ax[station+'diff'].tick_params(axis='x', which='both', labelbottom=False)
            ax[station+'diff'].tick_params(axis='y', which='both', labelleft=False)
            fig.subplots_adjust(top=0.8, left=0.1, right=0.9)
            fig.suptitle(vn+' '+units, fontsize = ts)

            # make colorbar
            # No-loading
            if gtagex == 'cas7_t1noDIN_x11ab':
                bbox = ax['mainbasin'].get_position()
                pad = 0.02
                h = 0.02
                cax1 = fig.add_axes([bbox.x0, bbox.y1 + pad, bbox.width, h])
                cb1 = fig.colorbar(cs_noload, cax=cax1, orientation='horizontal')
                cb1.outline.set_visible(False)
                cb1.ax.xaxis.set_ticks_position('top')
                cb1.ax.xaxis.set_label_position('top')
                cb1.set_label('No-loading')
            # diff
            else:
                bbox2 = ax['mainbasindiff'].get_position()
                cax2 = fig.add_axes([bbox2.x0, bbox2.y1 + pad, bbox2.width, h])
                cb2 = fig.colorbar(cs_diff, cax=cax2, orientation='horizontal')
                cb2.outline.set_visible(False)
                cb2.ax.xaxis.set_ticks_position('top')
                cb2.ax.xaxis.set_label_position('top')
                cb2.set_label('Difference')

##############################################################
##                    HELPER FUNCTIONS                      ##
##############################################################

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = '../../../LO_data/trapsD00/LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_LO = repeatrivs_df.loc[repeatrivs_df['SSM_rname'] == rname, 'LO_rname'].values[0]
    return rname_LO

def LO2SSM_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = Ldir['data'] / 'trapsD00' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_SSM = repeatrivs_df.loc[repeatrivs_df['LO_rname'] == rname, 'SSM_rname'].values[0]
    return rname_SSM


if WWTP_loc == True:
    # set up the time index for the record
    Ldir = Lfun.Lstart()
    dsf = Ldir['ds_fmt']
    dt0 = datetime.strptime('2020.01.01',dsf)
    dt1 = datetime.strptime('2020.12.31',dsf)
    days = (dt0, dt1)
        
    # pandas Index objects
    dt_ind = pd.date_range(start=dt0, end=dt1)
    yd_ind = pd.Index(dt_ind.dayofyear)

    # Get LiveOcean grid info --------------------------------------------------

    # get the grid data
    ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
    lon_wwtp = ds.lon_rho.values
    lat_wwtp = ds.lat_rho.values
    X = lon_wwtp[0,:] # grid cell X values
    Y = lat_wwtp[:,0] # grid cell Y values

    # get flow, nitrate, and ammonium values
    fp_wwtps = '../../../LO_output/pre/trapsP01/moh20_wwtps/lo_base/Data_historical/'
    moh20_flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow.p')    # m3/s
    moh20_no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3.p')      # mmol/m3
    moh20_nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4.p')      # mmol/m3

    fp_wwtps = '../../../LO_output/pre/trapsP01/was24_wwtps/lo_base/Data_historical/'
    was24_flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow.p')    # m3/s
    was24_no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3.p')      # mmol/m3
    was24_nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4.p')      # mmol/m3

    # calculate total DIN concentration in mg/L
    moh20_dindf_wwtps = (moh20_no3df_wwtps + moh20_nh4df_wwtps)/71.4    # mg/L
    was24_dindf_wwtps = (was24_no3df_wwtps + was24_nh4df_wwtps)/71.4    # mg/L

    # calculate daily loading timeseries in kg/d
    moh20_dailyloaddf_wwtps = 86.4*moh20_dindf_wwtps*moh20_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s
    was24_dailyloaddf_wwtps = 86.4*was24_dindf_wwtps*was24_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s

    # calculate average daily load over the year (kg/d)
    moh20_avgload_wwtps = moh20_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')
    was24_avgload_wwtps = was24_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

    # add row and col index for plotting on LiveOcean grid
    griddf0_wwtps = pd.read_csv('../../../LO_data/grids/cas7/moh20_wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    moh20_avgload_wwtps = moh20_avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    moh20_avgload_wwtps = moh20_avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    griddf0_wwtps = pd.read_csv('../../../LO_data/grids/cas7/was24_wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    was24_avgload_wwtps = was24_avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    was24_avgload_wwtps = was24_avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    # get point source lat and lon
    moh20_lon_wwtps = [X[int(col)] for col in moh20_avgload_wwtps['col_py']]
    moh20_lat_wwtps = [Y[int(row)] for row in moh20_avgload_wwtps['row_py']]
    was24_lon_wwtps = [X[int(col)] for col in was24_avgload_wwtps['col_py']]
    was24_lat_wwtps = [Y[int(row)] for row in was24_avgload_wwtps['row_py']]
    
    # define marker sizes (minimum size is 10 so dots don't get too small)
    moh20_sizes_wwtps = [max(0.05*load,5) for load in moh20_avgload_wwtps['avg-daily-load(kg/d)']]
    was24_sizes_wwtps = [max(0.05*load,5) for load in was24_avgload_wwtps['avg-daily-load(kg/d)']]

##########################################################
##                   PECS  Plotting                     ##
##########################################################

def generate_axes(fig):
    # 12 cols total:
    # - map uses first 3 cols
    # - center uses next 4 cols
    # - right uses last 5 cols
    # Make center and right equal total width:
    #   4 * r_center = 5 * r_right  -> choose r_center=1.25, r_right=1.0
    width_ratios = [1,1,1,  1.25,1.25,1.25,1.25,  1,1,1,1,1]
    gridspec = fig.add_gridspec(nrows=4, ncols=12,width_ratios=width_ratios,wspace=0.02, hspace=0.04)
    axes = {}
    axes['map'] = fig.add_subplot(gridspec[0:4, 0:3])
    axes['phytoplankton'] = fig.add_subplot(gridspec[0:1, 3:7])
    axes['SdetritusN'] = fig.add_subplot(gridspec[1:2, 3:7])
    axes['oxygen'] = fig.add_subplot(gridspec[2:3, 3:7])
    axes['NO3'] = fig.add_subplot(gridspec[3:4, 3:7])
    axes['phytoplanktondiff'] = fig.add_subplot(gridspec[0:1, 7:12])
    axes['SdetritusNdiff'] = fig.add_subplot(gridspec[1:2, 7:12])
    axes['oxygendiff'] = fig.add_subplot(gridspec[2:3, 7:12])
    axes['NO3diff'] = fig.add_subplot(gridspec[3:4, 7:12])
    return axes

# Initialize figure
fig = plt.figure(figsize=(12,8), constrained_layout=True)
ax = generate_axes(fig)

# calculate lat/lon for station
station = 'mainbasin'
lon = sta_dict[station][0]
lat = sta_dict[station][1]

# loop through different state variables
for i,vn in enumerate(vn_list_pecs):

    # plot map with station locations
    xmin = -123.29
    xmax = -122.1 # to make room for legend key
    ymin = 46.95 
    ymax = 48.5#48.93
    ax['map'].pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho), vmin=0, vmax=1.2, cmap='bone' )
    ax['map'].scatter(lon,lat,marker='D',color='crimson',s=150,zorder=5,edgecolor='pink')
    ax['map'].tick_params(axis='both', labelbottom=False, labelleft=False)
    ax['map'].set_xlim([xmin,xmax])
    ax['map'].set_ylim([ymin,ymax])
    ax['map'].set_xticks([])
    ax['map'].set_yticks([])
    pfun.dar(ax['map'])

    # add wwtp locations
    if WWTP_loc == True:
        ax['map'].scatter(moh20_lon_wwtps,moh20_lat_wwtps,color='none', edgecolors='k', linewidth=1, s=moh20_sizes_wwtps, label='WWTPs')
        ax['map'].scatter(was24_lon_wwtps,was24_lat_wwtps,color='none', edgecolors='k', linewidth=1, s=was24_sizes_wwtps)
        leg_szs = [100, 1000, 10000]
        szs = [0.05*(leg_sz) for leg_sz in leg_szs]
        l0 = plt.scatter([],[], s=szs[0], color='none', edgecolors='k', linewidth=1)
        l1 = plt.scatter([],[], s=szs[1], color='none', edgecolors='k', linewidth=1)
        l2 = plt.scatter([],[], s=szs[2], color='none', edgecolors='k', linewidth=1)
        labels = ['< 100', '1,000', '10,000']
        legend = ax['map'].legend([l0, l1, l2], labels, fontsize = 10, markerfirst=False,
            title='WWTP loading \n'+r' (kg N d$^{-1}$)',loc='upper left', labelspacing=1, borderpad=0.8)
        plt.setp(legend.get_title(),fontsize=9)

    # loop through two model runss
    for gtagex in gtagexes:

        # no loading run
        if gtagex == 'cas7_t1noDIN_x11ab':
            title = 'No-loading'
            # download .nc files
            fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
            ds = xr.open_dataset(fn)
            # get depth values
            z_rho = ds['z_rho'].transpose() # depth of u and v-velocities
            z_w   = ds['z_w'].transpose()   # depth of w-velocities
            z_min = np.min(z_w.values)
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

        elif gtagex == 'cas7_t1_x11ab':
            title = r'Loading $-$ No-loading'
            # download .nc files
            fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
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
            # save NO3 and oxygen differences to look at their molar ratios
            if vn == 'oxygen':
                DO_diff_mmol_m3 = list(ds_withLoading['t0_minus_t0noN'].values.flatten())
            elif vn == 'NO3':
                NO3_diff_mmol_m3 = list(ds_withLoading['t0_minus_t0noN'].values.flatten())
            val = ds_withLoading['t0_minus_t0noN'].transpose() * scale
            
            # get min and max for plotting
            # get min and max for plotting
            vlims = pinfo.vlims_dict[vn]
            if vn == 'NH4':
                scale = 1
            elif vn == 'NO3':
                scale = 0.2
            else:
                scale = 0.05
            vmin = vlims[1] * scale * -1
            vmax = vlims[1] * scale
            cmap = cmocean.cm.balance_r

        # create time vector
        dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
        dates_local = [pfun.get_dt_local(x) for x in dates]

        # Plot pcolormesh
        if gtagex == 'cas7_t1noDIN_x11ab':
            cs_noload = ax[vn].pcolormesh(dates_local, z_rho, val, vmin=vmin, vmax=vmax, cmap=cmap)
            cb1 = fig.colorbar(cs_noload, fraction=0.08, aspect=10, pad=0.02, ax=ax[vn], extend='both')
            cb1.outline.set_visible(False)
            cb1.ax.tick_params(labelsize=12)  
        else:
            cs_diff = ax[vn+'diff'].pcolormesh(dates_local, z_rho, val, vmin=vmin, vmax=vmax, cmap=cmap)
            cb1 = fig.colorbar(cs_diff, fraction=0.08, aspect=10, pad=0.02, ax=ax[vn+'diff'], extend='both')
            cb1.outline.set_visible(False)
            cb1.ax.tick_params(labelsize=12)  

        # formatting
        # ax[vn].text(0.03, 0.08, vn, fontweight='bold', transform=ax[vn].transAxes, fontsize=ls, color = 'k')
        if vn == 'phytoplankton':
            label = 'Phytoplankton' + units
            labelcolor = 'black'
        elif vn == 'NO3':
            label = 'NO3' + units
            labelcolor = 'white'
        elif vn == 'SdetritusN':
            label = 'Small Detritus' + units
            labelcolor = 'white'
        elif vn == 'oxygen':
            label = 'Oxygen' + units
            labelcolor = 'white'
        ax[vn].text(0.03, 0.08, label, fontweight='bold', transform=ax[vn].transAxes, fontsize=15, color = labelcolor)
        # gridlines at the month
        month_locator = mdates.MonthLocator(interval=1)
        month_fmt = mdates.DateFormatter('%b')
        ax[vn].xaxis.set_major_locator(month_locator)
        ax[vn].xaxis.set_major_formatter(month_fmt)
        ax[vn+'diff'].xaxis.set_major_locator(month_locator)
        ax[vn+'diff'].xaxis.set_major_formatter(month_fmt)
        for a in [ax[vn], ax[vn+'diff']]:
            a.grid(True, axis='x', which='major', linestyle=':', color='gray')
        # tick labels
        if vn == 'NO3':
            ax[vn].tick_params(axis='x', which='major', labelsize=12, rotation=30)
            ax[vn+'diff'].tick_params(axis='x', which='major', labelsize=12, rotation=30)
        else:
            ax[vn].tick_params(axis='x', which='both', labelbottom=False)
            ax[vn+'diff'].tick_params(axis='x', which='both', labelbottom=False)
        ax[vn].set_ylabel('z [m]', fontsize=12)
        ax[vn].tick_params(axis='y', which='both', labelsize=12)
        ax[vn+'diff'].tick_params(axis='y', which='both', labelleft=False)
        fig.subplots_adjust(top=0.8, left=0.1, right=0.9)

##########################################################
##               Bottom DO time series                  ##
##########################################################

# plt.close('all')

nwin = 10

# Initialize figure
fig,ax = plt.subplots(1,1,figsize=(6,3))

# create time vector
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

# colors
colors = ['black', 'gray', 'black', 'turquoise']
# linestyles
linestyles = ['-', '--', '-', '-']
# linewidths
linewidths = [3, 1.5, 1.3, 3]

for s,station in enumerate(list(sta_dict.keys())):
# for s,station in enumerate(['mainbasin']):
    # no-loading run only
    gtagex = 'cas7_t1noDIN_x11ab'
    # look at DO
    vn = 'oxygen'
    # download .nc files
    fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # get scale and units
    scale =  pinfo.fac_dict[vn]
    units = pinfo.units_dict[vn]
    vlims = pinfo.vlims_dict[vn]
    vmin = vlims[0]
    vmax = vlims[1]
    # get dataset
    val = ds[vn][:,0] * scale # dimension = t,0 where 0 is the bottom

    # plot
    ax.plot(dates_local, zfun.lowpass(val.values,n=nwin), color = colors[s],
             linestyle=linestyles[s], linewidth=linewidths[s], label=station)
    # ax.plot(dates_local,val, color = colors[s],alpha=0.5,
    #          linestyle=linestyles[s], linewidth=linewidths[s], label=station)

# format plot
ax.legend(loc='upper right', fontsize=12, ncol=4)
ax.set_ylim([vmin,vmax])
ax.set_xlim([dates_local[0],dates_local[-1]])
ax.tick_params(axis='both', which='major', labelsize=12, rotation=30)
ax.grid(True, linestyle=':', color='silver')
# gridlines at the month
month_locator = mdates.MonthLocator(interval=1)
month_fmt = mdates.DateFormatter('%b')
ax.xaxis.set_major_locator(month_locator)
ax.xaxis.set_major_formatter(month_fmt)
ax.xaxis.set_major_locator(month_locator)
ax.xaxis.set_major_formatter(month_fmt)
plt.tight_layout()


##########################################################
##                NO3 and DO diff ratio                 ##
##########################################################

# plt.close('all')
# Initialize figure
fig,ax = plt.subplots(figsize=(5,5))

# calculate lat/lon for station
station = 'mainbasin'
lon = sta_dict[station][0]
lat = sta_dict[station][1]

# loop through two model runss
for gtagex in gtagexes:

    # no loading run
    if gtagex == 'cas7_t1noDIN_x11ab':
        # download .nc files
        fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds = xr.open_dataset(fn)
        # get NO3 and DO
        DO_noloading_mmol_m3  = list(ds['oxygen'].values.flatten())
        NO3_noloading_mmol_m3 = list(ds['NO3'].values.flatten())


    elif gtagex == 'cas7_t1_x11ab':
        # download .nc files
        fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds_withLoading = xr.open_dataset(fn)
        # get NO3 and DO
        DO_loading_mmol_m3  = list(ds_withLoading['oxygen'].values.flatten())
        NO3_loading_mmol_m3 = list(ds_withLoading['NO3'].values.flatten())

# calculate diffs
DO_diff_mmol_m3   = np.array(DO_loading_mmol_m3) - np.array(DO_noloading_mmol_m3)
NO3_diff_mmol_m3  = np.array(NO3_loading_mmol_m3) - np.array(NO3_noloading_mmol_m3)

# plot
# ax.hist2d(NO3_diff_mmol_m3, DO_diff_mmol_m3,bins=100,cmap=cmc.tokyo_r)
ax.scatter(NO3_diff_mmol_m3, DO_diff_mmol_m3,alpha=0.2,edgecolor='None',s=2,color='crimson')
# ax.scatter(NO3_diff_mmol_m3, DO_diff_mmol_m3,alpha=0.03,s=20,color='cadetblue',
#            edgecolor='none',zorder=5)

# add slope of 2 line
x = np.linspace(2,3,10)
y = -2*x + 1
ax.plot(x,y,color='black',linestyle='--',linewidth=2,label='slope = -2')

# format plot
ax.grid(True, axis='both', which='major', linestyle=':', color='silver')
ax.tick_params(axis='both', which='major', labelsize=12)
ax.set_xlabel(r'$\Delta NO3$ [mmol/m3]', fontsize=12)
ax.set_ylabel(r'$\Delta DO$ [mmol/m3]', fontsize=12)
# ax.set_aspect('equal', adjustable='box')
plt.tight_layout()

# ax.set_ylim([-8,0])