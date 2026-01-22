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
    'elliot': 'ELB015',
    'lynchcove': 'HCB007',
    'lynchcove2': 'HCB004',
    'commencement': 'CMB003',
    # 'hammersley': 'OAK004',
    # 'totten': 'TOT002',
    # 'budd': 'BUD005',
    'carr': 'CRR001',
    'sinclair': 'SIN001',
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

# dictionary of shallow inlets
shallow_stns = {
    'hammersley': 'OAK004',
    'budd': 'BUD005',
}

# dictionary of deep inlets
deep_stns = {
    'elliot': 'ELB015',
    'lynchcove': 'HCB007',
    'lynchcove2': 'HCB004',
    'commencement': 'CMB003',
    'carr': 'CRR001',
    'sinclair': 'SIN001',
}

# get Puget Sound stations that aren't in terminal inlets
PS_stns = ['ADM002','PTH005','ADM001',
           'ADM003','PSB003','EAP001']


stations = ['ADM002','PTH005','ADM001',
           'ADM003','PSB003','EAP001',
           'ELB015','HCB007','HCB004',
           'CMB003','CRR001','SIN001']
# only look at inlets of interest
df_obs = df_obs.loc[df_obs.name.isin(stations),:]
df_model = df_model.loc[df_model.name.isin(stations),:]

# only keep columns of interest
df_obs = df_obs[['time', 'lat', 'lon', 'name', 'z','DO (uM)']]
df_model = df_model[['time', 'lat', 'lon', 'name', 'z','DO (uM)']]

# save as csv files
df_obs.to_csv('CTD_cast_observations.csv', index=False)
df_model.to_csv('CTD_cast_model.csv', index=False)

#####################################################################################################################################

##########################################################
##       Plotting all Ecology CTD locations             ##
##########################################################

letter = ['(a)','(b)']

# initialize figure
plt.figure(figsize=(9.5, 6))
ax_map = plt.subplot(1,2,1)
ax_do = plt.subplot(1,2,2)

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
ax_map.add_patch(Rectangle((lon_low, lat_low), lon_high-lon_low,lat_high-lat_low, facecolor='white'))#, facecolor='#EEEEEE'))
ax_map.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

# format
ax_map.tick_params(axis='x', labelrotation=30)
ax_map.set_ylabel('Latitude',fontsize=14)
ax_map.set_xlabel('Longitude',fontsize=14)
ax_map.set_xlim(-123.29, -122.1) # Puget Sound
ax_map.set_ylim(lat_low, lat_high) # Puget Sound
ax_map.tick_params(axis='both', labelsize=12)
pfun.dar(ax_map)
    

# add terminal inlets cast locations
for sta in terminl_stns:
    sta_lon = sta_dict[sta][0]
    sta_lat = sta_dict[sta][1]
    color = stn_color[sta]
    ecol = '\n('+ terminl_stns[sta] +')'
    ax_map.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=8,
    color='deeppink',markeredgecolor='white',markeredgewidth=0.8)
    # add inlet label
    lon_off = 0.14
    lat_off = -0.03
    ha = 'center'
    ax_map.text(sta_lon+lon_off,sta_lat+lat_off,terminl_stns[sta],va='center',
                ha=ha,color='deeppink',fontsize=9, fontweight='bold')
    
# add the remaining cast locations
for sta in PS_stns:
    # sta_lon = sta_dict[sta][0]
    # sta_lat = sta_dict[sta][1]
    sta_lat = df_obs.loc[df_obs['name'] == sta, 'lat'].values[0]
    sta_lon = df_obs.loc[df_obs['name'] == sta, 'lon'].values[0]
    ax_map.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=8,
    color='navy',markeredgecolor='white',markeredgewidth=0.8)
    # add inlet label
    lon_off = 0.14
    lat_off = 0.03
    ha = 'center'
    ax_map.text(sta_lon+lon_off,sta_lat+lat_off,sta,va='center',
                ha=ha,color='navy',fontsize=9, fontweight='bold')

# types of stations
t1 = ax_map.text(0.03, 0.98, 'Terminal Inlet stations',
            verticalalignment='top', horizontalalignment='left',fontweight='bold',
            transform=ax_map.transAxes, fontsize=10, color = 'deeppink')
t2 = ax_map.text(0.03, 0.94, 'Main Basin stations',
            verticalalignment='top', horizontalalignment='left',fontweight='bold',
            transform=ax_map.transAxes, fontsize=10, color = 'navy')
t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

# ax_map.set_title('(a) Inlet obs comparison locations',fontsize=12,loc='left')
ax_map.set_title('(a) Locations',fontsize=14,loc='left', fontweight='bold')

##########################################################
##         Puget Sound Property-property plots          ##
##########################################################

axes = [ax_do]
vars = [ r'DO [mg L$^{-1}$]']
vns = ['DO (uM)']

# format grid
for i,ax in enumerate(axes):
    ax.tick_params(axis='x', labelrotation=30)
    ax.grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax.tick_params(axis='both', labelsize=12)
    ax.set_title(letter[i+1] + ' ' + vars[i],loc='left',fontsize=14, fontweight='bold')

# plot property-property plots for terminal inlets
for i,vn in enumerate(vns):
    # initialize arrays to save values
    terminl_obs = np.array([])
    terminl_mod = np.array([])
    for stn,station in enumerate(terminl_stns):
        # get observational and model information
        df_ob_stn = df_obs.loc[df_obs.name==terminl_stns[station],:]
        df_mo_stn = df_model.loc[df_model.name==terminl_stns[station],:]
        if vn == 'DO (uM)':
            obsvals = df_ob_stn[vn].values * 32/1000
            modvals = df_mo_stn[vn].values * 32/1000
        else:
            obsvals = df_ob_stn[vn].values 
            modvals = df_mo_stn[vn].values
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
        axes[i].scatter(obsvals,modvals, color='deeppink', alpha=0.3, s=10, zorder=10, edgecolor='none')
        axes[i].plot([0,maxval*1.01], [0,maxval*1.01], color='grey', linestyle='-', zorder=5)
        axes[i].set_ylim([0,maxval])
        axes[i].set_xlim([0,maxval])
        axes[i].set_aspect('equal', adjustable='box')
        axes[i].xaxis.set_major_locator(plt.MaxNLocator(5))
        axes[i].yaxis.set_major_locator(plt.MaxNLocator(5))
        # add data to arrays
        terminl_obs = np.concatenate((terminl_obs,obsvals))
        terminl_mod = np.concatenate((terminl_mod,modvals))
    # # add bias and rmse label
    # calculate bias and rmse
    bias = np.nanmean(terminl_mod-terminl_obs)
    rmse = np.sqrt(np.nanmean((terminl_mod-terminl_obs)**2))
    nse = 1 - (
        np.nansum( (terminl_mod-terminl_obs)**2 )/
        np.nansum( (terminl_obs- np.nanmean(terminl_obs) )**2 )
        )
    t01 = axes[i].text(0.17, 0.9, 'Bias', fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        transform=axes[i].transAxes, fontsize=12, color = 'black')
    t02 = axes[i].text(0.38, 0.9, 'RMSE', fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        transform=axes[i].transAxes, fontsize=12, color = 'black')
    t1 = axes[i].text(0.17, 0.85, '{:.2f}'.format(round(bias,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=12, color = 'deeppink')
    t2 = axes[i].text(0.37, 0.85, '{:.2f}'.format(round(rmse,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=12, color = 'deeppink')
    # t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    # t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    # t3.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))

    axes[i].set_ylabel('Modeled',fontsize=14)
    axes[i].set_xlabel('Observed',fontsize=14)


    # initialize arrays to save values
    PSstns_obs = np.array([])
    PSstns_mod = np.array([])
    for stn,station in enumerate(PS_stns):
        # get observational and model information
        df_ob_stn = df_obs.loc[df_obs.name==station,:]
        df_mo_stn = df_model.loc[df_model.name==station,:]
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
        minval = 0
        maxval = 17
        # plot data and y = x line
        axes[i].scatter(obsvals,modvals, color='navy', alpha=0.1, s=10, zorder=10, edgecolor='none')
        axes[i].plot([minval,maxval*1.01], [minval,maxval*1.01], color='grey', linestyle='-', zorder=5)
        axes[i].set_ylim([minval,maxval])
        axes[i].set_xlim([minval,maxval])
        axes[i].set_aspect('equal', adjustable='box')
        axes[i].xaxis.set_major_locator(plt.MaxNLocator(5))
        axes[i].yaxis.set_major_locator(plt.MaxNLocator(5))
    # add data to arrays
        PSstns_obs = np.concatenate((PSstns_obs,obsvals))
        PSstns_mod = np.concatenate((PSstns_mod,modvals))
    # # add bias and rmse label
    # calculate bias and rmse and nash-sutcliffe model efficiency index
    bias = np.nanmean(PSstns_mod-PSstns_obs)
    rmse = np.sqrt(np.nanmean((PSstns_mod-PSstns_obs)**2))
    
    t1 = axes[i].text(0.17, 0.81, '{:.2f}'.format(round(bias,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=12, color = 'navy')
    t2 = axes[i].text(0.37, 0.81, '{:.2f}'.format(round(rmse,2)),
        verticalalignment='top', horizontalalignment='right',fontweight='bold',
        transform=axes[i].transAxes, fontsize=12, color = 'navy')
    # t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    # t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))


# plt.suptitle('Model vs. Observations property-property plots (Puget Sound 2017)', fontsize=15, x=0.65)
plt.subplots_adjust(wspace=-0.2)      
plt.tight_layout()


############################################################################################################################
## Selected time series

# create time vector
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
# dates_local = [pfun.get_dt_local(x) for x in dates]
dates_local = dates

# dictionary of selected stations in terminal inlets for time series
selected_stns = {
    # 'elliot': 'ELB015',
    # 'lynchcove': 'HCB007',
    # # 'hammersley': 'OAK004',
    'elliot': 'ELB015',
    'commencement': 'CMB003',
    'lynchcove': 'HCB007',
    'lynchcove2': 'HCB004',
    # 'hammersley': 'OAK004',
    # 'totten': 'TOT002',
    # 'budd': 'BUD005',
    'carr': 'CRR001',
    'sinclair': 'SIN001',
}

letter = ['(a)','(b)','(c)',
          '(d)','(e)','(f)',
          '(g)','(h)']

# initialize figure
fig, ax = plt.subplots(3,2,figsize = (10,9),sharey='row',sharex='col')
axes = ax.ravel()

# add data
for i,vn in enumerate(vns):
    for stn, station in enumerate(selected_stns):
        # get axis
        axis = axes[stn]
        # add a title
        if station == 'elliot':
            stn_name = 'elliott'
        else:
            stn_name = station
        axis.set_title(stn_name + ' ('+selected_stns[station]+')',fontsize=14, fontweight='bold')
        

        # calculate lat/lon for station
        lon = sta_dict[station][0]
        lat = sta_dict[station][1]

        # get observational information
        df_ob_stn = df_obs.loc[df_obs.name==selected_stns[station],:]
        df_mo_stn = df_model.loc[df_model.name==selected_stns[station],:]
        # get depth and time of observations
        z = df_ob_stn['z']
        time = [pd.Timestamp(x) for x in df_ob_stn['time']]


        # get station depth
        fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds_moor = xr.open_dataset(fn)
        h = ds_moor.h.values
        # get depth values
        z_w   = ds_moor['z_w'].transpose()   # depth of w-velocities
        z_min = np.min(z_w.values)
        
        # get all variables and convert to correct units
        if vn == 'DO (uM)': # already in uM
            ds_moor['DO (mg/L)'] = ds_moor['oxygen'] * 32/1000
            val = ds_moor['DO (mg/L)'].transpose()

        # if stn == 4:
        #     print(h)
            # print(df_ob_stn.to_string())

        # if deeper than 10 m, split into top 5 m and bottom 5 m layer
        d = 10
        if h > d:

            # model
            surf_z_ds = ds_moor.where((ds_moor.z_w >= -(d/2)))
            bott_z_ds = ds_moor.where((ds_moor.z_w <= (-1*h) + (d/2)))
            surf_mod_ds = ds_moor.where((ds_moor.z_rho >= -(d/2)))
            bott_mod_ds = ds_moor.where((ds_moor.z_rho <= (-1*h) + (d/2)))
            # get depth of each vertical layer
            z_thick_surf = np.diff(surf_z_ds['z_w'],axis=1) # 365,30
            z_thick_bott = np.diff(bott_z_ds['z_w'],axis=1) # 365,30
            # get model data 
            if vn == 'DO (uM)':
                surf_mod = surf_mod_ds['DO (mg/L)'].values
                bott_mod = bott_mod_ds['DO (mg/L)'].values
            else:
                surf_mod = surf_mod_ds[vn].values
                bott_mod = bott_mod_ds[vn].values
            # check if thickness array is all nan
            # (this occurs if the first z_w is already greater than the threshold, so we don't have two z_w values to diff)
            # in which case, average value is just the single value
            if np.isnan(z_thick_surf).all():
                surf_mod_avg = np.nansum(surf_mod,axis=1)
            # take weighted average given depth
            # first, multiply by thickness of each layer and sum in z, then divide by layer thickness
            else:
                surf_mod_avg = np.nansum(surf_mod * z_thick_surf, axis=1)/np.nansum(z_thick_surf,axis=1)
            # repeat check for surface
            if np.isnan(z_thick_bott).all():
                bott_mod_avg = np.nansum(bott_mod,axis=1)
            else:
                bott_mod_avg = np.nansum(bott_mod * z_thick_bott, axis=1)/np.nansum(z_thick_bott,axis=1)
            # plot model output
            axis.plot(dates_local, surf_mod_avg, color='lightseagreen', linewidth=2, alpha=0.5, zorder=5)
            axis.plot(dates_local, bott_mod_avg, color='black', linewidth=2, alpha=0.5, zorder=5,label='model')
            if stn == 0:
                print(dates_local)

            # observation
            surf_obs_df = df_ob_stn.where(z >= -(d/2))
            bott_obs_df = df_ob_stn.where(z <= (-1*h) + (d/2))
            # get average value
            if vn == 'DO (uM)':
                surf_obs_avg = surf_obs_df.groupby('time')[vn].mean() * 32/1000
                bott_obs_avg = bott_obs_df.groupby('time')[vn].mean() * 32/1000
            else:
                surf_obs_avg = surf_obs_df.groupby('time')[vn].mean()
                bott_obs_avg = bott_obs_df.groupby('time')[vn].mean()
            # plot observations
            unique_time = [pd.Timestamp(x) for x in df_ob_stn['time'].unique()] # one point per timestampe
            axis.scatter(unique_time, surf_obs_avg, color='lightseagreen', s=20,zorder=10)
            # add value for Carr Inlet
            if stn == 4:
                bott_obs_avg.loc[pd.Timestamp("2017-08-10 20:18:08")] = np.nan
                bott_obs_avg = bott_obs_avg.sort_index()    
                # print(unique_time)
                # print(bott_obs_avg)
            axis.scatter(unique_time, bott_obs_avg, color='navy', s=20,zorder=10, label='obs')
            # if stn == 4:
            #     axis.scatter(unique_time[1::], bott_obs_avg, color='black', s=20,zorder=10, label='obs')
            # else:
            #     axis.scatter(unique_time, bott_obs_avg, color='black', s=20,zorder=10, label='obs')


            # label
            if i == 0 and stn == 0:
                axis.text(0.7, 0.9, 'surface {} m'.format(str(round(d/2))),
                    verticalalignment='top', horizontalalignment='left',
                    transform=axis.transAxes, fontsize=12, color = 'lightseagreen',fontweight='bold')
                axis.text(0.7, 0.8, 'bottom {} m'.format(str(round(d/2))),
                    verticalalignment='top', horizontalalignment='left',
                    transform=axis.transAxes, fontsize=12, color = 'black',fontweight='bold')
                axis.legend(loc='lower right',fontsize=11, frameon=False, handletextpad=0.1,
                            handlelength=1, markerfirst = False)
                

        # otherwise just one layer
        else:

            # get depth of each vertical layer
            z_thick = np.diff(z_w,axis=0).transpose() # 365,30

            # model
            mod_ds = ds_moor
            # get model data 
            if vn == 'DO (uM)':
                mod = mod_ds['DO (mg/L)']
            else:
                mod = mod_ds[vn]
            # take weighted average given depth
            # first, multiply by thickness of each layer and sum in z, then divide by water column depth
            mod_avg = np.nansum(mod * z_thick, axis=1)/np.nansum(z_thick,axis=1)
            # plot model output
            axis.plot(dates_local, mod_avg, color='mediumorchid', linewidth=2, alpha=0.5, zorder=5, label='model')

            # observation
            obs_df = df_ob_stn
            # get average value
            if vn == 'DO (uM)':
                obs_avg = obs_df.groupby('time')[vn].mean() * 32/1000
            else:
                obs_avg = obs_df.groupby('time')[vn].mean()
            # plot observations
            unique_time = [pd.Timestamp(x) for x in df_ob_stn['time'].unique()]
            # axis.scatter(time, obs_df[var], color='k',s=3)
            axis.scatter(unique_time, obs_avg, color='mediumorchid', s=20, zorder=10, label='obs')

        # set y-axis title and limits
        # if i == 0:
        #     axis.set_ylim([0,25])
        # if i == 1:
        #     axis.set_ylim([0,36])
        # if i == 2:
        #     axis.set_ylim([0,25])
        # if i == 3:
        axis.set_ylim([0,16])
        if np.mod(stn,2) == 0:
            axis.set_ylabel(vars[i],fontsize=14)


# format grids and add titles
for i,axis in enumerate(axes):
    axis.set_xlim([dates_local[0],dates_local[-1]])
    axis.grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    axis.tick_params(axis='both', labelrotation=30,  labelsize=12)
    axis.set_yticks(np.arange(0, 18, 2))
    # axis.tick_params(axis='both', labelsize=12)
    axis.text(0.03, 0.95, letter[i],
            verticalalignment='top', horizontalalignment='left',
            transform=axis.transAxes, fontsize=14, color = 'k',fontweight='bold')

plt.subplots_adjust(wspace=-0.2)      
plt.tight_layout()