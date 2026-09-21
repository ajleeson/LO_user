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

plt.close('all')

gtagex = 'cas7_t0_x4b'
jobname = 'terminlet_EcolCTD'
startdate = '2015.01.01'
enddate = '2020.12.31'
years = ['2015','2016','2017','2018','2019','2020']

vn = 'DO (uM)'

jobname = 'terminlet_EcolCTD'
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
fig, ax = plt.subplots(6,1,figsize = (10,9),sharey='row',sharex='col')
axes = ax.ravel()

# Initialize empty dataframe with specific dtypes
new_obs_df = pd.DataFrame({
    'time': pd.Series(dtype='datetime64[ns]'),
    'inlet': pd.Series(dtype='str'),
    'layer': pd.Series(dtype='str'),
    'DO_mgL': pd.Series(dtype='float64')
})
new_mod_df = pd.DataFrame({
    'time': pd.Series(dtype='datetime64[ns]'),
    'inlet': pd.Series(dtype='str'),
    'layer': pd.Series(dtype='str'),
    'DO_mgL': pd.Series(dtype='float64')
})

# add data
for y,year in enumerate(years):
    # get observations
    in_dir = Ldir['parent'] / 'LO_output' / 'obsmod'
    # in_fn = in_dir / ('multi_ctd_' + year + '.p')
    in_fn = in_dir / ('combined_bottle_' + year + '_cas7_t1_x11ab.p')
    df_dict = pickle.load(open(in_fn, 'rb'))
    # only look at ecology stations
    source = 'ecology_nc'
    df_obs = df_dict['obs'].loc[df_dict['obs'].source==source,:]

    # loop through stations
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
        lon = sta_dict[selected_stns[station]][0]
        lat = sta_dict[selected_stns[station]][1]

        # get observational information
        df_ob_stn = df_obs.loc[df_obs.name==selected_stns[station],:]
        # get depth and time of observations
        z = df_ob_stn['z']
        time = [pd.Timestamp(x) for x in df_ob_stn['time']]


        # get station depth
        fn2 = '../../../../LO_output/extract/' + 'cas7_t1_x11ab' + '/moor/' + jobname + '/' + selected_stns[station] + '_' + startdate + '_' + enddate + '.nc'
        ds_moor2 = xr.open_dataset(fn2)
        h = ds_moor2.h.values
        # get depth values
        z_w   = ds_moor2['z_w'].transpose()   # depth of w-velocities
        z_min = np.min(z_w.values)
        
        # get all variables and convert to correct units
        if vn == 'DO (uM)': # already in uM
            ds_moor2['DO (mg/L)'] = ds_moor2['oxygen'] * 32/1000
            # val = ds_moor['DO (mg/L)'].transpose()

        # if stn == 4:
        #     print(h)
            # print(df_ob_stn.to_string())

        # if deeper than 10 m, split into top 5 m and bottom 5 m layer
        d = 10

        # model

        surf_z_ds2 = ds_moor2.where((ds_moor2.z_w >= -(d/2)))
        bott_z_ds2 = ds_moor2.where((ds_moor2.z_w <= (-1*h) + (d/2)))
        surf_mod_ds2 = ds_moor2.where((ds_moor2.z_rho >= -(d/2)))
        bott_mod_ds2 = ds_moor2.where((ds_moor2.z_rho <= (-1*h) + (d/2)))

        # get depth of each vertical layer
        z_thick_surf2 = np.diff(surf_z_ds2['z_w'],axis=1) # 365,30
        z_thick_bott2 = np.diff(bott_z_ds2['z_w'],axis=1) # 365,30
        # get model data 
        if vn == 'DO (uM)':
            surf_mod2 = surf_mod_ds2['DO (mg/L)'].values
            bott_mod2 = bott_mod_ds2['DO (mg/L)'].values
        else:
            surf_mod2 = surf_mod_ds2[vn].values
            bott_mod2 = bott_mod_ds2[vn].values
        # check if thickness array is all nan
        # (this occurs if the first z_w is already greater than the threshold, so we don't have two z_w values to diff)
        # in which case, average value is just the single value
        if np.isnan(z_thick_surf2).all():
            surf_mod_avg2 = np.nansum(surf_mod2,axis=1)
        # take weighted average given depth
        # first, multiply by thickness of each layer and sum in z, then divide by layer thickness
        else:
            surf_mod_avg2 = np.nansum(surf_mod2 * z_thick_surf2, axis=1)/np.nansum(z_thick_surf2,axis=1)

        # repeat check for bottom
        if np.isnan(z_thick_bott2).all():
            bott_mod_avg2 = np.nansum(bott_mod2,axis=1)
        else:
            bott_mod_avg2 = np.nansum(bott_mod2 * z_thick_bott2, axis=1)/np.nansum(z_thick_bott2,axis=1)



        # plot model output
        if y == 0:
            axis.plot(dates_local, surf_mod_avg2, color='lightseagreen', linewidth=2.5, alpha=0.3, zorder=5)
            axis.plot(dates_local, bott_mod_avg2, color='navy', linewidth=2.5, alpha=0.3, zorder=3,label='model')

        # add surface values to dataframe
        temp_df = pd.DataFrame({
            'time': pd.to_datetime(dates_local).values,
            'inlet': stn_name + ' ('+selected_stns[station]+')',  
            'layer': 'surface',
            'DO_mgL': surf_mod_avg2})
        new_mod_df = pd.concat([new_mod_df,temp_df], ignore_index=True)
        # add bottom values to dataframe
        temp_df = pd.DataFrame({
            'time': pd.to_datetime(dates_local).values,
            'inlet': stn_name + ' ('+selected_stns[station]+')',  
            'layer': 'bottom',
            'DO_mgL': bott_mod_avg2})
        new_mod_df = pd.concat([new_mod_df,temp_df], ignore_index=True)

        if stn == 0:
            pass
            # print(dates_local)


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
        # unique_time = [pd.Timestamp(x) for x in df_ob_stn['time'].unique()] # one point per timestampe
        axis.scatter(surf_obs_avg.index, surf_obs_avg, color='lightseagreen', s=20,zorder=10)

        # add surface values to dataframe
        temp_df = pd.DataFrame({
            # 'time': unique_time,
            'time': surf_obs_avg.index,
            'inlet': stn_name + ' ('+selected_stns[station]+')',  
            'layer': 'surface',
            'DO_mgL': surf_obs_avg})
        new_obs_df = pd.concat([new_obs_df,temp_df], ignore_index=True)


        # add value for Carr Inlet
        if stn == 4:
            bott_obs_avg.loc[pd.Timestamp("2017-08-10 20:18:08")] = np.nan
            bott_obs_avg = bott_obs_avg.sort_index()    
            # print(unique_time)
            # print(bott_obs_avg)
        axis.scatter(bott_obs_avg.index, bott_obs_avg, color='navy', s=20,zorder=10, label='obs')

        # add bottom values to dataframe
        temp_df = pd.DataFrame({
            # 'time': unique_time,
            'time': bott_obs_avg.index,
            'inlet': stn_name + ' ('+selected_stns[station]+')',  
            'layer': 'bottom',
            'DO_mgL': bott_obs_avg})
        new_obs_df = pd.concat([new_obs_df,temp_df], ignore_index=True)

        # print(bott_obs_avg)
        # if stn == 4:
        #     axis.scatter(unique_time[1::], bott_obs_avg, color='black', s=20,zorder=10, label='obs')
        # else:
        #     axis.scatter(unique_time, bott_obs_avg, color='black', s=20,zorder=10, label='obs')


        # label
        if y == 0 and stn == 0:
            axis.text(0.7, 0.9, 'surface {} m'.format(str(round(d/2))),
                verticalalignment='top', horizontalalignment='left',
                transform=axis.transAxes, fontsize=12, color = 'lightseagreen',fontweight='bold')
            axis.text(0.7, 0.8, 'bottom {} m'.format(str(round(d/2))),
                verticalalignment='top', horizontalalignment='left',
                transform=axis.transAxes, fontsize=12, color = 'navy',fontweight='bold')
            axis.legend(loc='lower right',fontsize=11, frameon=False, handletextpad=0.1,
                        handlelength=1, markerfirst = False,ncol=3)


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
            axis.set_ylabel('hi',fontsize=14)


# format grids and add titles
for i,axis in enumerate(axes):
    axis.set_xlim([dates_local[0],dates_local[-1]])
    axis.grid(True,color='gainsboro',linewidth=0.5,linestyle='--',axis='both')
    axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    axis.tick_params(axis='both', labelrotation=30,  labelsize=12)
    axis.set_yticks(np.arange(0, 18, 2))
    # axis.tick_params(axis='both', labelsize=12)
    axis.text(0.03, 0.95, letter[i],
            verticalalignment='top', horizontalalignment='left',
            transform=axis.transAxes, fontsize=14, color = 'k',fontweight='bold')

plt.subplots_adjust(wspace=-0.2)      
plt.tight_layout()

print(new_mod_df)
new_obs_df.to_csv('obs_DO_timeseries.csv', index=False)
new_mod_df.to_csv('mod_DO_timeseries.csv', index=False)