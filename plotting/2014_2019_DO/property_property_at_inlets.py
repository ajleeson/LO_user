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
startdate = '2014.01.01'
enddate = '2014.12.31'
year = '2014' # for making a date label

vn_list = ['temp','salt','phytoplankton','oxygen']
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
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures' / 'twentyone' / 'inlet_obs_comparison'
Lfun.make_dir(out_dir)

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

##########################################################
##                      Plotting                        ##
##########################################################

letter = ['(a)','(b)','(c)','(d)']

# Loop through all of the mooring stations
for i,station in enumerate(ecol_stn): # enumerate(['commencement']): #
    plt.close('all')

    print(station)

    # calculate lat/lon for station
    lon = sta_dict[station][0]
    lat = sta_dict[station][1]

    # get observational information
    df_ob_stn = df_obs.loc[df_obs.name==ecol_stn[station],:]
    df_mo_stn = df_model.loc[df_model.name==ecol_stn[station],:]
    # get depth and time of observations
    z = df_ob_stn['z']
    time = [pd.Timestamp(x) for x in df_ob_stn['time']]

    # Initialize Figure
    width = 11
    fig, ax = plt.subplots(rows,3,figsize = (width,(4/5)*width), gridspec_kw={'width_ratios': [2, 2, 1]}) # adjust subplot sizes
    fig.suptitle(station + ' ' + year, fontsize = 18)
    
    # loop through different state variables
    for i,vn in enumerate(vn_list):

        # PLOT MODEL DEPTH VS TIME PCOLORMESH -------------------------------------------------------------
        # download .nc files
        col = 0
        axis = ax[i,col]
        fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds = xr.open_dataset(fn)
        h = ds.h.values
        # get depth values
        z_rho = ds['z_rho'].transpose() # depth of u and v-velocities
        z_rho = z_rho.values
        z_w   = ds['z_w'].transpose()   # depth of w-velocities
        z_min = np.min(z_w.values)

        # set upper vertical limit (to see top of sea surface)
        z_max = -0.3*z_min
        if z_max < 3:
            z_max = 3

        # get scale and units
        scale =  pinfo.fac_dict[vn]
        vlims = pinfo.vlims_dict[vn]
        vmin = vlims[0] / scale
        vmax = vlims[1] /scale
        cmap = pinfo.cmap_dict[vn]

        # convert to correct units
        if vn == 'oxygen': # already in uM
            ds['DO (uM)'] = ds[vn]
            val = ds['DO (uM)'].transpose()
        if vn == 'phytoplankton': # need to multiply by 2.5
            ds['Chl (mg m-3)'] = ds[vn] * scale
            val = ds['Chl (mg m-3)'].transpose()
            vmin = vmin * scale
            vmax = vmax * scale
        if vn == 'temp': # convert to conservative temprature
            ds['SA'] = gsw.conversions.SA_from_SP(ds['salt'], ds['z_rho'], lon, lat)
            ds['CT'] = gsw.conversions.CT_from_pt(ds['SA'], ds['temp'])
            val = ds['CT'].transpose()
        if vn == 'salt': # convert to absolute salinity
            ds['SA'] = gsw.conversions.SA_from_SP(ds['salt'], ds['z_rho'], lon, lat)
            val = ds['SA'].transpose()
            # print(val.values[0,:])

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

        if vn == 'temp':
            var = 'CT'
            name = r'Cons Temp ($\degree C$)'
        if vn == 'salt':
            var = 'SA'
            name = r'Abs Salinity ($g\ kg^{-1}$)'
        if vn == 'phytoplankton':
            var = 'Chl (mg m-3)'
            name = r'Chl ($mg\ m^{-3}$)'
        if vn == 'oxygen':
            var = 'DO (uM)'
            name = r'DO ($\mu M$)'

        axis.text(0.02, 0.96, letter[i] + ' ' + name, fontweight='bold',
                verticalalignment='top', horizontalalignment='left',
                transform=axis.transAxes, fontsize=13, color = font_color)
        axis.set_ylim((z_min,z_max))
        axis.set_xlim((dates_local[0],dates_local[-1]))

        # title
        if i == 0:
            axis.set_title('Model')

        # add observation locations
        # get current station
        if i == 0:
            ax[i,0].plot(time,z,linestyle='none',marker='o', markersize=3,
                    markeredgecolor='none',markerfacecolor='black',label='obs')
            ax[i,0].legend(loc='upper right',fontsize=11, frameon=False, handletextpad=0.1, handlelength=1)


        # PLOT PROPERTY PROPERTY SCATTER -------------------------------------------------------------
        col = 2
        axis = ax[i,col]
        # get obs and model values
        obsvals = df_ob_stn[var].values
        modvals = df_mo_stn[var].values
        # calculate bias and rmse
        bias = np.nanmean(modvals-obsvals)
        rmse = np.sqrt(np.nanmean((modvals-obsvals)**2))
        # get min and max limits
        minval = np.nanmin([np.nanmin(obsvals),np.nanmin(modvals)])
        maxval = np.nanmax([np.nanmax(obsvals),np.nanmax(modvals)])
        # plot data and y = x line
        axis.scatter(obsvals,modvals, color='mediumorchid', alpha=0.6, s=10, zorder=10)
                     #c=z, cmap=cmocean.tools.crop_by_percent(cmocean.cm.thermal, 15, which='both'), 
        axis.plot([0,maxval*1.2], [0,maxval*1.2], 'k-', zorder=5)
        axis.set_ylim([0,maxval*1.2])
        axis.set_xlim([0,maxval*1.2])
        axis.set_aspect('equal')
        if i ==0:
            axis.set_title('Modeled vs. Observed')
        # add bias and rmse label
        t1 = axis.text(0.05, 0.95, 'bias: {}'.format(str(round(bias,2))),
            verticalalignment='top', horizontalalignment='left',
            transform=axis.transAxes, fontsize=ls, color = 'k')
        t2 = axis.text(0.05, 0.8, 'rmse: {}'.format(str(round(rmse,2))),
            verticalalignment='top', horizontalalignment='left',
            transform=axis.transAxes, fontsize=ls, color = 'k')
        t1.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
        t2.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
        


        # PLOT TIME SERIES -------------------------------------------------------------
        col = 1
        axis = ax[i,col]

        # if deeper than 10 m, split into top 5 m and bottom 5 m layer
        d = 10
        if h > d:

            # model
            surf_z_ds = ds.where((ds.z_w >= -(d/2)))
            bott_z_ds = ds.where((ds.z_w <= (-1*h) + (d/2)))
            surf_mod_ds = ds.where((ds.z_rho >= -(d/2)))
            bott_mod_ds = ds.where((ds.z_rho <= (-1*h) + (d/2)))
            # get depth of each vertical layer
            z_thick_surf = np.diff(surf_z_ds['z_w'],axis=1) # 365,30
            z_thick_bott = np.diff(bott_z_ds['z_w'],axis=1) # 365,30
            # get model data 
            surf_mod = surf_mod_ds[var].values
            bott_mod = bott_mod_ds[var].values
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
                print(bott_mod)
                bott_mod_avg = np.nansum(bott_mod,axis=1)
            else:
                bott_mod_avg = np.nansum(bott_mod * z_thick_bott, axis=1)/np.nansum(z_thick_bott,axis=1)
            # plot model output
            axis.plot(dates_local, surf_mod_avg, color='deeppink', linewidth=2, alpha=0.3, zorder=5)
            axis.plot(dates_local, bott_mod_avg, color='royalblue', linewidth=2, alpha=0.3, zorder=5,label='model')

            # observation
            surf_obs_df = df_ob_stn.where(z >= -(d/2))
            bott_obs_df = df_ob_stn.where(z <= (-1*h) + (d/2))
            # get average value
            surf_obs_avg = surf_obs_df.groupby('time')[var].mean()
            bott_obs_avg = bott_obs_df.groupby('time')[var].mean()
            # plot observations
            unique_time = [pd.Timestamp(x) for x in df_ob_stn['time'].unique()] # one point per timestampe
            axis.scatter(unique_time, surf_obs_avg, color='deeppink', s=20,zorder=10)
            axis.scatter(unique_time, bott_obs_avg, color='royalblue', s=20,zorder=10, label='obs')


            # label
            if i == 0:
                axis.text(0.05, 0.95, 'surface {} m'.format(str(round(d/2))),
                    verticalalignment='top', horizontalalignment='left',
                    transform=axis.transAxes, fontsize=ls, color = 'deeppink',fontweight='bold')
                axis.text(0.05, 0.8, 'bottom {} m'.format(str(round(d/2))),
                    verticalalignment='top', horizontalalignment='left',
                    transform=axis.transAxes, fontsize=ls, color = 'royalblue',fontweight='bold')
                

        # otherwise just one layer
        else:

            # get depth of each vertical layer
            z_thick = np.diff(z_w,axis=0).transpose() # 365,30
            print(z_thick.shape) 

            # model
            mod_ds = ds
            # get model data 
            mod = mod_ds[var]
            # take weighted average given depth
            # first, multiply by thickness of each layer and sum in z, then divide by water column depth
            mod_avg = np.nansum(mod * z_thick, axis=1)/np.nansum(z_thick,axis=1)
            # plot model output
            axis.plot(dates_local, mod_avg, color='mediumorchid', linewidth=2, alpha=0.3, zorder=5, label='model')

            # observation
            obs_df = df_ob_stn
            # get average value
            obs_avg = obs_df.groupby('time')[var].mean()
            # plot observations
            unique_time = [pd.Timestamp(x) for x in df_ob_stn['time'].unique()]
            # axis.scatter(time, obs_df[var], color='k',s=3)
            axis.scatter(unique_time, obs_avg, color='mediumorchid', s=20, zorder=10, label='obs')

        axis.set_xlim((dates_local[0],dates_local[-1]))
        axis.set_ylabel(name)
        if i ==0:
            axis.set_title('Time series comparison')
            axis.legend(loc='upper right',fontsize=11, frameon=False, handletextpad=0.1,
                        handlelength=1, markerfirst = False)



        # FORMATTING ---------------------------------------------------------

        def month_initial(x, pos):
                return mdates.num2date(x).strftime('%b')[0]

        # add bottom axis
        for col in [0,1,2]:

            axis = ax[i,col]
            # format grid
            axis.tick_params(axis='both', which='major', labelsize=ls)
            if col > 0 :
                axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')

            # format axis labels
            if i == rows-1 and col < 2:
                axis.set_xlabel(year, fontsize = fs)
                axis.tick_params(axis='both', which='major', labelsize=ls)
                axis.xaxis.set_major_formatter(ticker.FuncFormatter(month_initial))
                # axis.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
                axis.tick_params(axis='x', labelrotation=0, labelsize=ls)
            elif i < rows-1 and col != 2:
                axis.set_xticklabels([])
                

            # format colors
            axis.set_facecolor('#EEEEEE')
            for border in ['top','right','bottom','left']:
                axis.spines[border].set_visible(False)

    #         # place colorbar
    #         cbar = fig.colorbar(cs, ax=axis, aspect=10)
    #         cbar.outline.set_visible(False)

    plt.tight_layout
    plt.subplots_adjust(hspace = 0.2, wspace=0.3, top=0.9, left=0.08, right=0.98, bottom = 0.08) # wspace=0.02, 
    plt.savefig(out_dir / (station + '(' + str(round(lon,2)) + ',' + str(round(lat,2)) + ').png'))