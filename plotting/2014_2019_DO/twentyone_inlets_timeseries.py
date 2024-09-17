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
enddate_hrly = '2015.01.01 00:00:00'
year = '2014' # for making a date label

vn_list = ['oxygen','Urms','delta_rho','NO3','NH4','phytoplankton','zooplankton','SdetritusN','LdetritusN']
rows = len(vn_list)

# figure settings
fs = 14 # figure font size
ls = 14 # label size
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
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures' / 'twentyone' / 'timeseries'
Lfun.make_dir(out_dir)

# helper functions

def get_surf_average(vn,d):
    '''
    calculate average value of vn in surface d meters
    '''
    z_ds = ds.where((ds.z_w >= -d))
    var_ds = ds.where((ds.z_rho >= -d))
    # get thickness of each vertical layer
    z_thick = np.diff(z_ds['z_w'],axis=1) # 365,30
    # get variable data 
    var = var_ds[vn].values
    # check if thickness array is all nan
    # (this occurs if the first z_w is already greater than the threshold, so we don't have two z_w values to diff)
    # in which case, average value is just the single value
    if np.isnan(z_thick).all():
        val = np.nansum(var,axis=1)
    # take weighted average given depth
    # first, multiply by thickness of each layer and sum in z, then divide by layer thickness
    else:
        val = np.nansum(var * z_thick, axis=1)/np.nansum(z_thick,axis=1)
    return val

def get_bott_average(vn,d,h):
    '''
    calculate average value of vn in bottom d meters
    '''
    z_ds = ds.where((ds.z_w <= (-1*h) + d))
    var_ds = ds.where((ds.z_rho <= (-1*h) + d))
    # get thickness of each vertical layer
    z_thick = np.diff(z_ds['z_w'],axis=1) # 365,30
    # get variable data 
    var = var_ds[vn].values
    # check if thickness array is all nan
    # (this occurs if the first z_w is already greater than the threshold, so we don't have two z_w values to diff)
    # in which case, average value is just the single value
    if np.isnan(z_thick).all():
        val = np.nansum(var,axis=1)
    # take weighted average given depth
    # first, multiply by thickness of each layer and sum in z, then divide by layer thickness
    else:
        val = np.nansum(var * z_thick, axis=1)/np.nansum(z_thick,axis=1)
    return val

##########################################################
##                      Plotting                        ##
##########################################################

letter = ['(a)','(b)','c']

# create time vector
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local_hrly = [pfun.get_dt_local(x) for x in dates_hrly]

# Loop through all of the mooring stations
for i,station in enumerate(sta_dict): # enumerate(['commencement']): #
    plt.close('all')

    print(station)

    # calculate lat/lon for station
    lon = sta_dict[station][0]
    lat = sta_dict[station][1]

    # Initialize Figure
    fig, ax = plt.subplots(3,1,figsize = (18,6),sharex=True,
                           gridspec_kw={'height_ratios': [2,1,1]})
    
    # loop through different state variables
    for i,vn in enumerate(vn_list):
            
        if vn == 'oxygen':
                axis = ax[2]
                plotlabel = r'(c) $DO_{bottom}$'
        elif vn in ['Urms','delta_rho']:
                axis = ax[1]
                plotlabel = r'(b) $U_{rms}$ and $\Delta\rho$'
        else:
                axis = ax[0]
                plotlabel = '(a) Pools of nitrogen'

        # get data
        fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds = xr.open_dataset(fn)
        # hourly data for Urms
        fn = '../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '_hourly/' + station + '_' + startdate + '_' + enddate + '.nc'
        ds_hrly = xr.open_dataset(fn)

        # handle each state variable differently
        # vn_list = ['oxygen','delta_rho','NO3','NH4','phytoplankton','zooplankton','SdetritusN','LdetritusN']

        # BOTTOM DO
        if vn == 'oxygen':
            color = 'crimson'
            bott_ox_raw = ds['oxygen'].values[:,0] # 0 = bottom layer
            val = bott_ox_raw *  pinfo.fac_dict['oxygen'] # convert to mg/L
            # plot
            axis.plot(dates_local,val,linewidth=5,color=color,
                      alpha=0.8,label='Bottom DO')
            # format
            axis.set_ylim([0,13])
            axis.set_ylabel('Bottom DO\n'+r'[mg $L^{-1}$]',fontsize=ls+2,color=color)
            axis.tick_params(axis='y', labelcolor=color)
            axis.set_yticks(np.arange(0, 12, 2))

        # U_RMS (RMS tidal currents)
        if vn == 'Urms':
            color = 'mediumorchid'
            # get rms tidal current
            ubar = ds_hrly['ubar'].values
            vbar = ds_hrly['vbar'].values
            # use equation from Parker to get U_rms
            val = np.sqrt( zfun.lowpass(ubar**2 + vbar**2, f='godin')) 
            # plot
            axis.plot(dates_local_hrly,val,linewidth=3,color=color,
                      alpha=0.6,label='RMS Tidal Currents')
            # format
            ymax = np.nanmax(val) * 1.2 # 0.4
            ticks = 4
            axis.set_ylim([0,ymax])
            axis.set_ylabel(r'$U_{rms}$'+'\n'+r'[m $s^{-1}$]',fontsize=ls+2,color=color)
            axis.tick_params(axis='y', labelcolor=color)
            axis.set_yticks(np.arange(0, ymax, ymax/ticks))

        # DELTA RHO
        if vn == 'delta_rho':
            color = 'black'
            # press_top = [gsw.p_from_z(z,lat) for z in ds['z_rho'].values[:,-1]]
            # press_bott = [gsw.p_from_z(z,lat) for z in ds['z_rho'].values[:,0]]
            # calculate absolute salinity from practical salinity
            salt_abs_top_all = gsw.conversions.SA_from_SP(ds['salt'].values[:,-1], ds['z_rho'].values[:,-1], lon, lat)
            salt_abs_bott_all = gsw.conversions.SA_from_SP(ds['salt'].values[:,0], ds['z_rho'].values[:,0], lon, lat)
            # calculate conservative temperature from potential temperature
            cons_temp_top_all = gsw.conversions.CT_from_pt(salt_abs_top_all, ds['temp'].values[:,-1])
            cons_temp_bott_all = gsw.conversions.CT_from_pt(salt_abs_bott_all, ds['temp'].values[:,0])
            # calculate density
            # rho_top_all = gsw.rho(salt_abs_top_all,cons_temp_top_all,press_top)
            # rho_bott_all = gsw.rho(salt_abs_bott_all,cons_temp_bott_all,press_bott)
            rho_top_all = gsw.density.sigma0(salt_abs_top_all,cons_temp_top_all)
            rho_bott_all = gsw.density.sigma0(salt_abs_bott_all,cons_temp_bott_all)
            # calculate density difference
            rho_diff_all = rho_bott_all - rho_top_all
            # save value
            val = rho_diff_all
            # plot
            ax_rho = ax[1].twinx()
            axis = ax_rho
            axis.plot(dates_local,val,linewidth=1.5,color=color,
                      alpha=1,label=r'$\Delta\rho$')
            # format
            # axis.set_ylim([0,24])
            ymax = np.nanmax(val) * 1.2 # 24
            axis.set_ylim([0,ymax])
            axis.set_ylabel(r'$\rho_{bottom}-\rho_{surface}$'+'\n'+r'[kg m$^{-3}$]',fontsize=ls+2,color=color)
            axis.tick_params(axis='y', labelcolor=color)
            # axis.set_yticks(np.arange(0, 24, 4))
            axis.set_yticks(np.arange(0, ymax, ymax/ticks))

        # NO3 (average surface 20 m)
        if vn == 'NO3':
            color = 'navy'
            # only look at surface 20 m
            d=20
            val = get_surf_average('NO3',d)
            # plot
            axis.plot(dates_local,val,linewidth=3,color=color,
                      alpha=0.5,label='NO3 (Surf {}m)'.format(str(round(d))))
            # format
            # axis.set_ylim([10e-9,10e2])
            # axis.set_yscale('log')
            axis.set_ylim([0,45])
            axis.set_ylabel('Concentration\n'+r'[mmol N $m^{-3}$]',fontsize=ls+2)

        # NH4 (average surface 20 m)
        if vn == 'NH4':
            color = 'mediumvioletred'
            # only look at surface 20 m
            d=20
            val = get_surf_average('NH4',d)
            # plot
            axis.plot(dates_local,val,linewidth=2,color=color,
                      alpha=0.8,label='NH4  (Surf {}m)'.format(str(round(d))))
            
        # phyto (average surface 20 m)
        if vn == 'phytoplankton':
            color = 'olivedrab' #'#79c000'
            # only look at surface 20 m
            d=20
            val = get_surf_average('phytoplankton',d)
            # plot
            axis.plot(dates_local,val,linewidth=3,color=color,
                      alpha=0.8,label='phytoplankton (Surf {}m)'.format(str(round(d))))
            
        # zoop (average surface 20 m)
        if vn == 'zooplankton':
            color = 'deepskyblue'
            # only look at surface 20 m
            d=20
            val = get_surf_average('zooplankton',d)
            # plot
            axis.plot(dates_local,val,linewidth=3,color=color,
                      alpha=1,label='zooplankton (Surf 20m)')
            
        # small detritus (average surface 80 m)
        if vn == 'SdetritusN':
            color = 'black'
            # only look at surface 80 m
            d=80
            val = get_surf_average('SdetritusN',d)
            # plot
            axis.plot(dates_local,val,linewidth=1.5,color=color,
                      alpha=1,label='small detritus (Surf {}m)'.format(str(round(d))))
            
        # large detritus (average bottom 5 m)
        if vn == 'LdetritusN':
            color = 'saddlebrown'
            # only look at bottom 5m
            d=5
            val = get_bott_average('LdetritusN',d,ds['h'].values)
            # plot
            axis.plot(dates_local,val,linewidth=3,color=color,
                      alpha=0.5,label='large detritus (Bott {}m)'.format(str(round(d))))


        # FORMATTING ---------------------------------------------------------
        depth = round(int(ds.h.values))
        fig.suptitle(station + ' ' + year + ' ({} m deep)'.format(str(depth)), fontsize = 18)
        axis.set_xlim([dates_local[0],dates_local[-1]])
        ax[0].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        ax[1].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        ax[2].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        axis.tick_params(axis='both', which='major', labelsize=ls)
        axis.xaxis.set_major_formatter(mdates.DateFormatter("%b-%d"))
        axis.tick_params(axis='x', labelrotation=30, labelsize=ls)
        axis.text(0.02, 0.96, plotlabel,
                verticalalignment='top', horizontalalignment='left',
                transform=axis.transAxes, fontsize=15)
        if vn == 'oxygen':
                axis.set_xlabel(year, fontsize = fs)

        # format colors
        axis.set_facecolor('#EEEEEE')
        for border in ['top','right','bottom','left']:
            axis.spines[border].set_visible(False)

    ax[0].legend(loc='lower center',fontsize=ls-1,ncol=3, handletextpad=0.2, bbox_to_anchor=(0.5,0.67))
    # ax[1].legend(loc='upper center',fontsize=ls-1, frameon=False, handletextpad=0.2)
    # ax_rho.legend(loc='upper center',fontsize=ls-1,frameon=False, handletextpad=0.2,
    #               bbox_to_anchor=(0.7,1))

    # plt.tight_layout
    plt.subplots_adjust(hspace=0.1,top=0.95, bottom=0.15)
    plt.savefig(out_dir / (station + '(' + str(round(lon,2)) + ',' + str(round(lat,2)) + ').png'))
