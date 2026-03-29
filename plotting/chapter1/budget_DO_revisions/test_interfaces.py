"""
Test different interface depths
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import get_two_layer
import matplotlib.dates as mdates
from matplotlib import colormaps
from matplotlib.colors import ListedColormap

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t1_x11b'
jobname = 'twentyoneinlets'
year = '2017'

stations = 'all'

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# set up dates
startdate = year + '.01.01'
enddate = year + '.12.31'
enddate_hrly = str(int(year)+1)+'.01.01 00:00:00'

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get stations:
if stations == 'all':
    sta_dict = job_lists.get_sta_dict(jobname)
    # remove lynchcove2
    del sta_dict['lynchcove2']
    # remove shallow inlets (< 10 m deep)
    del sta_dict['hammersley']
    del sta_dict['henderson']
    del sta_dict['oak']
    del sta_dict['totten']
    del sta_dict['similk']
    del sta_dict['budd']
    del sta_dict['eld']
    del sta_dict['killsut']
    # del sta_dict['dabob']
else:
    sta_dict = stations

# # where to put output figures
# out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
# Lfun.make_dir(out_dir)

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')[2::]
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]
# crop time vector (because we only have jan 2 - dec 30)
dates_no_crop = dates_local_daily
dates_local_daily = dates_local_daily

print('\n')

##########################################################
##            Get all variables for analysis            ##
##########################################################

print('Getting all data for analysis\n')

# create dictionaries with interface depths
interface_dict = dict()

# create dictionaries to save dataframe of budget terms for each station
bottomlay_dict = {}
surfacelay_dict = {}
DOconcen_dict = {}
dimensions_dict = {}
for i,station in enumerate(sta_dict):
    bottomlay_dict[station] = pd.DataFrame()
    surfacelay_dict[station] = pd.DataFrame()
    DOconcen_dict[station] = pd.DataFrame()
    dimensions_dict[station] = pd.DataFrame()

####################################################################################
#                     COMPARE INTERFACE DEPTH METHODS (d/dt(DO))                  ##
####################################################################################

# initialize figure
fig,axes = plt.subplots(5,1, figsize=(12,8), sharex=True, sharey=True)
ax = axes.ravel()

# inlets roughly sorted by DOdeep
inlets = ['dyes','sinclair','quartermaster','case','crescent','carr',
          'elliot','commencement','penn','portsusan','holmes','dabob','lynchcove']

for i,station in enumerate(inlets):#enumerate(sta_dict):

# ---------------------------------- get BGC terms --------------------------------------------
        bgc_dir = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_bgc' / station
        # get months
        months = [year+'.01.01_'+year+'.01.31',
                    year+'.02.01_'+year+'.02.28',
                    year+'.03.01_'+year+'.03.31',
                    year+'.04.01_'+year+'.04.30',
                    year+'.05.01_'+year+'.05.31',
                    year+'.06.01_'+year+'.06.30',
                    year+'.07.01_'+year+'.07.31',
                    year+'.08.01_'+year+'.08.31',
                    year+'.09.01_'+year+'.09.30',
                    year+'.10.01_'+year+'.10.31',
                    year+'.11.01_'+year+'.11.30',
                    year+'.12.01_'+year+'.12.31',]
        # months = [year+'.07.01_'+year+'.07.31',
        #           year+'.08.01_'+year+'.08.31',]
        #         #   year+'.09.01_'+year+'.09.30',
        #         #   year+'.10.01_'+year+'.10.31',
        #         #   year+'.11.01_'+year+'.11.30',
        #         #   year+'.12.01_'+year+'.12.31',]
        

        interface_types = ['og','drdz','tef','halocline','oxycline'] # how to define dividing depth
        # loop through different interface types
        for t,type in enumerate(interface_types):
             
            # initialize arrays to save values
            photo_surf_unfiltered = []
            photo_deep_unfiltered = []
            cons_surf_unfiltered = [] # nitrification, respiration
            cons_deep_unfiltered = [] # nitrification, respiration, sediment oxygen demand
            airsea_surf_unfiltered = []
            o2vol_surf_unfiltered = []
            o2vol_deep_unfiltered = []
            vol_surf_unfiltered = []
            vol_deep_unfiltered = []
            DO_surf_unfiltered = []
            DO_deep_unfiltered = []

            # combine all months
            for month in months:
                fn = type + '_' + station + '_' + month + '.p'
                df_bgc = pd.read_pickle(bgc_dir/fn)
                # conversion factor to go from mmol O2/hr to kmol O2/s
                conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
                # get photosynthesis
                photo_surf_unfiltered = np.concatenate((photo_surf_unfiltered, df_bgc['surf photo [mmol/hr]'].values * conv)) # kmol/s
                photo_deep_unfiltered = np.concatenate((photo_deep_unfiltered, df_bgc['deep photo [mmol/hr]'].values * conv)) # kmol/s
                # get consumption
                surf_cons_terms = df_bgc['surf nitri [mmol/hr]'].values + df_bgc['surf respi [mmol/hr]'].values
                deep_cons_terms = df_bgc['deep nitri [mmol/hr]'].values + df_bgc['deep respi [mmol/hr]'].values + df_bgc['SOD [mmol/hr]'].values
                cons_surf_unfiltered = np.concatenate((cons_surf_unfiltered, surf_cons_terms * conv * -1)) # kmol/s; multiply by -1 b/c loss term
                cons_deep_unfiltered = np.concatenate((cons_deep_unfiltered, deep_cons_terms * conv * -1)) # kmol/s; multiply by -1 b/c loss term
                # get air-sea gas exchange
                airsea_surf_unfiltered = np.concatenate((airsea_surf_unfiltered, df_bgc['airsea [mmol/hr]'].values * conv)) # kmol/s
                # get (DO*V)
                o2vol_surf_unfiltered = np.concatenate((o2vol_surf_unfiltered, df_bgc['surf DO*V [mmol]'].values)) # mmol
                o2vol_deep_unfiltered = np.concatenate((o2vol_deep_unfiltered, df_bgc['deep DO*V [mmol]'].values)) # mmol
                # get volume
                vol_surf_unfiltered = np.concatenate((vol_surf_unfiltered, df_bgc['surf vol [m3]'].replace(0, np.nan).values)) # m3, and replace zeros with nan
                vol_deep_unfiltered = np.concatenate((vol_deep_unfiltered, df_bgc['deep vol [m3]'].values)) # m3
            # get DO concentration [mg/L]
            DO_surf_unfiltered = o2vol_surf_unfiltered/vol_surf_unfiltered * 32/1000 # mg/L
            DO_deep_unfiltered = o2vol_deep_unfiltered/vol_deep_unfiltered * 32/1000 # mg/L


            # take time derivative of (DO*V) to get d/dt (DO*V)
            ddtDOV_surf_unfiltered = np.diff(o2vol_surf_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
            ddtDOV_deep_unfiltered = np.diff(o2vol_deep_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s

            # apply Godin filter
            photo_surf  = zfun.lowpass(photo_surf_unfiltered, f='godin')[36:-34:24]
            photo_deep  = zfun.lowpass(photo_deep_unfiltered, f='godin')[36:-34:24]
            cons_surf   = zfun.lowpass(cons_surf_unfiltered, f='godin')[36:-34:24]
            cons_deep   = zfun.lowpass(cons_deep_unfiltered, f='godin')[36:-34:24]
            airsea_surf = zfun.lowpass(airsea_surf_unfiltered, f='godin')[36:-34:24]
            ddtDOV_surf = zfun.lowpass(ddtDOV_surf_unfiltered, f='godin')[36:-34:24]
            ddtDOV_deep = zfun.lowpass(ddtDOV_deep_unfiltered, f='godin')[36:-34:24]
            DO_surf = zfun.lowpass(DO_surf_unfiltered, f='godin')[36:-34:24]
            DO_deep = zfun.lowpass(DO_deep_unfiltered, f='godin')[36:-34:24]
            vol_surf = zfun.lowpass(vol_surf_unfiltered, f='godin')[36:-34:24]
            vol_deep = zfun.lowpass(vol_deep_unfiltered, f='godin')[36:-34:24]
            # get volume-normalized d/dt(DO) and mg/L/day units (godin filtered)
            ddtDOV_surf_mgL_day = ddtDOV_surf/vol_surf * 1000 * 32 * 60 * 60 * 24
            ddtDOV_deep_mgL_day = ddtDOV_deep/vol_deep * 1000 * 32 * 60 * 60 * 24

            # plot
            minday = 194
            maxday = 256
            # get colormaps
            cmap_temp = colormaps['rainbow_r'].resampled(256)
            cmap_oxy = ListedColormap(cmap_temp(np.linspace(0, 1, 256)))# get range of colormap
            # get mean and standard deviation
            mean_ddt_DO = np.nanmean(ddtDOV_deep_mgL_day[minday:maxday])
            std_ddt_DO = np.nanstd(ddtDOV_deep_mgL_day[minday:maxday])
            # plot standard deviation
            ax[t].errorbar(i,mean_ddt_DO, yerr=std_ddt_DO, fmt="o", color='black')
            # plot mean, colored by DO
            cs = ax[t].scatter(i,mean_ddt_DO,s=100, zorder=5, c=np.nanmean(DO_deep), cmap=cmap_oxy,
                               vmin=2,vmax=10)
            # create colorbarlegend
            if i == 0 and t == 0:
                cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
                cbar = fig.colorbar(cs, cax=cbar_ax)
                cbar.ax.tick_params(labelsize=12)
                cbar.set_label(r'Mean DO$_{deep}$ [mg/L]', fontsize=12)
                cbar.outline.set_visible(False)
            # # label mean
            ax[t].text(i-0.3,mean_ddt_DO, str(round(mean_ddt_DO,3)),
                        rotation=90,color='black', va='center', fontweight='bold')
            # format x axis
            ax[4].set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12], inlets, rotation=45, fontsize=12)

            # label y-axis
            if i == 0:
                 ax[t].set_ylabel(type + '\ninterface', fontsize=12)
                 ax[t].tick_params(axis='y', labelsize=12)
                 ax[t].axhline(0,-1,15,linestyle=':', color='silver')

        plt.suptitle(r'Mean mid-Jul to mid-Sep d/dt(DO$_{deep}$)',fontsize=14,fontweight='bold')
        plt.subplots_adjust(right=0.9,hspace=0,bottom=0.15)
        plt.show()


####################################################################################
#         COMPARE INTERFACE DEPTH METHODS (DOdeep vs. DOin colored by Tflush)     ##
####################################################################################

# initialize figure
fig,axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
ax = axes.ravel()

# inlets roughly sorted by DOdeep
inlets = ['dyes','sinclair','quartermaster','case','crescent','carr',
          'elliot','commencement','penn','portsusan','holmes','dabob','lynchcove']

months = [year+'.01.01_'+year+'.01.31',
            year+'.02.01_'+year+'.02.28',
            year+'.03.01_'+year+'.03.31',
            year+'.04.01_'+year+'.04.30',
            year+'.05.01_'+year+'.05.31',
            year+'.06.01_'+year+'.06.30',
            year+'.07.01_'+year+'.07.31',
            year+'.08.01_'+year+'.08.31',
            year+'.09.01_'+year+'.09.30',
            year+'.10.01_'+year+'.10.31',
            year+'.11.01_'+year+'.11.30',
            year+'.12.01_'+year+'.12.31',]

for i,station in enumerate(inlets):#enumerate(sta_dict):
        
# --------------------------- get TEF exchange flow terms ----------------------------------------
    in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'c21' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
    bulk = xr.open_dataset(in_dir)
    tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
    Qin = tef_df['q_p'] # Qin [m3/s]
    DOin = tef_df['oxygen_p'] # DOin [mmol/m3]

    for m,month in enumerate(months):

        if m == 0:
            minday = 0 #1
            maxday = 30 #32
        elif m == 1:
            minday = 30# 32
            maxday = 58# 60
        elif m == 2:
            minday = 58# 60
            maxday = 89# 91
        elif m == 3:
            minday = 89# 91
            maxday = 119# 121
        elif m == 4:
            minday = 119# 121
            maxday = 150# 152
        elif m == 5:
            minday = 150# 152
            maxday = 180# 182
        elif m == 6:
            minday = 180# 182
            maxday = 211# 213
        elif m == 7:
            minday = 211# 213
            maxday = 242# 244
        elif m == 8:
            minday = 242# 244
            maxday = 272# 274
        elif m == 9:
            minday = 272# 274
            maxday = 303# 305
        elif m == 10:
            minday = 303# 305
            maxday = 333# 335
        elif m == 11:
            minday = 333# 335
            maxday = 363

        DOin_monthlymean = np.nanmean(DOin[minday:maxday]) * 32/1000 # [mg/L]
        Qin_daily = Qin[minday:maxday] # [m3/s]

# ---------------------------------- get BGC terms --------------------------------------------
        bgc_dir = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_bgc' / station

        interface_types = ['og','drdz','tef','halocline','oxycline'] # how to define dividing depth
        # loop through different interface types
        for t,type in enumerate(interface_types):

            fn = type + '_' + station + '_' + month + '.p'
            df_bgc = pd.read_pickle(bgc_dir/fn)
            # conversion factor to go from mmol O2/hr to kmol O2/s
            conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
            # get (DO*V)
            o2vol_deep_unfiltered = df_bgc['deep DO*V [mmol]'].values # mmol
            # get volume
            vol_deep_unfiltered = df_bgc['deep vol [m3]'].values # m3
            vol_inlet_unfiltered = df_bgc['surf vol [m3]'].values + df_bgc['deep vol [m3]'].values
            # get DO concentration [mg/L]
            DO_deep_unfiltered = o2vol_deep_unfiltered/vol_deep_unfiltered * 32/1000 # mg/L
            # apply Godin filter
            DO_deep = zfun.lowpass(DO_deep_unfiltered, f='godin')[36:-34:24]
            vol_inlet = zfun.lowpass(vol_inlet_unfiltered, f='godin')[36:-34:24]

            # get monthly mean
            DOdeep_monthlymean = np.nanmean(DO_deep)
            # subsample Qin because we applied godin filter to the month of vol
            if m == 0:
                Qin_daily_resized = Qin_daily[:-1]
            elif m == 11:
                Qin_daily_resized = Qin_daily[1:]
            else:
                Qin_daily_resized = Qin_daily[1:-1]
            Tflush_monthlymean = np.nanmean(vol_inlet/Qin_daily_resized) / (60*60*24) # [days]

            # plot scatter plot
            # define colormap for Tflush and DO concentration
            cmap_temp = colormaps['gnuplot2_r'].resampled(256)
            cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.12, 0.8, 256)))# get range of colormap
            cs = ax[t].scatter(DOin_monthlymean,DOdeep_monthlymean, zorder=5,
                          c =Tflush_monthlymean, cmap=cmap_tflush, vmin=0, vmax=82)

            # format figure
            if i == 0 and m == 0:
                ax[t].text(0.3,11,type, size=14, ha='left', fontweight='bold')
                ax[t].tick_params(axis='x', labelrotation=30)
                ax[t].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
                ax[t].tick_params(axis='both', labelsize=12)
                if t > 2:
                    ax[t].set_xlabel(r'Monthly mean DO$_{in}$ [mg L$^{-1}$]', fontsize=12)
                if t in [0,3]:
                    ax[t].set_ylabel(r'Monthly mean DO$_{deep}$ [mg L$^{-1}$]', fontsize=12)
                # plot
                ax[t].plot([0,12],[0,12],color='gray')
                ax[t].text(1,1,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=9)
                ax[t].set_xlim([0,12])
                ax[t].set_ylim([0,12])
                ax[t].xaxis.set_ticks(np.arange(0, 13, 2))
                ax[t].yaxis.set_ticks(np.arange(0, 13, 2))
            # create colorbarlegend
            if i == 0 and t == 0:
                cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
                cbar = fig.colorbar(cs, cax=cbar_ax)
                cbar.ax.tick_params(labelsize=12)
                cbar.set_label(r'Monthly Mean T$_{flush}$ [days]', fontsize=12)
                cbar.outline.set_visible(False)

        plt.subplots_adjust(hspace=0.1,wspace=0.1)
        plt.show()

####################################################################################
#               COMPARE INTERFACE DEPTH METHODS (DOdeep time series)              ##
####################################################################################

# initialize figure
fig,axes = plt.subplots(3,2, figsize=(10,8.5), sharex=True, sharey=True)
ax = axes.ravel()

# inlets roughly sorted by DOdeep
inlets = ['dyes','sinclair','quartermaster','case','crescent','carr',
          'elliot','commencement','penn','portsusan','holmes','dabob','lynchcove']

for i,station in enumerate(inlets):#enumerate(sta_dict):

# ---------------------------------- get BGC terms --------------------------------------------
        bgc_dir = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_bgc' / station
        # get months
        months = [year+'.01.01_'+year+'.01.31',
                    year+'.02.01_'+year+'.02.28',
                    year+'.03.01_'+year+'.03.31',
                    year+'.04.01_'+year+'.04.30',
                    year+'.05.01_'+year+'.05.31',
                    year+'.06.01_'+year+'.06.30',
                    year+'.07.01_'+year+'.07.31',
                    year+'.08.01_'+year+'.08.31',
                    year+'.09.01_'+year+'.09.30',
                    year+'.10.01_'+year+'.10.31',
                    year+'.11.01_'+year+'.11.30',
                    year+'.12.01_'+year+'.12.31',]
        

        interface_types = ['og','drdz','tef','halocline','oxycline'] # how to define dividing depth
        # loop through different interface types
        for t,type in enumerate(interface_types):
             
            # initialize arrays to save values
            o2vol_surf_unfiltered = []
            o2vol_deep_unfiltered = []
            vol_surf_unfiltered = []
            vol_deep_unfiltered = []
            DO_surf_unfiltered = []
            DO_deep_unfiltered = []

            # combine all months
            for month in months:
                fn = type + '_' + station + '_' + month + '.p'
                df_bgc = pd.read_pickle(bgc_dir/fn)
                # conversion factor to go from mmol O2/hr to kmol O2/s
                conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
                # get (DO*V)
                o2vol_surf_unfiltered = np.concatenate((o2vol_surf_unfiltered, df_bgc['surf DO*V [mmol]'].values)) # mmol
                o2vol_deep_unfiltered = np.concatenate((o2vol_deep_unfiltered, df_bgc['deep DO*V [mmol]'].values)) # mmol
                # get volume
                vol_surf_unfiltered = np.concatenate((vol_surf_unfiltered, df_bgc['surf vol [m3]'].replace(0, np.nan).values)) # m3, and replace zeros with nan
                vol_deep_unfiltered = np.concatenate((vol_deep_unfiltered, df_bgc['deep vol [m3]'].values)) # m3
            # get DO concentration [mg/L]
            DO_surf_unfiltered = o2vol_surf_unfiltered/vol_surf_unfiltered * 32/1000 # mg/L
            DO_deep_unfiltered = o2vol_deep_unfiltered/vol_deep_unfiltered * 32/1000 # mg/L

            # apply Godin filter
            DO_surf = zfun.lowpass(DO_surf_unfiltered, f='godin')[36:-34:24]
            DO_deep = zfun.lowpass(DO_deep_unfiltered, f='godin')[36:-34:24]
            # apply hanning window
            DO_deep_hanning = zfun.lowpass(DO_deep,n=30)

            # plot
            if station == 'lynchcove':
                # lynchcove is a different color
                ax[t].plot(dates_local_daily,DO_deep_hanning,linewidth=3,color='white',alpha=0.5)
                ax[t].plot(dates_local_daily,DO_deep_hanning,linewidth=2,color='crimson')
                ax[t].text(0.78,0.08,'Lynch Cove',transform=ax[t].transAxes, rotation=25,
                        color='crimson', fontsize=10, fontweight='bold')
            else:
                ax[t].plot(dates_local_daily,DO_deep_hanning,linewidth=3,color='white',alpha=0.5)
                ax[t].plot(dates_local_daily,DO_deep_hanning,linewidth=2,color='black',alpha=0.5)

            # format figure
            if i == 0:
                ax[t].text(dates_local_daily[10],10,type, size=14, ha='left', fontweight='bold')
                ax[t].set_xlim(dates_hrly[0],dates_hrly[-25])
                ax[t].tick_params(axis='x', labelrotation=30)
                loc = mdates.MonthLocator(interval=1)
                ax[t].xaxis.set_major_locator(loc)
                ax[t].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
                ax[t].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
                ax[t].tick_params(axis='both', labelsize=12)
                if t in [0,2,4]:
                    ax[t].set_ylabel(r'DO$_{deep}$ [mg/L]', fontsize=12)
                # add decline period
                ax[t].axvline(dates_local_daily[194],0,12,color='teal',alpha=0.5)
                ax[t].axvline(dates_local_daily[256],0,12,color='teal',alpha=0.5)
                ax[t].axvspan(dates_local_daily[194],dates_local_daily[256],
                        alpha=0.3, color='paleturquoise')
                # add hypoxic period (Sep 1 - Oct31)
                ax[t].axvline(dates_local_daily[242],0,12,color='hotpink',alpha=0.5)
                ax[t].axvline(dates_local_daily[302],0,12,color='hotpink',alpha=0.5)
                ax[t].axvspan(dates_local_daily[242],dates_local_daily[302],
                        alpha=0.3, color='pink')

        plt.subplots_adjust(hspace=0.02,wspace=0.02)
        plt.tight_layout()
        plt.show()


####################################################################################
#                      SAVE DO_deep TIME SERIES FOR MANUSCRIPT                    ##
####################################################################################

# inlets roughly sorted by DOdeep
inlets = ['dyes','sinclair','quartermaster','case','crescent','carr',
          'elliot','commencement','penn','portsusan','holmes','dabob','lynchcove']

# initialize dictionary for results
DOdeep_timeseries = {}

for i,station in enumerate(inlets):#enumerate(sta_dict):

# ---------------------------------- get BGC terms --------------------------------------------
        bgc_dir = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_bgc' / station
        # get months
        months = [year+'.01.01_'+year+'.01.31',
                    year+'.02.01_'+year+'.02.28',
                    year+'.03.01_'+year+'.03.31',
                    year+'.04.01_'+year+'.04.30',
                    year+'.05.01_'+year+'.05.31',
                    year+'.06.01_'+year+'.06.30',
                    year+'.07.01_'+year+'.07.31',
                    year+'.08.01_'+year+'.08.31',
                    year+'.09.01_'+year+'.09.30',
                    year+'.10.01_'+year+'.10.31',
                    year+'.11.01_'+year+'.11.30',
                    year+'.12.01_'+year+'.12.31',]
        

        interface_types = ['tef'] # how to define dividing depth
        # loop through different interface types
        for t,type in enumerate(interface_types):
             
            # initialize arrays to save values
            o2vol_surf_unfiltered = []
            o2vol_deep_unfiltered = []
            vol_surf_unfiltered = []
            vol_deep_unfiltered = []
            DO_surf_unfiltered = []
            DO_deep_unfiltered = []

            # combine all months
            for month in months:
                fn = type + '_' + station + '_' + month + '.p'
                df_bgc = pd.read_pickle(bgc_dir/fn)
                # get (DO*V)
                o2vol_deep_unfiltered = np.concatenate((o2vol_deep_unfiltered, df_bgc['deep DO*V [mmol]'].values)) # mmol
                # get volume
                vol_deep_unfiltered = np.concatenate((vol_deep_unfiltered, df_bgc['deep vol [m3]'].values)) # m3
            # get DO concentration [mg/L]
            DO_deep_unfiltered = o2vol_deep_unfiltered/vol_deep_unfiltered * 32/1000 # mg/L

            # apply Godin filter
            DO_deep = zfun.lowpass(DO_deep_unfiltered, f='godin')[36:-34:24]
            # apply hanning window
            DO_deep_hanning = zfun.lowpass(DO_deep,n=30)

            # save to dictionary
            if station == 'elliot':
                station = 'elliott' # correct typo
            DOdeep_timeseries[station] = DO_deep_hanning

# save to csv file
hyp_vol_df = pd.DataFrame.from_dict(DOdeep_timeseries)
dates = pd.date_range(start='2017-01-02', end='2017-12-30', freq='D')
date_list = dates.strftime('%Y-%m-%d').tolist()
hyp_vol_df.insert(0, 'date', date_list)
hyp_vol_df.to_csv('../../../../terminal_inlet_DO_rev2/deeplayerDO_mgL_30dayHanning.csv', index=False)

####################################################################################
#                          SAVE MONTHLY MEANS FOR MANUSCRIPT                      ##
####################################################################################

# initialize dataframe to save values
monthly_mean_df = pd.DataFrame(columns=['month','inlet','DOdeep(mg/L)','DOin(mg/L)','Tflush(days)','Perc_hypoxic_volume(%)'])

# inlets roughly sorted by DOdeep
inlets = ['dyes','sinclair','quartermaster','case','crescent','carr',
          'elliot','commencement','penn','portsusan','holmes','dabob','lynchcove']

months = [year+'.01.01_'+year+'.01.31',
            year+'.02.01_'+year+'.02.28',
            year+'.03.01_'+year+'.03.31',
            year+'.04.01_'+year+'.04.30',
            year+'.05.01_'+year+'.05.31',
            year+'.06.01_'+year+'.06.30',
            year+'.07.01_'+year+'.07.31',
            year+'.08.01_'+year+'.08.31',
            year+'.09.01_'+year+'.09.30',
            year+'.10.01_'+year+'.10.31',
            year+'.11.01_'+year+'.11.30',
            year+'.12.01_'+year+'.12.31',]

month_name = ['Jan','Feb','Mar','Apr','May','Jun',
              'Jul','Aug','Sep','Oct','Nov','Dec']

# open files
box_fn = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
ds_box = xr.open_dataset(box_fn)
DX = (ds_box.pm.values)**-1
DY = (ds_box.pn.values)**-1
DA = DX*DY # get area of each grid cell in m^2
fn0 = xr.open_dataset(Ldir['roms_out'] / 'cas7_t1_x11b' / 'f2017.01.01' / 'ocean_his_0002.nc')
lonr = fn0.lon_rho.values
latr = fn0.lat_rho.values
seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
seg_df = pd.read_pickle(seg_name)
ds_oxy = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / ('cas7_t1_x11ab' + '_pugetsoundDO_' + year + '_DO_info.nc'))

for i,station in enumerate(inlets):#enumerate(sta_dict):

# --------------------------- get percent hypoxic volume ----------------------------------------

    ji_list = seg_df[station+'_p']['ji_list']
    jj_LO = [x[0] for x in ji_list] # y; lat; jj
    ii_LO = [x[1] for x in ji_list] # x; lon; ii
    # get lat and lon corresponding to ii and jj indices
    lat_LO = latr[jj_LO,0]
    lon_LO = lonr[0,ii_LO]
    # get corresponding ii and jj indices in box extraction
    lat_box_all = ds_box['lat_rho'].values[:,0]
    lon_box_all = ds_box['lon_rho'].values[0,:]
    jj = np.zeros(len(jj_LO))
    ii = np.zeros(len(ii_LO))
    for j_ind,lat in enumerate(lat_LO):
        jj[j_ind] = np.where(lat_box_all==lat)[0][0]
    for i_ind,lon in enumerate(lon_LO):
        ii[i_ind] = np.where(lon_box_all==lon)[0][0]
    # convert to array of ints
    jj = jj.astype(int)
    ii = ii.astype(int)

    # calculate percent hypoxic volume
    hyp_thick = ds_oxy['hyp_thick'].values # m

    # get hypoxic thickness and cell area corresponding to the inlet
    hyp_thick = hyp_thick[:,jj,ii] # m
    DA_inlet = DA[jj,ii]
    # calculate hypoxic volume
    hyp_vol = np.sum(hyp_thick * DA_inlet, axis=(1)) # m^3
    # crop to 363 days
    hyp_vol = hyp_vol[1:-1]
        
# --------------------------- get TEF exchange flow terms ----------------------------------------
    in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'c21' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
    bulk = xr.open_dataset(in_dir)
    tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
    Qin = tef_df['q_p'] # Qin [m3/s]
    DOin = tef_df['oxygen_p'] # DOin [mmol/m3]

    for m,month in enumerate(months):

        if m == 0:
            minday = 0 #1
            maxday = 30 #32
        elif m == 1:
            minday = 30# 32
            maxday = 58# 60
        elif m == 2:
            minday = 58# 60
            maxday = 89# 91
        elif m == 3:
            minday = 89# 91
            maxday = 119# 121
        elif m == 4:
            minday = 119# 121
            maxday = 150# 152
        elif m == 5:
            minday = 150# 152
            maxday = 180# 182
        elif m == 6:
            minday = 180# 182
            maxday = 211# 213
        elif m == 7:
            minday = 211# 213
            maxday = 242# 244
        elif m == 8:
            minday = 242# 244
            maxday = 272# 274
        elif m == 9:
            minday = 272# 274
            maxday = 303# 305
        elif m == 10:
            minday = 303# 305
            maxday = 333# 335
        elif m == 11:
            minday = 333# 335
            maxday = 363

        DOin_monthlymean = np.nanmean(DOin[minday:maxday]) * 32/1000 # [mg/L]
        Qin_daily = Qin[minday:maxday] # [m3/s]

        hyp_vol_daily = hyp_vol[minday:maxday]

# ---------------------------------- get BGC terms --------------------------------------------
        bgc_dir = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_bgc' / station

        interface_types = ['tef'] # how to define dividing depth
        # loop through different interface types
        for t,type in enumerate(interface_types):

            fn = type + '_' + station + '_' + month + '.p'
            df_bgc = pd.read_pickle(bgc_dir/fn)
            # conversion factor to go from mmol O2/hr to kmol O2/s
            conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
            # get (DO*V)
            o2vol_deep_unfiltered = df_bgc['deep DO*V [mmol]'].values # mmol
            # get volume
            vol_deep_unfiltered = df_bgc['deep vol [m3]'].values # m3
            vol_inlet_unfiltered = df_bgc['surf vol [m3]'].values + df_bgc['deep vol [m3]'].values
            # get DO concentration [mg/L]
            DO_deep_unfiltered = o2vol_deep_unfiltered/vol_deep_unfiltered * 32/1000 # mg/L
            # apply Godin filter
            DO_deep = zfun.lowpass(DO_deep_unfiltered, f='godin')[36:-34:24]
            vol_inlet = zfun.lowpass(vol_inlet_unfiltered, f='godin')[36:-34:24]

            # get monthly mean
            DOdeep_monthlymean = np.nanmean(DO_deep)
            # subsample Qin because we applied godin filter to the month of vol
            if m == 0:
                Qin_daily_resized = Qin_daily[:-1]
                hyp_vol_daily_resized = hyp_vol_daily[:-1]
            elif m == 11:
                Qin_daily_resized = Qin_daily[1:]
                hyp_vol_daily_resized = hyp_vol_daily[1:]
            else:
                Qin_daily_resized = Qin_daily[1:-1]
                hyp_vol_daily_resized = hyp_vol_daily[1:-1]
            Tflush_monthlymean = np.nanmean(vol_inlet/Qin_daily_resized) / (60*60*24) # [days]

            perc_hyp_vol_monthlymean = np.nanmean(hyp_vol_daily_resized/vol_inlet) * 100

            # add data to dataframe
            # save to dictionary
            if station == 'elliot':
                station_name = 'elliott' # correct typo
            else:
                station_name = station
            new_data = {'month': [month_name[m]], 
                        'inlet':   [station_name], 
                        'DOdeep(mg/L)':  [DOdeep_monthlymean],
                        'DOin(mg/L)': [DOin_monthlymean],
                        'Tflush(days)': [Tflush_monthlymean],
                        'Perc_hypoxic_volume(%)': [perc_hyp_vol_monthlymean]}
            df_new_rows = pd.DataFrame(new_data)
            monthly_mean_df = pd.concat([monthly_mean_df, df_new_rows],ignore_index=True)


# save to csv file
print(monthly_mean_df)
monthly_mean_df.to_csv('../../../../terminal_inlet_DO_rev2/inlet_monthly_means.csv', index=False)