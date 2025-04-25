"""
2025.04.21

Compare QinDOin and QinNin for the hourly data
and the artificially generated daily averages
"""

import numpy as np
import xarray as xr
import matplotlib.dates as mdates
from datetime import datetime, timedelta
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pylab as plt
import get_two_layer

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
jobname = 'twentyoneinlets'
year = '2017'


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
del sta_dict['dabob']

# get list of inlet names
inlets = list(sta_dict.keys())

    
# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
Lfun.make_dir(out_dir)

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')[2::]
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]
# crop time vector (because we only have jan 2 - dec 30)
dates_no_crop = dates_local_daily
dates_local_daily = dates_local_daily


print('\n')

###################################################################
##         Get hourly data from original tef2 extractions        ##
###################################################################

print('Getting all data for analysis\n')

# create dictionaries to save dataframe of budget terms for each station
hourly_dict = {}
for i,station in enumerate(inlets):
    hourly_dict[station] = pd.DataFrame()

# get lat and lon of grid
Ldir['ds0'] = startdate
in_dir = Ldir['roms_out'] / Ldir['gtagex']
# G, S, T = zrfun.get_basic_info(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
fn0 = xr.open_dataset(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
lonr = fn0.lon_rho.values
latr = fn0.lat_rho.values
# open box extraction
box_fn = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
ds_box = xr.open_dataset(box_fn)
DX = (ds_box.pm.values)**-1
DY = (ds_box.pn.values)**-1
DA = DX*DY # get area of each grid cell in m^2

for i,station in enumerate(inlets):

# --------------------------- get hourly TEF exchange flow terms ----------------------------------------
# Note: lowpass Godin filter already applied

        in_dir = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
        bulk = xr.open_dataset(in_dir)
        tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
        Q_p = tef_df['q_p'] # Qin [m3/s]
        Q_m = tef_df['q_m'] # Qout [m3/s]
        DO_p = tef_df['oxygen_p'] # DOin [mmol/m3]
        DO_m = tef_df['oxygen_m'] # DOout [mmol/m3]
        # convert from mmol/s to kmol/s
        QinDOin = (Q_p.values * DO_p.values) * (1/1000) * (1/1000) # Qin * DOin
        QoutDOout = (Q_m.values * DO_m.values) * (1/1000) * (1/1000) # Qout * DOout

        # nutrients (NO3 + NH4)
        DIN_p = tef_df['NO3_p']+tef_df['NH4_p'] # NH4in [mmol/m3] # NO3in [mmol/m3]
        DIN_m = tef_df['NO3_m']+tef_df['NH4_m'] # NH4out [mmol/m3] # NO3out [mmol/m3]
        # print(np.nanmean(DIN_p))
        # convert from mmol/s to kmol/s
        QoutDINout = (Q_m.values * DIN_m.values) * (1/1000) * (1/1000) # Qout * DINout
        QinDINin = (Q_p.values * DIN_p.values) * (1/1000) * (1/1000) # Qin * DINin

        
# ------------------------- save data in dataframe dict -----------------------------------
        # Note: everything is in units of kmol O2 /s

        hourly_dict[station]['QinDOin hourly (kmol/s)'] = QinDOin
        hourly_dict[station]['QoutDOout hourly (kmol/s)'] = QoutDOout
        hourly_dict[station]['QinDINin hourly (kmol/s)'] = QinDINin
        hourly_dict[station]['QoutDINout hourly (kmol/s)'] = QoutDINout
        hourly_dict[station]['Qin m3/s'] = Q_p.values 

        # print keys
        if i == 0:
            print(list(hourly_dict[station].keys()))


###################################################################
##       Get artificially created Eulerian daily averages        ##
################################################################### 

daily_dict = {}

for i,station in enumerate(inlets):
        
        daily_dict[station] = pd.DataFrame()

# --------------------------- get hourly Eulerian exchange flow terms ----------------------------------------

        in_dir = Ldir['LOo'] / 'loading_test' / 'artificial_daily_averages' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
        bulk = xr.open_dataset(in_dir)
        tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
        Q_p = tef_df['q_p'] # Qin [m3/s]
        Q_m = tef_df['q_m'] # Qout [m3/s]
        DO_p = tef_df['oxygen_p'] # DOin [mmol/m3]
        DO_m = tef_df['oxygen_m'] # DOout [mmol/m3]
        # convert from mmol/s to kmol/s
        QinDOin = (Q_p.values * DO_p.values) * (1/1000) * (1/1000) # Qin * DOin
        QoutDOout = (Q_m.values * DO_m.values) * (1/1000) * (1/1000) # Qout * DOout

        # nutrients (NO3 + NH4)
        DIN_p = tef_df['NO3_p']+tef_df['NH4_p'] # NH4in [mmol/m3] # NO3in [mmol/m3]
        DIN_m = tef_df['NO3_m']+tef_df['NH4_m'] # NH4out [mmol/m3] # NO3out [mmol/m3]
        # print(np.nanmean(DIN_p))
        # convert from mmol/s to kmol/s
        QoutDINout = (Q_m.values * DIN_m.values) * (1/1000) * (1/1000) # Qout * DINout
        QinDINin = (Q_p.values * DIN_p.values) * (1/1000) * (1/1000) # Qin * DINin

        
# ------------------------- save data in dataframe dict -----------------------------------
        # Note: everything is in units of kmol O2 /s

        daily_dict[station]['QinDOin hourly (kmol/s)'] = QinDOin
        daily_dict[station]['QoutDOout hourly (kmol/s)'] = QoutDOout
        daily_dict[station]['QinDINin hourly (kmol/s)'] = QinDINin
        daily_dict[station]['QoutDINout hourly (kmol/s)'] = QoutDINout
        daily_dict[station]['Qin m3/s'] = Q_p.values 



###################################################################
##          Plot comparison of hourly and daily values           ##
###################################################################

for i,station in enumerate(inlets):
     
    # initialize figure
    fig, ax = plt.subplots(3,1,figsize = (10,8), sharex=True)
    ax = ax.ravel()

    # format grid and tick label sizes
    for axis in ax:
        # axis.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')
        axis.tick_params(axis='both', labelsize=12)
        axis.set_xlim(dates_local[0],dates_local[-1])
        loc = mdates.MonthLocator(interval=1)
        axis.xaxis.set_major_locator(loc)
        axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        axis.tick_params(axis='x', labelrotation=30)
        axis.set_facecolor('#EEEEEE')
        axis.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
        axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
        axis.tick_params(axis='x', labelrotation=30)
        for border in ['top','right','bottom','left']:
            axis.spines[border].set_visible(False)

    # Qin
    ax[0].text(0.02,0.83,r'(a) Q$_{in}$ [m$^3$ s$^{-1}$]',
               transform=ax[0].transAxes, fontsize=12,
               fontweight='bold')
    ax[0].set_ylabel(r'm$^3$ s$^{-1}$',fontsize=12)
    # plot hourly data
    ax[0].plot(dates_local_daily,hourly_dict[station]['Qin m3/s'],
               color='hotpink',linewidth=2,alpha=0.8,
               label='Hourly (Godin-filtered TEF)')
    # plot daily data
    ax[0].plot(dates_local_daily,daily_dict[station]['Qin m3/s'][1:-1],
               color='royalblue',linewidth=1.5,alpha=0.65,
               label='Daily averages (Eulerian)')
    # set y-axis
    ax[0].set_ylim([0,1.1*np.nanmax(hourly_dict[station]['Qin m3/s'])])
    # add difference
    divider = make_axes_locatable(ax[0])
    ax2 = divider.append_axes("bottom", size='20%', pad=0.2)
    ax[0].figure.add_axes(ax2)
    ax2.plot(dates_local_daily,
             hourly_dict[station]['Qin m3/s'].values-daily_dict[station]['Qin m3/s'].values[1:-1],
             color='mediumturquoise')
    ax2.axhline(0,color='k', linestyle=':')
    ax2.set_xticks([])
    ax2.set_xlim(dates_local[0],dates_local[-1])
    # format difference
    ax2.set_facecolor('#EEEEEE')
    ax2.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax2.tick_params(axis='x', labelrotation=30)
    for border in ['top','right','bottom','left']:
        ax2.spines[border].set_visible(False)
    ax2.text(0.02,0.9,'Hourly - Daily',
               transform=ax2.transAxes, fontsize=10,
               fontweight='bold')

    # QinDOin
    ax[1].text(0.02,0.83,r'(b) Q$_{in}$DO$_{in}$ [kmol O$_2$ s$^{-1}$]',
               transform=ax[1].transAxes, fontsize=12,
               fontweight='bold')
    ax[1].set_ylabel(r'kmol O$_2$ s$^{-1}$',fontsize=12)
    # plot hourly data
    ax[1].plot(dates_local_daily,hourly_dict[station]['QinDOin hourly (kmol/s)'],
               color='hotpink',linewidth=2,alpha=0.8)
    # plot daily data
    ax[1].plot(dates_local_daily,daily_dict[station]['QinDOin hourly (kmol/s)'][1:-1],
               color='royalblue',linewidth=1.5,alpha=0.65)
    # set y-axis
    ax[1].set_ylim([0,1.1*np.nanmax(hourly_dict[station]['QinDOin hourly (kmol/s)'])])
    # add difference
    divider = make_axes_locatable(ax[1])
    ax2 = divider.append_axes("bottom", size='20%', pad=0.2)
    ax[1].figure.add_axes(ax2)
    ax2.plot(dates_local_daily,
             hourly_dict[station]['QinDOin hourly (kmol/s)'].values-daily_dict[station]['QinDOin hourly (kmol/s)'].values[1:-1],
             color='mediumturquoise')
    ax2.axhline(0,color='k', linestyle=':')
    ax2.set_xticks([])
    ax2.set_xlim(dates_local[0],dates_local[-1])
    # format difference
    ax2.set_facecolor('#EEEEEE')
    ax2.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax2.tick_params(axis='x', labelrotation=30)
    for border in ['top','right','bottom','left']:
        ax2.spines[border].set_visible(False)

    # QinDINin
    ax[2].text(0.02,0.83,r'(c) Q$_{in}$DIN$_{in}$ [kmol N s$^{-1}$]',
               transform=ax[2].transAxes, fontsize=12,
               fontweight='bold')
    ax[2].set_ylabel(r'kmol N s$^{-1}$',fontsize=12)
    # plot hourly data
    ax[2].plot(dates_local_daily,hourly_dict[station]['QinDINin hourly (kmol/s)'],
               color='hotpink',linewidth=2,alpha=0.8)
    # plot daily data
    ax[2].plot(dates_local_daily,daily_dict[station]['QinDINin hourly (kmol/s)'][1:-1],
               color='royalblue',linewidth=1.5,alpha=0.65)
    # set y-axis
    ax[2].set_ylim([0,1.1*np.nanmax(hourly_dict[station]['QinDINin hourly (kmol/s)'])])
    # remove x-axis
    ax[2].tick_params(axis='x', labelcolor='white')
    # add difference
    divider = make_axes_locatable(ax[2])
    ax2 = divider.append_axes("bottom", size='20%', pad=0.2)
    ax[2].figure.add_axes(ax2)
    ax2.plot(dates_local_daily,
             hourly_dict[station]['QinDINin hourly (kmol/s)'].values-daily_dict[station]['QinDINin hourly (kmol/s)'].values[1:-1],
             color='mediumturquoise')
    ax2.axhline(0,color='k', linestyle=':')
    ax2.set_xlim(dates_local[0],dates_local[-1])
    loc = mdates.MonthLocator(interval=1)
    ax2.xaxis.set_major_locator(loc)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax2.tick_params(axis='x', labelrotation=30)
    # format difference
    ax2.set_facecolor('#EEEEEE')
    ax2.grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax2.tick_params(axis='x', labelrotation=30)
    for border in ['top','right','bottom','left']:
        ax2.spines[border].set_visible(False)
    ax2.set_xlabel('2017', fontsize=12)

    # format figure
    ax[0].set_title(station,fontsize=14,fontweight='bold')
    ax[0].legend(loc='upper right', fontsize=12)

plt.subplots_adjust(hspace=3)      
plt.tight_layout()

