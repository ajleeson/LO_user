"""
Plot surface and deep budget for all 21 inlets
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import get_two_layer
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec

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

# initialize figure
fig,ax = plt.subplots(2,1, figsize=(8,6))

# plot DOin-DOout for all inlets
for i,station in enumerate(sta_dict):

    # Get DO concentrations
    in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'c21' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
    bulk = xr.open_dataset(in_dir)
    tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
    DO_in  = tef_df['oxygen_p'].values * 32/1000 # DOin [mg/L] = 32/1000 * [mmol/m3]
    DO_out = tef_df['oxygen_m'].values * 32/1000# DOout [mg/L] = 32/1000 * [mmol/m3]
    Qin  = tef_df['q_p'].values # [m3/s]
    Qout = tef_df['q_m'].values # [m3/s]
    # calculate DOin - DOout
    DOin_DOout = DO_in - DO_out
    DOin_DOout_filtered = zfun.lowpass(DOin_DOout,n=30)
    # calculate Qin - Qout
    Qin_Qout = Qin - Qout
    Qin_Qout_filtered = zfun.lowpass(Qin_Qout,n=30)

    # plot differences
    # DO
    ax[0].plot(dates_local_daily,DOin_DOout_filtered,linewidth=3,color='white',alpha=0.5)
    ax[0].plot(dates_local_daily,DOin_DOout_filtered,linewidth=2,color='black',alpha=0.5)
    ax[0].text(0.02,0.9,'(a) DOin - DOout',fontsize=14,fontweight='bold',transform=ax[0].transAxes)
    ax[0].set_ylabel('DOin - DOout' + r'[mg L$^{-1}$]',size=12)
    # Q
    ax[1].plot(dates_local_daily,Qin_Qout_filtered,linewidth=3,color='white',alpha=0.5)
    ax[1].plot(dates_local_daily,Qin_Qout_filtered,linewidth=2,color='black',alpha=0.5)
    ax[1].text(0.02,0.9,'(b) Qin - Qout',fontsize=14,fontweight='bold',transform=ax[1].transAxes)
    ax[1].set_ylabel('Qin - Qout' + r'[m3 s$^{-1}$]',size=12)

# format figure
minday = 194
maxday = 256
for axis in ax:
    axis.plot([dates_local_daily[0],dates_local_daily[-1]],[0,0],color='black',linewidth=2,linestyle=':')
    axis.set_xlim([dates_local[0],dates_local[-25]])
    axis.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')
    axis.tick_params(axis='x', labelrotation=30, labelsize=12)
    axis.tick_params(axis='y', labelsize=12)
    loc = mdates.MonthLocator(interval=1)
    axis.xaxis.set_major_locator(loc)
    axis.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    # add decline period
    axis.axvline(dates_local_daily[minday],0,12,color='teal',alpha=0.5)
    axis.axvline(dates_local_daily[maxday],0,12,color='teal',alpha=0.5)
    axis.axvspan(dates_local_daily[minday],dates_local_daily[maxday],
            alpha=0.3, color='paleturquoise')