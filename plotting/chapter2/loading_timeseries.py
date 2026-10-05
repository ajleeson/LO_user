

#################################################################################
#                              Import packages                                  #
#################################################################################

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
import pandas as pd
import xarray as xr
import numpy as np
from datetime import datetime, timedelta
from time import time
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import tef_fun

Ldir = Lfun.Lstart()


#################################################################################
#                                   Get data                                    #
#################################################################################

# start_date = '2017.01.01'
# end_date = '2017.12.31'
start_date = '2015.01.01'
end_date = '2020.12.31'
year = '2017'

# create time_vector
dates_daily = pd.date_range(start= start_date, end=end_date, freq= 'd')
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]

ds_loading = xr.open_dataset('../../../LO_output/chapter_2/data/loading_riverwwtp_extraction_'+start_date+'_'+end_date+'.nc')
ds_noloading = xr.open_dataset('../../../LO_output/chapter_2/data/noloading_riverwwtp_extraction_'+start_date+'_'+end_date+'.nc')

# get list of rivers and wwtps that discharge to Puget Sound
wwtps = ['LOTT', 'Puyallup', 'Suquamish', 'Tulalip', 'Fort Lewis', 'Swinomish',
         'STANWOOD STP', 'PORT ORCHARD WWTP', 'OAK HARBOR STP', 'LANGLEY STP',
        'ALDERWOOD STP', 'Lake Stevens Sewer District WWTP',
        'BAINBRIDGE ISLAND WWTP', 'MIDWAY SEWER DISTRICT WWTP',
        'Port Ludlow Wastewater Treatment Plant', 'Port Gamble WWTP',
        'LA CONNER STP', 'MARYSVILLE STP', 'King County Vashon WWTP',
        'LAKOTA WWTP', 'MILLER CREEK WWTP', 'SALMON CREEK WWTP', 'SHELTON STP',
        'MUKILTEO WATER AND WASTEWATER DISTRICT WWTP', 'REDONDO WWTP',
        'MESSENGER HOUSE CARE CENTER WWTP', 'Kitsap County Manchester WWTP',
        'GIG HARBOR STP', 'LYNNWOOD STP', 'EDMONDS STP', 'MT VERNON WWTP',
        'Gardner - Everett Water Pollution Control Facility',
        'Snohomish - Everett Water Pollution Control Facility',
        'King County West Point WWTP', 'BREMERTON STP', 'COUPEVILLE STP',
        'PENN COVE WWTP', 'SNOHOMISH STP', 'King County South WWTP',
        'WARM BEACH CAMPGROUND WWTP',
        'Kitsap County Sewer District #7 Water Reclamation Facility',
        'Kitsap County Central Kitsap WWTP',
        'SKAGIT COUNTY SEWER DIST 2 BIG LAKE WWTP',
        'Kitsap County Kingston WWTP', 'King County Brightwater WWTP',
        'TACOMA CENTRAL NO 1', 'TACOMA NORTH NO 3', 'SEASHORE VILLA STP',
        'TAMOSHAN STP', 'TAYLOR BAY STP', 'ALDERBROOK RESORT & SPA',
        'CARLYON BEACH STP', 'RUSTLEWOOD STP', 'HARTSTENE POINTE STP',
        'CHAMBERS CREEK STP', 'McNeil Island Special Commitment Center WWTP',
        'BOSTON HARBOR STP']
rivers = ['Agate East', 'Agate West', 'Anderson east', 'Anderson west', 'Artondale',
        'Bainbridge Island East', 'Bainbridge Island West', 'Blackjack Cr', 'Blake Island', 'Buenna',
        'Burley Cr+Purdy Cr', 'Butler Cr',  'Campbell Cr',  'Chambers Cr',  'Chico Cr',
        'Coulter Cr', 'Cranberry Cr', 'Curley Cr', 'Dabob Bay', 'Dana Passage North',
        'Dana Passage South', 'Deer Cr+Mable Taylor Cr', 'Des Moines Cr', 'Dutcher Cove', 'Dyes Inlet',
        'Ellis_Mission Cr', 'Ellisport', 'Federal Way', 'Filucy Bay', 'Fox Island',
        'Frye Cove', 'Gallagher Cove', 'Gig Harbor R', 'Glen Cove', 'Goldsborough Cr',
        'Gorst Cr','Grant East','Grant West','Green Cove','Green Valley Cr',
        'Gull Harbor','Hale Passage','Henderson Inlet','Herron','Hope Island',
        'Hylebos Cr','Jarrel Cove','Johns Cr','Judd Cr','Kennedy_Schneider',
        'Ketron','Ketron Island','Kitsap NE','Kitsap_Hood','Liberty Bay',
        'Lynch Cove','Magnolia Bch','Maury Island','Mayo Cove',
        'McAllister Cr','McCormick Cr','McLane Cove','Perry Cr+McLane Cr','McNeil Isl',
        'Mill Cr','Miller Bay','Miller Cr','Minter Cr','Moxlie Cr',
        'NW Hood','Olalla Cr','Peale Passage','Port Gamble R',
        'Port Townsend R','Quilcene','Rocky Cr','Rosedale',
        'Saltwater St Pk','Schneider Cr','Sequalitchew Cr','Sherwood Cr','Shingle Mill Cr',
        'Skookum Cr','Snodgrass Cr','South Snohomish','Squaxin Island East','Squaxin Island West',
        'Sun Pt','Tahlequah','Tahuya','Tolmie','University Place',
        'Van Gelden','Vaughn','Whidbey east','Whidbey west','Whitman Cove',
        'Wilson Pt','Woodard Cr','Woodland Cr','Young Cove','skagit', 'snohomish', 'stillaguamish', 'puyallup',
       'cedar', 'green', 'skokomish', 'dosewallips', 'hamma', 'duckabush', 'deschutes']

# crop ds to rivers and wwtps that discharge to Puget Sound
# Loading
ds_loading_wwtps = ds_loading.sel(riv=wwtps)
ds_loading_rivs  = ds_loading.sel(riv=rivers)
# No-loading
ds_noloading_wwtps = ds_noloading.sel(riv=wwtps)
ds_noloading_rivs  = ds_noloading.sel(riv=rivers)

# get flow and concentration data
# RIVER -------------------------------------
# Loading
Qr_loading_all = ds_loading_rivs['transport'].values # [m3/s] size = 365, nrivs
TNr_loading_all = (ds_loading_rivs['NO3'].values +
                   ds_loading_rivs['NH4'].values +
                   ds_loading_rivs['LDeN'].values +
                   ds_loading_rivs['SDeN'].values +
                   ds_loading_rivs['Phyt'].values +
                   ds_loading_rivs['Zoop'].values)  # mmol/m3
NO3r_loading_all = ds_loading_rivs['NO3'].values # mmol/m3
NH4r_loading_all = ds_loading_rivs['NH4'].values # mmol/m3
Qr_loading = np.sum(Qr_loading_all, axis=1)  # [m3/s], size = 365
TNr_loading = np.sum(Qr_loading_all * TNr_loading_all, axis=1) / Qr_loading  # mmol/m3, size = 365
NO3r_loading = np.sum(Qr_loading_all * NO3r_loading_all, axis=1) / Qr_loading  # mmol/m3, size = 365
NH4r_loading = np.sum(Qr_loading_all * NH4r_loading_all, axis=1) / Qr_loading  # mmol/m3, size = 365

# WWTP -------------------------------------
# Loading
# Loading
Qw_loading_all = ds_loading_wwtps['transport'].values # [m3/s] size = 365, nrivs
TNw_loading_all = (ds_loading_wwtps['NO3'].values +
                   ds_loading_wwtps['NH4'].values +
                   ds_loading_wwtps['LDeN'].values +
                   ds_loading_wwtps['SDeN'].values +
                   ds_loading_wwtps['Phyt'].values +
                   ds_loading_wwtps['Zoop'].values)  # mmol/m3
NO3w_loading_all = ds_loading_wwtps['NO3'].values # mmol/m3
NH4w_loading_all = ds_loading_wwtps['NH4'].values # mmol/m3
Qw_loading = np.sum(Qw_loading_all, axis=1)  # [m3/s], size = 365
TNw_loading = np.sum(Qw_loading_all * TNw_loading_all, axis=1) / Qw_loading  # mmol/m3, size = 365
NO3w_loading = np.sum(Qw_loading_all * NO3w_loading_all, axis=1) / Qw_loading  # mmol/m3, size = 365
NH4w_loading = np.sum(Qw_loading_all * NH4w_loading_all, axis=1) / Qw_loading  # mmol/m3, size = 365

# get loading in kg/day
# Loading
QrNO3r_loading_avg = (Qr_loading*NO3r_loading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QwNO3w_loading_avg = (Qw_loading*NO3w_loading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QrNH4r_loading_avg = (Qr_loading*NH4r_loading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QwNH4w_loading_avg = (Qw_loading*NH4w_loading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)


######################################################
# Get Ocean load from TEF

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t1_x11b'

section = 'ai' # Admiralty Inlet

# NO-LOADING RUN
in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1noDIN_x11b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31')
tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir,section)
# get inflowing values
Q_in = tef_df['q_p'] # Qin [m3/s]
NO3_in = tef_df['NO3_p'] # NO3in [mmol/m3]
NH4_in = tef_df['NH4_p'] # NH4in [mmol/m3]
# combine terms
DIN_in = NO3_in + NH4_in # DINin [mmol/m3]
# determine Qin*DINin
QinNO3in_mmol_s = Q_in * NO3_in # [mmol/m3]
QinNH4in_mmol_s = Q_in * NH4_in # [mmol/m3]
QinDINin_mmol_s = Q_in * DIN_in # [mmol/m3]
# convert to kg/d
QinNO3in_kg_d = QinNO3in_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QinNH4in_kg_d = QinNH4in_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QinDINin_kg_d = QinDINin_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)

# TEF dates
tef_start_date = '2017.01.02'
tef_end_date = '2017.12.30'
dates_daily_tef = pd.date_range(start=tef_start_date, end=tef_end_date, freq= 'd')

# pad ocean values for stacked time series
QinDINin_kg_d_padded = np.pad(QinDINin_kg_d,(732, 1097),constant_values=np.nan)

#######################################################
##               Plot time series                    ##
#######################################################

plt.close('all')
fig,ax = plt.subplots(4,1,figsize=(9,7),sharex=True)

# get colors
rivcolor = 'mediumpurple'
wwtpcolor = 'yellowgreen'
oceancolor = 'cornflowerblue'


# WWTP 
ax[0].plot(dates_local_daily,QwNO3w_loading_avg,linewidth=3,alpha=0.5,color=wwtpcolor,label='NO3')
ax[0].plot(dates_local_daily,QwNH4w_loading_avg,linewidth=1.5,color='olivedrab',label='NH4')
ax[0].legend(loc='upper right',fontsize=12)
ax[0].text(0.016,0.85,'(a) WWTP Loading',transform=ax[0].transAxes,fontsize=14,fontweight='bold')
ax[0].set_ylabel(r'Load [kg d$^{-1}$]',fontsize=12)

# River
ax[1].plot(dates_local_daily,QrNO3r_loading_avg,linewidth=3,alpha=0.5,color=rivcolor,label='NO3')
ax[1].plot(dates_local_daily,QrNH4r_loading_avg,linewidth=1.5,color='rebeccapurple',label='NH4')
ax[1].legend(loc='upper right',fontsize=12)
ax[1].text(0.016,0.85,'(b) River Loading',transform=ax[1].transAxes,fontsize=14,fontweight='bold')
ax[1].set_ylabel(r'Load [kg d$^{-1}$]',fontsize=12)

# Ocean
ax[2].plot(dates_daily_tef,QinNO3in_kg_d,linewidth=3,alpha=0.5,color=oceancolor,label='NO3')
ax[2].plot(dates_daily_tef,QinNH4in_kg_d,linewidth=1.5,color='royalblue',label='NH4')
ax[2].legend(loc='upper right',fontsize=12)
ax[2].text(0.016,0.85,'(c) Ocean Loading',transform=ax[2].transAxes,fontsize=14,fontweight='bold')
ax[2].set_ylabel(r'Load [kg d$^{-1}$]',fontsize=12)

# Stacked DIN
ax[3].stackplot(dates_local_daily,
             QwNO3w_loading_avg + QwNH4w_loading_avg, 
             QrNO3r_loading_avg + QrNH4r_loading_avg,
             QinDINin_kg_d_padded,
             labels=['WWTPs', 'Rivers', 'Ocean'],
             colors=[wwtpcolor, rivcolor, oceancolor],
             edgecolor='black',
             linewidth=0.4,
             alpha=0.5)
ax[3].text(0.016,0.85,'(d) Stacked DIN Loads',transform=ax[3].transAxes,fontsize=14,fontweight='bold')
ax[3].set_yscale('log') 
ax[3].legend(loc='upper right',fontsize=12)
ax[3].set_ylabel(r'Load [kg d$^{-1}$]',fontsize=12)

# format figure
ax[0].set_xlim([dates_local_daily[0],dates_local_daily[-1]])
ax[0].set_ylim([0,50000])
ax[1].set_ylim([0,200000])
ax[2].set_ylim([0,2500000])
ax[3].set_ylim([0,2500000])
ax[0].tick_params(axis='x', labelrotation=30)
ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[2].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[3].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[0].tick_params(axis='both', labelsize=12)
ax[1].tick_params(axis='both', labelsize=12)
ax[2].tick_params(axis='both', labelsize=12)
ax[3].tick_params(axis='both', labelsize=12)

plt.tight_layout()
plt.show()