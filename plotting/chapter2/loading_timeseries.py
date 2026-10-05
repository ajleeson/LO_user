

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

Ldir = Lfun.Lstart()


#################################################################################
#                                   Get data                                    #
#################################################################################

# start_date = '2017.01.01'
# end_date = '2017.12.31'
start_date = '2015.01.01'
end_date = '2020.12.31'

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


############# ##########################################
##               Plot time series                    ##
#######################################################

plt.close('all')
fig,ax = plt.subplots(2,1,figsize=(10,6),sharex=True)

# River
ax[0].plot(dates_local_daily,QrNO3r_loading_avg)
ax[0].plot(dates_local_daily,QrNH4r_loading_avg)

# WWTP
ax[1].plot(dates_local_daily,QwNO3w_loading_avg)
ax[1].plot(dates_local_daily,QwNH4w_loading_avg)

# format figure
ax[0].set_xlim([dates_local_daily[0],dates_local_daily[-1]])
ax[0].set_ylim([0,200000])
ax[1].set_ylim([0,40000])
ax[0].tick_params(axis='x', labelrotation=30)
ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[0].tick_params(axis='both', labelsize=12)
ax[1].tick_params(axis='both', labelsize=12)

plt.tight_layout()
plt.show()