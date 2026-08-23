

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

# No-loading
Qr_noloading_all = ds_noloading_rivs['transport'].values # [m3/s] size = 365, nrivs
TNr_noloading_all = (ds_noloading_rivs['NO3'].values +
                     ds_noloading_rivs['NH4'].values +
                     ds_noloading_rivs['LDeN'].values +
                     ds_noloading_rivs['SDeN'].values +
                     ds_noloading_rivs['Phyt'].values +
                     ds_noloading_rivs['Zoop'].values)  # mmol/m3
NO3r_noloading_all = ds_noloading_rivs['NO3'].values # mmol/m3
NH4r_noloading_all = ds_noloading_rivs['NH4'].values # mmol/m3
Qr_noloading = np.sum(Qr_noloading_all, axis=1)  # [m3/s], size = 365
TNr_noloading = np.sum(Qr_noloading_all * TNr_noloading_all, axis=1) / Qr_noloading  # mmol/m3, size = 365
NO3r_noloading = np.sum(Qr_noloading_all * NO3r_noloading_all, axis=1) / Qr_noloading  # mmol/m3, size = 365
NH4r_noloading = np.sum(Qr_noloading_all * NH4r_noloading_all, axis=1) / Qr_noloading  # mmol/m3, size = 365

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

# No-loading
Qw_noloading_all = ds_noloading_wwtps['transport'].values # [m3/s] size = 365, nrivs
TNw_noloading_all = (ds_noloading_wwtps['NO3'].values +
                     ds_noloading_wwtps['NH4'].values +
                     ds_noloading_wwtps['LDeN'].values +
                     ds_noloading_wwtps['SDeN'].values +
                     ds_noloading_wwtps['Phyt'].values +
                     ds_noloading_wwtps['Zoop'].values)  # mmol/m3
NO3w_noloading_all = ds_noloading_wwtps['NO3'].values # mmol/m3
NH4w_noloading_all = ds_noloading_wwtps['NH4'].values # mmol/m3
Qw_noloading = np.sum(Qw_noloading_all, axis=1)  # [m3/s], size = 365
TNw_noloading = np.sum(Qw_noloading_all * TNw_noloading_all, axis=1) / Qw_noloading  # mmol/m3, size = 365
NO3w_noloading = np.sum(Qw_noloading_all * NO3w_noloading_all, axis=1) / Qw_noloading  # mmol/m3, size = 365
NH4w_noloading = np.sum(Qw_noloading_all * NH4w_noloading_all, axis=1) / Qw_noloading  # mmol/m3, size = 365

# # get averages
# # Loading
# Qr_loading_avg = np.nanmean(Qr_loading)
# TNr_loading_avg = np.nanmean(TNr_loading)
# Qw_loading_avg = np.nanmean(Qw_loading)
# TNw_loading_avg = np.nanmean(TNw_loading)
# # No-loading
# Qr_noloading_avg = np.nanmean(Qr_noloading)
# TNr_noloading_avg = np.nanmean(TNr_noloading)
# Qw_noloading_avg = np.nanmean(Qw_noloading)
# TNw_noloading_avg = np.nanmean(TNw_noloading)

# print('Loading averages')
# print('River flow: ', Qr_loading_avg, 'm3/s')
# print('River TN: ', TNr_loading_avg, 'mmol/m3')
# print('WWTP flow: ', Qw_loading_avg, 'm3/s')
# print('WWTP TN: ', TNw_loading_avg, 'mmol/m3')
# print('--------------')
# print('No-loading averages')
# print('River flow: ', Qr_noloading_avg, 'm3/s')
# print('River TN: ', TNr_noloading_avg, 'mmol/m3')
# print('WWTP flow: ', Qw_noloading_avg, 'm3/s')
# print('WWTP TN: ', TNw_noloading_avg, 'mmol/m3')

# CALCULATE DO LOSS OVER THREE MOTHS, ASSUMING ALL WWTP NH4 GETS NITRIFIED
# print average WWTP NH4 loading
QwNH4w = np.nanmean(Qw_loading*NH4w_loading) # [mmol/s]
print('Avg. WWTP NH4 load: '+str(QwNH4w)+' mmol/s')
# convert to mmol/day
QwNH4w_mmol_day = QwNH4w * 86400 # [mmol/day]
# calculate load over three months, assuming three months = 91 days
QwNH4w_3months_mmol = QwNH4w_mmol_day * 91 # [mmol over three months]
# get corresponding amount of DO loss, where 2 moles of oxygen are used per mole of ammonium oxidized
DO_nitrification_mmol = QwNH4w_3months_mmol * 2 # [mmol over three months]
# Get nominal Puget Sound volume (non time-varying)
basin_mask_ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / 'basin_masks_from_pugetsoundDObox.nc')
mask_ps = basin_mask_ds.mask_pugetsound.values
fp = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
box_ds = xr.open_dataset(fp)
DX = (box_ds.pm.values)**-1
DY = (box_ds.pn.values)**-1
DA = DX*DY # grid cell area [m2]
PS_vol = np.nansum(basin_mask_ds['h'].values * DA * mask_ps) # [m^3]
print('Puget Sound volume: {} m3'.format(round(PS_vol,1)))
# Divide this amount of DO by the volume of Puget Sound to get DO concentration loss over three months
DO_loss_to_nitrification_mmol_m3 = DO_nitrification_mmol / (PS_vol) # [mmol/m3]
# convert to mg/L
DO_loss_to_nitrification_mg_L = DO_loss_to_nitrification_mmol_m3 * 32 / 1000 # [mg/L]
print('DO loss to nitrification over three months: {} mg/L'.format(round(DO_loss_to_nitrification_mg_L,3)))

# get averages
# Loading
QrTNr_loading_avg = np.nanmean(Qr_loading*TNr_loading)
QwTNw_loading_avg = np.nanmean(Qw_loading*TNw_loading)
QrNO3r_loading_avg = np.nanmean(Qr_loading*NO3r_loading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QwNO3w_loading_avg = np.nanmean(Qw_loading*NO3w_loading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QrNH4r_loading_avg = np.nanmean(Qr_loading*NH4r_loading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QwNH4w_loading_avg = np.nanmean(Qw_loading*NH4w_loading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
# No-loading
QrTNr_noloading_avg = np.nanmean(Qr_noloading*TNr_noloading)
QwTNw_noloading_avg = np.nanmean(Qw_noloading*TNw_noloading)
QrNO3r_noloading_avg = np.nanmean(Qr_noloading*NO3r_noloading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QwNO3w_noloading_avg = np.nanmean(Qw_noloading*NO3w_noloading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QrNH4r_noloading_avg = np.nanmean(Qr_noloading*NH4r_noloading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QwNH4w_noloading_avg = np.nanmean(Qw_noloading*NH4w_noloading) / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)

# print('Loading averages')
# print('River QrTNr: ', QrTNr_loading_avg, 'mmol/s')
# print('WWTP QwTNw: ', QwTNw_loading_avg, 'mmol/s')
# print('--------------')
# print('No-loading averages')
# print('River QrTNr: ', QrTNr_noloading_avg, 'mmol/s')
# print('WWTP QwTNw: ', QwTNw_noloading_avg, 'mmol/s')

# print WWTP nutrient concentrations
print('WWTP NH4 concentration: {} mmol/m3'.format(np.nanmean(NH4w_loading)))
print('WWTP NO3 concentration: {} mmol/m3'.format(np.nanmean(NO3w_loading)))


#######################################################
##                Stacked bar chart                  ##
#######################################################

plt.close('all')
fig,ax = plt.subplots(1,1,figsize=(6,4))

run = ('No-loading','Loading')
x = np.arange(len(run))

# Get loads
loads = {'River NO3': np.array([QrNO3r_noloading_avg,QrNO3r_loading_avg]),
        'River NH4': np.array([QrNH4r_noloading_avg,QrNH4r_loading_avg]),
        'WWTP NO3': np.array([QwNO3w_noloading_avg,QwNO3w_loading_avg]),
        'WWTP NH4': np.array([QwNH4w_noloading_avg,QwNH4w_loading_avg])}
# create colors
colors = {'River NO3': 'mediumpurple',
          'River NH4': 'mediumpurple',
          'WWTP NO3':  'yellowgreen',
          'WWTP NH4':  'yellowgreen'}
# alphas
alphas = {'River NO3': 0.4,
          'River NH4': 1,
          'WWTP NO3':  0.4,
          'WWTP NH4':  1}
# bar style
styles = {'River NO3': 'xx',
          'River NH4': '',
          'WWTP NO3':  'xx',
          'WWTP NH4':  ''}

width = 0.8

bottom = np.zeros(len(run))  # running stack base

# for load_type, values in loads.items():
#     ax.bar(x, values, width=width, bottom=bottom, label=load_type, color=colors[load_type],
#            edgecolor=colors[load_type], alpha=alphas[load_type],hatch=styles[load_type])
#     bottom += values

no_idx = run.index('No-loading')
for load_type, values in loads.items():
    for i, xi in enumerate(x):
        ec = 'none' if ('WWTP' in load_type and i == no_idx and np.isclose(values[i], 0)) else colors[load_type]
        ax.bar(xi, values[i], width=width, bottom=bottom[i], label=load_type if i == 0 else None, color=colors[load_type],
               alpha=alphas[load_type], hatch=styles[load_type], edgecolor=ec)
    bottom += values

ax.set_xticks(x)
ax.set_xticklabels(run, fontsize=14)
ax.tick_params(axis='y', labelsize=14)
ax.set_ylabel('Load (kg d$^{-1}$)', fontsize=14)
ax.set_title('Land-based loads to Puget Sound\n(mean of 2015-2020)',fontsize=16)
handles = [mpatches.Patch(facecolor=colors[k], edgecolor=colors[k], alpha=alphas[k], hatch=styles[k], label=k) for k in loads.keys()]
ax.legend(handles=handles, loc='upper left', fontsize=14)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.show()