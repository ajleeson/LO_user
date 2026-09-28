

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

ds_loading   = xr.open_dataset('../../../../LO_output/chapter_2/data/loading_riverwwtp_extraction_'+start_date+'_'+end_date+'.nc')

# get list of rivers that discharge to Puget Sound
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

# get flow and concentration data
# RIVER -------------------------------------
Qr_all = ds_loading_rivs['transport'].values # [m3/s] size = time, nrivs
DO_all = ds_loading_rivs['Oxyg'].values * 32/1000 # mg/L

# get flow-weighted average DO concentration
avg_DO = np.sum(Qr_all * DO_all, axis=1) / np.sum(Qr_all, axis=1) # mg/L, size = 365
