"""
Make climatologies for tiny rivers.
Discharge rate, temperature, and biogeochemisty variables.

Based on Ecology's timeseries, stored in LO_data/traps

This code shows how powerful pandas is for this kind of task.
Really just one line to make a climatology (the groupby call)
"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import matplotlib.dates as mdates

# Suppress warnings
pd.options.mode.chained_assignment = None  # default='warn'

# define year range to create climatologies
year0 = 1999
year1 = 2017

# define gridname
gridname = 'cas6'

# file with all traps names and ID numbers
traps_info_fn = Ldir['data'] / 'traps' / 'SSM_source_info.xlsx'
# location of historical data to process
riv_dir = Ldir['data'] / 'traps' / 'nonpoint_sources'
all_his_fns = os.listdir(riv_dir)

# Get all river names and river IDs
traps_info_df = pd.read_excel(traps_info_fn,usecols='D,E,F')
# Only interested in rivers
riv_all_df = traps_info_df.loc[traps_info_df['Inflow_Typ'] == 'River']
# remove duplicate river names (some are listed in two rows)
riv_singles_df = riv_all_df.loc[traps_info_df['Name'].str.contains('- 2') == False]
# rename rivers that have a ' - 1' at the end
riv_names_df = riv_singles_df['Name'].str.replace(' - 1', '')
# get river names and river ids
rivnames = riv_names_df.values
rivids = riv_singles_df['ID'].values

# # Union River -------------------------------------------------
# rivnames = rivnames[90:91]
# rivids = rivids[90:91]

# # Skokomish River -------------------------------------------------
# rivnames = rivnames[87:88]
# rivids = rivids[87:88]

# # Tahuya River -------------------------------------------------
# rivnames = rivnames[81:82]
# rivids = rivids[81:82]

# Cushman 2 River -------------------------------------------------
rivnames = rivnames[76:77]
rivids = rivids[76:77]

# initialize dataframes for all rivers
flow_clim_df = pd.DataFrame()
salt_clim_df = pd.DataFrame()
temp_clim_df = pd.DataFrame()
NO3_clim_df  = pd.DataFrame()
NH4_clim_df  = pd.DataFrame()
# phyto_clim_df = pd.DataFrame()
# chlo_clim_df = pd.DataFrame()
TIC_clim_df  = pd.DataFrame()
Talk_clim_df = pd.DataFrame()
DO_clim_df   = pd.DataFrame()

# loop through all rivers
for i,rname in enumerate(rivnames):

    print('{}: {}'.format(i,rname))

    # get river index
    rID = rivids[i]
    
    riv_fn = ''

    # find Ecology's timeseries based on river id
    for fn in all_his_fns:
        root, ext = os.path.splitext(fn)
        if root.startswith(str(rID)) and ext == '.xlsx':
            riv_fn = fn
    
    # Let user know if couldn't find timeseries for a given river
    if riv_fn == '0':
        print('No history file found for {}'.format(rname))
        continue
    # Otherwise, read history file as df
    else:
        riv_fp = str(riv_dir)  + '/' + riv_fn
        riv_df = pd.read_excel(riv_fp, skiprows=[0])

    # rename columns so that they are standardized
    # I have previously verified that Ecology's .xlsx files all have the same parameters
    riv_df = riv_df.set_axis(['Date', 'Year', 'Month', 'Day',
                            'Hour', 'Minute', 'Bin1', 'Flow(m3/s)',
                            'Temp(C)','Salt(ppt)','NH4(mg/L)',
                            'NO3+NO2(mg/L)', 'PO4(mg/L)', 'DO(mg/L)',
                            'pH', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
                            'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)',
                            'POCR(mg/L)', 'DOCS(mg/L)', 'DOCF(mg/L)',
                            'Diatoms', 'Dinoflag', 'Chl', 'DIC(mmol/m3)',
                            'Alk(mmol/m3)'], axis=1, inplace=False)

    # replace all zeros with nans, so zeros don't bias data
    riv_df = riv_df.replace(0, np.nan)

    # calculate averages (compress 1999-2017 timeseries to single day, with an average for each day)
    riv_avgs_df = riv_df.groupby(['Month','Day']).mean().reset_index()
    riv_avgs_df.index = riv_avgs_df.index + 1 # start day index from 1 instead of 0

    # Plot averages (this was written to test one river at a time, so only pass one river through for-loop)
    plotting = True
    vn = 'Flow(m3/s)'
    style=['-',':','-',':']
    colors=['lightcoral','indigo','lightcoral','indigo']
    if plotting == True:
        print('plotting...')
        fig, ax = plt.subplots(2,1, figsize=(11, 6), sharex = True)
        yrday = np.linspace(1,367,366)
        # Plot individual years
        for i,yr in enumerate([2009,2013,2010,2014]):
            riv_yr_df = riv_df.loc[riv_df['Year'] == yr]
            # Insert a nan on Feb 29 if not a leap year
            if np.mod(yr,4) != 0:
                # print('{} is not a leap year'.format(yr)) # debugging
                nans = [np.nan]*29
                riv_yr_df = riv_yr_df.reset_index(drop=True) # reset all dataframes to index from 0
                riv_yr_df.loc[58.5] = nans # leap year is 60th Julian day, so add a new 59th index since Python indexes from 0
                riv_yr_df = riv_yr_df.sort_index().reset_index(drop=True) # sort indices and renumber
                # print(riv_yr_df[58:61]) # debugging
            if yr in [2010,2014]:
                if yr == 2010:
                    ax[1].plot(yrday,riv_yr_df[vn], label=yr, linestyle=style[i], color=colors[i])
                else:
                    ax[1].plot(yrday+91,riv_yr_df[vn], label='2014 shifted', linestyle=style[i], color=colors[i])
            else:
                if yr == 2009:
                    ax[0].plot(yrday,riv_yr_df[vn], label=yr, linestyle=style[i], color=colors[i])
                else:
                    ax[0].plot(yrday+91,riv_yr_df[vn], label='2013 shifted', linestyle=style[i], color=colors[i])
        ax[0].legend(loc='best')
        ax[1].legend(loc='best')
        ax[0].set_ylabel(vn)
        ax[1].set_ylabel(vn)
        ax[1].set_xlabel('Julian Day')
        fig.suptitle('{} Flow'.format(rname))
        plt.show()
