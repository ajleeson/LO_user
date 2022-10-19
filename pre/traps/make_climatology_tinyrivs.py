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
year0 = 1980
year1 = 2020

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
riv_singles_df = riv_all_df[traps_info_df['Name'].str.contains('- 2') == False]
# rename rivers that have a ' - 1' at the end
riv_names_df = riv_singles_df['Name'].str.replace(' - 1', '')
# get river names and river ids
rivnames = riv_names_df.values
rivids = riv_singles_df['ID'].values

# just test two rivers for now (delete all of this later)
rivnames = rivnames[90:91]
rivids = rivids[90:91]

# initialize dataframes for all rivers
flow_clim_df = pd.DataFrame()

# loop through all rivers
for i,rname in enumerate(rivnames):
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
    riv_avgs_df = riv_df.groupby(['Month','Day']).mean()

    # plot to visualize
    plotting = True
    vn = 'DO(mg/L)'
    if plotting == True:
        fig, ax = plt.subplots(1,1, figsize=(11, 6))
        yrday = np.linspace(1,367,366)
        # Plot individual years
        for yr in range(1999,2018):
            riv_yr_df = riv_df.loc[riv_df['Year'] == yr]
            # Insert a nan on Feb 29 if not a leap year
            if np.mod(yr,4) != 0:
                # print('{} is not a leap year'.format(yr)) # debugging
                nans = [np.nan]*29
                riv_yr_df = riv_yr_df.reset_index(drop=True) # reset all dataframes to index from 0
                riv_yr_df.loc[58.5] = nans # leap year is 60th Julian day, so add a new 59th index since Python indexes from 0
                riv_yr_df = riv_yr_df.sort_index().reset_index(drop=True) # sort indices and renumber
                # print(riv_yr_df[58:61]) # debugging
            if yr == 2017:
                yrday_17 = np.linspace(1,214,213) # don't have full 2017 dataset
                plt.plot(yrday_17,riv_yr_df[vn],alpha=0.5, label=yr)
            else:
                plt.plot(yrday,riv_yr_df[vn],alpha=0.5, label=yr)
        # Plot average
        plt.plot(yrday,riv_avgs_df[vn].values, label='average', color='black', linewidth=2)
        plt.legend(loc='best', ncol = 4)
        plt.ylabel(vn)
        plt.xlabel('Julian Day')
        plt.title(rname)
        plt.show()

# flow_df = pd.read_pickle(riv_dir / ('ALL_flow_' + str(year0) + '_' + str(year1) + '.p'))
# temp_df = pd.read_pickle(riv_dir / ('ALL_temperature_' + str(year0) + '_' + str(year1) + '.p'))

# # make the climatologies (one line!)
# flow_clim_df = flow_df.groupby(flow_df.index.dayofyear).mean()
# temp_clim_df = temp_df.groupby(flow_df.index.dayofyear).mean()

# # drop temperature rivers with more than 50 missing yeardays
# temp_clim_df = temp_clim_df.loc[:, pd.isnull(temp_clim_df).sum() < 50]
# # fill missing temperature values
# temp_clim_df = temp_clim_df.interpolate()
# # make yearday 366 equal to yearday 365 (leap year is poorly sampled)
# temp_clim_df.loc[366,:] = temp_clim_df.loc[365,:]

# # check for missing values:
# if pd.isnull(flow_clim_df).sum().sum() != 0:
#     print('Warning, there are missing flow values!')
# if pd.isnull(temp_clim_df).sum().sum() != 0:
#     print('Warning, there are missing temperature values!')

# # save results
# flow_clim_df.to_pickle(riv_dir / ('CLIM_flow_' + str(year0) + '_' + str(year1) + '.p'))
# temp_clim_df.to_pickle(riv_dir / ('CLIM_temp_' + str(year0) + '_' + str(year1) + '.p'))

# # Plotting
# plt.close('all')

# fig = plt.figure(figsize=(16,10))
# rn_split = np.array_split(flow_clim_df.columns, 9)
# for ii in range(1,10):
#     ax = fig.add_subplot(3,3,ii)
#     flow_clim_df[rn_split[ii-1]].plot(ax=ax)
#     ax.set_xlim(0,366)
#     ax.set_ylim(bottom=0)
#     if ii >= 7:
#         ax.set_xlabel('Yearday')
#     if ii in [1, 4, 7]:
#         ax.set_ylabel(r'Flow [$m^{3}s^{-1}$]')
# fig.savefig(riv_dir / ('CLIM_flow_plot.png'))

# fig = plt.figure(figsize=(16,10))
# rn_split = np.array_split(temp_clim_df.columns, 9)
# for ii in range(1,10):
#     ax = fig.add_subplot(3,3,ii)
#     temp_clim_df[rn_split[ii-1]].plot(ax=ax)
#     ax.set_xlim(0,366)
#     ax.set_ylim(0, 25)
#     if ii >= 7:
#         ax.set_xlabel('Yearday')
#     if ii in [1, 4, 7]:
#         ax.set_ylabel(r'Temperature [$^{\circ}C$]')
# fig.savefig(riv_dir / ('CLIM_temp_plot.png'))
    
# plt.show()