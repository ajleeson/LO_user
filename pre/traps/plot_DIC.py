"""
Plot DIC, alkalinity, and pH for tiny rivers.

Based on Ecology's timeseries, stored in LO_data/traps

To run, from ipython:
run plot_DIC.py
"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import matplotlib.dates as mdates

# define year range to create climatologies
year0 = 1999
year1 = 2017

# location to save files
clim_dir = Ldir['LOo'] / 'pre' / 'traps' / 'tiny_rivers' /'Data_historical'

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

# # just Union River for now -------------------------------------------------
# rivnames = rivnames[90:91]
# rivids = rivids[90:91]

# # just Columbia River for now -------------------------------------------------
# rivnames = rivnames[93:94]
# rivids = rivids[93:94]

# # just Tsitsika River for now -------------------------------------------------
# rivnames = rivnames[29:30]
# rivids = rivids[29:30]

# initialize dataframes for all rivers
flow_clim_df = pd.DataFrame()
temp_clim_df = pd.DataFrame()
NO3_clim_df  = pd.DataFrame()
NH4_clim_df  = pd.DataFrame()
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

    # calculate averages
    # (compress 1999-2017 timeseries to single day, with an average for each day)
    riv_avgs_df = riv_df.groupby(['Month','Day']).mean().reset_index()
    # calculate standard deviation
    riv_sds_df = riv_df.groupby(['Month','Day']).std().reset_index()

    # replace all nans with zeros, so I'm no longer injecting nans
    riv_avgs_df = riv_avgs_df.replace(np.nan,0)
    riv_sds_df = riv_sds_df.replace(np.nan,0)

    # Set any negative TIC concentrations to zero
    riv_avgs_df['DIC(mmol/m3)'] = riv_avgs_df['DIC(mmol/m3)'].mask(riv_avgs_df['DIC(mmol/m3)'].lt(0),0)

    # Plot and save averages for each source
    vns = ['DIC(mmol/m3)','pH','Alk(mmol/m3)']
    fig, axes = plt.subplots(4,1, figsize=(10, 9), sharex=True)
    ax = axes.ravel()
    # create one-year date range for plotting
    yrday = pd.date_range(start ='1/1/2020', end ='12/31/2020', freq ='D')
    for j,vn in enumerate(vns):
        i = j+1
        # label subplot
        ax[i].set_title(vn,fontsize=14)
        # Plot individual years
        for yr in range(1999,2017):
            riv_yr_df = riv_df.loc[riv_df['Year'] == yr]
            # Insert a nan on Feb 29 if not a leap year
            if np.mod(yr,4) != 0:
                nans = [np.nan]*29
                riv_yr_df = riv_yr_df.reset_index(drop=True) # reset all dataframes to index from 0
                riv_yr_df.loc[58.5] = nans # leap year is 60th Julian day, so add a new 59th index since Python indexes from 0
                riv_yr_df = riv_yr_df.sort_index().reset_index(drop=True) # sort indices and renumber
            if yr == 2017:
                yrday_17 = pd.date_range(start ='1/1/2020', end ='8/02/2020', freq ='D') # don't have full 2017 dataset
                ax[i].plot(yrday_17,riv_yr_df[vn],alpha=0.5, label=yr, linewidth=1)
            else:
                ax[i].plot(yrday,riv_yr_df[vn],alpha=0.5, label=yr, linewidth=1)
        # fontsize of tick labels
        ax[i].tick_params(axis='both', which='major', labelsize=12)
        ax[i].tick_params(axis='x', which='major', rotation=30)
        ax[i].set_xlim([datetime.date(2020, 1, 1), datetime.date(2020, 12, 31)])
        # create legend
        if i == 3:
            handles, labels = ax[3].get_legend_handles_labels()
            ax[0].legend(handles, labels, loc='center', ncol = 6,fontsize=14)
            ax[0].axis('off')
        # Define the date format
        if i >= 3:
            date_form = mdates.DateFormatter("%b")
            ax[i].xaxis.set_major_formatter(date_form)
    # plot title is name of source
    ax[0].set_title(rname,fontsize=18)
    # Save figure
    figname = rname + '.png'
    save_path = clim_dir / 'DIC_plots' / figname
    fig.savefig(save_path)
    plt.close('all')
    # plt.show()