"""
Create scatter plot of max daily load
vs mean yearly load for all point sources

from ipython:
run ranked_pointsource_discharge.py

"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import matplotlib.dates as mdates
import datetime


# helper function 
def monthly2daily(df):
    '''
    turn a monthly dataframe into daily data
    '''
    # duplicate last row
    double_lr_df = pd.concat([df, df.iloc[-1:]], ignore_index=True)
    # picking arbitrary year to fill with daily data (but includes a leap year)
    start_date = datetime.date(1999, 1, 1)
    end_date = datetime.date(2017, 8, 1)
    dates = pd.date_range(start_date, end_date, freq='MS')
    # Replace month column with things that look like dates
    double_lr_df['Month'] = dates
    double_lr_df = double_lr_df.set_index('Month')
    # Change monthly to daily
    double_lr_daily_df = double_lr_df.resample('D').ffill()
    # delete last row (1/1 on the next year)
    daily_df = double_lr_daily_df[:-1]
    # make index start from 1 and go to 366
    daily_i_df = daily_df.reset_index(inplace=True)
    return daily_df
    

# define year range to create climatologies
year0 = 1999
year1 = 2017

# location to save file
clim_dir = Ldir['LOo'] / 'pre' / 'traps' / 'point_sources' /'Data_historical'

# file with all traps names and ID numbers
traps_info_fn = Ldir['data'] / 'traps' / 'SSM_source_info.xlsx'
# location of historical data to process
wwtp_dir = Ldir['data'] / 'traps' / 'point_sources'
all_his_fns = os.listdir(wwtp_dir)

# Get all wwtp names and wwtp IDs
traps_info_df = pd.read_excel(traps_info_fn,usecols='D,E,F')
# Only interested in wwtps
wwtp_all_df = traps_info_df.loc[traps_info_df['Inflow_Typ'] == 'Point Source']
#get river names
wwtp_names_df = wwtp_all_df['Name'].str.replace(' - 1', '')
# get river names and river ids
wwtpnames = wwtp_names_df.values
wwtpids = wwtp_all_df['ID'].values

# # just test Brightwater for now -------------------------------------------------
# wwtpnames = wwtpnames[46:47]
# wwtpids = wwtpids[46:47]

# # just test Birch Bay for now -------------------------------------------------
# wwtpnames = wwtpnames[26:27]
# wwtpids = wwtpids[26:27]

# # just test 5 WWTPs for now -------------------------------------------------
# wwtpnames = wwtpnames[12:17]
# wwtpids = wwtpids[12:17]

# initialize dataframes to save results
loads_stats_df = pd.DataFrame(index=wwtpnames, columns=['Max Daily Load (kg/day)', 'Mean Yearly Load (kg/yr)'])

# loop through all rivers
for i,wname in enumerate(wwtpnames):

    print('{}: {}'.format(i,wname))

    # get river index
    wID = wwtpids[i]
    
    wwtp_fn = ''

    # find Ecology's timeseries based on wwtp id
    for fn in all_his_fns:
        root, ext = os.path.splitext(fn)
        if root.startswith(str(wID)) and ext == '.xlsx':
            wwtp_fn = fn
    
    # Let user know if couldn't find timeseries for a given point source
    if wwtp_fn == '0':
        print('No history file found for {}'.format(wname))
    # Otherwise, read history file as df
    else:
        wwtp_fp = str(wwtp_dir)  + '/' + wwtp_fn
        wwtp_df = pd.read_excel(wwtp_fp, skiprows=[0])

    # rename columns so that they are standardized
    # I have previously verified that Ecology's .xlsx files all have the same parameters
    wwtp_df = wwtp_df.set_axis(['Date', 'Year', 'Month', 'Day',
                            'Hour', 'Minute', 'Bin1', 'Flow(m3/s)',
                            'Temp(C)','Salt(ppt)','NH4(mg/L)',
                            'NO3+NO2(mg/L)', 'PO4(mg/L)', 'DO(mg/L)',
                            'pH', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
                            'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)',
                            'POCR(mg/L)', 'DOCS(mg/L)', 'DOCF(mg/L)',
                            'Diatoms', 'Dinoflag', 'Chl', 'DIC(mmol/m3)',
                            'Alk(mmol/m3)'], axis=1, inplace=False)

    # replace all zeros with nans, so zeros don't bias data
    wwtp_df = wwtp_df.replace(0, np.nan)

    # converty monthly values to daily values
    wwtp_df = monthly2daily(wwtp_df)

    # calculate daily loads
    flow_daily = wwtp_df['Flow(m3/s)'].values    # [m3/s]
    no3_daily  = wwtp_df['NO3+NO2(mg/L)'].values # [mg/L]
    nh4_daily  = wwtp_df['NH4(mg/L)'].values     # [mg/L]
    daily_load_arr = 86.4 * (no3_daily + nh4_daily) * flow_daily # kg/d = 86.4 * mg/L * m3/s

    # add daily loads to dataframe
    # daily_load_df = pd.concat([daily_load_df, pd.Series(daily_load_arr, name=wname)], axis = 1)  # [kg/day] on monthly frequency

    max_daily_load = np.nanmax(daily_load_arr) # [kg/day]
    mean_yearly_load = 365*np.nanmean(daily_load_arr) # [kg/yr] = (q_1(kg/d)*1d + ... q_n(kg/d)*1d) / (n days * 1yr/365days)

    # add max loads to dataframe
    # max_loads_df = pd.concat([max_loads_df,pd.Series(max_daily_load,name=wname)], axis = 1)
    loads_stats_df.loc[wname, 'Max Daily Load (kg/day)'] = max_daily_load
    # add mean loads to dataframe
    # mean_loads_df = pd.concat([mean_loads_df,pd.Series(mean_yearly_load,name=wname)], axis = 1)
    loads_stats_df.loc[wname, 'Mean Yearly Load (kg/yr)'] = mean_yearly_load

# replace residual nans with zeros
loads_stats_df = loads_stats_df.replace(np.nan,0)

# separate birch bay wwtp from the rest
birchbay_stats_df = loads_stats_df.loc['Birch Bay']
loads_stats_df = loads_stats_df.drop('Birch Bay')

# Get max and mean values for plotting
max_daily_all = loads_stats_df['Max Daily Load (kg/day)']
mean_yearly_all = loads_stats_df['Mean Yearly Load (kg/yr)']
max_daily_birchbay = birchbay_stats_df['Max Daily Load (kg/day)']
mean_yearly_birchbay = birchbay_stats_df['Mean Yearly Load (kg/yr)']

plt.close('all')
fig,ax = plt.subplots(figsize=(8, 6))
ax.scatter(mean_yearly_all,max_daily_all,alpha=0.6,edgecolors='none',s=75)
ax.scatter(mean_yearly_birchbay,max_daily_birchbay,edgecolors='none',s=75,marker='D',color='coral',label='Birch Bay WWTP')
ax.set_yscale('log')
ax.set_xscale('log')
ax.grid(True)
ax.set_ylabel('Max Daily Load (kg/day)', fontsize = 14)
ax.set_xlabel('Mean Yearly Load (kg/yr)', fontsize = 14)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.set_title('Point Source Max Daily vs. Mean Yearly DIN Load', fontsize = 16)
ax.legend(loc='best',fontsize=14)
plt.show()

