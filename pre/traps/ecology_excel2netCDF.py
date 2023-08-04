"""
This script compiles all of Ecology's excel
loading data into two xarray files:
one for point sources
one for rivers

In theory, this script only needs to be run once.
Then, the xarray files can be referenced to generate climatologies.

Takes about 5 minutes to run on my local machine.

To run from ipython:
run ecology_excel2xarray.py
"""

#################################################################################
#                              Import packages                                  #
#################################################################################
from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import numpy as np
import os
import datetime
import xarray as xr

#################################################################################
#                              Helper functions                                 #
#################################################################################

def monthly2daily(df):
    '''
    turn a monthly dataframe into daily data (for the lenght of Ecology's timeseries)
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
    daily_df.reset_index(inplace=True)
    return daily_df

#################################################################################
#                              Get path to data                                 #
#################################################################################

# location of historical data to process
wwtp_dir = Ldir['data'] / 'traps' / 'point_sources'
wwtp_fns = os.listdir(wwtp_dir)
NWWTP = np.shape(wwtp_fns)[0]

# location of historical data to process
riv_dir = Ldir['data'] / 'traps' / 'nonpoint_sources'
riv_fns = os.listdir(riv_dir)
NTRIV = np.shape(riv_fns)[0]

# SSM metadata with lat/lon coordinates
trapsll_fn = Ldir['data'] / 'traps' / 'SSM_source_info.xlsx'
latlon_df = pd.read_excel(trapsll_fn,usecols='D,E,F,G,N,O')

#################################################################################
#                         Create point source dataset                           #
#################################################################################

# Start with one point source to get date information
# note that need to use river data because river data is monthly, wwtp is only monthly
riv_fp = str(riv_dir)  + '/' + riv_fns[0]
wwtp_df_example = pd.read_excel(riv_fp, skiprows=[0])
numdates = len(wwtp_df_example['Date'])

# Start Dataset (with empty data)
date = wwtp_df_example['Date']
source = np.arange(1,NWWTP+1)
pointsource_ds = xr.Dataset(data_vars=dict(ID=(['source'], np.ones((NWWTP,), dtype=int)),
        lon=(['source'], np.ones((NWWTP,))),
        lat=(['source'], np.ones((NWWTP,))),
        name=(['source'], ['placeholder placeholder']*NWWTP),
        flow=(['source', 'date'], np.zeros((NWWTP, numdates))),
        temp=(['source', 'date'], np.zeros((NWWTP, numdates))),
        NO3=(['source', 'date'], np.zeros((NWWTP, numdates))),
        NH4=(['source', 'date'], np.zeros((NWWTP, numdates))),
        TIC=(['source', 'date'], np.zeros((NWWTP, numdates))),
        Talk=(['source', 'date'], np.zeros((NWWTP, numdates))),
        DO=(['source', 'date'], np.zeros((NWWTP, numdates))),),
    coords=dict(source=source, date=date,),
    attrs=dict(description='Ecology data for point sources.'),)

# Add dataset metadata
pointsource_ds['ID'].attrs['long_name'] = 'source ID used in Salish Sea Model'
pointsource_ds['lon'].attrs['long_name'] = 'point source longitude'
pointsource_ds['lat'].attrs['long_name'] = 'point source latitude'
pointsource_ds['flow'].attrs['long_name'] = 'discharge rate'
pointsource_ds['flow'].attrs['units'] = 'm3/s'
pointsource_ds['temp'].attrs['long_name'] = 'discharge temperature'
pointsource_ds['temp'].attrs['units'] = 'C'
pointsource_ds['NO3'].attrs['long_name'] = 'nitrate+nitrite concentration'
pointsource_ds['NO3'].attrs['units'] = 'mmol/m3'
pointsource_ds['NH4'].attrs['long_name'] = 'ammonium concentration'
pointsource_ds['NH4'].attrs['units'] = 'mmol/m3'
pointsource_ds['TIC'].attrs['long_name'] = 'total inorganic carbon'
pointsource_ds['TIC'].attrs['units'] = 'mmol/m3'
pointsource_ds['Talk'].attrs['long_name'] = 'total alkalinity'
pointsource_ds['Talk'].attrs['units'] = 'meq/m3'
pointsource_ds['DO'].attrs['long_name'] = 'dissolved oxygen concentration'
pointsource_ds['DO'].attrs['units'] = 'mmol/m3'

print('Looping through point source files...')

# Loop through all WWTPs and add data to dataset
for i,fn in enumerate(wwtp_fns):

    # get ID and source name
    source_ID = int(fn.split('_', 1)[0])
    source_name_xlsx = fn.split('_', 1)[1]
    source_name = source_name_xlsx.split('.', 1)[0]
    print('{}/99: {}'.format(i+1,source_name))

    # load data as a dataframe
    wwtp_fp = str(wwtp_dir)  + '/' + fn
    wwtp_monthly_df = pd.read_excel(wwtp_fp, skiprows=[0]) 
    
    # rename columns so that they are standardized
    # I have previously verified that Ecology's .xlsx files all have the same parameters
    wwtp_monthly_df = wwtp_monthly_df.set_axis(['Date', 'Year', 'Month', 'Day',
                            'Hour', 'Minute', 'Bin1', 'Flow(m3/s)',
                            'Temp(C)','Salt(ppt)','NH4(mg/L)',
                            'NO3+NO2(mg/L)', 'PO4(mg/L)', 'DO(mg/L)',
                            'pH', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
                            'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)',
                            'POCR(mg/L)', 'DOCS(mg/L)', 'DOCF(mg/L)',
                            'Diatoms', 'Dinoflag', 'Chl', 'DIC(mmol/m3)',
                            'Alk(mmol/m3)'], axis=1, inplace=False)
    
    # point source data is monthly. Convert to daily
    wwtp_df = monthly2daily(wwtp_monthly_df)

    # Add source ID and name
    pointsource_ds['ID'][i] = source_ID
    pointsource_ds['name'][i] = source_name

    # Add source lat/lon
    pointsource_ds['lat'][i] = latlon_df.loc[latlon_df['ID'] == source_ID, 'Lat'].values[0]
    pointsource_ds['lon'][i] = latlon_df.loc[latlon_df['ID'] == source_ID, 'Lon'].values[0]
    
    # Add physics and biology data to dataset
    pointsource_ds.flow[i,:] = wwtp_df['Flow(m3/s)']
    pointsource_ds.temp[i,:] = wwtp_df['Temp(C)']
    pointsource_ds.NO3[i,:]  = wwtp_df['NO3+NO2(mg/L)'] * 71.4 # convert to mmol/m3
    pointsource_ds.NH4[i,:]  = wwtp_df['NH4(mg/L)']     * 71.4 # convert to mmol/m3
    pointsource_ds.TIC[i,:]  = wwtp_df['DIC(mmol/m3)']
    pointsource_ds.Talk[i,:] = wwtp_df['Alk(mmol/m3)']
    pointsource_ds.DO[i,:]   = wwtp_df['DO(mg/L)']      * 31.26 # convert to mmol/m3

# save dataset as .nc file in LO_data
out_fn = '../../../LO_data/traps/all_point_source_data.nc'
pointsource_ds.to_netcdf(out_fn)
pointsource_ds.close()
print('Point sources complete --------------------------------------------\n')

#################################################################################
#                       Create nonpoint source dataset                          #
#################################################################################

# Start with one nonpoint source to get date information
riv_fp = str(riv_dir)  + '/' + riv_fns[0]
riv_df_example = pd.read_excel(riv_fp, skiprows=[0])
numdates = len(riv_df_example['Date'])

# Start Dataset (with empty data)
date = riv_df_example['Date']
source = np.arange(1,NTRIV+1)
nonpointsource_ds = xr.Dataset(data_vars=dict(ID=(['source'], np.ones((NTRIV,), dtype=int)),
        lon=(['source'], np.ones((NTRIV,))),
        lat=(['source'], np.ones((NTRIV,))),
        name=(['source'], ['placeholder placeholder']*NTRIV),
        flow=(['source', 'date'], np.zeros((NTRIV, numdates))),
        temp=(['source', 'date'], np.zeros((NTRIV, numdates))),
        NO3=(['source', 'date'], np.zeros((NTRIV, numdates))),
        NH4=(['source', 'date'], np.zeros((NTRIV, numdates))),
        TIC=(['source', 'date'], np.zeros((NTRIV, numdates))),
        Talk=(['source', 'date'], np.zeros((NTRIV, numdates))),
        DO=(['source', 'date'], np.zeros((NTRIV, numdates))),),
    coords=dict(source=source, date=date,),
    attrs=dict(description='Ecology data for nonpoint sources.'),)

# Add dataset metadata
nonpointsource_ds['ID'].attrs['long_name'] = 'source ID used in Salish Sea Model'
nonpointsource_ds['lon'].attrs['long_name'] = 'point source longitude'
nonpointsource_ds['lat'].attrs['long_name'] = 'point source latitude'
nonpointsource_ds['flow'].attrs['long_name'] = 'discharge rate'
nonpointsource_ds['flow'].attrs['units'] = 'm3/s'
nonpointsource_ds['temp'].attrs['long_name'] = 'discharge temperature'
nonpointsource_ds['temp'].attrs['units'] = 'C'
nonpointsource_ds['NO3'].attrs['long_name'] = 'nitrate+nitrite concentration'
nonpointsource_ds['NO3'].attrs['units'] = 'mmol/m3'
nonpointsource_ds['NH4'].attrs['long_name'] = 'ammonium concentration'
nonpointsource_ds['NH4'].attrs['units'] = 'mmol/m3'
nonpointsource_ds['TIC'].attrs['long_name'] = 'total inorganic carbon'
nonpointsource_ds['TIC'].attrs['units'] = 'mmol/m3'
nonpointsource_ds['Talk'].attrs['long_name'] = 'total alkalinity'
nonpointsource_ds['Talk'].attrs['units'] = 'meq/m3'
nonpointsource_ds['DO'].attrs['long_name'] = 'dissolved oxygen concentration'
nonpointsource_ds['DO'].attrs['units'] = 'mmol/m3'

print('Looping through nonpoint source files...')

# Loop through all rivers and add data to dataset
for i,fn in enumerate(riv_fns):

    # get ID and source name
    source_ID = int(fn.split('_', 1)[0])
    source_name_xlsx = fn.split('_', 1)[1]
    source_name = source_name_xlsx.split('.', 1)[0]
    print('{}/161: {}'.format(i+1,source_name))

    # load data as a dataframe
    riv_fp = str(riv_dir)  + '/' + fn
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


    # Add source ID and name
    nonpointsource_ds['ID'][i] = source_ID
    nonpointsource_ds['name'][i] = source_name

    # Add source lat/lon
    # NOTE: we take the mean here because in SSM, large rivers are spread across two grid cell
    #       meaning that they have two lat/lon coordinates. We average to consolidate into one
    #       lat/lon coordinate. For rivers that are already in a single grid cell, the average of
    #       itself is itself.
    nonpointsource_ds['lat'][i] = np.mean(latlon_df.loc[latlon_df['ID'] == source_ID, 'Lat'].values)
    nonpointsource_ds['lon'][i] = np.mean(latlon_df.loc[latlon_df['ID'] == source_ID, 'Lat'].values)
    
    # Add physics and biology data to dataset
    nonpointsource_ds.flow[i,:] = riv_df['Flow(m3/s)']
    nonpointsource_ds.temp[i,:] = riv_df['Temp(C)']
    nonpointsource_ds.NO3[i,:]  = riv_df['NO3+NO2(mg/L)'] * 71.4 # convert to mmol/m3
    nonpointsource_ds.NH4[i,:]  = riv_df['NH4(mg/L)']     * 71.4 # convert to mmol/m3
    nonpointsource_ds.TIC[i,:]  = riv_df['DIC(mmol/m3)']
    nonpointsource_ds.Talk[i,:] = riv_df['Alk(mmol/m3)']
    nonpointsource_ds.DO[i,:]   = riv_df['DO(mg/L)']      * 31.26 # convert to mmol/m3

# save dataset as .nc file in LO_data
out_fn = '../../../LO_data/traps/all_nonpoint_source_data.nc'
nonpointsource_ds.to_netcdf(out_fn)
nonpointsource_ds.close()
print('Nonpoint sources complete ---------------------------------------')

print(pointsource_ds)
print(nonpointsource_ds)