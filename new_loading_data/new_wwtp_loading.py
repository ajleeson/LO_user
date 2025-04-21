"""
Script to plot and compare new WWTP data and my older climatologies.
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import datetime
import matplotlib.dates as mdates
import datetime
import os
from pathlib import Path

#################################################################################
#                               WWTP information                                #
#################################################################################

old_WWTP_name = 'West Point'
new_WWTP_name = 'King County West Point WWTP'


# old_WWTP_name = 'Stanwood' # discharges into Stillaguamish
# new_WWTP_name = 'STANWOOD STP'
    

#################################################################################
#                          Get old WWTP climatology data                        #
#################################################################################

# get flow [m3/s]
fn_oldWWTP_flow = Ldir['LOo'] / 'pre' / 'trapsP00' / 'point_sources' / 'lo_base' / 'Data_historical' / 'CLIM_flow.p'
df_oldWWTP_flow = pd.read_pickle(fn_oldWWTP_flow)
# Get WWTP flow
old_flow = df_oldWWTP_flow[old_WWTP_name].values # [m3/s]

# get NO3 and NH4 [mmol/m3]
fn_oldWWTP_NO3 = Ldir['LOo'] / 'pre' / 'trapsP00' / 'point_sources' / 'lo_base' / 'Data_historical' / 'CLIM_NO3.p'
df_oldWWTP_NO3 = pd.read_pickle(fn_oldWWTP_NO3)
fn_oldWWTP_NH4 = Ldir['LOo'] / 'pre' / 'trapsP00' / 'point_sources' / 'lo_base' / 'Data_historical' / 'CLIM_NH4.p'
df_oldWWTP_NH4 = pd.read_pickle(fn_oldWWTP_NH4)
# Get WWTP NO3 and NH4
old_NO3 = df_oldWWTP_NO3[old_WWTP_name].values # [mmol/m3]
old_NH4 = df_oldWWTP_NH4[old_WWTP_name].values # [mmol/m3]
# sum nutrients
old_nutrients = old_NO3 + old_NH4 # [mmol/m3]

# get total nutrient loading [kg/day]
old_nutrient_load_mmol_s = old_flow * old_nutrients # [m3/s * mmol/m3 = mmol/s]
# convert to kg/day: [1g = 14.01/1000^2mmol] [1day = 60*60*24sec]
old_nutrient_load_kg_day = old_nutrient_load_mmol_s * 14.01 * 60 * 60 * 24 / 1000 / 1000

# get one value every month
old_nutrient_load_subsampled = old_nutrient_load_kg_day[15::30]

# repeat for 15 years
old_nutrient_load_15years = np.tile(old_nutrient_load_subsampled,16)

#################################################################################
#                               Get new WWTP data                               #
#################################################################################

# read csv files
df_fac_attr = pd.read_csv('fac_attributes.csv')
df_nut_loads_new = pd.read_csv('nutrient_loads.csv')

# get WWTP ID
fac_ID = df_fac_attr.loc[df_fac_attr['FAC_NAME'] == new_WWTP_name, 'FAC_ID'].values[0]

# get flow corresponding to desired WWTP
new_flow_millGall_day = df_nut_loads_new.loc[df_nut_loads_new['FAC_ID'] == fac_ID, 'FLOW_MGD'].values.astype(float)
# convert to m3/s
new_flow_m3_s = new_flow_millGall_day * 0.0438126364

# get nutrient data [mg/L]
new_NO3 = df_nut_loads_new.loc[df_nut_loads_new['FAC_ID'] == fac_ID, 'NO2NO3N_MG_L'].values.astype(float)
new_NH4 = df_nut_loads_new.loc[df_nut_loads_new['FAC_ID'] == fac_ID, 'NH4N_MG_L'].values.astype(float)
# sum NO3 and NH4
new_nutrients = new_NO3 + new_NH4 # [mg/L]

# get total nutrient loading [kg/day]
new_nutrient_load_m3mg_sL = new_flow_m3_s * new_nutrients # [m3/s * mg/L]
# convert to kg/day: [1 kg = 1/1000 m3mg/L] [1day = 60*60*24sec]
new_nutrient_load_kg_day = new_nutrient_load_m3mg_sL * 60 * 60 * 24 / 1000

# rename for consistency
new_nutrient_load_15years = new_nutrient_load_kg_day


#################################################################################
#                                     Plot                                      #
#################################################################################

# get time array
t = pd.date_range(start='2005-01-01', end='2020-12-31', freq='MS').to_list()

plt.close('all')
fig, ax = plt.subplots(1,1,figsize = (10,6))

# plot new data
ax.plot(t,new_nutrient_load_15years, label='New Data (2005-2020 data)',
        linewidth=3,color='orchid', alpha=0.6)

# plot WWTP climatologies
ax.plot(t,old_nutrient_load_15years, label='Climatology (1999-2017 data)',
        linewidth=2,color='teal', alpha=0.8)

# format figure
ax.set_title('{} nutrient load comparison'.format(old_WWTP_name),
             fontsize=14,fontweight='bold')
ax.set_ylabel('Nutrient load [kg/day]',fontsize=12)
ax.legend(loc='best', fontsize=12)

ax.set_ylim([0,np.nanmax([np.nanmax(new_nutrient_load_15years),np.nanmax(old_nutrient_load_15years)])*1.1])
ax.set_xlim([pd.to_datetime('2005-01-01'),pd.to_datetime('2020-12-31')])
ax.xaxis.set_major_locator(mdates.YearLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax.tick_params(axis='x', labelrotation=30)
ax.tick_params(axis='both', labelsize=12)
ax.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')

