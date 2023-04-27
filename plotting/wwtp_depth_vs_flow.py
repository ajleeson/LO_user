"""
Code to plot water depth vs. average discharge for WWTPs
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as col
import pandas as pd
import math
import xarray as xr
import cmocean
from lo_tools import plotting_functions as pfun

# Get LiveOcean depths from grid info -----------------------------------------
ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc')
z = -ds.h.values

# Get wwtp average discharge rate
# get flowrate
fp_wwtps = '../../LO_output/pre/traps/point_sources/Data_historical/'
flowdf = pd.read_pickle(fp_wwtps+'CLIM_flow_1999_2017.p')    # m3/s
# calculate average flowrate
flow = flowdf.mean(axis=0).to_frame(name='avg-flow(m3/s)')
flow.index.name= 'rname'

# Get wwtp coordinates
fn_rowcol = '../../LO_data/grids/cas6/wwtp_info.csv'
df = pd.read_csv(fn_rowcol)
df = df.set_index('rname')

# add wwtp depth to the dataframe using apply function
df['depth'] = df.apply(lambda row: z[int(row.row_py),int(row.col_py)], axis = 1)

# add average flowrate to df
result = pd.concat([df, flow], axis=1)
print(result)

# plot
plt.close('all')
fig,ax = plt.subplots(figsize=(8, 6))
ax.scatter(result['avg-flow(m3/s)'],result['depth'],color='lightcoral',alpha=0.6,edgecolors='none',s=75)
ax.scatter(result.loc['Birch Bay']['avg-flow(m3/s)'],result.loc['Birch Bay']['depth'],
            edgecolors='k',s=75,marker='D',color='darkorchid', alpha=0.8, label='Birch Bay WWTP')
ax.scatter(result.loc['Oak Harbor Lagoon']['avg-flow(m3/s)'],result.loc['Oak Harbor Lagoon']['depth'],
            edgecolors='k',s=85,marker='^',color='yellowgreen', alpha=0.8, label='Oak Harbor Lagoon WWTP')
# ax.set_yscale('log')
ax.set_xscale('log')
ax.grid(True)
ax.set_xlabel('Average Flowrate (m$^3$ s$^{-1}$)', fontsize = 14)
ax.set_ylabel('Depth (m)', fontsize = 14)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.set_title('Point Source Depth vs. Discharge', fontsize = 16)
ax.legend(loc='best',fontsize=14)
ax.set_ylim([-200,0])
ax.set_xlim([1e-4,1e1])
plt.show()


#------------------------------------
# Calculate LwSrc flow statistics
total_flow = flowdf.sum(axis=1)
print('Annual Average of All Vertical Sources = {} m3/s'.format(round(total_flow.mean(),1)))
print('Annual Range of All Vertical Sources = {} m3/s'.format(round(total_flow.max()-total_flow.min(),1)))
print('Annual Max of All Vertical Sources = {} m3/s'.format(round(total_flow.max(),1)))
print('Annual Min of All Vertical Sources = {} m3/s'.format(round(total_flow.min(),1)))