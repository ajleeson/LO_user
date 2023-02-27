"""
Code to plot distance from wwtp current location to nearest coastal cell vs. average discharge for WWTPs

This is to visualize how much things would change if we were to implement wwtps as tiny rivers.
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
from lo_tools import zfun


# get the grid data
ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc')
z = -ds.h.values
mask_rho = np.transpose(ds.mask_rho.values)
lon = ds.lon_rho.values
lat = ds.lat_rho.values
X = lon[0,:] # grid cell X values
Y = lat[:,0] # grid cell Y values

# Get wwtp average discharge rate
# get flowrate
fp_wwtps = '../../LO_output/pre/traps/point_sources/Data_historical/'
flowdf = pd.read_pickle(fp_wwtps+'CLIM_flow_1999_2017.p')    # m3/s
# calculate average flowrate
df = flowdf.mean(axis=0).to_frame(name='avg-flow(m3/s)')
df.index.name= 'rname'

# Get original wwtp coordinates
fn_rowcol = '../../LO_data/grids/cas6/wwtp_info.csv'
rowcol_wwtp_df = pd.read_csv(fn_rowcol)
rowcol_wwtp_df = rowcol_wwtp_df.set_index('rname')

# get new wwtp cooridnates
fn_rowcol_triv = '../../LO_data/grids/cas6/wwtp_as_triv_info.csv'
rowcol_coast_df = pd.read_csv(fn_rowcol_triv)
rowcol_coast_df = rowcol_coast_df.set_index('rname')

# get distance from original wwtp to nearest coastal cell
df['dist2coast(km)'] = ''
for source in df.index:
    # get lat and lon coordinates of original wwtp location
    lon0 = X[int(rowcol_wwtp_df.loc[source]['col_py'])]
    lat0 = Y[int(rowcol_wwtp_df.loc[source]['row_py'])]
    # get lat and lon coordinates of coastal wwtp location
    lon = X[int(rowcol_coast_df.loc[source]['col_py'])]
    lat = Y[int(rowcol_coast_df.loc[source]['row_py'])]
    # get distance from original wwtp to nearest coastal cell
    xmeters, ymeters = zfun.ll2xy(lon, lat, lon0, lat0)
    distance = np.sqrt(xmeters**2 + ymeters**2)/1000
    df.loc[source,'dist2coast(km)'] = distance

print(df.sort_values(by=['avg-flow(m3/s)']))

# plot
plt.close('all')
fig,ax = plt.subplots(figsize=(8, 6))
ax.scatter(df['avg-flow(m3/s)'],df['dist2coast(km)'],color='lightcoral',alpha=0.6,edgecolors='none',s=75)
ax.scatter(df.loc['South King']['avg-flow(m3/s)'],df.loc['South King']['dist2coast(km)'],
            edgecolors='k',s=80,marker='D',color='darkorchid', alpha=0.8, label='South King')
ax.scatter(df.loc['West Point']['avg-flow(m3/s)'],df.loc['West Point']['dist2coast(km)'],
            edgecolors='k',s=140,marker='^',color='yellowgreen', alpha=0.8, label='West Point')
ax.scatter(df.loc['Annacis']['avg-flow(m3/s)'],df.loc['Annacis']['dist2coast(km)'],
            edgecolors='k',s=190,marker='*',color='deeppink', alpha=0.8, label='Annacis')
# ax.set_yscale('log')
# ax.set_xscale('log')
ax.grid(True)
ax.set_xlabel('Average Flowrate (m$^3$ s$^{-1}$)', fontsize = 14)
ax.set_ylabel('Distance to Coast (km)', fontsize = 14)
ax.tick_params(axis='both', which='major', labelsize=12)
ax.set_title('Distance from Point Source to Nearest Coastal Cell \n vs. Discharge', fontsize = 16)
ax.legend(loc='best',fontsize=14)
# ax.set_ylim([-200,0])
# ax.set_xlim([1e-4,1e1])
plt.show()
