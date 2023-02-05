"""
Code to plot DIN loads in Puget Sound
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

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = '../../LO_data/traps/LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_LO = repeatrivs_df.loc[repeatrivs_df['SSM_rname'] == rname, 'LO_rname'].values[0]

    return rname_LO

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc')
z = -ds.h.values
mask_rho = np.transpose(ds.mask_rho.values)
lon = ds.lon_rho.values
lat = ds.lat_rho.values
X = lon[0,:] # grid cell X values
Y = lat[:,0] # grid cell Y values
plon, plat = pfun.get_plon_plat(lon,lat)
# make a version of z with nans where masked
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
zm[np.transpose(mask_rho) != 0] = -1

# Point Sources -------------------------------------------------------------
# Prepare data for spatial summary plots

# get flow, nitrate, and ammonium values
fp_wwtps = '../../LO_output/pre/traps/point_sources/Data_historical/'
flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow_1999_2017.p')    # m3/s
no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3_1999_2017.p')      # mmol/m3
nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4_1999_2017.p')      # mmol/m3

# calculate total DIN concentration in mg/L
dindf_wwtps = (no3df_wwtps + nh4df_wwtps)/71.4    # mg/L

# calculate daily loading timeseries in kg/d
dailyloaddf_wwtps = 86.4*dindf_wwtps*flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s

# calculate average daily load over the year (kg/d)
avgload_wwtps = dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

# add row and col index for plotting on LiveOcean grid
griddf0_wwtps = pd.read_csv('../../LO_data/grids/cas6/wwtp_info.csv')
griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
avgload_wwtps = avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
avgload_wwtps = avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

# get point source lat and lon
lon_wwtps = [X[int(col)] for col in avgload_wwtps['col_py']]
lat_wwtps = [Y[int(row)] for row in avgload_wwtps['row_py']]

# Rivers -------------------------------------------------------------
# Prepare data for spatial summary plots

# get flow, nitrate, and ammonium values
fp_trivs = '../../LO_output/pre/traps/all_rivers/Data_historical/'
flowdf_rivs = pd.read_pickle(fp_trivs+'CLIM_flow_1999_2017.p')    # m3/s
no3df_rivs = pd.read_pickle(fp_trivs+'CLIM_NO3_1999_2017.p')      # mmol/m3
nh4df_rivs = pd.read_pickle(fp_trivs+'CLIM_NH4_1999_2017.p')      # mmol/m3

# calculate total DIN concentration in mg/L
dindf_rivs = (no3df_rivs + nh4df_rivs)/71.4    # mg/L

# calculate daily loading timeseries in kg/d
dailyloaddf_rivs = 86.4*dindf_rivs*flowdf_rivs # kg/d = 86.4 * mg/L * m3/s

# calculate average daily load over the year (kg/d)
avgload_rivs = dailyloaddf_rivs.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

# add row and col index for plotting on LiveOcean grid (tiny rivers)
griddf0_trivs = pd.read_csv('../../LO_data/grids/cas6/triv_info.csv')
griddf_trivs = griddf0_trivs.set_index('rname') # use river name as index
avgload_rivs = avgload_rivs.join(griddf_trivs['row_py']) # add row to avg load df (uses rname to index)
avgload_rivs = avgload_rivs.join(griddf_trivs['col_py']) # do the same for cols

# add row and col index for plotting on LiveOcean grid (pre-existing rivers)
griddf0_LOrivs = pd.read_csv('../../LO_data/grids/cas6/river_info.csv')
griddf_LOrivs = griddf0_LOrivs.set_index('rname') # use river name as index
# get naming conversion between Ecology and LiveOcean

# loop through and get river row and col for pre-existing LiveOcean rivers
for index, row in avgload_rivs.iterrows():
    # if missing row indices, then it is a pre-existing river
    if math.isnan(row['row_py']):
        # convert Ecology name to LO name
        avgload_rivs.loc[index]['row_py'] = griddf_LOrivs.loc[SSM2LO_name(index)]['row_py']
        avgload_rivs.loc[index]['col_py'] = griddf_LOrivs.loc[SSM2LO_name(index)]['col_py']

# get trivs lat and lon
lon_riv = [X[int(col)] for col in avgload_rivs['col_py']]
lat_riv = [Y[int(row)] for row in avgload_rivs['row_py']]

# Ocean -------------------------------------------------------------
# Data from Mackas et al. (1997)
avgload_SJdF = 2600 * 1000 # (tonnes per day * 1000 kg/ton)

# Just picking arbitrary point in SJdF
lat_SJdF = 48.4
lon_SJdF = -124.4

# Just picking arbitrary point in SoG

# PLOTTING ----------------------------------------------------

# define marker sizes
sizes_wwtps = [10*np.sqrt(load) for load in avgload_wwtps['avg-daily-load(kg/d)']]
sizes_rivs = [10*np.sqrt(load) for load in avgload_rivs['avg-daily-load(kg/d)']]
sizes_SJdF = 10*np.sqrt(avgload_SJdF)
sizes = [sizes_wwtps,sizes_rivs]

# pick colors
color_wwtps = '#AEDC3C' #'xkcd:pinkish red'
color_rivs = '#7148BC' #'mediumblue'
color_SJdF = '#70B0EA' #'gold'
colors = [color_wwtps, color_rivs]

# define labels
source_name = ['Point Source','River']

# lat and lon coords
lats = [lat_wwtps,lat_riv]
lons = [lon_wwtps,lon_riv]

add_ocean = False
# loop through all of the plots we need to make
for ii in range(3):

    if ii == 2:
        add_ocean = True
        i = 1
    else:
        i = ii

    # bathymetry
    fig = plt.figure(figsize=(8,9))
    ax = fig.add_subplot(111)
    # pfun.add_coast(ax,color='black')
    ax.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

    # plot DIN sources
    ax.scatter(lons[i],lats[i],
        color=colors[i], edgecolors='k', alpha=0.5, s=sizes[i])
    if add_ocean:
        ax.scatter(lon_SJdF,lat_SJdF,
            color=color_SJdF, edgecolors='k', alpha=0.5, s=sizes_SJdF)
        t = ax.text(lon_SJdF,lat_SJdF-0.7,'Ocean Load \n' r'$2.6\times10^6 \ kg \ d^{-1}$',
            horizontalalignment = 'center', fontsize = 16, color = 'k')
        t.set_bbox(dict(facecolor='white', alpha=0.6, edgecolor='none', boxstyle = 'Round'))

    # ax.set_xlim(-123.5,-122) # Puget Sound
    # ax.set_ylim(46.7,49.3) # Puget Sound
    ax.set_xlim(-125.5,-121.5) # Salish Sea
    ax.set_ylim(46.7,49.9) # Salish Sea
    # ax.set_xlim(-130,-121.5) # Full Grid
    # ax.set_ylim(42,52) # Full Grid

    # format
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.set_title('Average {} DIN Load'.format(source_name[i]), fontsize=20)
    pfun.dar(ax)

    # Create legend
    if i == 0:
        leg_szs = [10, 100, 1000, 10000]
        szs = [10*np.sqrt(leg_sz) for leg_sz in leg_szs]
        l0 = plt.scatter([],[], s=szs[0], color='grey', edgecolors='k', alpha=0.5)
        l1 = plt.scatter([],[], s=szs[1], color='grey', edgecolors='k', alpha=0.5)
        l2 = plt.scatter([],[], s=szs[2], color='grey', edgecolors='k', alpha=0.5)
        l3 = plt.scatter([],[], s=szs[3], color='grey', edgecolors='k', alpha=0.5)
        labels = ['10', '100', '1000', '10000']
        legend = ax.legend([l0, l1, l2, l3], labels, fontsize = 14,
            title=r'Loading (kg d$^{-1}$)', loc='lower left', labelspacing=1.5, borderpad=1.5)
        plt.setp(legend.get_title(),fontsize=16)

    plt.show()