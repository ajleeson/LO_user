"""
Code to plot DIN loads in Puget Sound
"""

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import xarray as xr
import cmocean
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun
Ldir = Lfun.Lstart()

# define grid indices to look at
j1 = 590
j2 = 1000# 1170 #1300
i1 = 400# 220 #160
i2 = 600# 652

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = Ldir['data'] / 'trapsD00' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_LO = repeatrivs_df.loc[repeatrivs_df['SSM_rname'] == rname, 'LO_rname'].values[0]

    return rname_LO

def LO2SSM_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = Ldir['data'] / 'trapsD00' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_SSM = repeatrivs_df.loc[repeatrivs_df['LO_rname'] == rname, 'SSM_rname'].values[0]

    return rname_SSM

# riv fun function to get biology
def get_bio_vec(vn, rn, yd_ind):
    ndt = len(yd_ind)
    yd = yd_ind.values
    ovec = np.ones(ndt)
    if vn == 'NO3':
        if rn == 'fraser':
            vv = 2 + (13/2) + (13/2)*np.cos(2*np.pi*((yd-30)/366))
        elif rn == 'columbia':
            vv = 5 + (35/2) + (35/2)*np.cos(2*np.pi*((yd)/366))
        else:
            vv = 5 * ovec
    elif vn == 'Oxyg':
        vv = 350 * ovec
    elif vn in ['TAlk', 'TIC']:
        if rn in ['columbia', 'deschutes', 'duwamish']:
            vv = 1000 * ovec
        else:
            vv = 300 * ovec
    else:
        vv = 0 * ovec # all others filled with zeros
    return vv

# set up the time index for the record
Ldir = Lfun.Lstart()
dsf = Ldir['ds_fmt']
dt0 = datetime.strptime('2020.01.01',dsf)
dt1 = datetime.strptime('2020.12.31',dsf)
days = (dt0, dt1)
    
# pandas Index objects
dt_ind = pd.date_range(start=dt0, end=dt1)
# we are essentially just getting indices from 1 - 366, to get the days of the year on a leap year
yd_ind = pd.Index(dt_ind.dayofyear)

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
fp = Ldir['data'] / 'grids' / 'cas7' / 'grid.nc'
ds = xr.open_dataset(fp)
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

# Rivers -------------------------------------------------------------
# Prepare data for spatial summary plots

# get flow, nitrate, and ammonium values
fp_trivs = Ldir['LOo'] / 'pre' / 'trapsP00' / 'tiny_rivers' / 'lo_base' / 'Data_historical'
fp_LOrivs = Ldir['LOo'] / 'pre'  / 'river1' / 'lo_base' / 'Data_historical'
fp_LObio = Ldir['LOo'] / 'pre' / 'trapsP00' / 'LO_rivbio' / 'lo_base' / 'Data_historical'
flowdf_rivs = pd.read_pickle(fp_trivs / 'CLIM_flow.p')    # m3/s
# pre-existing LO river flowrates
flowdf_LOrivs = pd.read_pickle(fp_LOrivs / 'CLIM_flow.p')    # m3/s

# calculate daily loading timeseries in kg/d
dailyflowdf_rivs = flowdf_rivs       # m3/s
dailyflowdf_LOrivs = flowdf_LOrivs   # m3/s

# calculate average daily load over the year (kg/d)
avgflow_trivs = dailyflowdf_rivs.mean(axis=0).to_frame(name='avg-daily-flow(m3/s)')
avgflow_LOrivs = dailyflowdf_LOrivs.mean(axis=0).to_frame(name='avg-daily-flow(m3/s)')
# avgflow_trivs = dailyflowdf_rivs.max(axis=0).to_frame(name='avg-daily-flow(m3/s)')
# avgflow_LOrivs = dailyflowdf_LOrivs.max(axis=0).to_frame(name='avg-daily-flow(m3/s)')

# add row and col index for plotting on LiveOcean grid (tiny rivers)
fp = Ldir['data'] / 'grids' / 'cas7' / 'triv_info.csv'
griddf0_trivs = pd.read_csv(fp)
griddf_trivs = griddf0_trivs.set_index('rname') # use river name as index
avgflow_trivs = avgflow_trivs.join(griddf_trivs['row_py']) # add row to avg load df (uses rname to index)
avgflow_trivs = avgflow_trivs.join(griddf_trivs['col_py']) # do the same for cols

# add row and col index for plotting on LiveOcean grid (pre-existing rivers)
fp = Ldir['data'] / 'grids' / 'cas7' / 'river_info.csv'
griddf0_LOrivs = pd.read_csv(fp)
griddf_LOrivs = griddf0_LOrivs.set_index('rname') # use river name as index
avgflow_LOrivs = avgflow_LOrivs.join(griddf_LOrivs['row_py']) # add row to avg load df (uses rname to index)
avgflow_LOrivs = avgflow_LOrivs.join(griddf_LOrivs['col_py']) # do the same for cols

# get average load of all rivers
avgflow_allrivs = pd.concat([avgflow_trivs, avgflow_LOrivs])
# drop nans
avgflow_allrivs = avgflow_allrivs.dropna()

# calculate total load
# sum up river load for all circles that are visible on the map
totload_rivs = avgflow_allrivs.loc[(avgflow_allrivs['col_py'] >= i1) &
                                   (avgflow_allrivs['col_py'] <= i2) &
                                   (avgflow_allrivs['row_py'] >= j1) &
                                   (avgflow_allrivs['row_py'] <= j2),
                                   'avg-daily-flow(m3/s)'].sum()
# totload_rivs = np.sum(avgload_allrivs['avg-daily-load(kg/d)']) # uncomment to sum up all rivers in model domain

# get trivs lat and lon
lon_riv = [X[int(col)] for col in avgflow_allrivs['col_py']]
lat_riv = [Y[int(row)] for row in avgflow_allrivs['row_py']]


# PLOTTING ----------------------------------------------------
plt.close('all')

# define marker sizes (minimum size is 10 so dots don't get too small)
sizes_rivs = avgflow_allrivs['avg-daily-flow(m3/s)'] #[max(4*flow,4*1) for flow in avgflow_allrivs['avg-daily-flow(m3/s)']]
sizes = [sizes_rivs]

source_name = ['River']

# pick colors
color_rivs = '#7148BC' # river purple
colors = [color_rivs]

# define subplot number
subplotnums = [111]

# lat and lon coords
lats = [lat_riv]
lons = [lon_riv]

for j in range(1):

    fig = plt.figure(figsize=(8,8))
    plt.tight_layout()
    # loop through all of the plots we need to make
    for i,sname in enumerate(source_name):

        # add water/land
        ax = fig.add_subplot(subplotnums[i])
        ax.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

        # plot DIN sources
        ax.scatter(lons[i],lats[i],
            color=colors[i], edgecolors='k', alpha=0.5, s=sizes[i])
        
        # add max river flow
        for j in range(len(lon_riv)):
            if sizes[i][j] > 5:
                ax.text(lons[i][j],lats[i][j]-0.02, str(round(sizes[i][j],1)),horizontalalignment='center')

        ax.set_xlim(X[i1],-121.4)#X[i2]) # Salish Sea
        ax.set_ylim(Y[j1],Y[j2]) # Salish Sea

        # format
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        ax.set_title('Mean {} Annual Flow'.format(sname), fontsize=20, loc='left')
        pfun.dar(ax)

        # Create legend
        if i == 0:
            leg_szs = [1, 10, 100]
            # szs = leg_szs
            szs = [4*(leg_sz) for leg_sz in leg_szs]
            # szs = [10*np.sqrt(leg_sz) for leg_sz in leg_szs]
            l0 = plt.scatter([],[], s=szs[0], color='grey', edgecolors='k', alpha=0.5)
            l1 = plt.scatter([],[], s=szs[1], color='grey', edgecolors='k', alpha=0.5)
            l2 = plt.scatter([],[], s=szs[2], color='grey', edgecolors='k', alpha=0.5)
            l3 = plt.scatter([],[], s=0, color='grey', edgecolors='k', alpha=0.5)
            labels = ['< 1', '10', '100', r'(m $^3$ s$^{-1}$)']
            legend = ax.legend([l0, l1, l2, l3], labels, fontsize = 14, markerfirst=False,
                title=r'Flow $\propto$ area', loc='lower right', labelspacing=1.5, borderpad=0.5)
            plt.setp(legend.get_title(),fontsize=16)

        # add label of total load
        tload = ax.text(0.77,0.87, r'$\Sigma$ ' + sname + 's: \n {:,}'.format(int(round(totload_rivs))) + r'$ \ m^3 \ s^{-1}$',
                horizontalalignment = 'center', fontsize = 16, color = 'k', transform=ax.transAxes)
        legcol = colors[i]
        alpha = 0.3
        tload.set_bbox(dict(facecolor=legcol, alpha=alpha, edgecolor='none', boxstyle = 'Round'))


        # reduce gap between subplots
        plt.subplots_adjust(wspace=0.05,left=0.05,right=0.95)
        
        plt.savefig(out_dir / ('river_flow_map.png'),transparent='True',bbox_inches='tight')


# # Print stuff for comparison with Ben

# print(avgload_wwtps.loc['West Point'])
# print('---------------------------')
# print(avgload_allrivs.loc['fraser'])
# print('---------------------------')
# print(avgload_allrivs.loc['puyallup'])

