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
j2 = 1170
i1 = 220
i2 = 652

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = '../../LO_data/trapsD00/LiveOcean_SSM_rivers.xlsx'
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
yd_ind = pd.Index(dt_ind.dayofyear)

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../LO_data/grids/cas7/grid.nc')
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
fp_wwtps = '../../LO_output/pre/trapsP00/point_sources/lo_base/Data_historical/'
flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow.p')    # m3/s
no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3.p')      # mmol/m3
nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4.p')      # mmol/m3

# calculate total DIN concentration in mg/L
dindf_wwtps = (no3df_wwtps + nh4df_wwtps)/71.4    # mg/L

# calculate daily loading timeseries in kg/d
dailyloaddf_wwtps = 86.4*dindf_wwtps*flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s

# calculate average daily load over the year (kg/d)
avgload_wwtps = dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

# add row and col index for plotting on LiveOcean grid
griddf0_wwtps = pd.read_csv('../../LO_data/grids/cas7/wwtp_info.csv')
griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
avgload_wwtps = avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
avgload_wwtps = avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

# get point source lat and lon
lon_wwtps = [X[int(col)] for col in avgload_wwtps['col_py']]
lat_wwtps = [Y[int(row)] for row in avgload_wwtps['row_py']]

# calculate total load
totload_wwtps = np.sum(avgload_wwtps['avg-daily-load(kg/d)'])

# Rivers -------------------------------------------------------------
# Prepare data for spatial summary plots

# get flow, nitrate, and ammonium values
# fp_trivs = '../../LO_output/pre/traps/all_rivers/Data_historical/'
fp_trivs = '../../LO_output/pre/trapsP00/tiny_rivers/lo_base/Data_historical/'
fp_LOrivs = '../../LO_output/pre/river1/lo_base/Data_historical/'
fp_LObio = '../../LO_output/pre/trapsP00/LO_rivbio/lo_base/Data_historical/'
flowdf_rivs = pd.read_pickle(fp_trivs+'CLIM_flow.p')    # m3/s
no3df_rivs = pd.read_pickle(fp_trivs+'CLIM_NO3.p')      # mmol/m3
nh4df_rivs = pd.read_pickle(fp_trivs+'CLIM_NH4.p')      # mmol/m3
# pre-existing LO river flowrates
flowdf_LOrivs = pd.read_pickle(fp_LOrivs+'CLIM_flow.p')    # m3/s
flowdf_LOrivs = flowdf_LOrivs.reset_index(drop=True)
# Ecology data for pre-existing LO rivers
no3df_LOrivs_ecol = pd.read_pickle(fp_LObio+'CLIM_NO3.p')      # mmol/m3
nh4df_LOrivs_ecol = pd.read_pickle(fp_LObio+'CLIM_NH4.p')      # mmol/m3
# get names of all pre-existing rivers for which Ecology has data (use LO naming convention)
LObio_names = [SSM2LO_name(riv) for riv in no3df_LOrivs_ecol.columns]

# get biology data for pre-existing rivers
no3df_LOrivs = pd.DataFrame()
nh4df_LOrivs = pd.DataFrame()
for rn in flowdf_LOrivs:
    # Use ecology data if there exists any
    if rn in LObio_names:
        no3df_LOrivs[rn] = no3df_LOrivs_ecol[LO2SSM_name(rn)]
        nh4df_LOrivs[rn] = nh4df_LOrivs_ecol[LO2SSM_name(rn)]
    # Otherwise use the rivfun method to guess
    else:
        no3df_LOrivs[rn] = get_bio_vec('NO3', rn, yd_ind)
        nh4df_LOrivs[rn] = get_bio_vec('NH4', rn, yd_ind)

# calculate total DIN concentration in mg/L
dindf_rivs = (no3df_rivs + nh4df_rivs)/71.4          # mg/L
dindf_LOrivs = (no3df_LOrivs + nh4df_LOrivs)/71.4    # mg/L

# calculate daily loading timeseries in kg/d
dailyloaddf_rivs = 86.4*dindf_rivs*flowdf_rivs       # kg/d = 86.4 * mg/L * m3/s
dailyloaddf_LOrivs = 86.4*dindf_LOrivs*flowdf_LOrivs # kg/d = 86.4 * mg/L * m3/s

# calculate average daily load over the year (kg/d)
avgload_trivs = dailyloaddf_rivs.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')
avgload_LOrivs = dailyloaddf_LOrivs.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

# add row and col index for plotting on LiveOcean grid (tiny rivers)
griddf0_trivs = pd.read_csv('../../LO_data/grids/cas7/triv_info.csv')
griddf_trivs = griddf0_trivs.set_index('rname') # use river name as index
avgload_trivs = avgload_trivs.join(griddf_trivs['row_py']) # add row to avg load df (uses rname to index)
avgload_trivs = avgload_trivs.join(griddf_trivs['col_py']) # do the same for cols

# add row and col index for plotting on LiveOcean grid (pre-existing rivers)
griddf0_LOrivs = pd.read_csv('../../LO_data/grids/cas7/river_info.csv')
griddf_LOrivs = griddf0_LOrivs.set_index('rname') # use river name as index
avgload_LOrivs = avgload_LOrivs.join(griddf_LOrivs['row_py']) # add row to avg load df (uses rname to index)
avgload_LOrivs = avgload_LOrivs.join(griddf_LOrivs['col_py']) # do the same for cols

# get average load of all rivers
avgload_allrivs = pd.concat([avgload_trivs, avgload_LOrivs])
# drop nans
avgload_allrivs = avgload_allrivs.dropna()

# calculate total load
totload_rivs = avgload_allrivs.loc[(avgload_allrivs['col_py'] >= i1) &
                                   (avgload_allrivs['col_py'] <= i2) &
                                   (avgload_allrivs['row_py'] >= j1) &
                                   (avgload_allrivs['row_py'] <= j2),
                                   'avg-daily-load(kg/d)'].sum()
# totload_rivs = np.sum(avgload_rivs['avg-daily-load(kg/d)'])

# get trivs lat and lon
lon_riv = [X[int(col)] for col in avgload_allrivs['col_py']]
lat_riv = [Y[int(row)] for row in avgload_allrivs['row_py']]

# Ocean -------------------------------------------------------------
# Data from Mackas et al. (1997)
avgload_SJdF = 2600 * 1000 # (tonnes per day * 1000 kg/ton)

# Just picking arbitrary point in SJdF
lat_SJdF = 48.4
lon_SJdF = -124.4

# Just picking arbitrary point in SoG

# PLOTTING ----------------------------------------------------

# define marker sizes (minimum size is 10 so dots don't get too small)
sizes_wwtps = [max(0.1*load,10) for load in avgload_wwtps['avg-daily-load(kg/d)']]
sizes_rivs = [max(0.1*load,10) for load in avgload_allrivs['avg-daily-load(kg/d)']]
sizes_SJdF = 0.1*avgload_SJdF
sizes = [sizes_wwtps,sizes_rivs,sizes_SJdF]

# pick colors
color_wwtps = '#AEDC3C' # wwtp green
color_rivs = '#7148BC' # river purple
color_SJdF = '#70B0EA' # ocean blue
colors = [color_wwtps, color_rivs,color_SJdF]

# define labels
source_name = ['WWTP','River','Ocean']

# define total loads
totalloads = [totload_wwtps,totload_rivs,avgload_SJdF]

# define percentages
netload = totload_wwtps+totload_rivs+avgload_SJdF
percent_wwtp = round(totload_wwtps/netload*100,1)
percent_rivs = round(totload_rivs/netload*100,1)
percent_ocean = round(avgload_SJdF/netload*100,1)
percentages = [percent_wwtp,percent_rivs,percent_ocean]

# define subplot number
subplotnums = [131,132,133]

# panel labels
letters = ['(a)','(b)','(c)']

# lat and lon coords
lats = [lat_wwtps,lat_riv,lat_SJdF]
lons = [lon_wwtps,lon_riv,lon_SJdF]

add_ocean = False

for j in range(1):
    # if j == 1:
    #     add_ocean = True

    fig = plt.figure(figsize=(21,8))
    plt.tight_layout()
    # loop through all of the plots we need to make
    for i,sname in enumerate(source_name):

        # add water/land
        ax = fig.add_subplot(subplotnums[i])
        # pfun.add_coast(ax,color='black')
        ax.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

        # plot DIN sources
        ax.scatter(lons[i],lats[i],
            color=colors[i], edgecolors='k', alpha=0.5, s=sizes[i])

        # ax.set_xlim(-123.5,-122) # Puget Sound
        # ax.set_ylim(46.7,49.3) # Puget Sound
        # ax.set_xlim(-125.5,-121.5) # Salish Sea
        # ax.set_ylim(46.7,49.9) # Salish Sea
        # ax.set_xlim(-130,-121.5) # Full Grid
        # ax.set_ylim(42,52) # Full Grid

        ax.set_xlim(X[i1],-121.4)#X[i2]) # Salish Sea
        ax.set_ylim(Y[j1],Y[j2]) # Salish Sea

        # format
        ax.axes.xaxis.set_visible(False)
        ax.axes.yaxis.set_visible(False)
        ax.set_title(letters[i]+' Average {} DIN Load'.format(sname), fontsize=20)
        pfun.dar(ax)

        # Create legend
        if i == 0:
            leg_szs = [100, 1000, 10000]
            # szs = leg_szs
            szs = [0.1*(leg_sz) for leg_sz in leg_szs]
            # szs = [10*np.sqrt(leg_sz) for leg_sz in leg_szs]
            l0 = plt.scatter([],[], s=szs[0], color='grey', edgecolors='k', alpha=0.5)
            l1 = plt.scatter([],[], s=szs[1], color='grey', edgecolors='k', alpha=0.5)
            l2 = plt.scatter([],[], s=szs[2], color='grey', edgecolors='k', alpha=0.5)
            l3 = plt.scatter([],[], s=0, color='grey', edgecolors='k', alpha=0.5)
            labels = ['< 100', '1,000', '10,000', r'(kg d$^{-1}$)']
            legend = ax.legend([l0, l1, l2, l3], labels, fontsize = 14, markerfirst=False,
                title=r'Loading $\propto$ area', loc='lower left', labelspacing=1.5, borderpad=0.5)
            plt.setp(legend.get_title(),fontsize=16)

        # add label of total load
        if sname == 'Ocean':
            tload = ax.text(0.40,0.87, '{:,}'.format(int(round(totalloads[i],-3))) + r'$ \ kg \ d^{-1}$' + '\n(Mackas & Harrison 1997)',
                horizontalalignment = 'center', fontsize = 16, color = 'k', transform=ax.transAxes)
        else:
            tload = ax.text(0.77,0.87, r'$\Sigma$ ' + sname + 's: \n {:,}'.format(int(round(totalloads[i],-3))) + r'$ \ kg \ d^{-1}$',
                    horizontalalignment = 'center', fontsize = 16, color = 'k', transform=ax.transAxes)
        legcol = colors[i]
        alpha = 0.3
        tload.set_bbox(dict(facecolor=legcol, alpha=alpha, edgecolor='none', boxstyle = 'Round'))

        # add percentages
        ax.text(0.82,0.04,'{}%'.format(percentages[i]),fontsize=18,color='k',transform=ax.transAxes)

        # reduce gap between subplots
        plt.subplots_adjust(wspace=0.05)
        
        plt.savefig(out_dir / ('loading_plot.png'))


# Print stuff for comparison with Ben

print(avgload_wwtps.loc['West Point'])
print('---------------------------')
print(avgload_allrivs.loc['fraser'])
print('---------------------------')
print(avgload_allrivs.loc['puyallup'])