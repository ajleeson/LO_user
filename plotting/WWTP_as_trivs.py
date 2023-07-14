"""
Plots a comparison of bottom DO for two identical grids with different oxygen value.

This is a custom function for a particular experiment, but code can be adapted for other use cases in the future.

From the terminal: python bottom_DO_nutrient_comparison.py

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import matplotlib.patches as mpatches
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

# ----------------------------------------------------------------------

# # Provide information about models to compare
# # basecase (has nutrients)
# bc_gtagex = 'cas6_traps2_x2b'
# # test condition (no nutrients)
# c1_gtagex = 'cas6_traps3_x2b'
# gtagexes = [bc_gtagex, c1_gtagex]

# hr = '0025'
# date = '2017.03.07'

#'Puget Sound','Whidbey Basin','North King County','Lynch Cove'
region = 'North King County'

# # Variables to compare
# vn_list = ['oxygen']#['NO3','phytoplankton','oxygen']

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# get plotting limits based on region
if region == 'North King County':
    xmin = -122.55
    xmax = -122.35
    ymin = 47.55
    ymax = 47.75

# Get grid data
G = zrfun.get_basic_info('/home/aleeson/LO_data/grids/cas6/grid.nc', only_G=True)
grid_ds = xr.open_dataset('/home/aleeson/LO_data/grids/cas6/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values
X = lon[0,:] # grid cell X values
Y = lat[:,0] # grid cell Y values
px, py = pfun.get_plon_plat(G['lon_rho'],G['lat_rho'])
z = grid_ds.h.values
mask_rho = np.transpose(grid_ds.mask_rho.values)
# get indices of min and max
imin_x = find_nearest(px[0,:],xmin)
imax_x = find_nearest(px[0,:],xmax)
imin_y = find_nearest(py[:,0],ymin)
imax_y = find_nearest(py[:,0],ymax)

# Initialize figure
fs = 10
pfun.start_plot(fs=fs, figsize=(8,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
# make a version of z with nans where masked
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
zm[np.transpose(mask_rho) != 0] = -1
ax.pcolormesh(px, py, zm, edgecolor='powderblue', linewidth=0.5, vmin=-1000, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])
pfun.add_coast(ax)
pfun.dar(ax)
ax.set_title('Locations of WWTPs as Tiny Rivers')
ax.set_ylabel('Lat')
ax.set_xlabel('Lon')

# add locations as wwtps
wwtplabel=True
og_loc_df = pd.read_csv('../../LO_data/grids/cas6/wwtp_info.csv')
for rn in og_loc_df.index:
    ii = int(og_loc_df.loc[rn,'col_py'])
    jj = int(og_loc_df.loc[rn,'row_py'])
    if wwtplabel==True:
        ax.scatter(lon[jj,ii], lat[jj,ii],marker='o', s=60, color='#AEDC3C', edgecolor = 'k',label='Original WWTP Location')
        wwtplabel=False
    else:
        ax.scatter(lon[jj,ii], lat[jj,ii],marker='o', s=60, color='#AEDC3C', edgecolor = 'k')

# ax.scatter(lon_wwtps,lat_wwtps,color='lawngreen', edgecolors='green', alpha=0.7)

# add locations as tiny rivers
triv_loc_df = pd.read_csv('../../LO_data/grids/cas6/wwtp_as_triv_info.csv')
new_loc_df = pd.DataFrame(columns=['rname','col_py', 'row_py'])
rivlabel=True
for rn in triv_loc_df.index:
# u or v grids.
    ii = int(triv_loc_df.loc[rn,'col_py'])
    jj = int(triv_loc_df.loc[rn,'row_py'])
    uv = triv_loc_df.loc[rn,'uv']
    isign = triv_loc_df.loc[rn,'isign']
    idir = triv_loc_df.loc[rn,'idir']
    if uv == 'u' and isign == 1:
        new_loc_df.loc[len(new_loc_df.index)] = [triv_loc_df.loc[rn,'rname'], ii+1, jj] 

    if uv == 'u' and isign == -1:
        new_loc_df.loc[len(new_loc_df.index)] = [triv_loc_df.loc[rn,'rname'], ii, jj]

    if uv == 'v' and isign == 1:
        new_loc_df.loc[len(new_loc_df.index)] = [triv_loc_df.loc[rn,'rname'], ii, jj+1]
                
    if uv == 'v' and isign == -1:
        new_loc_df.loc[len(new_loc_df.index)] = [triv_loc_df.loc[rn,'rname'], ii, jj]

for rn in new_loc_df.index:
    ii = int(new_loc_df.loc[rn,'col_py'])
    jj = int(new_loc_df.loc[rn,'row_py'])
    if rivlabel==True:
        ax.scatter(lon[jj,ii], lat[jj,ii],marker='*', s=60, color='#7148BC', label='Moved Location')
        rivlabel=False
    else:
        ax.scatter(lon[jj,ii], lat[jj,ii],marker='*', s=60, color='#7148BC')

plt.legend(loc='upper right')

# ==========================================================================================

# Assess changes

# add wwtp depth to the dataframe using apply function
og_loc_df['depth'] = og_loc_df.apply(lambda row: z[int(row.row_py),int(row.col_py)], axis = 1)
new_loc_df['depth'] = new_loc_df.apply(lambda row: z[int(row.row_py),int(row.col_py)], axis = 1)

# Get WWTP loading size
# set up the time index for the record
Ldir = Lfun.Lstart()
dsf = Ldir['ds_fmt']
dt0 = datetime.strptime('2020.01.01',dsf)
dt1 = datetime.strptime('2020.12.31',dsf)
days = (dt0, dt1)
    
# pandas Index objects
dt_ind = pd.date_range(start=dt0, end=dt1)
yd_ind = pd.Index(dt_ind.dayofyear)

# # Get LiveOcean grid info --------------------------------------------------

# # get the grid data
# ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc')
# z = -ds.h.values
# mask_rho = np.transpose(ds.mask_rho.values)
# lon = ds.lon_rho.values
# lat = ds.lat_rho.values
# X = lon[0,:] # grid cell X values
# Y = lat[:,0] # grid cell Y values
# plon, plat = pfun.get_plon_plat(lon,lat)
# # make a version of z with nans where masked
# zm = z.copy()
# zm[np.transpose(mask_rho) == 0] = np.nan
# zm[np.transpose(mask_rho) != 0] = -1

# Point Sources Loading -------------------------------------------------------------

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

# PLOTTING ----------------------------------------------------

# define marker sizes (minimum size is 10 so dots don't get too small)
sizes_wwtps = [max(0.1*load,10) for load in avgload_wwtps['avg-daily-load(kg/d)']]
print(avgload_wwtps)
print(og_loc_df)
print(new_loc_df)

# Plot the two different depths
fs = 14
pfun.start_plot(fs=fs, figsize=(6,6))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.yscale('log')  
plt.xscale('log')  
x_text = [0.9,0.6,1.0,1.0,0.7,0.2]
y_text = [0.3,0.25,0.4,0.85,0.9,0.5]
ind = 0
ax.plot([0,300],[0,300],':k')
# color based on whether the WWTP moved to a different regime (euphotic zone, outflow layer, inflow layer)
for i,size in enumerate(sizes_wwtps):
    og_depth = og_loc_df['depth'][i]
    new_depth = new_loc_df['depth'][i]
    if og_depth < 60 and og_depth > 25 and new_depth < 25:
        ax.scatter(og_depth,new_depth,color='orange',alpha=0.5,zorder=3,s=size)
        ax.annotate(og_loc_df['rname'][i], xy=(og_loc_df['depth'][i], new_loc_df['depth'][i]), xycoords='data',
            xytext=(x_text[ind],y_text[ind]), textcoords='axes fraction', va='center', ha='center', size=10,
            arrowprops=dict(facecolor='black', width=0.1, headwidth=5, headlength=5, edgecolor='black'))
        ind += 1
    elif og_depth > 60 and new_depth < 60:
        ax.scatter(og_depth,new_depth,color='olivedrab',alpha=0.5,zorder=3,s=size)
        ax.annotate(og_loc_df['rname'][i], xy=(og_loc_df['depth'][i], new_loc_df['depth'][i]), xycoords='data',
            xytext=(x_text[ind],y_text[ind]), textcoords='axes fraction', va='center', ha='center', size=10,
            arrowprops=dict(facecolor='black', width=0.1, headwidth=5, headlength=5, edgecolor='black'))
        ind += 1
    else:
        ax.scatter(og_depth,new_depth,color='purple',alpha=0.5,zorder=3,s=size)
# plot key depths
# ax.plot([0,300],[25,25],color='gray')
# ax.plot([0,300],[60,60],color='gray')
# ax.plot([25,25],[0,300],color='gray')
# ax.plot([60,60],[0,300],color='gray')
# format figure
ax.set_ylabel('New Depth [m]')
ax.set_xlim(0,300)
ax.set_ylim(0,300)
ax.set_xlabel('Original Depth [m]')
ax.set_facecolor('#EEEEEE')
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)
ax.grid(True,which='both', color='w')
plt.title('Change in depth of WWTPs\nimplemented as tiny rivers')
# color legend
purple = mpatches.Patch(color='purple', label='No change to depth regime', alpha=0.5)
yellow = mpatches.Patch(color='orange', label=r'Below photic zone $\rightarrow$ photic zone', alpha=0.5)
green = mpatches.Patch(color='olivedrab', label=r'Deep inflow $\rightarrow$ surface outflow', alpha=0.5)
second_legend = plt.legend(handles=[purple,yellow,green],loc='lower right',fontsize=10)
plt.gca().add_artist(second_legend)
# add wwtp discharge legend
leg_szs = [100, 1000, 10000]
szs = [0.1*(leg_sz) for leg_sz in leg_szs]
l0 = plt.scatter([],[], s=szs[0], color='purple', alpha=0.5)
l1 = plt.scatter([],[], s=szs[1], color='purple', alpha=0.5)
l2 = plt.scatter([],[], s=szs[2], color='purple', alpha=0.5)
labels = ['< 100', '1,000', '10,000']
legend1 = ax.legend([l0, l1, l2], labels, fontsize = 10, markerfirst=False,
    title='WWTP Discharge \n'+r' (kg N d$^{-1}$)',loc='upper left', labelspacing=1.6, borderpad=1.5)
plt.setp(legend1.get_title(),fontsize=12)

plt.show()
