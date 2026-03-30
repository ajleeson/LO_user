"""
This focuses on property-property plots and obs-mod plots.

It is focused on a single run and year.

In the "SET CHOICES" we pull in optional command line arguments that
allow more specific plotting choices, such as using a single obs source,
or only plotting Salish Sea stations. These are best used with
-test True when plotting by hand, as opposed to using the one_step_bottle_val_plot.py
driver.

For -test True it shows the plot on the screen instead of saving to a file.

Testing on mac:
run plot_val -gtx cas7_t1_x11ab -test True

"""
import sys
import pandas as pd
import numpy as np
import pickle
from lo_tools import plotting_functions as pfun
import cmocean
from lo_tools import Lfun, zfun, zrfun
import obsmod_functions as omfun
import xarray as xr

# command line arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-gtx', '--gtagex', type=str)   # e.g. cas7_t1_x11ab
parser.add_argument('-otype', type=str, default='bottle') # observation type, e.g. ctd, bottle, etc.
parser.add_argument('-year', type=int) # e.g. 2019
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
# more optional arguments
parser.add_argument('-single_source', type=str, default='all') # 'all' or a souce name like 'kc'
parser.add_argument('-dividing_depth', type=int, default=10) # depth [m] to delineate shallow
parser.add_argument('-do_arag', default=True, type=Lfun.boolean_string) # plot arag and pCO2
parser.add_argument('-coast', default=False, type=Lfun.boolean_string) # only plot coastal stations
parser.add_argument('-salish', default=False, type=Lfun.boolean_string) # only plot Salish stations
args = parser.parse_args()

Ldir = Lfun.Lstart()

if '_mac' in Ldir['lo_env']: # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt

in_dir = Ldir['LOo'] / 'obsmod'

# run info
# year = str(args.year)
gtx = args.gtagex
otype = args.otype

years = ['2015','2016','2017','2018','2019','2020']

# specify input
df0_dict = {}
for y,year in enumerate(years):
    in_fn = in_dir / ('combined_' + otype + '_' + year + '_' + gtx + '.p')
    df0_dict_temp = pickle.load(open(in_fn, 'rb'))
    if y == 0:
        df0_dict['obs'] = df0_dict_temp['obs']
        df0_dict[gtx] = df0_dict_temp[gtx]
    else:
        df0_dict['obs'] = pd.concat([df0_dict['obs'], df0_dict_temp['obs']])
        df0_dict[gtx] = pd.concat([df0_dict[gtx], df0_dict_temp[gtx]])


# add DIN field
for gtxo in df0_dict.keys():
    df0_dict[gtxo]['DIN (uM)'] = df0_dict[gtxo]['NO3 (uM)'] + df0_dict[gtxo]['NH4 (uM)']

# add DO in mg/L
for og in ['obs',gtx]:
    df0_dict[og]['DO (mg L-1)'] = (32 / 1000) * df0_dict[og]['DO (uM)']

###########################################################################
# PLOTTING - STRAIT OF GEORGIA

plt.close('all')
      
fig = plt.figure()
df_dict = df0_dict.copy()
        
# station map
# initialize figure
fig = plt.figure(figsize=(9.5, 6))
gs = fig.add_gridspec(nrows=4, ncols=4)
ax_map = fig.add_subplot(gs[:, 0])
ax_SA = fig.add_subplot(gs[0, 1:4])
ax_CT = fig.add_subplot(gs[1, 1:4])#,sharex=ax_SA)
ax_DO = fig.add_subplot(gs[2, 1:4])#,sharex=ax_SA)
ax_NO3 = fig.add_subplot(gs[3, 1:4])#,sharex=ax_SA)

# get just one region from dfo observations
obs_df = df_dict['obs']
lat_min, lat_max = 49.155, 49.17
lon_min, lon_max = -124, -123 #-123.5055, -125.5045
mask = ((obs_df['source'] == 'dfo1') &
        obs_df['lat'].between(lat_min, lat_max) &
        obs_df['lon'].between(lon_min, lon_max))
station_df = obs_df.loc[mask]
print(station_df)
ax_map.scatter(station_df['lon'], station_df['lat'], s=10, color='black', zorder=5)

# ax_map.scatter(-123.55, 49.1635, s=10, color='deeppink', zorder=5)


# get the grid data
grid_fn = Ldir['data'] / 'grids' / 'cas7' / 'grid.nc'
grid_ds = xr.open_dataset(grid_fn)
z = -grid_ds.h.values
mask_rho = np.transpose(grid_ds.mask_rho.values)
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
plon, plat = pfun.get_plon_plat(lon,lat)
# make a version of z with nans where masked
# this gives us a binary map of land and water cells
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
zm[np.transpose(mask_rho) != 0] = -1
ax_map.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-1.5, vmax=0, cmap='Greys')
# pfun.add_coast(ax)
ax_map.axis([-125,-122,45,52])
pfun.dar(ax_map)
ax_map.set_xlabel('')
ax_map.set_ylabel('')

############################
# plot timeseries

print(list(station_df.keys()))

# ax_time.scatter(station_df['time'],station_df['z'])

# separate into depth bins
deep_min, deep_max = -255, -245 # 250 m
midd_min, midd_max = -22, -17   # 20 m
shal_min, shal_max = -6, -4     # 5 m
deep_mask = (station_df['z'].between(deep_min, deep_max))
midd_mask = (station_df['z'].between(midd_min, midd_max))
shal_mask = (station_df['z'].between(shal_min, shal_max))
deep_df = station_df.loc[deep_mask]
midd_df = station_df.loc[midd_mask]
shal_df = station_df.loc[shal_mask]

axes = [ax_SA,ax_CT,ax_DO,ax_NO3]
vns = ['SA','CT','DO (mg L-1)','NO3 (uM)']
deep_color = 'navy'
midd_color = 'lightseagreen'
shal_color = 'yellowgreen'

for a,axis in enumerate(axes):
    axis.text(0.02,0.8,vns[a],transform=axis.transAxes,fontweight='bold',fontsize=12)
    axis.scatter(deep_df['time'],deep_df[vns[a]],color=deep_color,marker='o',s=15)
    axis.scatter(midd_df['time'],midd_df[vns[a]],color=midd_color,marker='D',s=15)
    axis.scatter(shal_df['time'],shal_df[vns[a]],color=shal_color,marker='*',s=18)
    if a < 3:
        axis.set_xticks([])
    axis.tick_params(axis='both',labelsize=12)





plt.subplots_adjust(hspace=0.01)
plt.tight_layout()
plt.show()


###########################################################################
# PLOTTING - COAST
      
fig = plt.figure()
df_dict = df0_dict.copy()
        
# station map
# initialize figure
fig = plt.figure(figsize=(9.5, 6))
gs = fig.add_gridspec(nrows=4, ncols=4)
ax_map = fig.add_subplot(gs[:, 0])
ax_SA = fig.add_subplot(gs[0, 1:4])
ax_CT = fig.add_subplot(gs[1, 1:4])#,sharex=ax_SA)
ax_DO = fig.add_subplot(gs[2, 1:4])#,sharex=ax_SA)
ax_NO3 = fig.add_subplot(gs[3, 1:4])#,sharex=ax_SA)

# get just one region from dfo observations
obs_df = df_dict['obs']
lat_min, lat_max = 48.36, 48.37
lon_min, lon_max = -125.60, -125.56
mask = ((obs_df['source'] == 'dfo1') &
        obs_df['lat'].between(lat_min, lat_max) &
        obs_df['lon'].between(lon_min, lon_max))
station_df = obs_df.loc[mask]
print(station_df)
ax_map.scatter(station_df['lon'], station_df['lat'], s=10, color='black', zorder=5)

# ax_map.scatter(-125.58, 48.367, s=10, color='deeppink', zorder=5)


# get the grid data
grid_fn = Ldir['data'] / 'grids' / 'cas7' / 'grid.nc'
grid_ds = xr.open_dataset(grid_fn)
z = -grid_ds.h.values
mask_rho = np.transpose(grid_ds.mask_rho.values)
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
plon, plat = pfun.get_plon_plat(lon,lat)
# make a version of z with nans where masked
# this gives us a binary map of land and water cells
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
zm[np.transpose(mask_rho) != 0] = -1
ax_map.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-1.5, vmax=0, cmap='Greys')
# pfun.add_coast(ax)
ax_map.axis([-126.5,-123.5,45,52])
pfun.dar(ax_map)
ax_map.set_xlabel('')
ax_map.set_ylabel('')

############################
# plot timeseries

print(list(station_df.keys()))

# ax_time.scatter(station_df['time'],station_df['z'])

# separate into depth bins
deep_min, deep_max = -145, -135 # 140 m
midd_min, midd_max = -22, -18   # 20 m
shal_min, shal_max = -6, -4     # 5 m
deep_mask = (station_df['z'].between(deep_min, deep_max))
midd_mask = (station_df['z'].between(midd_min, midd_max))
shal_mask = (station_df['z'].between(shal_min, shal_max))
deep_df = station_df.loc[deep_mask]
midd_df = station_df.loc[midd_mask]
shal_df = station_df.loc[shal_mask]

# ax_SA.scatter(station_df['time'],station_df['z'])

axes = [ax_SA,ax_CT,ax_DO,ax_NO3]
vns = ['SA','CT','DO (mg L-1)','NO3 (uM)']
deep_color = 'navy'
midd_color = 'lightseagreen'
shal_color = 'yellowgreen'

for a,axis in enumerate(axes):
    axis.text(0.02,0.8,vns[a],transform=axis.transAxes,fontweight='bold',fontsize=12)
    axis.scatter(deep_df['time'],deep_df[vns[a]],color=deep_color,marker='o',s=15)
    axis.scatter(midd_df['time'],midd_df[vns[a]],color=midd_color,marker='D',s=15)
    axis.scatter(shal_df['time'],shal_df[vns[a]],color=shal_color,marker='*',s=18)
    if a < 3:
        axis.set_xticks([])
    axis.tick_params(axis='both',labelsize=12)


plt.subplots_adjust(hspace=0.01)
plt.tight_layout()
plt.show()
