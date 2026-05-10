"""
Get CO2 uptake capacity from model history file
and freshwater content
and other useful values for understanding the drivers of CO2 uptake capacity

.nc files are saved in LO_output/chapter_3/data
"""

# import things
import numpy as np
import xarray as xr
import pandas as pd
import csv
import matplotlib.pylab as plt
from pathlib import Path
from datetime import datetime
import gsw
import PyCO2SYS as pyco2
import cmcrameri.cm as cmc
from lo_tools import Lfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

plt.close('all')

##############################################################
##                       USER INPUTS                        ##
##############################################################

years = ['2015','2016','2017','2018','2019','2020','2021','2022','2023','2024']

# which  model run to look at?
gtagex = 'cas7_t1_x11ab'

out_dir = Ldir['LOo'] / 'chapter_3' / 'data'




# # uncomment below to save data ---------------------------------------------------------------------------------------------------------

# # where to put output files
# Lfun.make_dir(out_dir)

# ##############################################################
# ##                GET MONTHLY MIN/MAX/MEAN                  ##
# ##############################################################

# # get all of the filepaths
# file_paths = [Ldir['LOo'] / 'chapter_3' / 'data' / f"{gtagex}_{year}_freshwatercontent_CO2uptake.nc" 
#     for year in years]

# # Open all of the data at once 
# ds_all = xr.open_mfdataset(file_paths, combine='by_coords')

# # get list of vns
# vns = list(ds_all.data_vars)

# # 3. Create a dictionary to map integer months (1-12) to strings ('Jan'-'Dec')
# month_map = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
#              7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}

# # Loop through months 1 to 12 explicitly
# for month_num in range(1, 13):
#     month_str = month_map[month_num]
#     print(f'--- Processing {month_str} ---')
    
#     month_list = []
    
#     # Loop through years and get one month at a time
#     for year in years:
#         filepath = Ldir['LOo'] / 'chapter_3' / 'data' / f"{gtagex}_{year}_freshwatercontent_CO2uptake.nc"
        
#         # Open one year of data
#         ds_year = xr.open_dataset(filepath)
        
#         # Get data for one month
#         ds_month_year = ds_year.sel(ocean_time=ds_year['ocean_time'].dt.month == month_num)
        
#         # Add to the list
#         month_list.append(ds_month_year.load())
#         ds_year.close()
        
#     # Concatenate the 10 years of this specific month
#     ds_month_all_years = xr.concat(month_list, dim='ocean_time')
    
#     print(f'Calculating stats for {month_str}...')
#     # Calculate min, max, mean (this will be instantly fast because it's in RAM)
#     ds_mean = ds_month_all_years.mean(dim='ocean_time', keep_attrs=True)
#     ds_min  = ds_month_all_years.min(dim='ocean_time', keep_attrs=True)
#     ds_max  = ds_month_all_years.max(dim='ocean_time', keep_attrs=True)
    
#     # Rename variables
#     ds_mean = ds_mean.rename({v: f"{v}_mean" for v in ds_mean.data_vars})
#     ds_min  = ds_min.rename({v: f"{v}_min" for v in ds_min.data_vars})
#     ds_max  = ds_max.rename({v: f"{v}_max" for v in ds_max.data_vars})

#     # Combine into a single dataset
#     ds_merged = xr.merge([ds_mean, ds_min, ds_max])
    
#     # Save to disk
#     output_filename = out_dir / f'monthly_climatology_freshwaterCO2_{month_str}.nc'
#     print(f'Saving {month_str} to disk...')
#     ds_merged.to_netcdf(output_filename)
    
#     print(f'Finished {month_str}!\n')
# # uncomment above to save new data ---------------------------------------------------------------------------------------------------------



#############################################################
#                        PLOT DATA                         ##
#############################################################

# vns = ['CO2_capacity_ideal', 'CO2_flux_actual', 'surf_temp', 'surf_salt',
# 'surf_alk', 'surf_TIC', 'wind2', 'delta_pCO2', 'Fs_31', 'Fs_30', 'Fs_29']
vns = ['CO2_flux_actual']#['CO2_capacity_ideal', 'CO2_flux_actual', 'Fs_30']

stats = ['mean', 'min', 'max']

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lons = lon[0,:]
lats = lat[:,0]
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values
px, py = pfun.get_plon_plat(lon,lat)

# lon/lat limits (Study Domain)
xmin = -126
xmax = -122
ymin = 45.5
ymax = 50.5

# loop through months
for month in ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']:

    # open dataset for this month
    ds_month = xr.open_dataset(out_dir / f'monthly_climatology_freshwaterCO2_{month}.nc')

    # loop through variables
    for vn in vns:

        # initialize figure
        fig,ax = plt.subplots(1,3, figsize=(14,8))
        ax = ax.ravel()

        # get colorbar limits and cmap and units
        if vn == 'CO2_capacity_ideal':
            vmin_mean = -20
            vmax_mean =  20
            vmin_lims = -40
            vmax_lims =  40
            cmap = cmc.vik
            title = r'Ideal CO$_2$ uptake capacity [mmol m$^{-3}$]'
        elif vn == 'CO2_flux_actual':
            vmin_mean = -15
            vmax_mean =  15
            vmin_lims = -60
            vmax_lims =  60
            cmap = cmc.vik
            title = r'CO$_2$ flux with wind [mmol m$^{-2}$ day$^{-1}$]'
        elif vn == 'Fs_30':
            vmin_mean = 0
            vmax_mean = 3
            vmin_lims = 0
            vmax_lims = 10
            cmap = cmc.batlowW_r
            title = r'Freshwater content (s$_0$ = 30) [m]'
        elif vn == 'wind2':
            vmin_mean = 0
            vmax_mean = 10
            vmin_lims = 0
            vmax_lims = 20
            cmap = cmc.tokyo
            title = r'Wind speed [m s$^{-1}$]'

        # loop through stats:
        for i,stat in enumerate(stats):

            if stat == 'mean':
                vmin = vmin_mean
                vmax = vmax_mean
            else:
                vmin = vmin_lims
                vmax = vmax_lims

            # get values to plot
            var = ds_month[vn + '_' + stat].values

            if vn == 'wind2':
                # take square root
                var = np.sqrt(var)

            # plot values
            cs = ax[i].pcolormesh(px,py,var,vmin=vmin,vmax=vmax,cmap=cmap)

            # Add Puget Sound Inset
            # [x0, y0, width, height]
            axins = ax[i].inset_axes([0.71, 0.0, 0.45, 0.6])
            # plot values in inset
            axins.pcolormesh(px, py, var,vmin=vmin, vmax=vmax,cmap=cmap)
            # Puget Sound limits
            axins.set_xlim(-123.2, -122.1)
            axins.set_ylim(46.95, 48.4)
            axins.tick_params(left=False, bottom=False)
            # format
            axins.set_xticklabels([])
            axins.set_yticklabels([])
            for spine in axins.spines.values():
                spine.set_edgecolor('grey')
                spine.set_linewidth(2)

            # create colorbarlegend
            if i == 0:
                cbar_ax = fig.add_axes([0.12, 0.08, 0.23, 0.03])
                cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal')
                cbar.ax.tick_params(labelsize=14, rotation=30)
                cbar.outline.set_visible(False)
            elif i == 1:
                cbar_ax = fig.add_axes([0.4, 0.08, 0.5, 0.03])
                cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal')
                cbar.ax.tick_params(labelsize=14, rotation=30)
                cbar.outline.set_visible(False)

            # # add testbed sites
            # testbed_color = 'crimson'
            # testbed_outline = 'pink'
            # testbed_lats = [48.129378, 48.078611, 48.5453, 47.584536, 47.60079806]
            # testbed_lons = [-123.457488, -123.045000, -123.0121, -122.342908, -122.43131423]
            # for lat, lon in zip(testbed_lats, testbed_lons):
            #     ax[i].scatter(   lon, lat, color='None', s=80, edgecolor=testbed_outline, linewidth=4, zorder=5, marker='D')
            #     ax[i].scatter(   lon, lat, color='None', s=80, edgecolor=testbed_color,   linewidth=2, zorder=5, marker='D')
            #     axins.scatter(lon, lat, color='None', s=80, edgecolor=testbed_outline, linewidth=4, zorder=5, marker='D')
            #     axins.scatter(lon, lat, color='None', s=80, edgecolor=testbed_color,   linewidth=2, zorder=5, marker='D')

            # add selection of natural sites
            natural_color = 'black'
            natural_outline = 'white'
            locations=['Columbia River Plume',
                    'Saratoga Passage',
                    'Hood Canal',
                    'Van Island Coast',
                    'Quadra Island',
                    'S. SoG (Fraser plume)',
                    'S. of San Juan Islands']
            natural_ys = [464,875,716,1047,1199,1075,898]
            natural_xs = [360,578,504,186, 233, 468, 513]
            for x, y, loc in zip(natural_xs, natural_ys, locations):
                if x in [578,504]: # only plot Whidbey and Hood Canal in inset
                    axins.scatter(lons[x], lats[y], color='None', s=180, edgecolor=natural_outline, linewidth=4, zorder=5)
                    axins.scatter(lons[x], lats[y], color='None', s=180, edgecolor=natural_color,   linewidth=2, zorder=5)
                else:
                    ax[i].scatter(   lons[x], lats[y], color='None', s=180, edgecolor=natural_outline, linewidth=4, zorder=5)
                    ax[i].scatter(   lons[x], lats[y], color='None', s=180, edgecolor=natural_color,   linewidth=2, zorder=5)

            # format figure
            ax[i].set_xlim([xmin,xmax])
            ax[i].set_ylim([ymin,ymax])
            ax[i].set_xticklabels([])
            ax[i].set_yticklabels([])
            ax[i].set_title(stat,fontsize=12, fontweight='bold',loc='left')
            ax[i].tick_params(left=False, bottom=False)
            pfun.add_coast(ax[i], color='silver')
            pfun.add_coast(axins, color='silver')
            pfun.dar(ax[i])
            pfun.dar(axins)
            plt.suptitle('(2015 - 2024) '+ month + ' - ' + title,
                         fontsize=14, fontweight='bold')

        # Generate plot
        plt.tight_layout
        plt.subplots_adjust(bottom=0.12, top=0.9)
        plt.show()