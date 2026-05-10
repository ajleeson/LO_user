# """
# Violin plots of bottom DO concentrations in Puget Sound basins
# Using all daily data from 2015-2020 (no climatological averaging)
# """

# # import things
# import numpy as np
# import xarray as xr
# import pandas as pd
# import matplotlib.pylab as plt
# from statsmodels.graphics.boxplots import violinplot

# from lo_tools import Lfun
# from lo_tools import plotting_functions as pfun

# import sys
# from pathlib import Path
# pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
# if str(pth) not in sys.path:
#     sys.path.append(str(pth))
# import gfun

# Gr = gfun.gstart()
# Ldir = Lfun.Lstart()

# plt.close('all')

# ##############################################################
# ##                       USER INPUTS                        ##
# ##############################################################

# years = ['2015','2016','2017','2018','2019','2020']

# # which model run to look at?
# gtagexes = ['cas7_t1noDIN_x11ab','cas7_t1_x11ab']

# ##############################################################
# ##                      PROCESS DATA                        ##
# ##############################################################

# # open datasets
# ds_loading_dict = {}
# ds_noloading_dict = {}
# for year in years:
#     # open datasets
#     ds_loading = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / ('cas7_t1_x11ab_pugetsoundDO_' + year + '_DO_info.nc'))
#     ds_noloading = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / ('cas7_t1noDIN_x11ab_pugetsoundDO_' + year + '_DO_info.nc'))
    
#     # drop leap day if 2016 or 2020
#     if year in ['2016','2020']:
#         ds_loading = ds_loading.sel(ocean_time=~((ds_loading.ocean_time.dt.month == 2) & (ds_loading.ocean_time.dt.day == 29)))
#         ds_noloading = ds_noloading.sel(ocean_time=~((ds_noloading.ocean_time.dt.month == 2) & (ds_noloading.ocean_time.dt.day == 29)))
    
#     # add ds to dictionary
#     ds_loading_dict[year] = ds_loading
#     ds_noloading_dict[year] = ds_noloading

# # get bottom DO arrays for all years
# DO_loading_arrays = [ds['DO_bot'].data for ds in ds_loading_dict.values()]  # [365, eta, xi]
# DO_noloading_arrays = [ds['DO_bot'].data for ds in ds_noloading_dict.values()]  # [365, eta, xi]

# # stack the arrays (Keep all years! Shape: [6 years, 365 days, eta, xi])
# DO_loading_stacked = np.stack(DO_loading_arrays, axis=0)  
# DO_noloading_stacked = np.stack(DO_noloading_arrays, axis=0) 

# # Load basin masks
# basin_mask_ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / 'basin_masks_from_pugetsoundDObox.nc')

# # define basins to loop through
# basins = {
#     'Hood Canal': basin_mask_ds['mask_hoodcanal'].data == 1,
#     'South Sound': basin_mask_ds['mask_southsound'].data == 1,
#     'Main Basin': basin_mask_ds['mask_mainbasin'].data == 1,
#     'Whidbey Basin': basin_mask_ds['mask_whidbeybasin'].data == 1
# }

# # crop by season indices
# jan = 0
# apr = 90
# jul = 181
# oct = 273

# ##############################################################
# ##                        PLOTTING                          ##
# ##############################################################

# # Initialize figure (2x2 layout for 4 basins)
# fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharey=True, sharex=True)
# ax = axes.ravel()

# # Loop through each basin and plot
# for i, (basin_name, mask) in enumerate(basins.items()):
    
#     # Extract seasonal data for this specific basin using the mask
#     # The ":" grabs all 6 years, then we slice by day, then apply the spatial mask
    
#     # Winter (Jan-Apr)
#     winter_loading   = DO_loading_stacked[:, jan:apr, mask].flatten()
#     winter_noloading = DO_noloading_stacked[:, jan:apr, mask].flatten()
    
#     # Spring (Apr-Jul)
#     spring_loading   = DO_loading_stacked[:, apr:jul, mask].flatten()
#     spring_noloading = DO_noloading_stacked[:, apr:jul, mask].flatten()
    
#     # Summer (Jul-Oct)
#     summer_loading   = DO_loading_stacked[:, jul:oct, mask].flatten()
#     summer_noloading = DO_noloading_stacked[:, jul:oct, mask].flatten()
    
#     # Fall (Oct-Dec)
#     fall_loading   = DO_loading_stacked[:, oct::, mask].flatten()
#     fall_noloading = DO_noloading_stacked[:, oct::, mask].flatten()

#     # Drop NaNs before plotting
#     winter_loading = winter_loading[~np.isnan(winter_loading)]
#     winter_noloading = winter_noloading[~np.isnan(winter_noloading)]
    
#     spring_loading = spring_loading[~np.isnan(spring_loading)]
#     spring_noloading = spring_noloading[~np.isnan(spring_noloading)]
    
#     summer_loading = summer_loading[~np.isnan(summer_loading)]
#     summer_noloading = summer_noloading[~np.isnan(summer_noloading)]
    
#     fall_loading = fall_loading[~np.isnan(fall_loading)]
#     fall_noloading = fall_noloading[~np.isnan(fall_noloading)]

#     # Plot No-Loading (Left Side - Black)
#     violinplot([winter_noloading, spring_noloading, summer_noloading, fall_noloading],
#                 positions=[1,2,3,4], show_boxplot=False,
#                 side='left', ax=ax[i], plot_opts={'violin_fc':'black'})
    
#     # Plot Loading (Right Side - Blue)
#     violinplot([winter_loading, spring_loading, summer_loading, fall_loading],
#                 positions=[1,2,3,4], show_boxplot=False,
#                 side='right', ax=ax[i], plot_opts={'violin_fc':'cornflowerblue'})

#     # format subplot
#     ax[i].grid(visible=True, axis='both', color='silver', linestyle='--')
#     ax[i].tick_params(axis='both', labelsize=12)
#     xticks = np.arange(1, 5, 1)
#     ax[i].set_xticks(xticks, labels=['Winter','Spring','Summer','Fall'])
    
#     # Basin Title
#     ax[i].set_title(f'{basin_name}', fontweight='bold', fontsize=12)
    
#     # Add legends/text only to the first subplot to avoid clutter
#     if i == 0:
#         ax[i].text(0.6, ax[i].get_ylim()[1]*0.9, 'No-Loading', fontweight='bold', color='black', alpha=0.6, fontsize=11)
#         ax[i].text(0.6, ax[i].get_ylim()[1]*0.8, 'Loading', fontweight='bold', color='cornflowerblue', alpha=0.6, fontsize=11)

# # Set Y labels for left column
# ax[0].set_ylabel('DO [mg/L]', fontsize=12)
# ax[2].set_ylabel('DO [mg/L]', fontsize=12)

# plt.tight_layout()
# plt.show()



"""
Violin plots of bottom DO differences (Loading - No-Loading) in Puget Sound basins
Using all daily data from 2015-2020
"""

# import things
import numpy as np
import xarray as xr
import matplotlib.pylab as plt

from lo_tools import Lfun
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

years = ['2015','2016','2017','2018','2019','2020']

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# open datasets
ds_loading_dict = {}
ds_noloading_dict = {}
for year in years:
    ds_loading = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / ('cas7_t1_x11ab_pugetsoundDO_' + year + '_DO_info.nc'))
    ds_noloading = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / ('cas7_t1noDIN_x11ab_pugetsoundDO_' + year + '_DO_info.nc'))
    
    # drop leap day if 2016 or 2020
    if year in ['2016','2020']:
        ds_loading = ds_loading.sel(ocean_time=~((ds_loading.ocean_time.dt.month == 2) & (ds_loading.ocean_time.dt.day == 29)))
        ds_noloading = ds_noloading.sel(ocean_time=~((ds_noloading.ocean_time.dt.month == 2) & (ds_noloading.ocean_time.dt.day == 29)))
    
    ds_loading_dict[year] = ds_loading
    ds_noloading_dict[year] = ds_noloading

# get bottom DO arrays for all years
DO_loading_arrays = [ds['DO_bot'].data for ds in ds_loading_dict.values()]
DO_noloading_arrays = [ds['DO_bot'].data for ds in ds_noloading_dict.values()]

# stack the arrays [6 years, 365 days, eta, xi]
DO_loading_stacked = np.stack(DO_loading_arrays, axis=0)  
DO_noloading_stacked = np.stack(DO_noloading_arrays, axis=0) 

# Calculate the difference directly! (Loading - No-Loading)
# A negative value means the anthropogenic loading decreased the DO.
Delta_DO = DO_loading_stacked - DO_noloading_stacked

# Load basin masks
basin_mask_ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / 'basin_masks_from_pugetsoundDObox.nc')

# define basins to loop through
basins = {
    'Hood Canal': basin_mask_ds['mask_hoodcanal'].data == 1,
    'South Sound': basin_mask_ds['mask_southsound'].data == 1,
    'Main Basin': basin_mask_ds['mask_mainbasin'].data == 1,
    'Whidbey Basin': basin_mask_ds['mask_whidbeybasin'].data == 1
}

# crop by season indices
jan = 0
apr = 90
jul = 181
oct = 273

##############################################################
##                        PLOTTING                          ##
##############################################################

# Initialize figure (2x2 layout for 4 basins)
fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharey=True, sharex=True)
ax = axes.ravel()

for i, (basin_name, mask) in enumerate(basins.items()):
    
    # Extract seasonal delta data for this specific basin using the mask
    winter_delta = Delta_DO[:, jan:apr, mask].flatten()
    spring_delta = Delta_DO[:, apr:jul, mask].flatten()
    summer_delta = Delta_DO[:, jul:oct, mask].flatten()
    fall_delta   = Delta_DO[:, oct::, mask].flatten()

    # Drop NaNs before plotting
    winter_delta = winter_delta[~np.isnan(winter_delta)]
    spring_delta = spring_delta[~np.isnan(spring_delta)]
    summer_delta = summer_delta[~np.isnan(summer_delta)]
    fall_delta   = fall_delta[~np.isnan(fall_delta)]

    data_to_plot = [winter_delta, spring_delta, summer_delta, fall_delta]

    # Plot using standard matplotlib violins (no split needed for 1 variable)
    parts = ax[i].violinplot(data_to_plot, positions=[1, 2, 3, 4], 
                             showmeans=True, showextrema=True)
    
    # Make the violins look nice (color them based on the fact it's a difference)
    for pc in parts['bodies']:
        pc.set_facecolor('crimson')
        pc.set_edgecolor('black')
        pc.set_alpha(0.6)
    
    parts['cmeans'].set_color('black')
    parts['cmaxes'].set_color('black')
    parts['cmins'].set_color('black')
    parts['cbars'].set_color('black')

    # Add a horizontal line at exactly 0.0 to act as a reference
    ax[i].axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.8)

    # format subplot
    ax[i].grid(visible=True, axis='both', color='silver', linestyle=':', alpha=0.7)
    ax[i].tick_params(axis='both', labelsize=12)
    ax[i].set_xticks([1, 2, 3, 4])
    ax[i].set_xticklabels(['Winter', 'Spring', 'Summer', 'Fall'])
    
    # Basin Title
    ax[i].set_title(f'{basin_name}', fontweight='bold', fontsize=12)

# Set Y labels for left column
ax[0].set_ylabel(r'$\Delta$ DO (Load - NoLoad) [mg/L]', fontsize=12)
ax[2].set_ylabel(r'$\Delta$ DO (Load - NoLoad) [mg/L]', fontsize=12)

plt.tight_layout()
plt.show()