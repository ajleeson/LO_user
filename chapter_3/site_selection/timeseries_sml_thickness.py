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
import matplotlib.dates as mdates
import csv
import matplotlib.pylab as plt
from dask.diagnostics import ProgressBar
from pathlib import Path
from datetime import datetime
import gsw
import PyCO2SYS as pyco2
import cmcrameri.cm as cmc
from lo_tools import Lfun, zfun
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

years = ['2015','2016','2017','2018','2019','2020','2021']#['2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025']

# which  model run to look at?
gtagex = 'cas7_t1_x11ab'

out_dir = Ldir['LOo'] / 'chapter_3' / 'data'


# where to put output files
Lfun.make_dir(out_dir)

##############################################################
##                    INITIALIZE STUFF                      ##
##############################################################

# get all of the filepaths
file_paths = [Path('/dat1') / 'kmhewett' / 'LO_output' /'extract'/ 'cas7_t1_x11ab' / 'sml_plus' / 'threshold_p125' / ('LO_domain_'+year+'.01.01_'+year+'.12.31') / ('LO_domain_sml_plus_'+year+'.01.01_'+year+'.12.31.nc') 
    for year in years]


# Open all of the data at once 
ds_all = xr.open_mfdataset(file_paths, combine='by_coords')

# get list of vns
vns = list(ds_all.data_vars)

# y and x indices of sites
natural_ys = [464,875,716,1047,1199,1075,898]
natural_xs = [360,578,504,186, 233, 468, 513]
locations=['Columbia River Plume',
            'Saratoga Passage',
            'Hood Canal',
            'Tofino',
            'Quadra Island',
            'S. SoG (Fraser plume)',
            'S. of San Juan Islands']

##############################################################
##                   CREATE CLIMATOLOGIES                   ##
##############################################################

# Put sites in data arrays AND assign the locations list as the coordinate
y_idx = xr.DataArray(natural_ys, dims="site", coords={"site": locations})
x_idx = xr.DataArray(natural_xs, dims="site", coords={"site": locations})

# Extract data at the seven sites (only variables of interest)
ds_sites = ds_all.isel(eta_rho=y_idx, xi_rho=x_idx)

out_file = out_dir / 'site_daily_climatologies_sml.nc'

if not out_file.exists():
    print("Calculating and saving climatology data (this will take a few minutes)...")
    
    # 1. Create a 'MM-DD' string coordinate for perfect calendar alignment
    month_day = ds_sites['ocean_time'].dt.strftime('%m-%d')
    ds_sites = ds_sites.assign_coords(month_day=month_day)
    # Only keep one variable, e.g., 'temp'
    ds_sites = ds_sites[['SML_thickness']]  # Replace 'temp' with your variable name
    
    # 2. Group by 'month_day'. 
    # Non-leap years simply won't contribute to '02-29' (equivalent to being NaN)
    clim_mean = ds_sites.groupby('month_day').mean(dim='ocean_time')
    clim_std = ds_sites.groupby('month_day').std(dim='ocean_time')
    
    # 3. Replace the 'MM-DD' strings with integers 1 to 366
    # There are exactly 366 unique MM-DD combinations.
    clim_mean['month_day'] = np.arange(1, 367)
    clim_std['month_day'] = np.arange(1, 367)
    
    # Rename the dimension to 'dayofyear' so your plotting code works perfectly
    clim_mean = clim_mean.rename({'month_day': 'dayofyear', **{v: f"{v}_mean" for v in clim_mean.data_vars}})
    clim_std = clim_std.rename({'month_day': 'dayofyear', **{v: f"{v}_std" for v in clim_std.data_vars}})
    
    # Merge and save
    clim_combined = xr.merge([clim_mean, clim_std])
    
    with ProgressBar():
        clim_combined.to_netcdf(out_file)
    print("Saved!")
else:
    print(f"File {out_file} already exists. Ready to plot!")

# #############################################################
# #                        PLOT DATA                         ##
# #############################################################

# # Load the saved data (this is instant!)
# clim_ds = xr.open_dataset(out_dir / 'site_daily_climatologies.nc')

# vns = ['CO2_flux_actual', 'Fs_30', 'CO2_capacity_ideal', 'wind2']

# colors = ['#4477AA','#66CCEE','#228833','#CCBB44','#EE6677','#AA3377','#BBBBBB']

    
# # loop through variables
# for vn in vns:

#     # initialize figure
#     fig, ax = plt.subplots(1, 1, figsize=(11, 5.5))

#     if vn == 'Fs_30':
#         vmin = 0
#         vmax = 5
#         title = r'Freshwater content (s$_0$ = 30) [m]'
#     elif vn == 'CO2_flux_actual':
#         vmin = -40
#         vmax = 20
#         title = r'CO$_2$ flux with wind [mmol m$^{-2}$ day$^{-1}$]'
#         ax.axhline(0, color='black', linewidth=1, linestyle=":") # add horizontal line at y=0 for reference
#     elif vn == 'CO2_capacity_ideal':
#         vmin = -50
#         vmax = 10
#         title = r'Ideal CO$_2$ uptake capacity [mmol m$^{-3}$]'
#         ax.axhline(0, color='black', linewidth=1, linestyle=":") # add horizontal line at y=0 for reference
#     elif vn == 'wind2':
#         vmin = 0
#         vmax = 10
#         title = r'Wind speed [m s$^{-1}$]'

#     # Loop through each site coordinate
#     for i,site_name in enumerate(clim_ds.site.values):      
#         # Get mean and standard deviation
#         mean_data = clim_ds[f"{vn}_mean"].sel(site=site_name)
#         std_data = clim_ds[f"{vn}_std"].sel(site=site_name)

#         if vn == 'wind2':
#                 # take square root
#                 mean_data = np.sqrt(mean_data)
#                 std_data = np.sqrt(std_data)

#         # apply lowpass filter
#         nwin = 14
#         mean_data_filtered = zfun.lowpass(mean_data.values, n=nwin)
#         # get time vector
#         date = pd.to_datetime(clim_ds.dayofyear.values, format='%j')
#         # Plot the mean line
#         ax.plot(clim_ds.dayofyear, mean_data_filtered, label=site_name, color=colors[i],
#                 linewidth=3)
#         # # Add shaded standard deviation around the line
#         # ax.fill_between(clim_ds.dayofyear, 
#         #                 mean_data - std_data, 
#         #                 mean_data + std_data, 
#         #                 color=colors[i], alpha=0.3)

#     # format figure
#     ax.set_title(f'Daily Climatology; {nwin}-day Hanning Window',fontsize=14,fontweight='bold')
#     ax.set_ylabel(title, fontsize=14)
#     ax.set_ylim([vmin,vmax])
#     ax.set_xlim([0, 365])
#     # format grid
#     ax.tick_params(axis='x', labelrotation=30)
#     loc = mdates.MonthLocator(interval=1)
#     ax.xaxis.set_major_locator(loc)
#     ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
#     # ax.grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#     ax.tick_params(axis='both', labelsize=12)
#     ax.legend(loc='best', fontsize=12,ncol=2)

#     # Generate plot
#     plt.tight_layout() # added () to properly execute tight_layout
#     plt.show()