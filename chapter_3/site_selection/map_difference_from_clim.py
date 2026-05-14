"""
Generates monthly climatology + year-to-climatology anomaly maps.
Climatology panel is large on left, 10 yearly anomaly panels in a 2x5 grid on the right.
- Climatology colorbar: horizontal below the map.
- Anomaly matrix: one shared horizontal colorbar below matrix.
- Each anomaly panel: just the year in the title, with RMSE printed on its map.
Auto-saves PNG for each month to the figures directory.
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import colors
from pathlib import Path
import matplotlib.patheffects as patheffects

import cmcrameri.cm as cmc

from lo_tools import Lfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

plt.close('all')

##############################################################
##                       USER INPUTS                        ##
##############################################################

vn = 'SML_thickness'   # Variable to plot
gtagex = 'cas7_t1_x11ab'
years = ['2015','2016','2017','2018','2019','2020','2021','2022','2023','2024']
month_dict = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
             7: 'Jul', 8: 'Aug', 9: 'Sep', 10: 'Oct', 11: 'Nov', 12: 'Dec'}

# Set up output directories
Ldir = Lfun.Lstart()
Gr = gfun.gstart()
data_dir = Ldir['LOo'] / 'chapter_3' / 'data'
out_dir = Ldir['LOo'] / 'chapter_3' / 'figures'
out_dir.mkdir(exist_ok=True, parents=True)

# Get grid data
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lons = lon[0, :]
lats = lat[:, 0]
px, py = pfun.get_plon_plat(lon, lat)
xmin, xmax, ymin, ymax = -126, -122, 45.5, 50.5


# Plotting parameters for variable
if vn == 'CO2_capacity_ideal':
    clim_cmap = cmc.vik
    clim_title = r'Ideal CO$_2$ uptake capacity [mmol m$^{-3}$]'
    clim_vmin, clim_vmax = -20, 20
    diff_cmap = cmc.vik
    diff_vmin, diff_vmax = -10, 10
elif vn == 'CO2_flux_actual':
    clim_cmap = cmc.vik
    clim_title = r'CO$_2$ flux with wind [mmol m$^{-2}$ day$^{-1}$]'
    clim_vmin, clim_vmax = -15, 15
    diff_cmap = cmc.vik
    diff_vmin, diff_vmax = -8, 8
elif vn == 'Fs_30':
    clim_cmap = cmc.batlowW_r
    clim_title = r'Freshwater content (s$_0$ = 30) [m]'
    clim_vmin, clim_vmax = 0, 3
    diff_cmap = cmc.vik
    diff_vmin, diff_vmax = -2, 2
elif vn == 'wind2':
    clim_cmap = cmc.tokyo
    clim_title = r'Wind speed [m s$^{-1}$]'
    clim_vmin, clim_vmax = 0, 10
    diff_cmap = cmc.vik
    diff_vmin, diff_vmax = -4, 4
elif vn == 'surf_temp':
    clim_cmap = cmc.roma_r
    clim_title = r'Surface Temp [degC]'
    clim_vmin, clim_vmax = 5, 20
    diff_cmap = cmc.vik
    diff_vmin, diff_vmax = -5, 5
elif vn == 'SML_thickness':
    clim_cmap = cmc.lapaz_r
    clim_title = 'Surface mixed layer thickness [m]'
    clim_vmin, clim_vmax = 0, 35
    diff_cmap = cmc.vik
    diff_vmin, diff_vmax = -10, 10
    
else:
    raise ValueError(f"No plotting settings for variable '{vn}'! Add options above.")

##############################################################
##                PLOT MONTHLY CLIMATOLOGIES                ##
##############################################################

for month_num, month_str in month_dict.items(): # [(7, "Jul")]: 
    print(f'Processing month: {month_str}')
    if vn == 'SML_thickness':
        clim_path = data_dir / f'monthly_climatology_SML_{month_str}.nc'
    else:
        clim_path = data_dir / f'monthly_climatology_freshwaterCO2_{month_str}.nc'
    clim_ds = xr.open_dataset(clim_path)
    clim_mean = clim_ds[f'{vn}_mean'].values

    # --- Set up spacing for subplots
    fig = plt.figure(figsize=(12, 7), constrained_layout=True)
    gs = fig.add_gridspec(3, 6, height_ratios=[14, 14, 2], width_ratios=[2.2, 1, 1, 1, 1, 1],
                         wspace=0.02, hspace=0.12)
    ax_clim = fig.add_subplot(gs[0:2, 0])
    anomaly_axes = []
    for yrow in range(2):
        for ycol in range(5):
            ax = fig.add_subplot(gs[yrow, 1 + ycol])
            anomaly_axes.append(ax)
    # Colorbars
    cax_clim = fig.add_subplot(gs[2, 0])
    cax_anom = fig.add_subplot(gs[2, 1:])

    if vn == 'wind2':
        clim_mean = np.sqrt(clim_mean)

    # --- Plot climatological mean ---
    pc_clim = ax_clim.pcolormesh(px, py, clim_mean, vmin=clim_vmin, vmax=clim_vmax, cmap=clim_cmap)
    ax_clim.set_title('Climatological mean', fontsize=12, fontweight='bold', loc='left')
    pfun.add_coast(ax_clim, color='grey')
    ax_clim.set_xlim([xmin, xmax])
    ax_clim.set_ylim([ymin, ymax])
    ax_clim.set_xticks([]); ax_clim.set_yticks([])

    pfun.dar(ax_clim)
    # Horizontal colorbar for climatology
    cbar_clim = fig.colorbar(pc_clim, cax=cax_clim, orientation='horizontal')
    cbar_clim.set_label(clim_title, fontsize=12)
    cbar_clim.outline.set_visible(False)
    cbar_clim.ax.tick_params(labelsize=12, rotation=30)

    # Calculate difference between each year and climatology
    diffs = []
    for year in years:
        if vn == 'SML_thickness':
            fpath = data_dir / f'LO_domain_sml_plus_{year}.01.01_{year}.12.31.nc'
        else:
            fpath = data_dir / f"{gtagex}_{year}_freshwatercontent_CO2uptake.nc"
        ds_year = xr.open_dataset(fpath)
        ds_month = ds_year.sel(ocean_time=ds_year['ocean_time'].dt.month == month_num)
        var = ds_month[vn].values
        if var.ndim == 3:
            year_mean = np.nanmean(var, axis=0)  # (y, x)
        else:
            year_mean = var
        if vn == 'wind2':
            year_mean = np.sqrt(year_mean)
        diff = year_mean - clim_mean
        diffs.append(diff)


    # Plot difference between year and anomaly
    for i, (year, ax, diff) in enumerate(zip(years, anomaly_axes, diffs)):
        pc = ax.pcolormesh(px, py, diff, vmin=diff_vmin, vmax=diff_vmax, cmap=diff_cmap)
        ax.set_title(year, fontsize=12, loc='left', fontweight='bold')
        pfun.add_coast(ax, color='grey')
        ax.set_xlim([xmin, xmax]); ax.set_ylim([ymin, ymax])
        ax.set_xticks([]); ax.set_yticks([])

        pfun.dar(ax)
        # Only calculate RMSE inside the region shown in the axis
        lon_rho = grid_ds.lon_rho.values
        lat_rho = grid_ds.lat_rho.values
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        # create mask
        mask = ((lon_rho >= x0) & (lon_rho <= x1) &
                (lat_rho >= y0) & (lat_rho <= y1) & ~np.isnan(diff))
        # calculate RMSE and add label to subplot
        if np.any(mask):
            rmse = np.sqrt(np.nanmean(diff[mask]**2))
            ax.text(0.98, 0.96, f'RMSE=\n{rmse:.2f}',
                transform=ax.transAxes, fontsize=10, color='k', va='top', ha='right')

    # Remove any unused axes (should not be any, but for robustness)
    for ax in anomaly_axes[len(years):]:
        ax.remove()

    # --- Single shared anomaly colorbar
    norm = colors.Normalize(vmin=diff_vmin, vmax=diff_vmax)
    sm = plt.cm.ScalarMappable(cmap=diff_cmap, norm=norm)
    sm.set_array([])
    cbar_anom = fig.colorbar(sm, cax=cax_anom, orientation="horizontal")
    cbar_anom.set_label('Anomaly', fontsize=12)
    cbar_anom.ax.tick_params(labelsize=12, rotation=30)
    cbar_anom.outline.set_visible(False)

    # format figure
    fig.suptitle(month_str + ' ' + clim_title,
        fontsize=14, fontweight='bold', y=1.05)

    # save figure
    plt.show()
    png_out = out_dir / f'{vn}_{month_str}_clim_and_anoms_grid.png'
    fig.savefig(png_out, dpi=180, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {png_out}')