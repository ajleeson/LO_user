"""
Plot model domain properties for OAE site selection
"""

###################################################################
##                       import packages                         ##  
###################################################################      

import copy
import numpy as np
import xarray as xr
import cmocean
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import cmcrameri.cm as cmc

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()
Ldir = Lfun.Lstart()

plt.close('all')

###################################################################
##                          User Inputs                          ##  
################################################################### 

date = '2020.05.22'
date_formatted = '2020-05-22'

###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

gtagex = 'cas7_t1d_x11ad'

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values

# lon/lat limits (Study Domain)
xmin = -126
xmax = -122
ymin = 45.5
ymax = 50.5

# get model output
fp = Ldir['roms_out'] / gtagex / ('f'+date) / 'ocean_avg_0001.nc'
ds_dye = xr.open_dataset(fp)

# get surface mixed layer depth
fp = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'sml_plus' / 'threshold_p125' / 'LO_domain_2020.01.01_2020.12.31' / 'LO_domain_sml_plus_2020.01.01_2020.12.31.nc'
ds_sml = xr.open_dataset(fp)
ds_sml_may22 = ds_sml.sel(ocean_time=date_formatted)
# print(ds_sml_may22)

# get wind speeds
fp = Ldir['LOo'] / 'forcing' / 'cas7' / ('f'+date) / 'atm00'
ds_Uwind = xr.open_dataset(fp/'Uwind.nc') # wind data is hourly
ds_Vwind = xr.open_dataset(fp/'Vwind.nc') # wind data is hourly
# calculate hourly wind speed
hourly_windspeed = np.sqrt(ds_Uwind['Uwind']**2 + ds_Vwind['Vwind']**2)
# average to get daily average wind speed
windspeed_daily = np.nanmean(hourly_windspeed, axis=0)
# apply land mask
mask_rho = grid_ds.mask_rho.values   

# get pH in surface mixed layer
fp = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'sml_plus' / 'threshold_p125' / 'pH' / 'pH_2020_sml_p125.nc'
ds_pH = xr.open_dataset(fp)


###################################################################
##                     Get values to plot                        ##  
################################################################### 

# get pcolormesh values
# surface dye
surf_dye_01 = ds_dye['dye_01'][0,-1,:,:].values
# surface mixed layer depth
sml = ds_sml_may22['SML_thickness'][0,:,:].values
# salinity in surface mixed layer
SA = ds_sml_may22['SA'][0,:,:].values
# temperature in surface mixed layer
CT = ds_sml_may22['CT'][0,:,:].values
# surface wind speed (daily average)
windspeed_daily = np.ma.masked_where(mask_rho==0, windspeed_daily)
# pH in surface mixed layer
pH_may22 = ds_pH['pH'].values[142,:,:]
# alkalinity in surface mixed layer
alk = ds_sml_may22['ALK'][0,:,:].values

# make list of vars
vars = [surf_dye_01, sml, SA, CT, windspeed_daily, pH_may22, alk]

# list of vmins and vmax
vmins = [0,     0,  25, 8,  0,  7.5, 1000]
vmaxs = [0.5,   35, 33, 16, 12, 8.5, 2100]

# get colormaps
cmaps = [cmc.devon_r, cmc.batlowW_r, cmc.nuuk, cmc.roma_r,
         cmc.davos, cmc.bam, cmc.acton]

# titles
titles = ['Surface dye\n'+r'concentration [kg m$^{-3}$]',
          'Surface mixed layer\n'+ r'depth [m]',
          'Surface mixed layer\n'+ r'abs. salinity [g kg$^{-1}$]',
          'Surface mixed layer\n'+ r'cons. temp [$\degree$C]',
          'Surface wind\n'+ r'speed [m s$^{-1}$]',
          'Surface mixed layer\npH',
          'Surface mixed layer\n'+ r'total alkalinity [meq m$^{-3}$]']

# # labels
# labels = ['(a) ','(b) ', '(c) ', '(d) ', '(e) ']


###################################################################
##                         Make plots                            ##  
################################################################### 

px, py = pfun.get_plon_plat(lon,lat)

lons = lon[0,:]
lats = lat[:,0]

# for testing: print index of grid cell closest to a lat/lon from google maps
# print(min(range(len(lats)), key=lambda i: abs(lats[i]-48.298809)))
# print(min(range(len(lons)), key=lambda i: abs(lons[i]+123.005667)))


for var, vmin, vmax, cmap, title in zip(vars[::-1], vmins[::-1], vmaxs[::-1], cmaps[::-1], titles[::-1]):


    # Initialize figure
    fig,ax = plt.subplots(1,1, figsize=(7,9))

    # plot values
    cs = ax.pcolormesh(px,py,var,vmin=vmin, vmax=vmax, cmap=cmap)

    # Add Puget Sound Inset
    # [x0, y0, width, height]
    axins = ax.inset_axes([0.71, 0.0, 0.45, 0.6])
    # plot values in inset
    axins.pcolormesh(px, py, var, vmin=vmin, vmax=vmax, cmap=cmap)
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

    # add colorbar
    cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
    cbar.ax.tick_params(labelsize=18, rotation=30)
    cbar.outline.set_visible(False)

    # format figure
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.tick_params(left=False, bottom=False)
    pfun.add_coast(ax, color='silver')
    pfun.add_coast(axins, color='silver')
    pfun.dar(ax)
    pfun.dar(axins)
    ax.set_title(title, fontsize=20,
                loc='Left', fontweight='bold')
    
    print('======================================\n'+title)
    # add testbed sites
    testbed_color = 'crimson'
    testbed_outline = 'pink'
    testbed_lats = [48.129378, 48.078611, 48.5453, 47.584536, 47.60079806]
    testbed_lons = [-123.457488, -123.045000, -123.0121, -122.342908, -122.43131423]
    for lat, lon in zip(testbed_lats, testbed_lons):
        ax.scatter(   lon, lat, color='None', s=80, edgecolor=testbed_outline, linewidth=4, zorder=5, marker='D')
        ax.scatter(   lon, lat, color='None', s=80, edgecolor=testbed_color,   linewidth=2, zorder=5, marker='D')
        axins.scatter(lon, lat, color='None', s=80, edgecolor=testbed_outline, linewidth=4, zorder=5, marker='D')
        axins.scatter(lon, lat, color='None', s=80, edgecolor=testbed_color,   linewidth=2, zorder=5, marker='D')

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
            ax.scatter(   lons[x], lats[y], color='None', s=180, edgecolor=natural_outline, linewidth=4, zorder=5)
            ax.scatter(   lons[x], lats[y], color='None', s=180, edgecolor=natural_color,   linewidth=2, zorder=5)

        
        # print values
        print('{} ({})'.format(np.round(var[y,x],2), loc))

 

    # Generate plot
    plt.tight_layout
    plt.subplots_adjust(bottom=0.001, top=0.9)
    plt.show()
