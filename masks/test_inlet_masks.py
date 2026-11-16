"""
Test to make sure inlet masks worked
"""

# import things
import numpy as np
import xarray as xr
import matplotlib.pylab as plt

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun


Ldir = Lfun.Lstart()

##############################################################
##                       USER INPUTS                        ##
##############################################################

# DO threshold for calculating days of bottom DO < threshold,
# and bottom area with DO < threshold
DO_thresh = 2 # [mg/L]

# Show WWTP locations?
WWTP_loc = False

remove_straits = True

vn = 'oxygen'

years =  ['2015','2016','2017','2018','2019','2020']

# which  model run to look at?
gtagex = 'cas7_t1_x11ab' # long hindcast (anthropogenic)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures'
Lfun.make_dir(out_dir)

region = 'Puget Sound'

plt.close('all')

##############################################################
##                      PROCESS DATA                        ##
##############################################################

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

# read in masks
basin_mask_ds = xr.open_dataset('../../LO_output/chapter_2/data/basin_masks_from_pugetsoundDObox.nc')
mask_rho = basin_mask_ds.mask_rho.values
mask_ps = basin_mask_ds.mask_pugetsound.values


# get plotting limits based on region
if region == 'Puget Sound':
    # box extracion limits: [-123.29, -122.1, 46.95, 48.93]
    xmin = -123.29
    xmax = -122.1
    ymin = 46.95
    # ymax = 48.93
    ymax=48.5

colors = ['#62B6CB','#A8C256','#96031A','#957FEF','#F9627D','#476ad1','darkorange']

##############################################################
##                HYPOXIC VOLUME TIMESERIES                 ##
##############################################################

# initialize figure
plt.close('all)')
f, ax0 = plt.subplots(1,1,figsize =(6,8))

# format figure
# ax0.set_xlim([xmin,xmax])
# ax0.set_ylim([ymin,ymax])
# ax0.set_yticklabels([])
# ax0.set_xticklabels([])
# ax0.axis('off')
ax0.set_ylabel('Latitude', fontsize=12)
ax0.set_xlabel('Longitude', fontsize=12)
ax0.tick_params(axis='both', labelsize=12)
# ax0.pcolormesh(plon, plat, zm, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
lon = basin_mask_ds.lon_rho.values
lat = basin_mask_ds.lat_rho.values
plon, plat = pfun.get_plon_plat(lon,lat)
ax0.pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho),
                vmin=0, vmax=10, cmap='Blues')


# loop through inlets
inlets = ['dyes','sinclair','quartermaster','case','crescent','carr',
          'elliot','commencement','penn','portsusan','holmes','dabob','lynchcove']
inlet_mask_ds = xr.open_dataset('../../LO_output/hypVol_for_intrmdl_cmprsn/inlet_masks_from_seg_file.nc')
for inlet in inlets:
    mask_inlet = inlet_mask_ds[inlet].values
    ax0.pcolormesh(plon, plat, np.where(mask_inlet == 0, np.nan, mask_inlet),
                vmin=0, vmax=2, cmap='Blues')

pfun.dar(ax0)