"""
Create map of terminal inlets
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cmocean
import pandas as pd
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun

Ldir = Lfun.Lstart()

year = '2017'
# set up dates
startdate = year + '.01.01'
enddate = year + '.12.31'

jobname = 'twentyoneinlets'
# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')
# Get stations:
sta_dict = job_lists.get_sta_dict(jobname)

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../../../LO_data/grids/cas7/grid.nc')
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

# Create map
plt.close('all')
fig = plt.figure(figsize=(8,9))
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

# get terminal inlet locations
for stn,station in enumerate(sta_dict): # stations: 
    # get segment information
    seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
    seg_df = pd.read_pickle(seg_name)
    ji_list = seg_df[station+'_p']['ji_list']
    jj = [x[0] for x in ji_list]
    ii = [x[1] for x in ji_list]
    # set everything to nan that is not in the inlet
    # first make everything a nan
    inlet_loc = np.full(zm.shape, np.nan) 
    # set values of 1 for everything that is in the inlet
    inlet_loc[jj,ii] = 40
    # add inlet locations
    plt.pcolormesh(plon, plat, inlet_loc, linewidth=0.5, vmin=0, vmax=100, cmap=plt.get_cmap('spring'))

# format
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)
ax.set_xlim(-123.5, -122.1) #-123.29, -122.1) # Puget Sound
ax.set_ylim(46.95, 48.46) # Puget Sound

pfun.dar(ax)

# remove border
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)


# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
Lfun.make_dir(out_dir)
plt.subplots_adjust(left=0.05, top=0.9, bottom=0.05, right=0.95)
plt.savefig('terminal_inlet_map.png',transparent=True)
