#%%
"""
Plots a timeseries of different variables at mooring stations to compare different model runs.

Inputs: list of gtagex, label names, moor job name (As defined in job_lists), start date, end date

The labels names should be the distinguishin factor between the gtagexes.

I need to define these inputs at line 38 in this document. 

Then, from the terminal, I can simply run: python mooring_timeseries.py

"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

#%% DEFINE INPUTS ----------------------------------------------------------------------------

gtagexes = ['alpe2_vUw0_atme0', 'alpe2_vUw6_atme0', 'alpe2_vUwneg6_atme0']
labels = ['No Wind','6 m/s Eastward', '6 m/s Westward']
jobname = 'wind_ana'
startdate = '2020.02.01'
enddate = '2020.02.15'
dayafterend = '2020.02.16' # for making date arrays. The next day after we have no more data

#%%  PLOT THE MOORING LOCATIONS -----------------------------------------

# parse gtagex
gtagexample = gtagexes[0]
gridname, tag, ex_name = gtagexample.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)

# load grid data
grid_nc = '../../LO_data/grids/' + gridname + '/grid.nc'
ds = xr.open_dataset(grid_nc)
z = -ds.h.values
mask_rho = ds.mask_rho.values

lon = ds.lon_rho.values
lat = ds.lat_rho.values

plon, plat = pfun.get_plon_plat(lon,lat)
pad = 0.05*(plat[-1,0]-plat[0,0])
ax_lims = (plon[0,0]-pad, plon[0,-1]+pad, plat[0,0]-pad, plat[-1,0]+pad)

# make a version of z with nans where masked
zm = z.copy()
zm[mask_rho == 0] = np.nan

# Initialize plots
plt.close('all')
pfun.start_plot(figsize=(12,10))

# plot bathymetry
fig = plt.figure(figsize=(8,7))
ax = fig.add_subplot(111)
cs = ax.pcolormesh(plon, plat, zm, vmin=-20, vmax=0, cmap=plt.get_cmap(cmocean.cm.deep_r))
cbar = fig.colorbar(cs, ax=ax)
cbar.ax.set_title('Depth (m)', fontsize = 12)
ax.axis(ax_lims)
ax.set_title('Mooring Locations')
ax.set_ylabel('Latitude', fontsize = 16)
ax.set_xlabel('Longitude', fontsize = 16)

# plot mooring locations
moor_lats = [x[1] for x in sta_dict.values()]
moor_lons= [x[0] for x in sta_dict.values()]
print(moor_lats)
print(moor_lons)
plt.scatter(moor_lons,moor_lats, s=75, c='hotpink', marker='*')

# label mooring locations
for i,station in enumerate(sta_dict.keys()):
    ax.text(moor_lons[i],moor_lats[i]-0.1,station, fontsize=12, c='pink', horizontalalignment='center')


#%% CREATE THE TIMESERIES --------------------------------------------------------------------

# create time vector
dates = pd.date_range(start= startdate, end= dayafterend, freq= '1H')
dates_local = [pfun.get_dt_local(x) for x in dates]

# Loop through all of the mooring stations
for i,station in enumerate(sta_dict.keys()):

    print('creating new plots')

    # initalize a new figure
    fig,ax = plt.subplots(4,1,figsize=(12,8), sharex = True)
    # plt.tight_layout()
    ax = ax.ravel()
    ax[0].set_title(station)

    # label variables
    ax[0].set_ylabel(r'E-W Velocity ($m s^{-1}$)')
    ax[1].set_ylabel('Salinity')
    ax[2].set_ylabel(r'$\zeta$ (m)')
    ax[3].set_ylabel('E-W Wind Stress (Pa)')

    # loop through all of the model scenarios
    for j,gtagex in enumerate(gtagexes):
        # print('{},{}'.format(gtagex,station))

        # download .nc files
        fn = '../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
        # print(fn)
        ds = xr.open_dataset(fn)

        # depth layer
        dl = 29

        # plot E-W velocities
        ax[0].plot(dates_local, ds.u.values[:,dl], label=labels[j])
        ax[0].legend(loc = 'best')

        # plot salinities
        ax[1].plot(dates_local, ds.salt.values[:,dl])

        # plot free surface
        ax[2].plot(dates_local, ds.zeta.values)

        # E-W wind stress
        ax[3].plot(dates_local, ds.sustr.values)

        ax[3].xaxis.set_major_formatter(mdates.DateFormatter("%b-%d"))


#%% GENERATE PLOTS ------------------------------------------------------------------------------
plt.show()