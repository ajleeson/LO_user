"""
Compare average bottom DO between multiple years
(Set up to run for 6 years)

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator
import matplotlib.patches as patches
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import pandas as pd
import cmocean
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

##############################################################
##                       USER INPUTS                        ##
##############################################################

remove_straits = True

vn = 'oxygen'

year =  '2014'

# which  model run to look at?
gtagexes = ['cas7_t1_x11ab','cas7_t1noDIN_x11ab'] 

# where to put output figures
out_dir = Ldir['LOo'] / 'chapter_2' / 'figures'
Lfun.make_dir(out_dir)

region = 'Puget Sound'

plt.close('all')

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# open datasets
if remove_straits:
    straits = 'noStraits'
else:
    straits = 'withStraits'

letters = ['(a) ','(b) ','(c) ','(d) ','(e) ','(f) ']
vars = ['NO3','phytoplankton','zooplankton','NH4','LdetritusN','SdetritusN']

# initialize empty dictionaries
NO3_vert_dict = {}
phyto_vert_dict = {}
zoop_vert_dict = {}
NH4_vert_dict = {}
Ldet_vert_dict = {}
Sdet_vert_dict = {}

for gtagex in gtagexes:
    # add ds to dictionary
    ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / (gtagex + '_' + year + '_NPZD_vert_ints_' + straits + '.nc'))
    NO3_vert_int = ds['NO3_vert_int'].values
    phyto_vert_int = ds['phyto_vert_int'].values
    zoop_vert_int = ds['zoop_vert_int'].values
    NH4_vert_int = ds['NH4_vert_int'].values
    LdetritusN_vert_int = ds['LdetritusN_vert_int'].values
    SdetritusN_vert_int = ds['SdetritusN_vert_int'].values
    # if not a leap year, add a nan on feb 29 (julian day 60 - 1 because indexing from 0)
    if np.mod(int(year),4) != 0: 
        NO3_vert_int = np.insert(NO3_vert_int,59,'nan',axis=0)
        phyto_vert_int = np.insert(phyto_vert_int,59,'nan',axis=0)
        zoop_vert_int = np.insert(zoop_vert_int,59,'nan',axis=0)
        NH4_vert_int = np.insert(NH4_vert_int,59,'nan',axis=0)
        LdetritusN_vert_int = np.insert(LdetritusN_vert_int,59,'nan',axis=0)
        SdetritusN_vert_int = np.insert(SdetritusN_vert_int,59,'nan',axis=0)
    if gtagex == 'cas7_t1_x11ab': # july 23, 2014 was missing
        NO3_vert_int = np.insert(NO3_vert_int,205,'nan',axis=0)
        phyto_vert_int = np.insert(phyto_vert_int,205,'nan',axis=0)
        zoop_vert_int = np.insert(zoop_vert_int,205,'nan',axis=0)
        NH4_vert_int = np.insert(NH4_vert_int,205,'nan',axis=0)
        LdetritusN_vert_int = np.insert(LdetritusN_vert_int,205,'nan',axis=0)
        SdetritusN_vert_int = np.insert(SdetritusN_vert_int,205,'nan',axis=0)
    NO3_vert_dict[gtagex] = NO3_vert_int
    phyto_vert_dict[gtagex] = phyto_vert_int
    zoop_vert_dict[gtagex] = zoop_vert_int
    NH4_vert_dict[gtagex] = NH4_vert_int
    Ldet_vert_dict[gtagex] = LdetritusN_vert_int
    Sdet_vert_dict[gtagex] = SdetritusN_vert_int

# get grid cell area
fp = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
PSbox_ds = xr.open_dataset(fp)
DX = (PSbox_ds.pm.values)**-1
DY = (PSbox_ds.pn.values)**-1
DA = DX*DY # get area in m2

# initialize dictionary for vertical integrals [mol]
NO3_vol = {}
phyto_vol = {}
zoop_vol = {}
NH4_vol = {}
Ldet_vol = {}
Sdet_vol = {}
for gtagex in gtagexes:

    NO3_vert_int = NO3_vert_dict[gtagex]
    NO3_vol_timeseries = np.sum(NO3_vert_int * DA, axis=(1, 2)) # [mol]
    NO3_vol[gtagex] = NO3_vol_timeseries

    phyto_vert_int = phyto_vert_dict[gtagex]
    phyto_vol_timeseries = np.sum(phyto_vert_int * DA, axis=(1, 2)) # [mol]
    phyto_vol[gtagex] = phyto_vol_timeseries

    zoop_vert_int = zoop_vert_dict[gtagex]
    zoop_vol_timeseries = np.sum(zoop_vert_int * DA, axis=(1, 2)) # [mol]
    zoop_vol[gtagex] = zoop_vol_timeseries

    NH4_vert_int = NH4_vert_dict[gtagex]
    NH4_vol_timeseries = np.sum(NH4_vert_int * DA, axis=(1, 2)) # [mol]
    NH4_vol[gtagex] = NH4_vol_timeseries

    Ldet_vert_int = Ldet_vert_dict[gtagex]
    Ldet_vol_timeseries = np.sum(Ldet_vert_int * DA, axis=(1, 2)) # [mol]
    Ldet_vol[gtagex] = Ldet_vol_timeseries

    Sdet_vert_int = Sdet_vert_dict[gtagex]
    Sdet_vol_timeseries = np.sum(Sdet_vert_int * DA, axis=(1, 2)) # [mol]
    Sdet_vol[gtagex] = Sdet_vol_timeseries

# get plotting limits based on region
if region == 'Puget Sound':
    # box extracion limits: [-123.29, -122.1, 46.95, 48.93]
    xmin = -123.29
    xmax = -122.1 + 0.1 # to make room for legend key
    ymin = 46.95 - 0.1 # to make room for legend key
    ymax = 48.93

# get grid data
grid_ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
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

##############################################################
##             Plot hypoxic volume time series              ##
##############################################################


# initialize figure
fig, axes = plt.subplots(2,3,figsize=(14,9),sharex=True)
ax = axes.ravel()

# create time vector
startdate = '2020.01.01'
enddate = '2020.12.31'
dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
dates_local = [pfun.get_dt_local(x) for x in dates]

# define linestyles and linewidths and alpha
# first is with loading, second is no-loading
linestyles = ['-','--']
linewidths = [3,1]
alphas = [0.3,1]

# plot timeseries
for j,var_vol in enumerate([NO3_vol,phyto_vol,zoop_vol,NH4_vol,Ldet_vol,Sdet_vol]):

    for i,gtagex in enumerate(gtagexes):
        ax[j].plot(dates_local,var_vol[gtagex],linestyle=linestyles[i],
                   color='black',linewidth=linewidths[i],alpha=alphas[i])


    # format figure
    if j == 0:
        ax[j].text(0.02,0.08,'Natural: dashed', fontsize = 12, ha='left', 
               transform=ax[j].transAxes)
        ax[j].text(0.02,0.03,'Anthropogenic: solid', fontsize = 12, ha='left', 
               transform=ax[j].transAxes,fontweight='bold',color='black',alpha=0.5)
    ax[j].grid(visible=True, axis='both', color='silver', linestyle='--')
    ax[j].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax[j].tick_params(axis='both', labelsize=12)
    if np.mod(j,3) == 0:
        ax[j].set_ylabel(r'Moles of Nitrogen [mol]', fontsize=12)
    ax[j].set_ylim([0,np.nanmax(var_vol[gtagexes[0]])*1.1])
    ax[j].text(0.02,0.92,letters[j] + vars[j], fontsize = 12, ha='left', 
               fontweight='bold', transform=ax[j].transAxes)
    ax[j].set_xlim([dates_local[0],dates_local[-1]])

plt.tight_layout()
