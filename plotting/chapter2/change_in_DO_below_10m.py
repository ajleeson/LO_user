"""

Analyze change in DO for depths deeper than 10 m

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
from matplotlib.ticker import FuncFormatter
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


# Hanning window length
nwin = 20

# years =  ['2014','2015']
# years =  ['2015','2016','2017','2018','2019','2020']
years =  ['2017']#['2015','2016','2017','2018','2019','2020']

# which  model run to look at?
gtagexes = ['cas7_t1_x11ab','cas7_t1noDIN_x11ab'] 

# where to put output figures
out_dir = Ldir['LOo'] / 'chapter_2' / 'figures'
Lfun.make_dir(out_dir)

# regions = ['All Puget Sound']
regions = ['Hood Canal', 'South Sound', 'Whidbey Basin', 'Main Basin', 'All Puget Sound']
# colors = ['hotpink','mediumpurple','dodgerblue','yellowgreen','black']

plt.close('all')

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# read in masks
basin_mask_ds = grid_ds = xr.open_dataset('../../../LO_output/chapter_2/data/basin_masks_from_pugetsoundDObox.nc')
mask_rho = basin_mask_ds.mask_rho.values
mask_hc = basin_mask_ds.mask_hoodcanal.values
mask_ss = basin_mask_ds.mask_southsound.values
mask_wb = basin_mask_ds.mask_whidbeybasin.values
mask_mb = basin_mask_ds.mask_mainbasin.values
mask_ps = basin_mask_ds.mask_pugetsound.values
lon = basin_mask_ds['lon_rho'].values
lat = basin_mask_ds['lat_rho'].values
h = basin_mask_ds['h'].values
plon, plat = pfun.get_plon_plat(lon,lat)

##############################################################
# get average concentration per basin

# initialize empty dictionaries and fill with vertical integrals
DO_vert_dict = {}

# _shallow10m_deep_DO.nc

for year in years:
    for gtagex in gtagexes:
        ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / (gtagex + '_pugetsoundDO_' + year + '_shallow10m_deep_DO.nc'))
        DO_vert_int = ds['deep_DO_mgL'].values
        # add data to dictionaries
        DO_vert_dict[gtagex+year] = DO_vert_int

# # grid cell areas
# fp = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
# box_ds = xr.open_dataset(fp)
# DX = (box_ds.pm.values)**-1
# DY = (box_ds.pn.values)**-1
# DA = DX*DY # get area in m2


# # initialize dictionary for average concentration (volume integrals [mol], normalized by volume)
# DO_vol_norm = {}

# for year in years:
#     for region in regions:

#         # get mask for the region
#         if region == 'Hood Canal':
#             mask = mask_hc
#         elif region == 'South Sound':
#             mask = mask_ss
#         elif region == 'Whidbey Basin':
#             mask = mask_wb
#         elif region == 'Main Basin':
#             mask = mask_mb
#         elif region == 'All Puget Sound':
#             mask = mask_ps

#         # basin volume
#         h_masked = h * mask
#         basin_vol = np.sum(h_masked * DA) # [m3]

#         for gtagex in gtagexes:
#             DO_vert_int = DO_vert_dict[gtagex+year]
#             DO_vert_int_masked = DO_vert_int * mask
#             DO_vol_timeseries = np.sum(DO_vert_int_masked * DA, axis=(1, 2)) # [mol]
#             DO_vol_norm[gtagex+region+year] = DO_vol_timeseries


# DO_timeseries_noloading = []
# DO_timeseries_loading = []

# # get full time series
# for year in years:
#     for region in regions:
#         DO_timeseries_noloading.extend(DO_vol_norm['cas7_t1noDIN_x11ab'+region+year])
#         DO_timeseries_loading.extend(DO_vol_norm['cas7_t1_x11ab'+region+year])

# # get average molar quantity
# DO_mols_daily_avg_noloading = np.nanmean(DO_timeseries_noloading)
# DO_mols_daily_avg_loading = np.nanmean(DO_timeseries_loading)


# print('----------')

# print('Vol-integrated DO in Puget Sound decreased by: {} perc'.format( round(
#     (DO_mols_daily_avg_noloading-DO_mols_daily_avg_loading)/DO_mols_daily_avg_noloading * 100 ,2)))


# # calculate total Puget Sound volume
# h_masked = h * mask_ps
# PugetSound_vol = np.sum(h_masked * DA) # [m3]
# print('Change in concentrations ----------------')

# print('DO loading {} mmol/m3'.format(DO_mols_daily_avg_loading/PugetSound_vol*1000))
# print('DO no-loading {} mmol/m3'.format(DO_mols_daily_avg_noloading/PugetSound_vol*1000))
