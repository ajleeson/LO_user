"""
Plot OAE time series
"""

# import things
import numpy as np
import xarray as xr
import csv
import matplotlib.pylab as plt
import matplotlib.colors as colors
from pathlib import Path
from datetime import datetime
import pandas as pd
import matplotlib.dates as mdates
import gsw
from matplotlib.patches import Rectangle
import PyCO2SYS as pyco2
import cmcrameri.cm as cmc
from lo_tools import Lfun, zrfun, zfun
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

# which  model runs to look at?
baseline = 'cas7_t1_x11ab'
perturbation = 'cas7_t1alkNOISE_x11ab'

# get lat and lon
date = '2020.05.31'
ds_base = xr.open_dataset(Ldir['roms_out'] / baseline / ('f' + date) / 'ocean_his_0001.nc')
lons = ds_base.coords['lon_rho'].values
lats = ds_base.coords['lat_rho'].values

##############################################################
##                  Get time series data                    ##
##############################################################

# get the deltas for the baseline and perturbed runs
fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'coastal_oae_patch_noisetest_ogmodel_2020.05.31_2020.06.30.nc'
ds_ogmodel = xr.open_dataset(fp)

time = ds_ogmodel.ocean_time.values
delta_DIC_ogmodel = ds_ogmodel.delta_DIC_full.values # kmol
delta_Alk_ogmodel = ds_ogmodel.delta_Alk_full.values # kmol

###################################################################
##                         Plotting                              ##  
################################################################### 

# get grid info
px, py = pfun.get_plon_plat(lons,lats)

# Initialize figure
fig, ax = plt.subplots(1,1, figsize=(8, 5))

# calculate total amount of alkalinity added -----------------------------
# get the initialization his file
ds_init = xr.open_dataset(Ldir['roms_out'] / perturbation / ('f' + date) / 'ocean_his_0001.nc')
# get difference in alkalinity between the perturbed and baseline runs
delta_Alk_init = ds_init['alkalinity'].values - ds_base['alkalinity'].values # [meq/m3]
# get grid cell volume...
DX = (ds_base.pm.values)**-1
DY = (ds_base.pn.values)**-1
DA = DX*DY # [m2]
# Get cel vertical thicknesses
# get S for the whole grid
Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
reader = csv.DictReader(open(Sfp))
S_dict = {}
for row in reader:
    S_dict[row['ITEMS']] = row['VALUES']
S = zrfun.get_S(S_dict)
# get cell thickness. Note, this is the same for both runs because hydrodynamics are identical
h = ds_base['h'].values # height of water column
zeta = ds_base['zeta'].values
z_rho, z_w = zrfun.get_z(h, zeta, S)
dzr = np.diff(z_w, axis=0) # sigma layer thickness [m]
# Calculate volume of each cell in the whole grid
vol_all = dzr * DA # [m3]
# multiply delta alk by volume to get total alk added
total_alk_added = np.nansum(delta_Alk_init * vol_all)/1000/1000 # [kmol eq]


# plot efficiency (delta DIC / total alkalinity added) over time
ax.plot(time,delta_DIC_ogmodel/total_alk_added, linewidth=3, color='royalblue',alpha=0.3, label='OG model')

# format figure
ax.set_ylim([0,1])
ax.set_ylabel(r'$\eta = \Delta$ DIC / $\Delta$ Alk', fontsize=14)
ax.grid(True, color='silver', linestyle=':')
ax.tick_params(axis='both', labelsize=12, rotation=30)  
loc = mdates.MonthLocator(interval=1)
ax.xaxis.set_major_locator(loc)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
ax.set_xlim([np.datetime64('2020-05-31'), np.datetime64('2020-06-30')])