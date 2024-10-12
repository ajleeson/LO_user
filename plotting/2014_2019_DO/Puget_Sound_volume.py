
# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
import csv
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

# which  model run to look at?
gtagex = 'cas7_t0_x4b' # long hindcast (anthropogenic)

##############################################################
##              CALCULATE PUGET SOUND VOLUME                ##
##############################################################

# get percent hypoxic volume
fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
ds = xr.open_dataset(fp)
print()
if remove_straits:
        print('    Removing Straits...')
        lat_threshold = 48.14
        lon_threshold = -122.76
        # Create a mask for latitudes and longitudes in the Straits
        mask = (ds['lat_rho'] > lat_threshold) & (ds['lon_rho'] < lon_threshold)
        # Apply the mask
        ds['h'] = xr.where(mask, np.nan, ds['h'])
        ds['pm'] = xr.where(mask, np.nan, ds['pm'])
        ds['pn'] = xr.where(mask, np.nan, ds['pn'])
# print('calculating vertical thickness')
# # get S for the whole grid
# Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
# reader = csv.DictReader(open(Sfp))
# S_dict = {}
# for row in reader:
#     S_dict[row['ITEMS']] = row['VALUES']
# S = zrfun.get_S(S_dict)
# # get cell thickness
# h = ds['h'].values # height of water column
# z_rho, z_w = zrfun.get_z(h, 0*h, S) 
# dzr = np.diff(z_w, axis=0) # vertical thickness of all cells [m]
# print(dzr.shape)

h = ds['h'].values # height of water column
print(h.shape)
# replace 4m water height with nan
h[h==4]=np.nan


Z = h/1000 # water depth in km
DX = (ds.pm.values)**-1
DY = (ds.pn.values)**-1
DA = DX*DY*(1/1000)*(1/1000) # get area, but convert from m^2 to km^2
print(DA.shape)

# calculate volume of Puget Sound
print('Puget Sound volume with straits omitted [km^3]')
vol = np.nansum(Z * DA)
print(vol)