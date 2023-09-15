"""
This is the main program for making the TIDE forcing file.

Test on mac in ipython:

run make_wwtp_forcing_main.py 

"""

from pathlib import Path
import sys, os
from datetime import datetime, timedelta

# from lo_tools import forcing_argfun2 as ffun

# Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
# result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

# date_string = Ldir['date_string']
out_dir = '../../../LO_user/forcing/upwelling'

import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd

out_path = Path(out_dir)
out_fn = out_path / 'botwwtp-10.nc'
out_fn.unlink(missing_ok=True)

# Make the time vector.
NT = 21 #len(ot_vec)

N = 16 # 16 vertical layers

# RIVERS -------------------------------------------------------------------------------------

# create a single point source in the center of the upwelling grid (21,40)
NWWTP = 1

# Start Dataset
ds = xr.Dataset()

# Add time coordinate
his_ds = xr.open_dataset('../../upwelling-tests/results/roms_his_og.nc')
time_vec = his_ds.ocean_time.values
# print(time_vec)
ds['river_time'] = (('river_time',), time_vec)
# ds['river_time'].attrs['units'] = Lfun.roms_time_units
ds['river_time'].attrs['long_name'] = 'river time'

# Add river coordinate
ds['river'] = (('river',), np.arange(1,NWWTP+1))
ds['river'].attrs['long_name'] = 'river runoff identification number'

# # Add river names
rname = ['WWTP']
ds['river_name'] = (('river',), rname)
ds['river_name'].attrs['long_name'] = 'river name'

# Add Vshape
vn = 'river_Vshape'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = ('s_rho', 'river')
# All discharge coming from the bottom layer
Vshape = np.zeros((N, NWWTP))
Vshape[0,:] = 1
ds[vn] = (dims, Vshape)
ds['river_Vshape'].attrs['long_name'] = vinfo['long_name']

# Add position and direction
for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    if vn == 'river_direction':
        direc = [2] # vertical sources have direction 2
        ds[vn] = (('river',), direc)
    elif vn == 'river_Xposition':
        X_pos = [21]
        ds[vn] = (('river',), X_pos)
    elif vn == 'river_Eposition':
        E_pos = [40]
        ds[vn] = (('river',), E_pos) #E_vec)
    ds[vn].attrs['long_name'] = vinfo['long_name']
        
# Add transport
vn = 'river_transport'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = (vinfo['time'],) + ('river',)
Q_mat = np.zeros((NT, NWWTP))
Q_mat[:,0] = -10 * np.ones(NT)
ds[vn] = (dims, Q_mat)
ds[vn].attrs['long_name'] = vinfo['long_name']
ds[vn].attrs['units'] = vinfo['units']

print(ds)

# Save to NetCDF
ds.to_netcdf(out_fn)
ds.close()

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

