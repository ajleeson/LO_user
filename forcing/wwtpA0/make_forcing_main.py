"""
This is the main program for making the WWTP forcing file.

Test on mac in ipython:

run make_forcing_main.py -g awwtp1 -t v0 -r backfill -s continuation -d 2020.01.01 -f wwtp0 -test True

"""

from pathlib import Path
import sys, os
from datetime import datetime, timedelta

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

date_string = Ldir['date_string']
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + date_string) / Ldir['frc']

import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd

if Ldir['testing']:
    from importlib import reload
    reload(zrfun)

out_fn = out_dir / 'wwtps.nc'
out_fn.unlink(missing_ok=True)

# Make the time vector.  Here I just have two time points, at the start
# and end of the day, but you could have more, e.g. hourly.  You would still
# want the total time to just be one day.
dt0 = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
dt1 = dt0 + timedelta(days=1)
ot_vec = np.array([Lfun.datetime_to_modtime(dt0), Lfun.datetime_to_modtime(dt1)])
NT = len(ot_vec)

S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)
N = S['N']

grid_fn = Ldir['grid'] / 'grid.nc'
G = zrfun.get_basic_info(grid_fn, only_G=True)

# get the list of wwtps and indices for this grid
gri_fn = Ldir['grid'] / 'wwtp_info.csv'
gri_df = pd.read_csv(gri_fn, index_col='wname')
NWWTP = len(gri_df)

# Start Dataset
ds = xr.Dataset()

# Add time coordinate
ds['river_time'] = (('river_time',), ot_vec)
ds['river_time'].attrs['units'] = Lfun.roms_time_units
ds['river_time'].attrs['long_name'] = 'river time'

# Add river coordinate
ds['river'] = (('river',), np.arange(1,NWWTP+1))
ds['river'].attrs['long_name'] = 'point source identification number'

# Add river names
ds['wwtp_name'] = (('river',), list(gri_df.index))
ds['wwtp_name'].attrs['long_name'] = 'wwtp name'

# Add Vshape
vn = 'river_Vshape'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = ('s_rho', 'river')
# For Vtransform = 2, even spacing is a good approximation, and
# we implement this by using 1/N as the fraction in each vertical cell.
Vshape = (1/N) * np.ones((N, NWWTP))
ds[vn] = (dims, Vshape)
ds[vn].attrs['long_name'] = vinfo['long_name'] 

# Add position and direction
for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    if vn == 'river_direction':
        wwtp_dir = 2 * np.ones(NWWTP) # set point source diretion to enter vertically (2)
        ds[vn] = (('river',), wwtp_dir)
    elif vn == 'river_Xposition':
        X_vec = np.nan * np.ones(NWWTP)
        ii = 0
        for wn in gri_df.index:
            # if gri_df.loc[wn, 'idir'] == 0:
            #     X_vec[ii] = gri_df.loc[wn, 'col_py'] + 1
            # elif gri_df.loc[wn, 'idir'] == 1:
            #     X_vec[ii] = gri_df.loc[wn, 'col_py']
            # Kind of just added this case, but need to double check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            X_vec[ii] = gri_df.loc[wn, 'col_py']
            ii += 1
        ds[vn] = (('river',), X_vec)
    elif vn == 'river_Eposition':
        E_vec = np.nan * np.ones(NWWTP)
        ii = 0
        for wn in gri_df.index:
            # if gri_df.loc[wn, 'idir'] == 0:
            #     E_vec[ii] = gri_df.loc[wn, 'row_py']
            # elif gri_df.loc[wn, 'idir'] == 1:
            #     E_vec[ii] = gri_df.loc[wn, 'row_py'] + 1
            # Kind of just added this case, but need to double check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            E_vec[ii] = gri_df.loc[wn, 'row_py']
            ii += 1
        ds[vn] = (('river',), E_vec)
    ds[vn].attrs['long_name'] = vinfo['long_name']
        
# Add transport
vn = 'river_transport'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = (vinfo['time'],) + ('river',)
Q_mat = np.zeros((NT, NWWTP))
ii = 0
for wn in gri_df.index:
    if wn == 'creek0':
        Q_mat[:,ii] = 1000 * np.ones(NT) * gri_df.loc[wn, 'isign']
        # You could make the transport a function of time, for example by making
        # dti = pd.DatetimeIndex([dt0, dt1]) and then using a function of
        # dti.dayofyear.
    else:
        # You could add other rivers here
        Q_mat[:,ii] = 1000 * np.ones(NT)
        pass
    ii += 1
ds[vn] = (dims, Q_mat)
ds[vn].attrs['long_name'] = vinfo['long_name']
ds[vn].attrs['units'] = vinfo['units']

# Add salinity and temperature
for vn in ['river_salt', 'river_temp']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = (vinfo['time'],) + ('s_rho', 'river')
    if vn == 'river_salt':
        TR_mat = np.zeros((NT, N, NWWTP))
    elif vn == 'river_temp':
        TR_mat = 10 * np.ones((NT, N, NWWTP))
    ds[vn] = (dims, TR_mat)
    ds[vn].attrs['long_name'] = vinfo['long_name']
    ds[vn].attrs['units'] = vinfo['units']
    
# Save to NetCDF
ds.to_netcdf(out_fn)
ds.close()

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)