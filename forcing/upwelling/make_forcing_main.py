"""
This is the main program for making the TIDE forcing file.

Test on mac in ipython:

run make_forcing_main.py -g alpe2 -t v0 -r backfill -s continuation -d 2020.01.01 -f rivps0 -test True

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

# if Ldir['testing']:
#     from importlib import reload
#     reload(zrfun)

out_path = Path(out_dir)
out_fn = out_path / 'wwtp.nc'
out_fn.unlink(missing_ok=True)

# Make the time vector.  Here I just have two time points, at the start
# and end of the day, but you could have more, e.g. hourly.  You would still
# want the total time to just be one day.
# dt0 = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
# dt1 = dt0 + timedelta(days=1)
# ot_vec = np.array([Lfun.datetime_to_modtime(dt0), Lfun.datetime_to_modtime(dt1)])
NT = 21 #len(ot_vec)

# S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
# S = zrfun.get_S(S_info_dict)
N = 16 # 16 vertical layers

# grid_fn = Ldir['grid'] / 'grid.nc'
# G = zrfun.get_basic_info(grid_fn, only_G=True)

# RIVERS -------------------------------------------------------------------------------------

# create a single point source in the center of the upwelling grid (21,40)
NWWTP = 1

# Start Dataset
ds = xr.Dataset()

# # Add time coordinate
# ds['river_time'] = (('river_time',), ot_vec)
# ds['river_time'].attrs['units'] = Lfun.roms_time_units
# ds['river_time'].attrs['long_name'] = 'river time'

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
# For Vtransform = 2, even spacing is a good approximation, and
# we implement this by using 1/N as the fraction in each vertical cell.
Vshape = (1/N) * np.ones((N, NWWTP))
ds[vn] = (dims, Vshape)
ds['river_Vshape'].attrs['long_name'] = vinfo['long_name']

# Add position and direction
for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    if vn == 'river_direction':
        direc = [2] # vertical sources have direction 2
        ds[vn] = (('river',), direc)
        # ri_ds[vn] = (('river',), gri_df.idir.to_numpy())
    elif vn == 'river_Xposition':
        X_pos = [21]
        # X_vec = np.nan * np.ones(NRIV)
        # ii = 0
        # for ii,rn in enumerate(gri_df.index):
        #     if gri_df.loc[rn, 'idir'] == 0:
        #         X_vec[ii] = gri_df.loc[rn, 'col_py'] + 1
        #     elif gri_df.loc[rn, 'idir'] == 1:
        #         X_vec[ii] = gri_df.loc[rn, 'col_py']
        #     # ii += 1
        ds[vn] = (('river',), X_pos) #X_vec)
    elif vn == 'river_Eposition':
        E_pos = [40]
        # E_vec = np.nan * np.ones(NRIV)
        # ii = 0
        # for ii,rn in enumerate(gri_df.index):
        #     if gri_df.loc[rn, 'idir'] == 0:
        #         E_vec[ii] = gri_df.loc[rn, 'row_py']
        #     elif gri_df.loc[rn, 'idir'] == 1:
        #         E_vec[ii] = gri_df.loc[rn, 'row_py'] + 1
            # ii += 1
        ds[vn] = (('river',), E_pos) #E_vec)
    ds[vn].attrs['long_name'] = vinfo['long_name']
        
# Add transport
vn = 'river_transport'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = (vinfo['time'],) + ('river',)
Q_mat = np.zeros((NT, NWWTP))
Q_mat[:,0] = 1 * np.ones(NT)
# ii = 0
# for ii,rn in enumerate(gri_df.index):
#     if rn == 'creek0':
#         Q_mat[:,ii] = 10 * np.ones(NT) * gri_df.loc[rn, 'isign']
#         # You could make the transport a function of time, for example by making
#         # dti = pd.DatetimeIndex([dt0, dt1]) and then using a function of
#         # dti.dayofyear.
#     else:
#         # You could add other rivers here
#         pass
    # ii += 1
ds[vn] = (dims, Q_mat)
ds[vn].attrs['long_name'] = vinfo['long_name']
ds[vn].attrs['units'] = vinfo['units']

# Add salinity and temperature
# for vn in ['river_salt', 'river_temp']:
#     vinfo = zrfun.get_varinfo(vn, vartype='climatology')
#     dims = (vinfo['time'],) + ('s_rho', 'river')
#     if vn == 'river_salt':
#         TR_mat = np.zeros((NT, N, NRIV))
#     elif vn == 'river_temp':
#         TR_mat = 10 * np.ones((NT, N, NRIV))
#     ri_ds[vn] = (dims, TR_mat)
#     ri_ds[vn].attrs['long_name'] = vinfo['long_name']
#     ri_ds[vn].attrs['units'] = vinfo['units']

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

# *******************************************************

# result_dict['end_dt'] = datetime.now()
# ffun.finale(Ldir, result_dict)
