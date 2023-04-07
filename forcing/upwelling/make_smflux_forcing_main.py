"""
This makes the surface momentum flux forcing files for an analytical run.

Designed to run only as backfill.

Testing:

run make_smflux_forcing_main.py

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
# out_fn = out_path / 'smflux.nc'
# out_fn.unlink(missing_ok=True)

# Make the time vector.
NT = 21

N = 16 # 16 vertical layers

NR = 82 # rows

NC = 42 # cols

# Create fields for the state variables.
vn_list = ['sustr','svstr']

# Add time coordinate
his_ds = xr.open_dataset('../../upwelling-tests/results/roms_his_og.nc')
time_vec = his_ds.ocean_time.values


for vn in vn_list:  
    out_fn = out_path / (vn + '.nc')
    out_fn.unlink(missing_ok=True)
    ds = xr.Dataset()
    vinfo = zrfun.get_varinfo(vn)
    tname =  vinfo['time_name']
    dims = (tname,) + vinfo['space_dims_tup']

    if vn == 'sustr':
        const = 0 # [N/m2] Pascal
        values = const*np.ones((NT, NR, NC)) # shape of u2dvar
        # values[mu2==0] = np.nan

    elif vn == 'svstr':
        const = 0 # [N/m2]
        values = const*np.ones((NT, NR-1, NC)) # shape of v2dvar
        # values[mv2==0] = np.nan

    ds[vn] = (dims, values)

    ds[vn].attrs['units'] = vinfo['units']
    ds[vn].attrs['long_name'] = vinfo['long_name']
    # time coordinate
    ds[tname] = ((tname,), time_vec)
    # ds[tname].attrs['units'] = Lfun.roms_time_units
    ds[tname].attrs['long_name'] = 'ocean time'

    print(ds)

    # and save to NetCDF
    ds.to_netcdf(out_fn)
    ds.close()

# def print_info(fn):
#     print('\n' + str(fn))
#     ds = xr.open_dataset(fn)#, decode_times=False)
#     print(ds)
#     ds.close()

# # Check results
# nc_list = [item + '.nc' for item in vn_list]
# if Ldir['testing']:
#     # print info about the files to the screen
#     for fn in nc_list:
#         print_info(out_dir / fn)
# result_dict['result'] = 'success'
# for fn in nc_list:
#     if (out_dir / fn).is_file():
#         pass
#     else:
#        result_dict['result'] = 'fail'

# *******************************************************

# result_dict['end_dt'] = datetime.now()
# ffun.finale(Ldir, result_dict)
