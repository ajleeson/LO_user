"""
This makes the surface momentum flux forcing files for an analytical run.

Designed to run only as backfill.

Testing:

run make_forcing_main.py -g fsg -t v0 -r backfill -s continuation -d 2020.01.01 -f upw -test True

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta
import matplotlib.pyplot as plt

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import xarray as xr
from time import time
import numpy as np
from lo_tools import Lfun, zfun, zrfun

if Ldir['testing']:
    from importlib import reload
    reload(zrfun)

# This directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']

# get grid and S info, and some sizes
G = zrfun.get_basic_info(Ldir['grid'] / 'grid.nc', only_G=True)
NR = G['M']; NC = G['L']

# Make the time vector.  Make hourly time points. This vector must span just one day. 
dt0 = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)
ot_vec = [Lfun.datetime_to_modtime(dt0)]
for i in range(0,24):
    dt1 = dt0 + timedelta(hours=1)
    ot_vec.append(Lfun.datetime_to_modtime(dt1))
    dt0 = dt1

NT = len(ot_vec)

# Create fields for the state variables.
vn_list = ['sustr','svstr']

# create a mask
# mu2 = np.ones((NT, NR, NC-1)) * G['mask_rho'].reshape((1, NR, NC-1)) # shape of u2dvar
# mv2 = np.ones((NT, NR-1, NC)) * G['mask_rho'].reshape((1, NR-1, NC)) # shape of v2dvar

for vn in vn_list:  
    out_fn = out_dir / (vn + '.nc')
    out_fn.unlink(missing_ok=True)
    ds = xr.Dataset()
    vinfo = zrfun.get_varinfo(vn)
    tname =  vinfo['time_name']
    dims = (tname,) + vinfo['space_dims_tup']

    if vn == 'sustr':
        const = 0.1 # [N/m2] Pascal
        values = const*np.ones((NT, NR, NC-1)) # shape of u2dvar
        # values[mu2==0] = np.nan

    elif vn == 'svstr':
        const = 0 # [N/m2]
        values = const*np.ones((NT, NR-1, NC)) # shape of v2dvar
        # values[mv2==0] = np.nan

    ds[vn] = (dims, values)

    ds[vn].attrs['units'] = vinfo['units']
    ds[vn].attrs['long_name'] = vinfo['long_name']
    # time coordinate
    ds[tname] = ((tname,), ot_vec)
    ds[tname].attrs['units'] = Lfun.roms_time_units
    ds[tname].attrs['long_name'] = 'ocean time'
    # and save to NetCDF
    Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
    ds.to_netcdf(out_fn, encoding=Enc_dict)
    ds.close()

def print_info(fn):
    print('\n' + str(fn))
    ds = xr.open_dataset(fn)#, decode_times=False)
    print(ds)
    ds.close()

# Check results
nc_list = [item + '.nc' for item in vn_list]
if Ldir['testing']:
    # print info about the files to the screen
    for fn in nc_list:
        print_info(out_dir / fn)
result_dict['result'] = 'success'
for fn in nc_list:
    if (out_dir / fn).is_file():
        pass
    else:
       result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
