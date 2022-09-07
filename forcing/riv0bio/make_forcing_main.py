"""
This is the main program for making the TIDE forcing file.

Test on mac in ipython:

run make_forcing_main.py -g alpe2 -t v0 -r backfill -s continuation -d 2020.01.01 -f riv0bio -test True

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

out_fn = out_dir / 'rivers.nc'
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

# RIVERS -------------------------------------------------------------------------------------

# get the list of rivers and indices for this grid
gri_fn = Ldir['grid'] / 'river_info.csv'
gri_df = pd.read_csv(gri_fn, index_col='rname')
NRIV = len(gri_df)

# Start Dataset
ri_ds = xr.Dataset()

# Add time coordinate
ri_ds['river_time'] = (('river_time',), ot_vec)
ri_ds['river_time'].attrs['units'] = Lfun.roms_time_units
ri_ds['river_time'].attrs['long_name'] = 'river time'

# Add river coordinate
ri_ds['river'] = (('river',), np.arange(1,NRIV+1))
ri_ds['river'].attrs['long_name'] = 'river runoff identification number'

# Add river names
ri_ds['river_name'] = (('river',), list(gri_df.index))
ri_ds['river_name'].attrs['long_name'] = 'river name'

# Add Vshape
vn = 'river_Vshape'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = ('s_rho', 'river')
# For Vtransform = 2, even spacing is a good approximation, and
# we implement this by using 1/N as the fraction in each vertical cell.
Vshape = (1/N) * np.ones((N, NRIV))
ri_ds[vn] = (dims, Vshape)
ri_ds['river_Vshape'].attrs['long_name'] = vinfo['long_name']

# Add position and direction
for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    if vn == 'river_direction':
        ri_ds[vn] = (('river',), gri_df.idir.to_numpy())
    elif vn == 'river_Xposition':
        X_vec = np.nan * np.ones(NRIV)
        # ii = 0
        for ii,rn in enumerate(gri_df.index):
            if gri_df.loc[rn, 'idir'] == 0:
                X_vec[ii] = gri_df.loc[rn, 'col_py'] + 1
            elif gri_df.loc[rn, 'idir'] == 1:
                X_vec[ii] = gri_df.loc[rn, 'col_py']
            # ii += 1
        ri_ds[vn] = (('river',), X_vec)
    elif vn == 'river_Eposition':
        E_vec = np.nan * np.ones(NRIV)
        # ii = 0
        for ii,rn in enumerate(gri_df.index):
            if gri_df.loc[rn, 'idir'] == 0:
                E_vec[ii] = gri_df.loc[rn, 'row_py']
            elif gri_df.loc[rn, 'idir'] == 1:
                E_vec[ii] = gri_df.loc[rn, 'row_py'] + 1
            # ii += 1
        ri_ds[vn] = (('river',), E_vec)
    ri_ds[vn].attrs['long_name'] = vinfo['long_name']
        
# Add transport
vn = 'river_transport'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = (vinfo['time'],) + ('river',)
Q_mat = np.zeros((NT, NRIV))
# ii = 0
for ii,rn in enumerate(gri_df.index):
    if rn == 'creek0':
        Q_mat[:,ii] = 1000 * np.ones(NT) * gri_df.loc[rn, 'isign']
        # You could make the transport a function of time, for example by making
        # dti = pd.DatetimeIndex([dt0, dt1]) and then using a function of
        # dti.dayofyear.
    else:
        # You could add other rivers here
        pass
    # ii += 1
ri_ds[vn] = (dims, Q_mat)
ri_ds[vn].attrs['long_name'] = vinfo['long_name']
ri_ds[vn].attrs['units'] = vinfo['units']

# Add salinity and temperature and biogeochemistry
for vn in ['river_salt', 'river_temp',
 'river_NO3','river_NH4','river_Chlo','river_Phyt',
 'river_Zoop','river_LDeN','river_SDeN','river_LDeC',
 'river_SDeC','river_TIC','river_TAlk','river_Oxyg']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = (vinfo['time'],) + ('s_rho', 'river')

    # values based on averages from Skagit River (411_Skagit R.xlsx)
    if vn == 'river_salt':
        TR_mat = np.zeros((NT, N, NRIV))
    elif vn == 'river_temp':
        TR_mat = 10 * np.ones((NT, N, NRIV))
    elif vn == 'river_NO3':
        TR_mat = 0.09*71.4 * np.ones((NT, N, NRIV))
    elif vn == 'river_NH4':
        TR_mat = 0 * np.ones((NT, N, NRIV))
    elif vn == 'river_Chlo':
        TR_mat = 0 * np.ones((NT, N, NRIV))
    elif vn == 'river_Phyt':
        TR_mat = 0 * np.ones((NT, N, NRIV))
    elif vn == 'river_Zoop':
        TR_mat = 0 * np.ones((NT, N, NRIV))
    elif vn == 'river_LDeN':
        TR_mat = 0 * np.ones((NT, N, NRIV))
    elif vn == 'river_SDeN':
        TR_mat = 0 * np.ones((NT, N, NRIV))
    elif vn == 'river_LDeC':
        TR_mat = 0 * np.ones((NT, N, NRIV))
    elif vn == 'river_SDeC':
        TR_mat = 0 * np.ones((NT, N, NRIV))
    elif vn == 'river_TIC':
        TR_mat = 454.67 * np.ones((NT, N, NRIV))
    elif vn == 'river_TAlk':
        TR_mat = 410.40 * np.ones((NT, N, NRIV))
    elif vn == 'river_Oxyg':
        TR_mat = 11.57*31.26 * np.ones((NT, N, NRIV))
    ri_ds[vn] = (dims, TR_mat)
    ri_ds[vn].attrs['long_name'] = vinfo['long_name']
    ri_ds[vn].attrs['units'] = vinfo['units']
    
# WWTPs --------------------------------------------------------------------------------------------------------------------------------------------------------

# get the list of wwtps and indices for this grid
gwwtp_fn = Ldir['grid'] / 'wwtp_info.csv'
gwwtp_df = pd.read_csv(gwwtp_fn, index_col='wname')
NWWTP = len(gwwtp_df)

# Start Dataset
wwtp_ds = xr.Dataset()

# Add time coordinate
wwtp_ds['river_time'] = (('river_time',), ot_vec)
wwtp_ds['river_time'].attrs['units'] = Lfun.roms_time_units
wwtp_ds['river_time'].attrs['long_name'] = 'river time'

# Add wwtp coordinate
wwtp_ds['river'] = (('river',), np.arange(NRIV+1,NRIV+1+NWWTP))
wwtp_ds['river'].attrs['long_name'] = 'point source identification number'

# Add river names
wwtp_ds['river_name'] = (('river',), list(gwwtp_df.index))
wwtp_ds['river_name'].attrs['long_name'] = 'river name'

# Add Vshape
vn = 'river_Vshape'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = ('s_rho', 'river')
# All discharge coming from the bottom layer
Vshape = np.zeros((N, NWWTP))
# Vshape[0,:] = 1
Vshape[20,:] = 1
wwtp_ds[vn] = (dims, Vshape)
wwtp_ds['river_Vshape'].attrs['long_name'] = vinfo['long_name']

# Add position and direction
for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    if vn == 'river_direction':
        # set point source diretion to enter vertically (2)
        wwtp_dir = 2 * np.ones(NWWTP) 
        wwtp_ds[vn] = (('river',), wwtp_dir)
    elif vn == 'river_Xposition':
        X_vec = np.nan * np.ones(NWWTP)
        for ii,wn in enumerate(gwwtp_df.index):
            X_vec[ii] = gwwtp_df.loc[wn, 'col_py']
        wwtp_ds[vn] = (('river',), X_vec)
    elif vn == 'river_Eposition':
        E_vec = np.nan * np.ones(NWWTP)
        # ii = 0
        for ii,wn in enumerate(gwwtp_df.index):
            E_vec[ii] = gwwtp_df.loc[wn, 'row_py']
            # ii += 1
        wwtp_ds[vn] = (('river',), E_vec)
            # ii += 1
        wwtp_ds[vn] = (('river',), E_vec)
    wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
        
# Add transport
vn = 'river_transport'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = (vinfo['time'],) + ('river',)
Q_mat = np.zeros((NT, NWWTP))
# ii = 0
for ii,rn in enumerate(gwwtp_df.index):
    if rn == 'user_specified_wwtp':
        Q_mat[:,ii] = 3000 * np.ones(NT)
    else:
        # wwtps all have the same flowrate
        # average from WestPoint
        Q_mat[:,ii] = 4.5e3 * np.ones(NT)# 4.5 * np.ones(NT)
    # ii += 1
wwtp_ds[vn] = (dims, Q_mat)
wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
wwtp_ds[vn].attrs['units'] = vinfo['units']

# Add salinity and temperature & biogeochemistry
for vn in ['river_salt', 'river_temp',
 'river_NO3','river_NH4','river_Chlo','river_Phyt',
 'river_Zoop','river_LDeN','river_SDeN','river_LDeC',
 'river_SDeC','river_TIC','river_TAlk','river_Oxyg']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = (vinfo['time'],) + ('s_rho', 'river')

    # values based on averages from West Point (539_West Point.xlsx)
    if vn == 'river_salt':
        TR_mat = np.zeros((NT, N, NWWTP))
    elif vn == 'river_temp':
        TR_mat = 10 * np.ones((NT, N, NWWTP))
    elif vn == 'river_NO3':
        TR_mat = 0 * np.ones((NT, N, NWWTP))# 4.82*71.4 * np.ones((NT, N, NWWTP))
    elif vn == 'river_NH4':
        TR_mat = 0 * np.ones((NT, N, NWWTP))# 4.48*71.4 * np.ones((NT, N, NWWTP))
    elif vn == 'river_Chlo':
        TR_mat = 0 * np.ones((NT, N, NWWTP))
    elif vn == 'river_Phyt':
        TR_mat = 0 * np.ones((NT, N, NWWTP))
    elif vn == 'river_Zoop':
        TR_mat = 0 * np.ones((NT, N, NWWTP))
    elif vn == 'river_LDeN':
        TR_mat = 0 * np.ones((NT, N, NWWTP))
    elif vn == 'river_SDeN':
        TR_mat = 0 * np.ones((NT, N, NWWTP))
    elif vn == 'river_LDeC':
        TR_mat = 0 * np.ones((NT, N, NWWTP))
    elif vn == 'river_SDeC':
        TR_mat = 0 * np.ones((NT, N, NWWTP))
    elif vn == 'river_TIC':
        TR_mat = 2638.45 * np.ones((NT, N, NWWTP))
    elif vn == 'river_TAlk':
        TR_mat = 2000 * np.ones((NT, N, NWWTP))
    elif vn == 'river_Oxyg':
        TR_mat = 5.9*31.26 * np.ones((NT, N, NWWTP))
    wwtp_ds[vn] = (dims, TR_mat)
    wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
    wwtp_ds[vn].attrs['units'] = vinfo['units']

# --------------------------------------------------------

# Merge the river and wwtp datasets
ds = xr.merge([ri_ds, wwtp_ds])
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

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
