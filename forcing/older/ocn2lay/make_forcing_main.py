"""
This makes the ocn forcing files for an analytical run.

Designed to run only as backfill.

Testing:

run make_forcing_main.py -g fsg -t vnow -r backfill -s continuation -d 2020.01.01 -f ocn2lay -test True

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import gsw

from lo_tools import forcing_argfun as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

import xarray as xr
from time import time
import numpy as np
from lo_tools import Lfun, zfun, zrfun, Ofun_nc

if Ldir['testing']:
    verbose = True
    from importlib import reload
    reload(Ofun_nc)
    reload(zrfun)
else:
    verbose = False

# This directory is created, along with Info and Data subdirectories, by ffun.intro()
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gtag'] / ('f' + Ldir['date_string']) / Ldir['frc']

# Datetime of the day we are working on
this_dt = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# get grid and S info, and some sizes
G = zrfun.get_basic_info(Ldir['grid'] / 'grid.nc', only_G=True)
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)
NZ = S['N']; NR = G['M']; NC = G['L']

# Make the time vector.  Here I just have two time points, at the start
# and end of the day, but you could have more, e.g. hourly.  You would still
# want the total time to just be one day.
dt0 = datetime.strptime(Ldir['date_string'], Lfun.ds_fmt)

# Pick a date to start first day initial conditions if different than other days
diffDay1 = False
day1 = datetime.strptime('2020.01.01', Lfun.ds_fmt)

dt1 = dt0 + timedelta(days=1)
ot_vec = np.array([Lfun.datetime_to_modtime(dt0), Lfun.datetime_to_modtime(dt1)])
NT = len(ot_vec)

# Create fields for the state variables.
# This would be the place to create more complex fields, e.g. salt(t,z)
V = dict()
V['zeta'] = np.zeros((NT, NR, NC))
V['ubar'] = np.zeros((NT, NR, NC-1))
V['vbar'] = np.zeros((NT, NR-1, NC))
V['temp'] = 10 * np.ones((NT, NZ, NR, NC))

# Create different first day initial conditions
if diffDay1:
    if dt0 == day1:
        fresh_lat = 45.12 # DEFINE LATITUDE ABOVE WHICH WATER IS FRESH AT T=0 ON THE FIRST RUN DAY
        # difference from lats to desired lat
        lats = G['lat_rho'][:,0]
        lats_diff = abs(lats - fresh_lat*np.ones(np.shape(lats)))
        # lat_index
        lat_index = int(np.where(lats_diff == np.min(lats_diff))[0])
        V['salt'][0,:,lat_index::,:] = 0 * V['salt'][0,:,lat_index::,:]
    
V['salt'] = np.ones((NT, NZ, NR, NC))
# make two layers with salinity. Bottom layer to be more saline
V['salt'][:,round(NZ/2)::,:,:] = 30 * V['salt'][:,round(NZ/2)::,:,:]
V['salt'][:,0:round(NZ/2),:,:] = 33 * V['salt'][:,0:round(NZ/2),:,:]

# Plotting salinity depth profile -------------------------------------------------------
fig,ax1 = plt.subplots()
# get depth values
depths = zrfun.get_z(100 * np.ones((1,1)), np.zeros((1,1)), S)
rho_depths = depths[0]
# get salinity
sal = V['salt'][0,:,0,0]
# get pressure
p = [gsw.p_from_z(depth,44.5) for depth in rho_depths]
# calculate density from pressure and salinity
rho = [gsw.rho(sal[i],10,p[i]) for i in range(len(p))]
# plot everything
col1 = 'tab:purple'
col2 = 'tab:cyan'
ax1.plot(sal,rho_depths, color = col1, linestyle='dashed')
ax1.set_xlabel(r'Salinity ($g \: kg^{-1}$)', color = col1, fontsize = 12)
ax1.tick_params(axis='x', labelcolor = col1, labelsize = 12)
ax1.set_ylabel('Depth (m)', fontsize = 12)
ax1.tick_params(axis = 'y', labelsize = 12)
ax2 = ax1.twiny()
ax2.plot(rho,rho_depths, color = col2)
ax2.set_xlabel(r'Density ($kg \: m^{-3}$)', color = col2, fontsize = 12)
ax2.tick_params(axis = 'x', labelcolor = col2, labelsize = 12)
plt.tight_layout
plt.title('Ocean Salinity and Density Profile',fontsize = 14)
plt.show()
# -------------------------------------------------------------------------------------

V['u'] = np.zeros((NT, NZ, NR, NC-1))
V['v'] = np.zeros((NT, NZ, NR-1, NC))

# Create masks
mr2 = np.ones((NT, NR, NC)) * G['mask_rho'].reshape((1, NR, NC))
mr3 = np.ones((NT, NZ, NR, NC)) * G['mask_rho'].reshape((1, 1, NR, NC))
mu2 = np.ones((NT, NR, NC-1)) * G['mask_u'].reshape((1, NR, NC-1))
mu3 = np.ones((NT, NZ, NR, NC-1)) * G['mask_u'].reshape((1, 1, NR, NC-1))
mv2 = np.ones((NT, NR-1, NC)) * G['mask_v'].reshape((1, NR-1, NC))
mv3 = np.ones((NT, NZ, NR-1, NC)) * G['mask_v'].reshape((1, 1, NR-1, NC))

# Apply masks
V['zeta'][mr2==0] = np.nan
V['ubar'][mu2==0] = np.nan
V['vbar'][mv2==0] = np.nan
V['salt'][mr3==0] = np.nan
V['temp'][mr3==0] = np.nan
V['u'][mu3==0] = np.nan
V['v'][mv3==0] = np.nan

# Write climatology file: first use of zrfun.get_varinfo().
tt0 = time()
out_fn = out_dir / 'ocean_clm.nc'
out_fn.unlink(missing_ok=True)
ds = xr.Dataset()
for vn in V.keys():
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    tname = vinfo['time_name']
    dims = (vinfo['time_name'],) + vinfo['space_dims_tup']
    ds[vn] = (dims, V[vn])
    ds[vn].attrs['units'] = vinfo['units']
    ds[vn].attrs['long_name'] = vinfo['long_name']
    # time coordinate
    ds[tname] = (('ocean_time',), ot_vec)
    ds[tname].attrs['units'] = Lfun.roms_time_units
# and save to NetCDF
Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
ds.to_netcdf(out_fn, encoding=Enc_dict)
ds.close()
print('- Write clm file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# Write initial condition file
tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_ini.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc.make_ini_file(in_fn, out_fn)
print('- Write ini file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# Write boundary file
tt0 = time()
in_fn = out_dir / 'ocean_clm.nc'
out_fn = out_dir / 'ocean_bry.nc'
out_fn.unlink(missing_ok=True)
Ofun_nc.make_bry_file(in_fn, out_fn)
print('- Write bry file: %0.2f sec' % (time()-tt0))
sys.stdout.flush()

def print_info(fn):
    #print('\n' + str(fn))
    ds = xr.open_dataset(fn)#, decode_times=False)
    # print(ds)
    # print("salt = {}".format(V['salt'][0,0,-1,0]))
    ds.close()

# Check results
nc_list = ['ocean_clm.nc', 'ocean_ini.nc', 'ocean_bry.nc']
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
