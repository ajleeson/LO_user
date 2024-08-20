# Volume terms, separated into 2 layers

# Original file from Jilian: https://github.com/Jilian0717/LO_user/blob/main/tracer_budget/two_layer/get_DO_bgc_air_sea_shallow_deep.py
# modified by me

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from lo_tools import zfun, zrfun, Lfun
from datetime import datetime, timedelta
import scipy
import pickle, sys
import pandas as pd
from time import time
tt0 = time()

#%------------------------------------------------
Ldir = Lfun.Lstart()
# Ldir['roms_out'] = Ldir['roms_out2']
# Ldir['roms_out'] = Ldir['roms_out1']
Ldir['roms_out'] = Ldir['roms_out5'] # for apogee
Ldir['gtagex'] = 'cas7_t0_x4b'

ds0 = '2014.01.01'
ds1 = '2014.01.31'
Ldir['ds0'] = ds0
in_dir = Ldir['roms_out'] / Ldir['gtagex']
G, S, T = zrfun.get_basic_info(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')

fn0 = xr.open_dataset(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
dx = 1/fn0.pm.values
dy = 1/fn0.pn.values
lonr = fn0.lon_rho.values
latr = fn0.lat_rho.values
area = dx * dy
NX, NY = dx.shape

dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
dt00 = dt0

#% load salish sea j,i
seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
seg_df = pd.read_pickle(seg_name)
ji_list = seg_df['lynchcove_p']['ji_list']
jj = [x[0] for x in ji_list]
ii = [x[1] for x in ji_list]

inDomain = np.zeros([NX, NY])
inDomain[jj,ii] = 1 # inside domain, the index=1, outside domian, the index =0, this step is to reduce the size of spatial variable (hopefully)

Vol_sum_shallow = []
Vol_sum_deep = []
t = []


cnt = 0
#%%
while dt00 <= dt1:  # loop each day and every history file
    print(dt00)
    sys.stdout.flush()
    ds00 = dt00.strftime(Lfun.ds_fmt)
    fn_list = Lfun.get_fn_list('hourly', Ldir, ds00, ds00)
    #%%
    for fn in fn_list[0:-1]: 
        print(fn)
        ds = xr.open_dataset(fn)
        zeta = ds.zeta.values.squeeze()
        h = ds.h.values
        zrho = zrfun.get_z(h, zeta, S, only_rho=True)
        
        z_w = zrfun.get_z(h, zeta, S, only_rho=False, only_w=True)
        vol = np.diff(z_w,axis=0) * area # grid cell volume
        
        tmp_zrho = zrho[:,jj,ii] # in domain
        ix_shallow = tmp_zrho>=-6 # shallower than 20 m
        ix_deep = tmp_zrho<-6  # deeper than 20m
        
        #Oxy_vol = Oxy * vol * stat # only account for Salish Sea
        #Oxy_vol_sum.append(np.nansum(Oxy_vol))
        tmp_DOV = vol[:,jj,ii]
        Vol_sum_shallow.append(np.nansum(tmp_DOV[ix_shallow]))
        Vol_sum_deep.append(   np.nansum(tmp_DOV[ix_deep]))
        
        cnt += 1
        t.append(ds.ocean_time.values)
        ds.close()       
    dt00 = dt00 + timedelta(days=1)
        
# save netcdf
from netCDF4 import Dataset
nc = Dataset('O2_bgc_shallow_deep_'+ds0+'_'+ds1+'.nc','w')
time = nc.createDimension('time', len(t))

times = nc.createVariable('time','f8',('time',))
times.units = 'seconds*1e9 since 1970-01-01 00:00:00'
Vol_sum_shallow_tmp = nc.createVariable('Vol_sum_shallow','f4', ('time',),compression='zlib',complevel=9)
Vol_sum_shallow_tmp.units = 'm3'
Vol_sum_deep_tmp = nc.createVariable('Vol_sum_deep','f4', ('time',),compression='zlib',complevel=9)
Vol_sum_deep_tmp.units = 'm3'

times[:] = t

Vol_sum_shallow_tmp[:] = Vol_sum_shallow
Vol_sum_deep_tmp[:] = Vol_sum_deep

nc.close()
