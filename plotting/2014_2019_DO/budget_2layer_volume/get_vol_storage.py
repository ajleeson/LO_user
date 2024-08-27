# Get volume of terminal inlets using previous lowpass box extraction I have created for Puget Sound

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
Ldir['gtagex'] = 'cas7_t0_x4b'

ds0 = '2014.01.01'
ds1 = '2014.12.31'
year = '2014'
Ldir['ds0'] = ds0
in_dir = Ldir['roms_out'] / Ldir['gtagex']
# G, S, T = zrfun.get_basic_info(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')

fn0 = xr.open_dataset(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
lonr = fn0.lon_rho.values
latr = fn0.lat_rho.values
# NX, NY = dx.shape

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_'+ds0+'_'+ds1) / '2layer_volume_storage'
Lfun.make_dir(out_dir)

dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
dt00 = dt0

#########################################################
# Get lowpass Godin filter box extraction for pugetsoundDO
#########################################################

# open box extraction
box_fn = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'box' / ('pugetsoundDO_'+year+'.01.01_'+year+'.12.31.nc')
ds_box = xr.open_dataset(box_fn)

print('Extracting volume throughout Puget Sound')

# # get grid cell area
dx = 1/ds_box.pm.values
dy = 1/ds_box.pn.values
area = dx * dy
z_w = ds_box['z_w'].values
# get grid cell volumes
volume_all = area * np.diff(z_w,axis=1) # m3 (time,s_rho,y,x)

# get grid cell depths
zrho = ds_box['z_rho'].values

#########################################################
#          Calculate shallow and deep volume           ##
#########################################################

stations = ['lynchcove','penn','budd','case','carr']
# create dictionaries with interface depths
interface_dict = dict()


for i,station in enumerate(stations): # enumerate(sta_dict):
    # print status
    print('({}/{}) Working on {}...'.format(i+1,len(stations),station))

    # get interface depth from csv file
    with open('interface_depths.csv', 'r') as f:
        for line in f:
            inlet, interface_depth = line.strip().split(',')
            interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
    z_interface = float(interface_dict[station])

    # initialize empty dataframe for saving
    df = pd.DataFrame()

    #% Get indices of inlet
    seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
    seg_df = pd.read_pickle(seg_name)
    ji_list = seg_df[station+'_p']['ji_list']
    jj_LO = [x[0] for x in ji_list] # y; lat; jj
    ii_LO = [x[1] for x in ji_list] # x; lon; ii
    # get lat and lon corresponding to ii and jj indices
    lat_LO = latr[jj_LO,0]
    lon_LO = lonr[0,ii_LO]
    # get corresponding ii and jj indices in box extraction
    lat_box_all = ds_box['lat_rho'].values[:,0]
    lon_box_all = ds_box['lon_rho'].values[0,:]
    jj = np.zeros(len(jj_LO))
    ii = np.zeros(len(ii_LO))
    for j,lat in enumerate(lat_LO):
        jj[j] = np.where(lat_box_all==lat)[0][0]
    for i,lon in enumerate(lon_LO):
        ii[i] = np.where(lon_box_all==lon)[0][0]
    # convert to array of ints
    jj = jj.astype(int)
    ii = ii.astype(int)
    
    # create boolean arrays for shallow and deep layers
    tmp_zrho = zrho[:,:,jj,ii] # in domain
    ix_shallow = tmp_zrho >= z_interface # shallower than interface
    ix_deep = tmp_zrho < z_interface  # deeper than interface

    # convert from boolean array to int
    ix_shallow = ix_shallow.astype(int)
    ix_deep = ix_deep.astype(int)
    
    # crop volume to shallow and deep layer within inlet
    tmp_V = volume_all[:,:,jj,ii]
    surf_vol = np.nansum(np.nansum(tmp_V*ix_shallow,axis=1),axis=1)
    deep_vol = np.nansum(np.nansum(tmp_V*ix_deep,axis=1),axis=1)

    # create dataframe of values
    df['surface [m3]'] = surf_vol
    df['deep [m3]'] = deep_vol
    # save to pickle file
    df.to_pickle(out_dir / (station + '.p'))
