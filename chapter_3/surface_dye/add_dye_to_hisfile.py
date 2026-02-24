'''
Modify ocean history file to add 
1 kg/m3 of dye in the surface 5 m of the water column

Based on Parker's driver_roms00oae script
and my modify_ocn_forcing script

run add_dye_to_hisfile -gtx cas7_t1_x11ab -0 2020.05.30 -ro 5

Note that this script can only be run for one day at a time

'''

import sys
import shutil
import argparse
from datetime import datetime, timedelta
from time import time
import xarray as xr
import numpy as np
import matplotlib.pylab as plt
import csv
import sys
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

####################################################
# argument parsing

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-gtx', '--gtagex', default='cas7_t1_x11ab', type=str) # e.g. cas7_t1_x11b
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.

args = parser.parse_args()

gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

###################################################

gtagex_new = 'cas7_t1d_x11ad'

##################################################

# get all information from arguments
argsd = args.__dict__
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]


ds0 = args.ds0
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)

date = 'f' + dt0.strftime(Lfun.ds_fmt)

# set output location
out_dir = ('../../../LO_roms/' + gtagex_new + '/' + date)
Lfun.make_dir(out_dir)

# get original history file (use the prior day's history file)

roms_out_dir = Ldir['roms_out'] / Ldir['gtagex'] / date
ds_og_his = xr.open_dataset(roms_out_dir / 'ocean_his_0002.nc')

# make a copy of the original dataset to modify
ds_new = ds_og_his.copy()

# remove biology variables
bio_vars = ['NO3', 'NH4', 'chlorophyll',
            'phytoplankton', 'zooplankton',
            'LdetritusN', 'SdetritusN',
            'LdetritusC', 'SdetritusC',
            'TIC', 'alkalinity', 'oxygen']
for bio_var in bio_vars:
    ds_new = ds_new.drop_vars([bio_var])

# add dye variable to dataset, using a copy of temp as a reference
# and set all values to be zero
dye_conc = xr.zeros_like(ds_new['temp'])
ds_new['dye_01'] = dye_conc
ds_new['dye_01'].attrs['long_name'] = 'Passive tracer dye concentration'
ds_new['dye_01'].attrs['standard_name'] = 'mass_concentration_of_dye_'
ds_new['dye_01'].attrs['units'] = 'kg m-3' 
ds_new['dye_01'].attrs['field'] = 'dye_' 

# determine which cells are in the top 5 m of the water column
# (or roughly 5 m)

# get cell thickness
# get S for the whole grid
Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
reader = csv.DictReader(open(Sfp))
S_dict = {}
for row in reader:
    S_dict[row['ITEMS']] = row['VALUES']
S = zrfun.get_S(S_dict)
# get cell thickness
h = ds_og_his['h'].values # height of water column
z_rho, z_w = zrfun.get_z(h, ds_new.zeta[0,:,:].to_numpy(), S) # depth of rho and w points. zero is the surface
# get vertical thickness of all cells [m] 
dzr = np.diff(z_w, axis=0) # [z,y,x]

# print(dzr[:,600,300])
# print(z_w[:,600,300])

# # test with water column depth > 5 m
# test_y = 600
# test_x = 300

# # test with water column depth < 5 m
# test_y = 600
# test_x = 353

# # test with first layer = 5 m
# test_y = 600
# test_x = 50

# # test land cell
# test_y = 10
# test_x = 600

# # test outflow of skagit river
# test_lat = 48.30
# test_lon = -122.42
# lons = ds_new['lon_rho'].values[0,:]
# lats = ds_new['lat_rho'].values[:,0]
# test_y = min(range(len(lats)), key=lambda i: abs(lats[i]-test_lat)) + 1
# test_x = min(range(len(lons)), key=lambda i: abs(lons[i]-test_lon)) + 1
# print(z_w[:,test_y,test_x])


# testing
print('Adding dye....')
for x in ds_new['xi_rho'].values: #[test_x]: # loop through lon
    print(x)
    for y in ds_new['eta_rho'].values: #[test_y]: # loop through lat
        # get ssh
        zeta = ds_new.zeta.values[0,y,x]
        # determine z_w depth that is nearest to -5 m
        int_num_of_layers = 0
        min_dist_to_5 = 1000 # arbitratily large starting value [m]
        depth_nearest_5 = 0 # [m]
        # loop through all depths
        # reverse so we check from surface to bottom
        # and add ssh, so we aren't indexing from zero
        # (which may not be where the water level is at)
        for z_w_depth in reversed(z_w[:,y,x]-zeta):
            # check how close z_w_depth is to 5
            dist_to_5 = np.abs(5-np.abs(z_w_depth))
            # if we are getting closer to 5, then increment
            if dist_to_5 < min_dist_to_5:
                int_num_of_layers += 1
                min_dist_to_5 = dist_to_5
                depth_nearest_5 = z_w_depth
            # otherwise, we have overshot, and we can exit the loop
            else:
                continue
        int_num_of_layers -= 1 # subtract one because we go from z_w to z_rho indexing
        # scale dye concentration by depth 
        scaled_concentration = 5/np.abs(depth_nearest_5)
        # add dye to the top 5 m of the water column
        ds_new['dye_01'][:,-int_num_of_layers::,y,x] = scaled_concentration
# # testing (only uncomment if testing one lat/lon grid cell at a time)
# print('\nNumber of layers from surface: {}'.format(int_num_of_layers))
# print('\nMinimum distance to 5 m: {} m'.format(min_dist_to_5))
# print('\nDepth of nearest layer boundary to 5 m: {} m'.format(depth_nearest_5))
# print('\n dye_01 values at this location:')
# print(ds_new['dye_01'][:,:,test_y,test_x].values)
# print('\nVertical integral of dye (kg/m2)')
# print(np.nansum(ds_new['dye_01'][:,:,test_y,test_x].values*dzr[:,test_y,test_x]))

# # apply land mask
# print('Applying land mask...')
# ds_new['dye_01'] = ds_new['dye_01'].where(ds_og_his['mask_rho'])

# print(ds_new['dye_01'][:,:,test_y,test_x].values)

# pfun.start_plot(fs=12, figsize=(8,8))
# fig,ax = plt.subplots(1,1)
# x = ds_new['lon_rho'].values
# y = ds_new['lat_rho'].values
# px, py = pfun.get_plon_plat(x,y)
# cs = ax.pcolormesh(px, py, ds_new['dye_01'][0,-1,:,:])
# # add colorbar
# cbar = plt.colorbar(cs,ax=ax, location='bottom')
# plt.show()


# check output
# get dz
G, S, T = zrfun.get_basic_info(roms_out_dir / 'ocean_his_0002.nc')
zr, zw = zrfun.get_z(G['h'],ds_new.zeta[0,:,:].to_numpy(),S)
dz = np.diff(zw, axis=0) # [m]
# get dye concentrations
dye_conc = ds_new['dye_01'][0,:,:,:].to_numpy() # [kg/m3]
# vertical integrals
dye_vert_int = (dye_conc * dz).sum(axis=0) # [kg/m2]
print('Vertical integral of dye should be nominally 5 kg/m2 +/- floating point error')
print('    Min vertical integral: {} kg/m2'.format(np.nanmin(dye_vert_int)))
print('    Max vertical integral: {} kg/m2'.format(np.nanmax(dye_vert_int)))

print('\n=================')
print('Time check:')
print(ds_new.ocean_time)


################################################################
# Save .nc files
print('Saving {}'.format(date))
ds_new.to_netcdf(str(out_dir) + '/ocean_his_0002.nc')

print('Done')
