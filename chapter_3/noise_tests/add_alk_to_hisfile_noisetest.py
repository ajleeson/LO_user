'''
Modify ocean history file to add 
1 kg/m3 of alkalinity in the surface 5 m of the water column

Based on Parker's driver_roms00oae script
and my modify_ocn_forcing script

run add_alk_to_hisfile_noisetest -gtx cas7_t1_x11ab -0 2020.05.30 -ro 5

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

gtagex_new = 'cas7_t1alkNOISE_x11ab'

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
ds_new = ds_og_his.copy(deep=True)


# add alkalinity as a blob ---------------------------------------


# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
mask_rho = grid_ds.mask_rho.values
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values
plon, plat = pfun.get_plon_plat(lon,lat)
lons = lon[0,:]
lats = lat[:,0]


# offshore location closer to model boundaries
site_lon_test = -127.25
site_lat_test = 45.5
# get nearest lat/lon indices to random point selected on google maps
site_y = min(range(len(lats)), key=lambda i: abs(lats[i]-site_lat_test))
site_x = min(range(len(lons)), key=lambda i: abs(lons[i]-site_lon_test))
# get actual lat/lon
site_lon = lons[site_x]
site_lat = lats[site_y]
# # plot to visualize site location
# plt.close('all')
# fig, ax = plt.subplots(1,1,figsize = (6,8))
# ax.pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho), vmin=0, vmax=5, cmap='Blues' )
# ax.scatter(site_lon, site_lat, color='red', s=100, marker='*')
# pfun.dar(ax)
# plt.show()

# Amount to increase alkalinity concentration
dalk = 200 # mmol m-3, same as ROMS units
# find cell location
G, S, T = zrfun.get_basic_info(roms_out_dir / 'ocean_his_0002.nc')
Lon = G['lon_rho'][0,:]
Lat = G['lat_rho'][:,0]
# error checking
if (site_lon < Lon[0]) or (site_lon > Lon[-1]):
    print('ERROR: lon out of bounds')
    sys.exit()
if (site_lat < Lat[0]) or (site_lat > Lat[-1]):
    print('ERROR: lat out of bounds')
    sys.exit()
ix = zfun.find_nearest_ind(Lon, site_lon)
iy = zfun.find_nearest_ind(Lat, site_lat)
# error checking
if G['mask_rho'][iy,ix] == 0:
    print('ERROR: point on land mask. Exiting.')
    sys.exit()

# get cell thickness
# get S for the whole grid
Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
reader = csv.DictReader(open(Sfp))
S_dict = {}
for row in reader:
    S_dict[row['ITEMS']] = row['VALUES']
S = zrfun.get_S(S_dict)
# get cell thickness
h = ds_new['h'].values # height of water column
z_rho, z_w = zrfun.get_z(h, ds_new.zeta[0,:,:].to_numpy(), S) # depth of rho and w points
# get vertical thickness of all cells [m] 
dzr = np.diff(z_w, axis=0) # [z,y,x]
# get thicknesses at test cite
dzr_site = dzr[:,iy,ix]

# # scale alkalinity by thickness of top 2 sigma layers
# nominal_thickness = 10 # [m] nominal thickness of sigma layers
# actual_thickness = np.nansum(dzr_site[-2:]) # thickness of top 2 sigma layers
# dalk_scaled = dalk * nominal_thickness/actual_thickness   

# add alkalinity
# In this case we are adding to the top 2 sigma layers in a single grid cell
# also, span multiple grid cells (15 x 15) to make a blob
pad = 7 # (15-1)/2
dalk_nominal = 100 # mmol m-3, same as ROMS units
ds_new.alkalinity.values[0,-1:,iy-pad:iy+pad+1,ix-pad:ix+pad+1] += dalk_nominal

# plot to visualize site alkalinity dosing
plt.close('all')
fig, ax = plt.subplots(1,1,figsize = (6,8))
diff = ds_new.alkalinity[0,-1,:,:] - ds_og_his.alkalinity[0,-1,:,:]
cs = ax.pcolormesh(plon, plat, diff, cmap='plasma', vmin=0, vmax=120)
cbar = fig.colorbar(cs)
pfun.dar(ax)
plt.show()

################################################################
# Save .nc files
print('Saving {}'.format(date))
ds_new.to_netcdf(str(out_dir) + '/ocean_his_0002.nc')

print('Done')
