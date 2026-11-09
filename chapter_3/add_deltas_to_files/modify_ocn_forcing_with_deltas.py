'''
Modify already-generated ocean forcing files
to add a dTIC and dTA variables with values of zero.

Based on Parker's driver_forcing00 script
and Jilian's script that modified nudging to climatology

run modify_ocn_forcing_with_deltas -g cas7 -f ocnG00 -0 2020.05.31 -1 2020.06.30
'''

import sys
import argparse
from datetime import datetime, timedelta
import xarray as xr
# from subprocess import Popen as Po
# from subprocess import PIPE as Pi

from lo_tools import Lfun, zfun

parser = argparse.ArgumentParser()
# # arguments without defaults are required
parser.add_argument('-g', '--gridname', type=str)   # e.g. cas7
parser.add_argument('-f', '--frc', type=str, default='ocnG00')        # ocean forcing only!
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-1', '--ds1', type=str, default='') # is set to ds0 if omitted
parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 0 = Ldir['roms_out1'], etc.


args = parser.parse_args()
Ldir = Lfun.Lstart(gridname=args.gridname)

ds0 = args.ds0
if len(args.ds1) == 0:
    ds1 = ds0
else:
    ds1 = args.ds1
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)

# loop over all days
dt = dt0
while dt <= dt1:

    date = 'f' + dt.strftime(Lfun.ds_fmt)

    # make output directory
    out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / \
        date / (args.frc + 'delta')
    Lfun.make_dir(out_dir)
    
    # get forcing files from Parker's apogee account
    # local pc for testing : ../../../LO_output/forcing/cas7/
    # parker's apogee: ../../../../parker/LO_output/forcing/cas7/
    ds_bry = xr.open_dataset('../../../../parker/LO_output/forcing/cas7/' + date + '/' + args.frc + '/ocean_bry.nc')
    ds_clm = xr.open_dataset('../../../../parker/LO_output/forcing/cas7/' + date + '/' + args.frc + '/ocean_clm.nc')
    # ds_bry = xr.open_dataset('../../../LO_output/forcing/cas7/' + date + '/' + args.frc + '/ocean_bry.nc')
    # ds_clm = xr.open_dataset('../../../LO_output/forcing/cas7/' + date + '/' + args.frc + '/ocean_clm.nc')

    # add dye to the different datasets, using a copy of temp as a reference

    # Ocean boundary conditions
    # North
    dTIC_time_north = xr.zeros_like(ds_bry['temp_north']).rename({'temp_time': 'bry_time'})
    ds_bry['dTIC_north'] = dTIC_time_north
    ds_bry['dTIC_north'].attrs['long_name'] = 'delta TIC northern boundary conditions'
    ds_bry['dTIC_north'].attrs['units'] = 'mmol m-3' 

    dTA_time_north = xr.zeros_like(ds_bry['temp_north']).rename({'temp_time': 'bry_time'})
    ds_bry['dTA_north'] = dTA_time_north
    ds_bry['dTA_north'].attrs['long_name'] = 'delta alkalinity northern boundary conditions'
    ds_bry['dTA_north'].attrs['units'] = 'meq m-3' 

    # South
    dTIC_time_south = xr.zeros_like(ds_bry['temp_south']).rename({'temp_time': 'bry_time'})
    ds_bry['dTIC_south'] = dTIC_time_south
    ds_bry['dTIC_south'].attrs['long_name'] = 'delta TIC southern boundary conditions'
    ds_bry['dTIC_south'].attrs['units'] = 'mmol m-3' 

    dTA_time_south = xr.zeros_like(ds_bry['temp_south']).rename({'temp_time': 'bry_time'})
    ds_bry['dTA_south'] = dTA_time_south
    ds_bry['dTA_south'].attrs['long_name'] = 'delta alkalinity southern boundary conditions'
    ds_bry['dTA_south'].attrs['units'] = 'meq m-3' 

    # East
    dTIC_time_east = xr.zeros_like(ds_bry['temp_east']).rename({'temp_time': 'bry_time'})
    ds_bry['dTIC_east'] = dTIC_time_east
    ds_bry['dTIC_east'].attrs['long_name'] = 'delta TIC eastern boundary conditions'
    ds_bry['dTIC_east'].attrs['units'] = 'mmol m-3' 

    dTA_time_east = xr.zeros_like(ds_bry['temp_east']).rename({'temp_time': 'bry_time'})
    ds_bry['dTA_east'] = dTA_time_east
    ds_bry['dTA_east'].attrs['long_name'] = 'delta alkalinity eastern boundary conditions'
    ds_bry['dTA_east'].attrs['units'] = 'meq m-3'

    # West
    dTIC_time_west = xr.zeros_like(ds_bry['temp_west']).rename({'temp_time': 'bry_time'})
    ds_bry['dTIC_west'] = dTIC_time_west
    ds_bry['dTIC_west'].attrs['long_name'] = 'delta TIC western boundary conditions'
    ds_bry['dTIC_west'].attrs['units'] = 'mmol m-3' 

    dTA_time_west = xr.zeros_like(ds_bry['temp_west']).rename({'temp_time': 'bry_time'})
    ds_bry['dTA_west'] = dTA_time_west
    ds_bry['dTA_west'].attrs['long_name'] = 'delta alkalinity western boundary conditions'
    ds_bry['dTA_west'].attrs['units'] = 'meq m-3'

    # Add dTIC_time coordinate
    ds_bry['bry_time'].attrs = ds_bry['temp_time'].attrs.copy()
    ds_bry['bry_time'].attrs = ds_bry['temp_time'].attrs.copy()

    # Ocean climatology
    dTIC_time = xr.zeros_like(ds_clm['temp']).rename({'temp_time': 'clm_time'})
    ds_clm['dTIC'] = dTIC_time
    ds_clm['dTIC'].attrs['long_name'] = 'delta TIC concentration'
    ds_clm['dTIC'].attrs['units'] = 'mmol m-3' 

    dTA_time = xr.zeros_like(ds_clm['temp']).rename({'temp_time': 'clm_time'})
    ds_clm['dTA'] = dTA_time
    ds_clm['dTA'].attrs['long_name'] = 'delta alkalinity concentration'
    ds_clm['dTA'].attrs['units'] = 'meq m-3'

    # Add dTIC_time coordinate
    ds_clm['clm_time'].attrs = ds_clm['temp_time'].attrs.copy()
    ds_clm['clm_time'].attrs = ds_clm['temp_time'].attrs.copy()


    ################################################################
    # Save .nc files
    print('Saving {}'.format(date))
    ds_bry.to_netcdf(str(out_dir) + '/ocean_bry.nc')
    ds_clm.to_netcdf(str(out_dir) + '/ocean_clm.nc')

    dt += timedelta(days=1)

print('Done')




