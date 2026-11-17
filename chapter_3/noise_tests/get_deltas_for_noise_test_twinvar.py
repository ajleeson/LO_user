"""
To be run on apogee
to get time series of delta DIC, delta alkalinity, 
over the whole domain
"""

# import things
import numpy as np
import xarray as xr
import pandas as pd
import csv
import matplotlib.pylab as plt
from pathlib import Path
from datetime import datetime
import gsw
import PyCO2SYS as pyco2
import cmcrameri.cm as cmc
from lo_tools import Lfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

##############################################################
##                       USER INPUTS                        ##
##############################################################

ds0 = '2020.05.31'
ds1 = '2020.06.30'

# which  model run to look at? twin variables!
gtagex = 'cas7_t1twin_x11aboaetwin'

out_dir = Ldir['LOo'] / 'chapter_3' / 'data'

# where to put output files
Lfun.make_dir(out_dir)

##############################################################
##                    HELPER FUNCTIONS                      ##
##############################################################

def start_ds(ocean_time):
    '''
    Initialize dataset to store processed values
    ocean_time = ocean time vector
    '''
    Ndays = len(ocean_time.values)

    ds = xr.Dataset(data_vars=dict(

        # difference in DIC between perturbation and baseline
        delta_DIC_full    = (['ocean_time'], np.zeros((Ndays))),
        # difference in total alkalinity between perturbation and baseline
        delta_Alk_full    = (['ocean_time'], np.zeros((Ndays))),
        # difference in CO2 flux between perturbation and baseline
        delta_CO2flux_full = (['ocean_time'], np.zeros((Ndays))),
            
        ),
    coords=dict(ocean_time=ocean_time))
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for processed CO2 uptake capacity
    '''
    ds['delta_DIC_full'].attrs['long_name'] = 'Full domain: Difference in DIC between perturbation and baseline'
    ds['delta_DIC_full'].attrs['units'] = 'kmol C'

    ds['delta_Alk_full'].attrs['long_name'] = 'Full domain: Difference in total alkalinity between perturbation and baseline'
    ds['delta_Alk_full'].attrs['units'] = 'kmol'

    ds['delta_CO2flux_full'].attrs['long_name'] = 'Full domain: Difference in CO2 flux between perturbation and baseline'
    ds['delta_CO2flux_full'].attrs['units'] = 'mol m-2 s-1'

    return ds

##############################################################
##                      PROCESS DATA                        ##
##############################################################

print('Processing started...\n')

# get info to find history files
gridname, tag, ex = gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex)
# add more entries to Ldir
Ldir['roms_out'] = Ldir['roms_out']
Ldir['ds0'] = ds0
Ldir['ds1'] = ds1
Ldir['list_type'] = 'average'
# print(Ldir.keys())
# get history files
fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])

# Initialize empty list of datasets (meant for daily datasets)
ds_list = []

# loop through history files for the year
for i,fn in enumerate(fn_list):

    # get data
    date_str = Path(fn).parent.name
    print('    ' + date_str) 
    # get data
    ds = xr.open_dataset(fn)

    # initalize dataset to store daily data for this day
    daily_ds = start_ds(ds['ocean_time'])

    # if this is the first time step  get grid cell areas (time-invariant)
    if i == 0:
        DX = (ds.pm.values)**-1
        DY = (ds.pn.values)**-1
        DA = DX*DY # [m2]
    
    # Get cel vertical thicknesses
    # get S for the whole grid
    Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
    reader = csv.DictReader(open(Sfp))
    S_dict = {}
    for row in reader:
        S_dict[row['ITEMS']] = row['VALUES']
    S = zrfun.get_S(S_dict)
    # get cell thickness. Note, this is the same for both runs because hydrodynamics are identical
    h = ds['h'].values # height of water column
    zeta = ds['zeta'].values
    z_rho, z_w = zrfun.get_z(h, zeta, S)
    dzr = np.diff(z_w, axis=0) # sigma layer thickness [m]

    # Calculate volume of each cell in the whole grid
    vol_all = dzr * DA # [m3]
    # Calculate volume of just the surface layer
    vol_surf = dzr[-1,:,:] * DA # [m3]

    # Get delta alkalinity and delta DIC in each cell
    dTA = ds['dTA'] # [mmol m-3]
    dTIC = ds['dTIC'] # [mmol m-3]
    dTA_kmol = dTA/1000/1000 # [kmol m-3]
    dTIC_kmol = dTIC/1000/1000 # [kmol m-3]

    # Multiply all variables by cell volume to get all final values in units of kmol
    delta_DIC_individualcells = dTIC_kmol * vol_all # [kmol] dims are (t,z,y,x)
    delta_Alk_individualcells = dTA_kmol  * vol_all # [kmol] dims are (t,z,y,x)

    # sum up values to get total volume integral in the full domain
    delta_DIC_full = np.nansum(delta_DIC_individualcells[0,:,:,:]) #[kmol]
    delta_Alk_full = np.nansum(delta_Alk_individualcells[0,:,:,:]) #[kmol]

    # add data to daily dataset
    daily_ds['delta_DIC_full'] = xr.DataArray(delta_DIC_full,
                                coords={'ocean_time': ds['ocean_time'].values},
                                dims=['ocean_time'])
    
    daily_ds['delta_Alk_full'] = xr.DataArray(delta_Alk_full,
                                coords={'ocean_time': ds['ocean_time'].values},
                                dims=['ocean_time'])
    
    # Cast to float32 HERE to keep RAM usage extremely low while building the list
    daily_ds = daily_ds.astype(np.float32)

    # Append to list
    ds_list.append(daily_ds)
    
    # Close the raw dataset to free up memory and prevent file-handle limits
    ds.close()

# PREPARE DATA FOR SAVING
# Combine all of the daily datasets into one dataset
ds_out = xr.concat(ds_list, dim='ocean_time')
# Add your metadata once at the end
ds_out = add_metadata(ds_out)

print('    Saving dataset')
# Create a compression dictionary
comp = dict(zlib=True, complevel=4)
encoding = {var: comp for var in ds_out.data_vars}
# Save with encoding
ds_out.to_netcdf(out_dir / ('coastal_oae_patch_noisetest_twinvar_'+ds0+'_'+ds1+'.nc'), encoding=encoding)

print(ds_out)

print('Done')