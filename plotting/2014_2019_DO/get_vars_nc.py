"""
This script gets the following variables from box extraction:
surface T
surface S
bottom T
bottom S
depth-integrated LdetN
depth-integrated SdetN
depth-integrated phytoplankton
depth-integrated zooplankton
depth-integrated NO3
depth-integrated NH4
surface NO3
surface NH4

Data are saved in a new .nc file

This script searches for yearly box extractions in LO_output, for the
region "pugetsoundDO"

It also crops out data from the Straits, so as to not bias the results
in Puget Sound. (optional using flag remove_straits)

.nc files are saved in LO_output/pugetsound_DO/data
"""

# import things
import numpy as np
import xarray as xr
import csv
import pinfo
from lo_tools import Lfun, zrfun

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

remove_straits = False

years = ['2014','2015','2016','2017','2018','2019']

# which  model run to look at?
gtagex = 'cas7_t0_x4b' # long hindcast (anthropogenic)

# where to put output files
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'data'
Lfun.make_dir(out_dir)

##############################################################
##                    HELPER FUNCTIONS                      ##
##############################################################

def start_ds(ocean_time,eta_rho,xi_rho):
    '''
    Initialize dataset to store processed DO data
    ocean_time = ocean time vector
    eta_rho = eta rho vector
    xi_rho = xi_rho vector
    '''
    Ndays = len(ocean_time.values)
    Neta = len(eta_rho.values)
    Nxi = len(xi_rho.values)

    ds = xr.Dataset(data_vars=dict(
        # depth of water column
        depth_bot   = (['eta_rho','xi_rho'], np.zeros((Neta,Nxi))),
        # surface temp
        surfT       = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # surface salinity
        surfS       = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # bottom temp
        botT       = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # bottom salinity
        botS       = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # depth-integrated LdetN
        intLdetN    = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # depth-integrated SdetN
        intSdetN    = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # depth-integrated phytoplankton
        intphyto    = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # depth-integrated zooplankton
        intzoop     = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # depth-integrated NO3
        intNO3      = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # depth-integrated NH4
        intNH4     = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # surface NO3
        surfNO3      = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # surface NH4
        surfNH4     = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),),

    coords=dict(ocean_time=ocean_time, eta_rho=eta_rho, xi_rho=xi_rho,),)
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for processed DO data
    '''

    ds['depth_bot'].attrs['long_name'] = 'watercolumn depth'
    ds['depth_bot'].attrs['units'] = 'm'

    ds['surfT'].attrs['long_name'] = 'surface temperature'
    ds['surfT'].attrs['units'] = 'C'

    ds['botT'].attrs['long_name'] = 'bottom temperature'
    ds['botT'].attrs['units'] = 'C'

    ds['surfS'].attrs['long_name'] = 'surface salinity'
    ds['surfS'].attrs['units'] = 'C'

    ds['botS'].attrs['long_name'] = 'bottom salinity'
    ds['botS'].attrs['units'] = 'C'

    ds['intLdetN'].attrs['long_name'] = 'depth-integrated large detritus'
    ds['intLdetN'].attrs['units'] = 'mmol N / m2'

    ds['intSdetN'].attrs['long_name'] = 'depth-integrated small detritus'
    ds['intSdetN'].attrs['units'] = 'mmol N / m2'

    ds['intphyto'].attrs['long_name'] = 'depth-integrated phytoplankton'
    ds['intphyto'].attrs['units'] = 'mmol N / m2'

    ds['intzoop'].attrs['long_name'] = 'depth-integrated zooplankton'
    ds['intzoop'].attrs['units'] = 'mmol N / m2'

    ds['intNO3'].attrs['long_name'] = 'depth-integrated nitrate'
    ds['intNO3'].attrs['units'] = 'mmol N / m2'

    ds['intNH4'].attrs['long_name'] = 'depth-integrated ammonium'
    ds['intNH4'].attrs['units'] = 'mmol N / m2'

    ds['surfNO3'].attrs['long_name'] = 'surface nitrate'
    ds['surfNO3'].attrs['units'] = 'mmol N / m3'

    ds['surfNH4'].attrs['long_name'] = 'surface ammonium'
    ds['surfNH4'].attrs['units'] = 'mmol N / m3'

    return ds


##############################################################
##                      PROCESS DATA                        ##
##############################################################

print('Processing started...\n')

for year in years:
    print(year)

    # get data
    fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_'+year+'.01.01_'+year+'.12.31.nc')
    ds_raw = xr.open_dataset(fp)

    # set values in strait of juan de fuca and strait of georgia to nan 
    # (so they don't interfere with analysis)
    if remove_straits:
        print('    Removing Straits...')
        lat_threshold = 48.14
        lon_threshold = -122.76
        # Create a mask for latitudes and longitudes in the Straits
        mask = (ds_raw['lat_rho'] > lat_threshold) & (ds_raw['lon_rho'] < lon_threshold)
        # Expand mask dimensions to match 'oxygen' dimensions
        expanded_mask = mask.expand_dims(ocean_time=len(ds_raw['ocean_time']), s_rho=len(ds_raw['s_rho']))
        # Apply the mask to the variable
        ds_raw['temp'] = xr.where(expanded_mask, np.nan, ds_raw['temp'])
        ds_raw['salt'] = xr.where(expanded_mask, np.nan, ds_raw['salt'])
        ds_raw['LdetritusN'] = xr.where(expanded_mask, np.nan, ds_raw['LdetritusN'])
        ds_raw['SdetritusN'] = xr.where(expanded_mask, np.nan, ds_raw['SdetritusN'])
        ds_raw['phytoplankton'] = xr.where(expanded_mask, np.nan, ds_raw['phytoplankton'])
        ds_raw['zooplankton'] = xr.where(expanded_mask, np.nan, ds_raw['zooplankton'])
        ds_raw['NO3'] = xr.where(expanded_mask, np.nan, ds_raw['NO3'])
        ds_raw['NH4'] = xr.where(expanded_mask, np.nan, ds_raw['NH4'])

    # initialize dataset
    ds = start_ds(ds_raw['ocean_time'],
                  ds_raw['eta_rho'],
                  ds_raw['xi_rho'],)
    # add metadata
    ds = add_metadata(ds)

    print('    Calculating depth-integrated bio vars')
    # get thickness of hypoxic layer in watercolumn at ever lat/lon cell (ocean_time: 365, eta_rho: 441, xi_rho: 177)
    # units are in mmol N / m2
    Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
    reader = csv.DictReader(open(Sfp))
    S_dict = {}
    for row in reader:
        S_dict[row['ITEMS']] = row['VALUES']
    S = zrfun.get_S(S_dict)
    # get cell thickness
    h = ds_raw['h'].values # height of water column
    z_rho, z_w = zrfun.get_z(h, 0*h, S) 
    dzr = np.diff(z_w, axis=0) # vertical thickness of all cells [m]  
    # Multiple cell height array by variable (m * mmol N / m3 = mmol N / m2)
    # then sum along z to depth-integrate
    print('         LdetN')
    height_times_LdetN = dzr * ds_raw['LdetritusN'].values
    intLdetN = np.nansum(height_times_LdetN,axis=1)
    print('         SdetN')
    height_times_SdetN = dzr * ds_raw['SdetritusN'].values
    intSdetN = np.nansum(height_times_SdetN,axis=1)
    print('         phytoplankton')
    height_times_phyto = dzr * ds_raw['phytoplankton'].values
    intphyto = np.nansum(height_times_phyto,axis=1)
    print('         zooplankton')
    height_times_zoop = dzr * ds_raw['zooplankton'].values
    intzoop = np.nansum(height_times_zoop,axis=1)
    print('         NO3')
    height_times_NO3 = dzr * ds_raw['NO3'].values
    intNO3 = np.nansum(height_times_NO3,axis=1)
    print('         NH4')
    height_times_NH4 = dzr * ds_raw['NH4'].values
    intNH4 = np.nansum(height_times_NH4,axis=1)

    # add data to ds
    print('    Adding data to dataset')

    # depth of water column
    ds['depth_bot'] = xr.DataArray(ds_raw['h'].values,
                                coords={'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['eta_rho', 'xi_rho'])
    
    # surfT
    ds['surfT'] = xr.DataArray(ds_raw['temp'].values[:,-1,:,:],
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # bottT
    ds['botT'] = xr.DataArray(ds_raw['temp'].values[:,0,:,:],
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # surfS
    ds['surfS'] = xr.DataArray(ds_raw['salt'].values[:,-1,:,:],
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # bottS
    ds['botS'] = xr.DataArray(ds_raw['salt'].values[:,0,:,:],
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # depth-integrated LdetN
    ds['intLdetN'] = xr.DataArray(intLdetN,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # depth-integrated LdetN
    ds['intSdetN'] = xr.DataArray(intSdetN,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # depth-integrated LdetN
    ds['intphyto'] = xr.DataArray(intphyto,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # depth-integrated LdetN
    ds['intzoop'] = xr.DataArray(intzoop,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # depth-integrated NO3
    ds['intNO3'] = xr.DataArray(intNO3,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # depth-integrated NH4
    ds['intNH4'] = xr.DataArray(intNH4,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # surface NO3
    ds['surfT'] = xr.DataArray(ds_raw['NO3'].values[:,-1,:,:],
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    
    # surface NH4
    ds['surfT'] = xr.DataArray(ds_raw['NH4'].values[:,-1,:,:],
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])

    print('    Saving dataset')
    # save dataset
    if remove_straits:
        straits = 'noStraits'
    else:
        straits = 'withStraits'
    ds.to_netcdf(out_dir / (year + '_vars_' + straits + '.nc'))

print('Done')

print(list(ds.keys()))