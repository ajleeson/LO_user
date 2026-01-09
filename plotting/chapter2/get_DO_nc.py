"""
This script determines the s-level, depth, and concentration of the 
DO minima in the watercolumn, and saves the data in a new .nc file.

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

region = 'pugetsoundDO'
# region = 'SS_and_HC_low'

remove_straits = True

years = ['2015']

# which  model run to look at?
gtagex = 'cas7_t1_x11ab'  
# gtagex = 'cas7_t1noDIN_x11ab'  

# where to put output files
out_dir = Ldir['LOo'] / 'chapter_2' / 'data'
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
        # depth of DO minima
        depth_min   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # slevel of DO minima
        slev_min    = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # concentration of DO minima
        DO_min      = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # depth of water column
        depth_bot   = (['eta_rho','xi_rho'], np.zeros((Neta,Nxi))),
        # DO concentration at bottom
        DO_bot      = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # thickness of hypoxic layer
        hyp_thick   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),),
    coords=dict(ocean_time=ocean_time, eta_rho=eta_rho, xi_rho=xi_rho,),)
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for processed DO data
    '''

    ds['depth_min'].attrs['long_name'] = 'depth of watercolumn DO minima'
    ds['depth_min'].attrs['units'] = 'm'

    ds['slev_min'].attrs['long_name'] = 's-level of watercolumn DO minima'
    ds['slev_min'].attrs['units'] = 'unitless'

    ds['DO_min'].attrs['long_name'] = 'concentration of watercolumn DO minima'
    ds['DO_min'].attrs['units'] = 'mg/L'

    ds['depth_bot'].attrs['long_name'] = 'watercolumn depth'
    ds['depth_bot'].attrs['units'] = 'm'

    ds['DO_bot'].attrs['long_name'] = 'DO concentration at bottom'
    ds['DO_bot'].attrs['units'] = 'mg/L'

    ds['hyp_thick'].attrs['long_name'] = 'thickness of hypoxic layer'
    ds['hyp_thick'].attrs['units'] = 'm'

    return ds


##############################################################
##                      PROCESS DATA                        ##
##############################################################

print('Processing started...\n')

for year in years:
    print(year)

    # get data
    fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / (region+'_'+year+'.01.01_'+year+'.12.31.nc')
    ds_raw = xr.open_dataset(fp)

    # set values in strait of juan de fuca and strait of georgia to nan 
    # (so they don't interfere with analysis)
    # only need to do this over the whole Puget Sound region
    if region == 'pugetsoundDO':
        if remove_straits:
            print('    Removing Straits...')
            lat_threshold = 48.14
            lon_threshold = -122.76
            # Create a mask for latitudes and longitudes in the Straits
            mask = (ds_raw['lat_rho'] > lat_threshold) & (ds_raw['lon_rho'] < lon_threshold)
            # Expand mask dimensions to match 'oxygen' dimensions
            expanded_mask = mask.expand_dims(ocean_time=len(ds_raw['ocean_time']), s_rho=len(ds_raw['s_rho']))
            # Apply the mask to the 'oxygen' variable
            ds_raw['oxygen'] = xr.where(expanded_mask, np.nan, ds_raw['oxygen'])

    # initialize dataset
    ds = start_ds(ds_raw['ocean_time'],
                  ds_raw['eta_rho'],
                  ds_raw['xi_rho'],)
    # add metadata
    ds = add_metadata(ds)

    print('    Calculating hypoxic thickness')
    # get thickness of hypoxic layer in watercolumn at ever lat/lon cell (ocean_time: 365, eta_rho: 441, xi_rho: 177)
    # units are in m (thickness of hypoxic layer)
    # get S for the whole grid
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
    # Now get oxygen values at every grid cell and convert to mg/L
    oxy_mgL = pinfo.fac_dict['oxygen'] * ds_raw['oxygen'].values
    # remove all non-hypoxic values (greater than 2 mg/L)
    hypoxic = np.where(oxy_mgL <= 2, 1, np.nan) # array of nans and ones. one means hypoxic, nan means nonhypoxic
    # Multiple cell height array by hypoxic array boolean array
    hyp_cell_thick = dzr * hypoxic
    # Sum along z to get thickness of hypoxic layer
    hyp_thick = np.nansum(hyp_cell_thick,axis=1)

    print('    Calculating depth of DO minima')
    # get s-rho of the lowest DO (array with dimensions of (ocean_time: 365, eta_rho: eta_size, xi_rho: xi_size))
    eta_size = ds_raw['h'].sizes['eta_rho']
    xi_size  = ds_raw['h'].sizes['xi_rho']
    srho_min = ds_raw['oxygen'].idxmin(dim='s_rho', skipna=True)#.values
    # get depths, but also flatten the time dimension
    depths = ds_raw['h'].values
    # reshape
    depths_reshape = depths.reshape((1,eta_size,xi_size))
    # convert srho to depths
    depth_min = depths_reshape * srho_min

    print('    Calculating s-level of DO minima')
    # get s-level of the lowest DO (array with dimensions of (ocean_time: 365, eta_rho: 441, xi_rho: 177))
    # add new dimension with s-levels
    s_level = np.linspace(0,29,30)
    # natural
    ds_raw['oxygen'] = ds_raw['oxygen'].assign_coords({'s_level': ('s_rho',s_level)})
    ds_raw['oxygen'] = ds_raw['oxygen'].swap_dims({'s_rho': 's_level'})
    # calculate slevel corresponding to DO min
    slev_min = ds_raw['oxygen'].idxmin(dim='s_level', skipna=True).values

    print('    Calculating concentration of DO minima')
    # get corresponding DO minima concentration (mg/L)
    DO_min = pinfo.fac_dict['oxygen'] * ds_raw['oxygen'].min(dim='s_level', skipna=True).values

    # get bottom DO concentration
    DO_bot = pinfo.fac_dict['oxygen'] * ds_raw['oxygen'][:,0,:,:].values

    # add data to ds
    print('    Adding data to dataset')
    # depth of DO minima
    ds['depth_min'] = xr.DataArray(depth_min,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    # slevel
    ds['slev_min'] = xr.DataArray(slev_min,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    # concentration of DO minima
    ds['DO_min'] = xr.DataArray(DO_min,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    # depth of water column
    ds['depth_bot'] = xr.DataArray(depths_reshape.reshape((eta_size,xi_size)),
                                coords={'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['eta_rho', 'xi_rho'])
    # DO concentration at bottom
    ds['DO_bot'] = xr.DataArray(DO_bot,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])
    # hypoxic layer thickness
    ds['hyp_thick'] = xr.DataArray(hyp_thick,
                                coords={'ocean_time': ds_raw['ocean_time'].values,
                                        'eta_rho': ds_raw['eta_rho'].values,
                                        'xi_rho': ds_raw['xi_rho'].values},
                                dims=['ocean_time','eta_rho', 'xi_rho'])

    print('    Saving dataset')
    # save dataset
    if region == 'pugetsoundDO':
        if remove_straits:
            straits = 'noStraits'
        else:
            straits = 'withStraits'
        ds.to_netcdf(out_dir / (gtagex + '_' + region + '_' + year + '_DO_info_' + straits + '.nc'))
    else:
        ds.to_netcdf(out_dir / (gtagex + '_' + region + '_' + year + '_DO_info.nc'))

print('Done')