"""
To be run on apogee
to get time series of delta DIC, delta alkalinity, 
surface dye, and total dye over the whole domain
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

# ds0 = '2020.06.01'
# ds0 = '2020.08.31' # testing ------------------------------ APOGEE CHANGE!!!
# ds1 = '2020.08.31'
ds0 = '2020.07.01'
ds1 = '2020.08.30'

# which  model runs to look at?
gtagex_base = 'cas7_t1_x11ab'
gtagex_pert  = 'cas7_t1dgeWB_x11abd' #'cas7_t1dgeWB_x11abd3monthscont' # 'cas7_t1dgeWB_x11abd' ------------------------------ APOGEE CHANGE!!!

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
        delta_DIC    = (['ocean_time'], np.zeros((Ndays))),
        # difference in total alkalinity between perturbation and baseline
        delta_Alk    = (['ocean_time'], np.zeros((Ndays))),
        # total amount of dye in domain
        total_dye   = (['ocean_time'], np.zeros((Ndays))),
        # amount of dye in surface layer only
        surf_dye    = (['ocean_time'], np.zeros((Ndays))),
            
        ),
    coords=dict(ocean_time=ocean_time))
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for processed CO2 uptake capacity
    '''
    ds['delta_DIC'].attrs['long_name'] = 'Difference in DIC between perturbation and baseline'
    ds['delta_DIC'].attrs['units'] = 'kmol C'

    ds['delta_Alk'].attrs['long_name'] = 'Difference in total alkalinity between perturbation and baseline'
    ds['delta_Alk'].attrs['units'] = 'kmol'

    ds['total_dye'].attrs['long_name'] = 'Total amount of dye in domain'
    ds['total_dye'].attrs['units'] = 'kmol (assuming OH-)'

    ds['surf_dye'].attrs['long_name'] = 'Total amount of dye in surface layer only'
    ds['surf_dye'].attrs['units'] = 'kmol (assuming OH-)'

    return ds

##############################################################
##                      PROCESS DATA                        ##
##############################################################

print('Processing started...\n')

# BASELINE: get info to find history files
gridname_base, tag_base, ex_base = gtagex_base.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname_base, tag=tag_base, ex_name=ex_base)
# add more entries to Ldir
Ldir['roms_out'] =  Ldir['roms_out5'] # Ldir['roms_out']------------------------------ APOGEE CHANGE!!!
Ldir['ds0'] = ds0
Ldir['ds1'] = ds1
Ldir['list_type'] = 'average'
# print(Ldir.keys())
# get history files
fn_list_base = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])

# PERTURBATION: get info to find history files
gridname_pert, tag_pert, ex_pert = gtagex_pert.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname_pert, tag=tag_pert, ex_name=ex_pert)
# add more entries to Ldir
Ldir['roms_out'] = Ldir['roms_out']
Ldir['ds0'] = ds0
Ldir['ds1'] = ds1
Ldir['list_type'] = 'average'
# print(Ldir.keys())
# get history files
fn_list_pert = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])

# Initialize empty list of datasets (meant for daily datasets)
ds_list = []

# loop through history files for the year
for i,fn_base in enumerate(fn_list_base):

    # get the corresponding perturbation file
    fn_pert = fn_list_pert[i]

    # get data
    date_str = Path(fn_base).parent.name
    print('    ' + date_str) 
    # get data
    ds_base = xr.open_dataset(fn_base)
    ds_pert = xr.open_dataset(fn_pert)

    # initalize dataset to store daily data for this day
    daily_ds = start_ds(ds_base['ocean_time'])

    # if this is the first time step  get grid cell areas (time-invariant)
    if i == 0:
        DX = (ds_base.pm.values)**-1
        DY = (ds_base.pn.values)**-1
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
    h = ds_base['h'].values # height of water column
    zeta = ds_base['zeta'].values
    z_rho, z_w = zrfun.get_z(h, zeta, S)
    dzr = np.diff(z_w, axis=0) # sigma layer thickness [m]

    # Calculate volume of each cell in the whole grid
    vol_all = dzr * DA # [m3]
    # Calculate volume of just the surface layer
    vol_surf = dzr[-1,:,:] * DA # [m3]

    # Get delta alkalinity and delta DIC in each cell
    alk_base = ds_base['alkalinity'] # [meq/m3]
    alk_pert = ds_pert['alkalinity'] # [meq/m3]
    alk_pert_minus_base = alk_pert - alk_base # [meq/m3]
    alk_pert_minus_base_kmol = alk_pert_minus_base/1000/1000 # [kmol/m3]

    TIC_base = ds_base['TIC'] # [mmol/m3]
    TIC_pert = ds_pert['TIC'] # [mmol/m3]
    TIC_pert_minus_base = TIC_pert - TIC_base # [mmol/m3]
    TIC_pert_minus_base_kmol = TIC_pert_minus_base/1000/1000 # [kmol/m3]

    # get total dye and surface dye
    dye_all = ds_pert['dye_01']          # [kg/m3] (t,z,y,x)
    dye_surf = ds_pert['dye_01'][0,-1,:,:] # [kg/m3] (y,x)
    # convert units to kmol/m3
    # where the 1.7e-5 converts dye in kg to mmol (assuming dye is a proxy for OH-)
    dye_all_kmol  = dye_all  / 1.7e-5  / 1000 / 1000 # [kmol/m3]
    dye_surf_kmol = dye_surf / 1.7e-5  / 1000 / 1000 # [kmol/m3]

    # Multiply all variables by cell volume to get all final values in units of kmol
    delta_DIC_individualcells = TIC_pert_minus_base_kmol * vol_all # [kmol] dims are (t,z,y,x)
    delta_Alk_individualcells = alk_pert_minus_base_kmol * vol_all # [kmol] dims are (t,z,y,x)
    total_dye_individualcells = dye_all_kmol * vol_all   # [kmol] dims are (t,z,y,x)
    surf_dye_individualcells  = dye_surf_kmol * vol_surf # [kmol] dims are (y,x)

    # Crop to sub-domain of model to eliminate boundary noise!!!
    # lon/lat limits (Study Domain)
    xmin = -126
    xmax = -122
    ymin = 45.5
    ymax = 50.5
    # get eta and xi indices corresponding to these limits
    lon = ds_base.lon_rho.values
    lat = ds_base.lat_rho.values
    lons = lon[0,:]
    lats = lat[:,0]
    eta_min = min(range(len(lats)), key=lambda i: abs(lats[i]-ymin))
    eta_max = min(range(len(lats)), key=lambda i: abs(lats[i]-ymax))
    xi_min  = min(range(len(lons)), key=lambda i: abs(lons[i]-xmin))
    xi_max  = min(range(len(lons)), key=lambda i: abs(lons[i]-xmax))
    # crop to sub-domain
    delta_DIC_individualcells_cropped = delta_DIC_individualcells[0,:, eta_min:eta_max, xi_min:xi_max]
    delta_Alk_individualcells_cropped = delta_Alk_individualcells[0,:, eta_min:eta_max, xi_min:xi_max]
    total_dye_individualcells_cropped = total_dye_individualcells[0,:, eta_min:eta_max, xi_min:xi_max]
    surf_dye_individualcells_cropped  = surf_dye_individualcells[eta_min:eta_max, xi_min:xi_max]

    # sum up values to get total volume integral in the sub-domain
    delta_DIC = np.nansum(delta_DIC_individualcells_cropped) #[kmol]
    delta_Alk = np.nansum(delta_Alk_individualcells_cropped) #[kmol]
    total_dye = np.nansum(total_dye_individualcells_cropped) #[kmol]
    surf_dye  = np.nansum(surf_dye_individualcells_cropped)  #[kmol]

    # add data to daily dataset
    daily_ds['delta_DIC'] = xr.DataArray(delta_DIC,
                                coords={'ocean_time': ds_base['ocean_time'].values},
                                dims=['ocean_time'])
    daily_ds['delta_Alk'] = xr.DataArray(delta_Alk,
                                coords={'ocean_time': ds_base['ocean_time'].values},
                                dims=['ocean_time'])
    daily_ds['total_dye'] = xr.DataArray(total_dye,
                                coords={'ocean_time': ds_base['ocean_time'].values},
                                dims=['ocean_time'])
    daily_ds['surf_dye'] = xr.DataArray(surf_dye,
                                coords={'ocean_time': ds_base['ocean_time'].values},
                                dims=['ocean_time'])

    # Cast to float32 HERE to keep RAM usage extremely low while building the list
    daily_ds = daily_ds.astype(np.float32)

    # Append to list
    ds_list.append(daily_ds)
    
    # Close the raw dataset to free up memory and prevent file-handle limits
    ds_base.close()
    ds_pert.close()

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
ds_out.to_netcdf(out_dir / ('onemonimpulse_oae_deltas_SUBDOMAIN_'+ds0+'_'+ds1+'.nc'), encoding=encoding)

print(ds_out)

print('Done')