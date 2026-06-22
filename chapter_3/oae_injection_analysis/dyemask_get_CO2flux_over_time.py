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
# ds1 = '2020.06.30'

ds0 = '2020.07.01'
ds1 = '2020.10.31'

# which  model runs to look at?
gtagex_base = 'cas7_t1_x11ab'
gtagex_pert  = 'cas7_t1dgeWB_x11abd'#'cas7_t1dgeWB_x11abd3monthscont' # 'cas7_t1dgeWB_x11abd' ------------------------------ APOGEE CHANGE!!!

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
        # total DIC in baseline
        total_DIC_base_full    = (['ocean_time'], np.zeros((Ndays))),
        # difference in CO2 flux between perturbation and baseline
        total_CO2_flux_full    = (['ocean_time'], np.zeros((Ndays))),

        
        # difference in total alkalinity between perturbation and baseline
        delta_Alk_full    = (['ocean_time'], np.zeros((Ndays))),
        # total alkalinity in baseline
        total_Alk_base_full    = (['ocean_time'], np.zeros((Ndays))),

        # total amount of dye in domain
        total_dye_full   = (['ocean_time'], np.zeros((Ndays))),
            
        ),
    coords=dict(ocean_time=ocean_time))
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for processed CO2 uptake capacity
    '''
    ds['delta_DIC_full'].attrs['long_name'] = 'Full domain: Difference in DIC between perturbation and baseline'
    ds['delta_DIC_full'].attrs['units'] = 'kmol C'

    ds['total_DIC_base_full'].attrs['long_name'] = 'Full domain: Total DIC in baseline'
    ds['total_DIC_base_full'].attrs['units'] = 'kmol C'

    ds['total_CO2_flux_full'].attrs['long_name'] = 'Full domain: Difference in CO2 flux between perturbation and baseline'
    ds['total_CO2_flux_full'].attrs['units'] = 'kmol C'

    ds['delta_Alk_full'].attrs['long_name'] = 'Full domain:Difference in total alkalinity between perturbation and baseline'
    ds['delta_Alk_full'].attrs['units'] = 'kmol'

    ds['total_Alk_base_full'].attrs['long_name'] = 'Full domain: Total alkalinity in baseline'
    ds['total_Alk_base_full'].attrs['units'] = 'kmol'

    ds['total_dye_full'].attrs['long_name'] = 'Full domain:Total amount of dye in domain'
    ds['total_dye_full'].attrs['units'] = 'kmol (assuming OH-)'


    return ds

##############################################################
##                 Important coefficients                   ##
##############################################################


# Schmidt number transfer coefficients (Wanninkhof, 1992)
A_CO2 = 2073.1
B_CO2 = 125.62
C_CO2 = 3.6276
D_CO2 = 0.043219
E_CO2 = 0.0

# Surface CO2 solubility coefficients (Weiss, 1974)
A1 = -60.2409
A2 = 93.4517
A3 = 23.3585
B1 = 0.023517
B2 = -0.023656
B3 = 0.0047036

# Coefficients to calculate secular trent in atmospehric pCO2 (citation???)
D0 = 282.6
D1 = 0.125
D2 =-7.18
D3 = 0.862
D4 =-0.99
D5 = 0.28
D6 =-0.80
D7 = 0.06

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
    
    # Get cell vertical thicknesses
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
    alk_base_kmol            = alk_base/1000/1000 # [kmol/m3]

    TIC_base = ds_base['TIC'] # [mmol/m3]
    TIC_pert = ds_pert['TIC'] # [mmol/m3]
    TIC_pert_minus_base = TIC_pert - TIC_base # [mmol/m3]
    TIC_pert_minus_base_kmol = TIC_pert_minus_base/1000/1000 # [kmol/m3]
    TIC_base_kmol            = TIC_base/1000/1000 # [kmol/m3]

    # get total dye and surface dye
    dye_all = ds_pert['dye_01']          # [kg/m3] (t,z,y,x)

    # convert units to kmol/m3
    # where the 1.7e-5 converts dye in kg to mmol (assuming dye is a proxy for OH-)
    dye_all_kmol  = dye_all  / 1.7e-5  / 1000 / 1000 # [kmol/m3]

    # apply mask to only include cells where dye concentration is greater than 5e-4 mmol/m3
    # 5e-4 mmol/m3 = 5e-10 kmol/m3
    dye_threshold = 5e-10 # [kmol/m3]
    # binary mask: 1 where dye >= threshold, 0 otherwise
    dye_mask = xr.where(dye_all_kmol >= dye_threshold, 1, 0)


    # Multiply all variables by cell volume to get all final values in units of kmol
    # and apply the dye mask
    delta_DIC_individualcells = TIC_pert_minus_base_kmol * vol_all * dye_mask # [kmol] dims are (t,z,y,x)
    delta_Alk_individualcells = alk_pert_minus_base_kmol * vol_all * dye_mask # [kmol] dims are (t,z,y,x)
    total_dye_individualcells = dye_all_kmol * vol_all * dye_mask   # [kmol] dims are (t,z,y,x)

    total_DIC_individualcells_base = TIC_base_kmol * vol_all * dye_mask # [kmol] dims are (t,z,y,x)
    total_Alk_individualcells_base = alk_base_kmol * vol_all * dye_mask # [kmol] dims are (t,z,y,x)

    # sum up values to get total volume integral in the full domain
    delta_DIC_full = np.nansum(delta_DIC_individualcells[0,:,:,:]) #[kmol]
    delta_Alk_full = np.nansum(delta_Alk_individualcells[0,:,:,:]) #[kmol]
    total_dye_full = np.nansum(total_dye_individualcells[0,:,:,:]) #[kmol]
    total_DIC_base_full = np.nansum(total_DIC_individualcells_base[0,:,:,:]) #[kmol]
    total_Alk_base_full = np.nansum(total_Alk_individualcells_base[0,:,:,:]) #[kmol]

    # Get change in CO2 flux     
    # CALCULATE SCHMIDT NUMBER (same for both runs because hydrodynamics identical)
    tempC = ds_base['temp'].values[:,-1,:,:] # surface temperature in Celsius
    Sc = A_CO2 + B_CO2*tempC + C_CO2*tempC**2 + D_CO2*tempC**3 + E_CO2*tempC**4
    # DETERMINE GAS TRANSFER DEPTH [m] (again, same for both runs)
    windspeed2 = ds_base['Uwind'].values**2 + ds_base['Vwind'].values**2 # wind speed squared
    cff3 = 24/100 * 0.31 * windspeed2 * (Sc/660)**(-1/2) # [m] (over one day)
    # GET SURFACE CO2 SOLUBILITY [moles/(kg atm)] (same for both runs, again)
    tempK0p01 = (ds_base['temp'].values[:,-1,:,:] + 273.15)/100 # temp in K divided by 100
    S = ds_base['salt'].values[:,-1,:,:] # surface salinity
    CO2_sol = np.exp(A1 + A2/tempK0p01 + A3*np.log(tempK0p01) + 
                        S*(B1 + B2*tempK0p01 + B3*tempK0p01**2))
    # GET PCO2 OF ATMOSPHERE [uatm]
    date_obj = datetime.strptime(date_str[1::], "%Y.%m.%d")
    yearday = int(date_obj.strftime("%j")) # yearday from date
    pmonth = int(2020) - 1951 + yearday/365 # months since January 1951
    pi2 = 2 * np.pi
    PCO2air_secular = D0 + D1*pmonth*12 + D2*np.sin(pi2*pmonth+D3) + D4*np.sin(pi2*pmonth+D5) + D6*np.sin(pi2*pmonth+D7)
    # GET PCO2 OF SURFACE OCEAN [uatm]
    P = gsw.p_from_z(z_rho, ds_base['lon_rho'].values, ds_base['lat_rho'].values)
    # in-situ density
    z = z_rho
    lons = ds_base['lon_rho'].values
    lats = ds_base['lat_rho'].values
    PT = ds_base['temp'].values # potential temperature
    SP = ds_base['salt'].values # practical salinity
    SA = gsw.SA_from_SP(SP, P, lons, lats) # absolute salinity
    CT = gsw.CT_from_pt(SA, PT) # conservative temperature
    rho = gsw.rho(SA, CT, P) # in situ density [kg m-3]
    # surface density
    rho_surf = rho[:,-1,:,:] # [kg m-3]
    # Up until here, everything should be the same between runs!
    # Now we get out of hydrodynamics and into carbonate chem, and
    # need to be careful about the two different model conditions
    # Convert from micromol/L to micromol/kg using in situ dentity because these are the
    # units expected by pyco2.
    alk_surf_umolkg_base = 1000 * ds_base['alkalinity'].values[:,-1,:,:] / rho_surf
    TIC_surf_umolkg_base = 1000 * ds_base['TIC'].values[:,-1,:,:] / rho_surf
    alk_surf_umolkg_pert = 1000 * ds_pert['alkalinity'].values[:,-1,:,:] / rho_surf
    TIC_surf_umolkg_pert = 1000 * ds_pert['TIC'].values[:,-1,:,:] / rho_surf
    # get surface values and apply mask to get only water values (so script runs faster)
    mask_rho = ds_base['mask_rho'].values 
    aPres =  gsw.p_from_z(z_rho[-1,:,:], ds_base['lon_rho'].values,
                ds_base['lat_rho'].values).squeeze()[mask_rho==1]  # [dbar] # pressure is the same in both runs
    # alkalnity and TIC are different between runs, so get them separately
    aALK_base = alk_surf_umolkg_base.squeeze()[mask_rho==1]                 # [umol/kg]
    aTIC_base = TIC_surf_umolkg_base.squeeze()[mask_rho==1]                 # [umol/kg]
    aALK_pert = alk_surf_umolkg_pert.squeeze()[mask_rho==1]                 # [umol/kg]
    aTIC_pert = TIC_surf_umolkg_pert.squeeze()[mask_rho==1]                 # [umol/kg]
    # temp and salt are bit-reproducible
    aTemp = ds_base['temp'].values[:,-1,:,:].squeeze()[mask_rho==1]# [C]
    aSalt = ds_base['salt'].values[:,-1,:,:].squeeze()[mask_rho==1]# [psu]
    # Then calculate pCO2 in the two different runs
    # https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/
    CO2dict_base = pyco2.sys(par1=aALK_base, par1_type=1, par2=aTIC_base, par2_type=2,
            salinity=aSalt, temperature=aTemp, pressure=aPres,
            total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
    aPCO2_base = CO2dict_base['pCO2']
    CO2dict_pert = pyco2.sys(par1=aALK_pert, par1_type=1, par2=aTIC_pert, par2_type=2,
            salinity=aSalt, temperature=aTemp, pressure=aPres,
            total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
    aPCO2_pert = CO2dict_pert['pCO2']
    # Create a 2D array of NaNs with the same shape as mask_rho
    surf_pCO2_base = np.full(mask_rho.shape, np.nan)
    surf_pCO2_pert = np.full(mask_rho.shape, np.nan)
    # Map the 1D pCO2 results back into the 2D array, indexing with [mask_rho==1]
    surf_pCO2_base[mask_rho==1] = aPCO2_base
    surf_pCO2_pert[mask_rho==1] = aPCO2_pert
    # CALCULATE IDEAL CO2 CAPACITY [MMOL/M3]
    CO2_capacity_ideal_base = CO2_sol * (PCO2air_secular - surf_pCO2_base) # [mmol/m3]
    CO2_capacity_ideal_pert = CO2_sol * (PCO2air_secular - surf_pCO2_pert) # [mmol/m3]
    # GET ACTUAL CO2 TRANSFER CAPACITY BY INCLUDING WIND SPEED [MMOL/M2]
    CO2_flux_actual_base = CO2_capacity_ideal_base * cff3# [mmol/m3] * [m] = [mmol/m2/day] 
    CO2_flux_actual_pert = CO2_capacity_ideal_pert * cff3# [mmol/m3] * [m] = [mmol/m2/day] 
    # flatten first dimension
    CO2_flux_actual_base = CO2_flux_actual_base[0,:,:] # [mmol/m2/day]
    CO2_flux_actual_pert = CO2_flux_actual_pert[0,:,:] # [mmol/m2/day]
    # Get difference in flux
    delta_CO2_flux = CO2_flux_actual_pert - CO2_flux_actual_base # [mmol/m2/day]
    # Integrate over the surface to get total CO2 flux in moles per day
    delta_CO2_flux_overarea = delta_CO2_flux* DA # [mmol/m2/day] * [m2] = [mmol/day]
    # convert to kmol
    delta_CO2_flux_overarea_kmol = delta_CO2_flux_overarea / 1000 / 1000 # [kmol/day]: kmol/day * 1 day
    # crop using dye mask (surface only)
    # binary mask: 1 where dye >= threshold, 0 otherwise
    dye_mask_2D = xr.where(dye_all_kmol[:,-1,:,:] >= dye_threshold, 1, 0)
    # apply mask
    delta_CO2_flux_overarea_kmol_masked = delta_CO2_flux_overarea_kmol * dye_mask_2D
    # get toal flux
    total_CO2_flux_full = np.nansum(delta_CO2_flux_overarea_kmol_masked) # [kmol]

    # add data to daily dataset
    daily_ds['delta_DIC_full'] = xr.DataArray(delta_DIC_full,
                                coords={'ocean_time': ds_base['ocean_time'].values},
                                dims=['ocean_time'])
    daily_ds['total_DIC_base_full'] = xr.DataArray(total_DIC_base_full,
                                coords={'ocean_time': ds_base['ocean_time'].values},
                                dims=['ocean_time'])
    daily_ds['total_CO2_flux_full'] = xr.DataArray(total_CO2_flux_full,
                                coords={'ocean_time': ds_base['ocean_time'].values},
                                dims=['ocean_time'])

    daily_ds['delta_Alk_full'] = xr.DataArray(delta_Alk_full,
                                coords={'ocean_time': ds_base['ocean_time'].values},
                                dims=['ocean_time'])
    daily_ds['total_Alk_base_full'] = xr.DataArray(total_Alk_base_full,
                                coords={'ocean_time': ds_base['ocean_time'].values},
                                dims=['ocean_time'])
    
    daily_ds['total_dye_full'] = xr.DataArray(total_dye_full,
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
ds_out.to_netcdf(out_dir / ('dyemask_CO2flux_'+ds0+'_'+ds1+'.nc'), encoding=encoding)

print(ds_out)

print('Done')