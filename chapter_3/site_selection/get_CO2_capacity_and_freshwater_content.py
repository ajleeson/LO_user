"""
Get CO2 uptake capacity from model history file
and freshwater content
and other useful values for understanding the drivers of CO2 uptake capacity

.nc files are saved in LO_output/chapter_3/data
"""

# import things
import numpy as np
import xarray as xr
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

plt.close('all')

##############################################################
##                       USER INPUTS                        ##
##############################################################

years = ['2021','2022','2023','2024','2025'] #['2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025']

# which  model run to look at?
gtagex = 'cas7_t1_x11ab'

# where to put output files
out_dir = Ldir['LOo'] / 'chapter_3' / 'data'
Lfun.make_dir(out_dir)

##############################################################
##                    HELPER FUNCTIONS                      ##
##############################################################

def start_ds(ocean_time,eta_rho,xi_rho):
    '''
    Initialize dataset to store processed CO2 uptake capcity
    ocean_time = ocean time vector
    eta_rho = eta rho vector
    xi_rho = xi_rho vector
    '''
    Ndays = len(ocean_time.values)
    Neta = len(eta_rho.values)
    Nxi = len(xi_rho.values)

    ds = xr.Dataset(data_vars=dict(
        # Capacity for CO2 uptake (just based on temperature and concentration difference)
        CO2_capacity_ideal    = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # Actual CO2 uptake capacity (which includes wind speed)
        CO2_flux_actual   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        
        # Surface values for sensitivty testing later
        surf_temp   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        surf_salt   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        surf_alk    = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        surf_TIC    = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        wind2       = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        delta_pCO2  = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),

        # Freshwater content [m]
        Fs_31        = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        Fs_30        = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        Fs_29      = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
            
        ),
    coords=dict(ocean_time=ocean_time, eta_rho=eta_rho, xi_rho=xi_rho,),)
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for processed CO2 uptake capacity
    '''
    ds['CO2_capacity_ideal'].attrs['long_name'] = 'Ideal capacity of CO2 uptake based on just solubility and Delta pCO2-- given as a concentration'
    ds['CO2_capacity_ideal'].attrs['units'] = 'mmol CO2 m-3'

    ds['CO2_flux_actual'].attrs['long_name'] = 'Actual flux of CO2 over one day, taking into account winds'
    ds['CO2_flux_actual'].attrs['units'] = 'mmol CO2 m-2 day-1'

    ds['surf_temp'].attrs['long_name'] = 'Surface potential temperature'
    ds['surf_temp'].attrs['units'] = 'deg C'

    ds['surf_salt'].attrs['long_name'] = 'Surface practical salinity'
    ds['surf_salt'].attrs['units'] = 'psu'

    ds['surf_alk'].attrs['long_name'] = 'Surface alkalinity'
    ds['surf_alk'].attrs['units'] = 'meq/m3'

    ds['surf_TIC'].attrs['long_name'] = 'Surface TIC'
    ds['surf_TIC'].attrs['units'] = 'mmol/m3'

    ds['wind2'].attrs['long_name'] = 'wind speed squared'
    ds['wind2'].attrs['units'] = 'm2/s2'

    ds['delta_pCO2'].attrs['long_name'] = 'difference in pCO2 between atmosphere and surface ocean'
    ds['delta_pCO2'].attrs['units'] = 'uatm'

    ds['Fs_31'].attrs['long_name'] = 'freshwater content (s0 = 31)'
    ds['Fs_31'].attrs['units'] = 'm'

    ds['Fs_30'].attrs['long_name'] = 'freshwater content (s0 = 30)'
    ds['Fs_30'].attrs['units'] = 'm'

    ds['Fs_29'].attrs['long_name'] = 'freshwater content (s0 = 29)'
    ds['Fs_29'].attrs['units'] = 'm'

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

for year in years:

    print('{}'.format(year))

    # Initialize empty list of datasets (meant for daily datasets)
    ds_list = []

    # get info to find history files
    gridname, tag, ex_name = gtagex.split('_')
    # get the dict Ldir
    Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
    # add more entries to Ldir
    Ldir['roms_out'] = Ldir['roms_out5']
    Ldir['ds0'] = year + '.01.01'
    Ldir['ds1'] = year + '.12.31'
    Ldir['list_type'] = 'average'
    # print(Ldir.keys())
    # get history files
    fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])

    # loop through history files for the year
    for i,fn in enumerate(fn_list):

        date_str = Path(fn).parent.name
        print('    ' + date_str) 
        # get data
        ds_raw = xr.open_dataset(fn)
        # initalize dataset to store daily data
        daily_ds = start_ds(ds_raw['ocean_time'],
                        ds_raw['eta_rho'],
                        ds_raw['xi_rho'])
        
        # GET Z_RHO
        # get S for the whole grid
        Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
        reader = csv.DictReader(open(Sfp))
        S_dict = {}
        for row in reader:
            S_dict[row['ITEMS']] = row['VALUES']
        S = zrfun.get_S(S_dict)
        # get cell thickness
        h = ds_raw['h'].values # height of water column
        zeta = ds_raw['zeta'].values
        z_rho, z_w = zrfun.get_z(h, zeta[0, :, :], S)
        
        # CALCULATE SCHMIDT NUMBER
        tempC = ds_raw['temp'].values[:,-1,:,:] # surface temperature in Celsius
        Sc = A_CO2 + B_CO2*tempC + C_CO2*tempC**2 + D_CO2*tempC**3 + E_CO2*tempC**4

        # DETERMINE GAS TRANSFER DEPTH [m]
        windspeed2 = ds_raw['Uwind'].values**2 + ds_raw['Vwind'].values**2 # wind speed squared
        cff3 = 24/100 * 0.31 * windspeed2 * (Sc/660)**(-1/2) # [m] (over one day)

        # GET SURFACE CO2 SOLUBILITY [moles/(kg atm)]
        tempK0p01 = (ds_raw['temp'].values[:,-1,:,:] + 273.15)/100 # temp in K divided by 100
        S = ds_raw['salt'].values[:,-1,:,:] # surface salinity
        CO2_sol = np.exp(A1 + A2/tempK0p01 + A3*np.log(tempK0p01) + 
                         S*(B1 + B2*tempK0p01 + B3*tempK0p01**2))
        
        # GET PCO2 OF ATMOSPHERE [uatm]
        date_obj = datetime.strptime(date_str[1::], "%Y.%m.%d")
        yearday = int(date_obj.strftime("%j")) # yearday from date
        pmonth = int(year) - 1951 + yearday/365 # months since January 1951
        pi2 = 2 * np.pi
        PCO2air_secular = D0 + D1*pmonth*12 + D2*np.sin(pi2*pmonth+D3) + D4*np.sin(pi2*pmonth+D5) + D6*np.sin(pi2*pmonth+D7)

        # GET PCO2 OF SURFACE OCEAN [uatm]
        P = gsw.p_from_z(z_rho, ds_raw['lon_rho'].values, ds_raw['lat_rho'].values)
        # in-situ density
        z = z_rho
        lons = ds_raw['lon_rho'].values
        lats = ds_raw['lat_rho'].values
        PT = ds_raw['temp'].values # potential temperature
        SP = ds_raw['salt'].values # practical salinity
        SA = gsw.SA_from_SP(SP, P, lons, lats) # absolute salinity
        CT = gsw.CT_from_pt(SA, PT) # conservative temperature
        rho = gsw.rho(SA, CT, P) # in situ density [kg m-3]
        # surface density
        rho_surf = rho[:,-1,:,:] # [kg m-3]
        # Convert from micromol/L to micromol/kg using in situ dentity because these are the
        # units expected by pyco2.
        alk_surf_umolkg = 1000 * ds_raw['alkalinity'].values[:,-1,:,:] / rho_surf
        TIC_surf_umolkg = 1000 * ds_raw['TIC'].values[:,-1,:,:] / rho_surf
        # get surface values and apply mask to get only water values (so script runs faster)
        mask_rho = ds_raw['mask_rho'].values 
        aPres =  gsw.p_from_z(z_rho[-1,:,:], ds_raw['lon_rho'].values,
                    ds_raw['lat_rho'].values).squeeze()[mask_rho==1]  # [dbar]
        aALK = alk_surf_umolkg.squeeze()[mask_rho==1]                 # [umol/kg]
        aTIC = TIC_surf_umolkg.squeeze()[mask_rho==1]                 # [umol/kg]
        aTemp = ds_raw['temp'].values[:,-1,:,:].squeeze()[mask_rho==1]# [C]
        aSalt = ds_raw['salt'].values[:,-1,:,:].squeeze()[mask_rho==1]# [psu]
        # Note: here is where to get info on the inputs, outputs, and units:
        # https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/
        CO2dict = pyco2.sys(par1=aALK, par1_type=1, par2=aTIC, par2_type=2,
                salinity=aSalt, temperature=aTemp, pressure=aPres,
                total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
        aPCO2 = CO2dict['pCO2']
        # Create a 2D array of NaNs with the same shape as mask_rho
        surf_pCO2 = np.full(mask_rho.shape, np.nan)
        # Map the 1D pCO2 results back into the 2D array, indexing with [mask_rho==1]
        surf_pCO2[mask_rho==1] = aPCO2

        # CALCULATE IDEAL CO2 CAPACITY [MMOL/M3]
        CO2_capacity_ideal = CO2_sol * (PCO2air_secular - surf_pCO2) # [mmol/m3]

        # GET ACTUAL CO2 TRANSFER CAPACITY BY INCLUDING WIND SPEED [MMOL/M2]
        CO2_flux_actual = CO2_capacity_ideal * cff3# [mmol/m3] * [m] = [mmol/m2] 
        
        # -----------------------------------------------------------------------
        # GET THE CELL THICKNESS FOR VERTICAL INTEGRALS
        # get S for the whole grid
        Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
        reader = csv.DictReader(open(Sfp))
        S_dict = {}
        for row in reader:
            S_dict[row['ITEMS']] = row['VALUES']
        S = zrfun.get_S(S_dict)
        # get cell thickness
        h = ds_raw['h'].values # height of water column
        # loop over time to get z_rho and z_w:
        zeta = ds_raw['zeta'].values
        Nt = zeta.shape[0]
        Nz = S['N']
        dzr_all = np.empty((Nt, Nz, *h.shape))
        # loop over time -- but in this case, there is only one time step, so this loop will only run once
        for t in range(Nt):
            # make sure to use zeta as an input to account for SSH variability!!!
            z_rho, z_w = zrfun.get_z(h, zeta[t, :, :], S)
            dzr_all[t, :, :, :] = np.diff(z_w, axis=0) # sigma layer thickness at one time

        # CALCULATE SO - S(x,z) / S0 (use values of s0 = 31, 30, and 29)
        deltaS_31   = (31   - SP) / 31
        deltaS_30   = (30   - SP) / 30
        deltaS_29   = (29 - SP) / 29
        # set values < 0 to 0, so that we only sum where s < s0
        deltaS_31[deltaS_31 < 0]     = 0
        deltaS_30[deltaS_30 < 0]     = 0
        deltaS_29[deltaS_29 < 0]     = 0

        # vertically integrate to get freshwater content
        freshwater_content_29   = np.nansum(deltaS_29   * dzr_all,axis=1)
        freshwater_content_31   = np.nansum(deltaS_31   * dzr_all,axis=1)
        freshwater_content_30   = np.nansum(deltaS_30   * dzr_all,axis=1)
        # apply land mask
        land_mask_2d = ds_raw['mask_rho'].values == 0
        land_mask_3d = land_mask_2d[np.newaxis, :, :]
        freshwater_content_31   = np.ma.masked_where(land_mask_3d, freshwater_content_31)
        freshwater_content_30   = np.ma.masked_where(land_mask_3d, freshwater_content_30)
        freshwater_content_29   = np.ma.masked_where(land_mask_3d, freshwater_content_29)

        # ------------------------------------------------------------


        # add data to daily dataset
        daily_ds['CO2_capacity_ideal'] = xr.DataArray(CO2_capacity_ideal,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['CO2_flux_actual'] = xr.DataArray(CO2_flux_actual,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['surf_temp'] = xr.DataArray(ds_raw['temp'].values[:,-1,:,:],
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['surf_salt'] = xr.DataArray(ds_raw['salt'].values[:,-1,:,:],
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['surf_alk'] = xr.DataArray(ds_raw['alkalinity'].values[:,-1,:,:],
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['surf_TIC'] = xr.DataArray(ds_raw['TIC'].values[:,-1,:,:],
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['wind2'] = xr.DataArray(windspeed2,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['delta_pCO2'] = xr.DataArray(np.expand_dims(PCO2air_secular - surf_pCO2, axis=0),
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        
        daily_ds['Fs_31'] = xr.DataArray(freshwater_content_31,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['Fs_30'] = xr.DataArray(freshwater_content_30,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['Fs_29'] = xr.DataArray(freshwater_content_29,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])

        # Cast to float32 HERE to keep RAM usage extremely low while building the list
        daily_ds = daily_ds.astype(np.float32)

        # Append to list
        ds_list.append(daily_ds)
        
        # Close the raw dataset to free up memory and prevent file-handle limits
        ds_raw.close()

    # PREPARE DATA FOR SAVING
    # Combine all of the daily datasets into one dataset
    ds = xr.concat(ds_list, dim='ocean_time')
    # Add your metadata once at the end
    ds = add_metadata(ds)

    print('    Saving dataset')
    # Create a compression dictionary
    comp = dict(zlib=True, complevel=4)
    encoding = {var: comp for var in ds.data_vars}
    # Save with encoding
    ds.to_netcdf(out_dir / (gtagex + '_' + year + '_freshwatercontent_CO2uptake.nc'), encoding=encoding)


    # print(ds)

print('Done')


# ##############################################
# # TEST PLOTTING

# # get pcolormesh values
# CO2_capacity_ideal = ds['CO2_capacity_ideal'][-1,:,:].values
# CO2_flux_actual = ds['CO2_flux_actual'][-1,:,:].values
# surf_salt = ds['surf_salt'][-1,:,:].values
# Fs_30 = ds['Fs_30'][-1,:,:].values

# # make list of vars
# vars = [CO2_capacity_ideal,CO2_flux_actual, surf_salt,Fs_30]

# # list of vmins and vmax
# vmins = [-15,-50,25,0]
# vmaxs = [15,50,32,5]

# # get colormaps
# cmaps = [cmc.vik,cmc.vik,cmc.imola,cmc.batlowW_r]

# # titles
# titles = ['CO2 uptake capacity based on\n'+'solubility and $\Delta$ pCO2 [mmol CO2 m$^{-3}$]',
#           'CO2 flux over one day\n'+r'accounting for wind [mmol CO2 m$^{-2}$ day$^{-1}$]',
#           'Surface practical salinity',
#           'Freshwater content\n'+r'(s$_0$ = 30) [m]']

# # Get grid data
# G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
# grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
# lon = grid_ds.lon_rho.values
# lat = grid_ds.lat_rho.values
# lon_u = grid_ds.lon_u.values
# lat_u = grid_ds.lat_u.values
# lon_v = grid_ds.lon_v.values
# lat_v = grid_ds.lat_v.values
# px, py = pfun.get_plon_plat(lon,lat)

# # lon/lat limits (Study Domain)
# xmin = -126
# xmax = -122
# ymin = 45.5
# ymax = 50.5

# for var, vmin, vmax, cmap, title in zip(vars, vmins, vmaxs, cmaps, titles):

#     # Initialize figure
#     fig,ax = plt.subplots(1,1, figsize=(7,9))

#     # plot values
#     cs = ax.pcolormesh(px,py,var,vmin=vmin, vmax=vmax, cmap=cmap)

#     # Add Puget Sound Inset
#     # [x0, y0, width, height]
#     axins = ax.inset_axes([0.71, 0.0, 0.45, 0.6])
#     # plot values in inset
#     axins.pcolormesh(px, py, var, vmin=vmin, vmax=vmax, cmap=cmap)
#     # Puget Sound limits
#     axins.set_xlim(-123.2, -122.1)
#     axins.set_ylim(46.95, 48.4)
#     axins.tick_params(left=False, bottom=False)
#     # format
#     axins.set_xticklabels([])
#     axins.set_yticklabels([])
#     for spine in axins.spines.values():
#         spine.set_edgecolor('grey')
#         spine.set_linewidth(2)

#     # add colorbar
#     cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
#     cbar.ax.tick_params(labelsize=14, rotation=30)
#     cbar.outline.set_visible(False)

#     # format figure
#     ax.set_xlim([xmin,xmax])
#     ax.set_ylim([ymin,ymax])
#     ax.set_yticklabels([])
#     ax.set_xticklabels([])
#     ax.tick_params(left=False, bottom=False)
#     pfun.add_coast(ax, color='silver')
#     pfun.add_coast(axins, color='silver')
#     pfun.dar(ax)
#     pfun.dar(axins)
#     ax.set_title(title, fontsize=14,
#                 loc='Left', fontweight='bold')

#     # Generate plot
#     plt.tight_layout
#     plt.subplots_adjust(bottom=0.001, top=0.9)
#     plt.show()