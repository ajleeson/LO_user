"""
WARNING: THIS TAKES A LONG TIME TO RUN! 
USEFUL FOR ONE-DAY CASE STUDIES, BUT NOT TO BE RUN ON LONGER TIME PERIODS ON APOGEE
Get CO2 uptake capacity from model history file
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

years = ['2020']#['2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025']

# which  model run to look at?
gtagex = 'cas7_t1noDIN_x11ab' #'cas7_t1_x11ab'

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
        
        # CO2 flux simulated-- when we set all values constant and vary one at a time
        CO2_flux_varied_temp   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        CO2_flux_varied_alk   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        CO2_flux_varied_TIC   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        CO2_flux_varied_wind   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        CO2_flux_varied_salt   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        CO2_flux_varied_dpco2  = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
            
        ),
    coords=dict(ocean_time=ocean_time, eta_rho=eta_rho, xi_rho=xi_rho,),)
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for processed CO2 uptake capacity
    '''
    ds['CO2_capacity_ideal'].attrs['long_name'] = 'Ideal capacity of CO2 uptake based on just solubility and Delta pCO2-- given as a concentration'
    ds['CO2_capacity_ideal'].attrs['units'] = 'mmol CO2 m-3'

    ds['CO2_flux_actual'].attrs['long_name'] = 'Total flux of CO2 over one day, taking into account winds'
    ds['CO2_flux_actual'].attrs['units'] = 'mmol CO2 m-2 day-1'

    # CO2 flux simulated-- when we set all values constant and vary one at a time
    ds['CO2_flux_varied_temp'].attrs['long_name'] = 'CO2 flux over one day if all variables are constant except temp'
    ds['CO2_flux_varied_temp'].attrs['units'] = 'mmol CO2 m-2 day-1'

    ds['CO2_flux_varied_alk'].attrs['long_name'] = 'CO2 flux over one day if all variables are constant except alkalinity'
    ds['CO2_flux_varied_alk'].attrs['units'] = 'mmol CO2 m-2 day-1'

    ds['CO2_flux_varied_TIC'].attrs['long_name'] = 'CO2 flux over one day if all variables are constant except TIC'
    ds['CO2_flux_varied_TIC'].attrs['units'] = 'mmol CO2 m-2 day-1'

    ds['CO2_flux_varied_wind'].attrs['long_name'] = 'CO2 flux over one day if all variables are constant except wind'
    ds['CO2_flux_varied_wind'].attrs['units'] = 'mmol CO2 m-2 day-1'

    ds['CO2_flux_varied_salt'].attrs['long_name'] = 'CO2 flux over one day if all variables are constant except salt'
    ds['CO2_flux_varied_salt'].attrs['units'] = 'mmol CO2 m-2 day-1'
    
    ds['CO2_flux_varied_dpco2'].attrs['long_name'] = 'CO2 flux over one day varying only Delta pCO2 (constant solubility and wind)'
    ds['CO2_flux_varied_dpco2'].attrs['units'] = 'mmol CO2 m-2 day-1'

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
    Ldir['roms_out'] = Ldir['roms_out']#Ldir['roms_out5']
    Ldir['ds0'] = year + '.05.22' #'.01.01'
    Ldir['ds1'] = year + '.05.22' #'.12.31'
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
        
        # GET Z_RHO ---------------------------------------------------------------------------------------
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
        
        # 1. Get basic 2D variables and mask
        mask_rho = ds_raw['mask_rho'].values 
        water_idx = (mask_rho == 1)
        
        # 2. Extract 1D Arrays (only water cells)
        aLon = ds_raw['lon_rho'].values[water_idx]
        aLat = ds_raw['lat_rho'].values[water_idx]
        aPres = gsw.p_from_z(z_rho[-1,:,:], ds_raw['lon_rho'].values, ds_raw['lat_rho'].values).squeeze()[water_idx]

        aTemp = ds_raw['temp'].values[:,-1,:,:].squeeze()[water_idx]
        aSalt = ds_raw['salt'].values[:,-1,:,:].squeeze()[water_idx]
        aALK  = ds_raw['alkalinity'].values[:,-1,:,:].squeeze()[water_idx]
        aTIC  = ds_raw['TIC'].values[:,-1,:,:].squeeze()[water_idx]
        
        windspeed2 = (ds_raw['Uwind'].values**2 + ds_raw['Vwind'].values**2).squeeze()
        aWind2 = windspeed2[water_idx]

        # 3. Calculate True Domain Averages (Water only!)
        mTemp = aTemp.mean()
        mSalt = aSalt.mean()
        mALK  = aALK.mean()
        mTIC  = aTIC.mean()
        mWind2= aWind2.mean()
        
        # GET PCO2 OF ATMOSPHERE [uatm] -------------------------------------------------------------------
        date_obj = datetime.strptime(date_str[1::], "%Y.%m.%d")
        yearday = int(date_obj.strftime("%j")) # yearday from date
        pmonth = int(year) - 1951 + yearday/365 # months since January 1951
        pi2 = 2 * np.pi
        PCO2air_secular = D0 + D1*pmonth*12 + D2*np.sin(pi2*pmonth+D3) + D4*np.sin(pi2*pmonth+D5) + D6*np.sin(pi2*pmonth+D7)

        # 4. DEFINE HELPER FUNCTION TO DO ALL MATH AT ONCE ------------------------------------------------
        def calc_CO2_flux(temp, salt, alk_vol, tic_vol, wind2, pres, lons, lats, pco2_air):
            """ Computes 1D array of CO2 flux. Accepts scalars or 1D arrays """
            # Density & Units
            SA = gsw.SA_from_SP(salt, pres, lons, lats)
            CT = gsw.CT_from_pt(SA, temp)
            rho = gsw.rho(SA, CT, pres)
            alk_kg = 1000 * alk_vol / rho
            tic_kg = 1000 * tic_vol / rho
            
            # PyCO2SYS
            CO2dict = pyco2.sys(par1=alk_kg, par1_type=1, par2=tic_kg, par2_type=2,
                    salinity=salt, temperature=temp, pressure=pres,
                    total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
            pCO2_ocean = CO2dict['pCO2']
            
            # Transfer & Solubility
            Sc = A_CO2 + B_CO2*temp + C_CO2*temp**2 + D_CO2*temp**3 + E_CO2*temp**4
            cff3 = 24/100 * 0.31 * wind2 * (Sc/660)**(-1/2) # [m] (over one day)
            
            tempK0p01 = (temp + 273.15)/100
            CO2_sol = np.exp(A1 + A2/tempK0p01 + A3*np.log(tempK0p01) + 
                             salt*(B1 + B2*tempK0p01 + B3*tempK0p01**2))
            
            # Capacity and Flux calculations
            capacity_ideal = CO2_sol * (pco2_air - pCO2_ocean)
            flux = capacity_ideal * cff3
            
            return capacity_ideal, flux, pCO2_ocean

        # 5. RUN SCENARIOS --------------------------------------------------------------------------------
        cap_ideal, flux_actual, pco2_ocean_actual = calc_CO2_flux(aTemp, aSalt, aALK, aTIC, aWind2, aPres, aLon, aLat, PCO2air_secular)
        _, flux_var_temp, _ = calc_CO2_flux(aTemp, mSalt, mALK, mTIC, mWind2, aPres, aLon, aLat, PCO2air_secular)
        _, flux_var_alk, _  = calc_CO2_flux(mTemp, mSalt, aALK, mTIC, mWind2, aPres, aLon, aLat, PCO2air_secular)
        _, flux_var_tic, _  = calc_CO2_flux(mTemp, mSalt, mALK, aTIC, mWind2, aPres, aLon, aLat, PCO2air_secular)
        _, flux_var_salt, _ = calc_CO2_flux(mTemp, aSalt, mALK, mTIC, mWind2, aPres, aLon, aLat, PCO2air_secular)
        _, flux_var_wind, _ = calc_CO2_flux(mTemp, mSalt, mALK, mTIC, aWind2, aPres, aLon, aLat, PCO2air_secular)

        # 6. RUN DPCO2 SCENARIO ---------------------------------------------------------------------------
        # Calculate mean solubility and mean gas transfer velocity
        mTempK0p01 = (mTemp + 273.15) / 100
        mCO2_sol = np.exp(A1 + A2/mTempK0p01 + A3*np.log(mTempK0p01) + 
                          mSalt*(B1 + B2*mTempK0p01 + B3*mTempK0p01**2))
        
        mSc = A_CO2 + B_CO2*mTemp + C_CO2*mTemp**2 + D_CO2*mTemp**3 + E_CO2*mTemp**4
        mcff3 = 24/100 * 0.31 * mWind2 * (mSc/660)**(-1/2)
        
        # Calculate flux varying ONLY the actual gradient (pCO2air - pCO2ocean), holding solubility & wind constant
        flux_var_dpco2 = mCO2_sol * (PCO2air_secular - pco2_ocean_actual) * mcff3


        # Helper to map 1D back to 3D (adding the ocean_time dimension back)
        def map_to_3d(arr_1d, mask):
            arr_2d = np.full(mask.shape, np.nan)
            arr_2d[mask == 1] = arr_1d
            # Add the time dimension back to the front: shape becomes (1, eta_rho, xi_rho)
            return np.expand_dims(arr_2d, axis=0)

        # add data to daily dataset -----------------------------------------------------------------------
        daily_ds['CO2_capacity_ideal'] = xr.DataArray(map_to_3d(cap_ideal, mask_rho),
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        
        daily_ds['CO2_flux_actual'] = xr.DataArray(map_to_3d(flux_actual, mask_rho),
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        
        daily_ds['CO2_flux_varied_temp'] = xr.DataArray(map_to_3d(flux_var_temp, mask_rho),
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        
        daily_ds['CO2_flux_varied_alk'] = xr.DataArray(map_to_3d(flux_var_alk, mask_rho),
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        
        daily_ds['CO2_flux_varied_TIC'] = xr.DataArray(map_to_3d(flux_var_tic, mask_rho),
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        
        daily_ds['CO2_flux_varied_salt'] = xr.DataArray(map_to_3d(flux_var_salt, mask_rho),
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        
        daily_ds['CO2_flux_varied_wind'] = xr.DataArray(map_to_3d(flux_var_wind, mask_rho),
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
                                    
        daily_ds['CO2_flux_varied_dpco2'] = xr.DataArray(map_to_3d(flux_var_dpco2, mask_rho),
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])

        # Append to list
        ds_list.append(daily_ds)

    # PREPARE DATA FOR SAVING
    # Combine all of the daily datasets into one dataset
    ds = xr.concat(ds_list, dim='ocean_time')
    # Add your metadata once at the end
    ds = add_metadata(ds)
    # print('    Saving dataset')
    # ds.to_netcdf(out_dir / (gtagex + '_' + year + '_CO2uptake_capacity.nc'))

    print(ds)

print('Done')


##############################################
# TEST PLOTTING

# get pcolormesh values
val_cap_ideal   = ds['CO2_capacity_ideal'][-1,:,:].values
val_flux_actual = ds['CO2_flux_actual'][-1,:,:].values
val_flux_temp   = ds['CO2_flux_varied_temp'][-1,:,:].values
val_flux_alk    = ds['CO2_flux_varied_alk'][-1,:,:].values
val_flux_tic    = ds['CO2_flux_varied_TIC'][-1,:,:].values
val_flux_salt   = ds['CO2_flux_varied_salt'][-1,:,:].values
val_flux_wind   = ds['CO2_flux_varied_wind'][-1,:,:].values
val_flux_dpco2  = ds['CO2_flux_varied_dpco2'][-1,:,:].values

# make list of vars
vars = [
    val_cap_ideal, 
    val_flux_actual, 
    val_flux_temp, 
    val_flux_alk, 
    val_flux_tic, 
    val_flux_salt, 
    val_flux_wind,
    val_flux_dpco2
]

# list of vmins and vmax (Capacity gets 15, all fluxes get 50)
vmins = [-15, -50, -50, -50, -50, -50, -50, -50]
vmaxs = [ 15,  50,  50,  50,  50,  50,  50,  50]

# get colormaps
cmaps = [cmc.vik] * 8

# titles
flux_unit = r'[mmol CO2 m$^{-2}$ day$^{-1}$]'
titles = [
    'CO2 uptake capacity\n' + r'based on solubility and $\Delta$ pCO2 [mmol CO2 m$^{-3}$]',
    'Total daily CO2 flux\n' + flux_unit,
    'CO2 flux (Varying Temp only)\n' + flux_unit,
    'CO2 flux (Varying Alkalinity only)\n' + flux_unit,
    'CO2 flux (Varying TIC only)\n' + flux_unit,
    'CO2 flux (Varying Salt only)\n' + flux_unit,
    'CO2 flux (Varying Wind only)\n' + flux_unit,
    'CO2 flux (Varying $\Delta$ pCO2 only)\n' + flux_unit
]

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values
px, py = pfun.get_plon_plat(lon,lat)

# lon/lat limits (Study Domain)
xmin = -126
xmax = -122
ymin = 45.5
ymax = 50.5

for var, vmin, vmax, cmap, title in zip(vars, vmins, vmaxs, cmaps, titles):

    # Initialize figure
    fig,ax = plt.subplots(1,1, figsize=(7,9))

    # plot values
    cs = ax.pcolormesh(px,py,var,vmin=vmin, vmax=vmax, cmap=cmap)

    # Add Puget Sound Inset
    # [x0, y0, width, height]
    axins = ax.inset_axes([0.71, 0.0, 0.45, 0.6])
    # plot values in inset
    axins.pcolormesh(px, py, var, vmin=vmin, vmax=vmax, cmap=cmap)
    # Puget Sound limits
    axins.set_xlim(-123.2, -122.1)
    axins.set_ylim(46.95, 48.4)
    axins.tick_params(left=False, bottom=False)
    # format
    axins.set_xticklabels([])
    axins.set_yticklabels([])
    for spine in axins.spines.values():
        spine.set_edgecolor('grey')
        spine.set_linewidth(2)

    # add colorbar
    cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
    cbar.ax.tick_params(labelsize=14, rotation=30)
    cbar.outline.set_visible(False)

    # format figure
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.tick_params(left=False, bottom=False)
    pfun.add_coast(ax, color='silver')
    pfun.add_coast(axins, color='silver')
    pfun.dar(ax)
    pfun.dar(axins)
    ax.set_title(title, fontsize=14,
                loc='Left', fontweight='bold')

    # Generate plot
    plt.tight_layout
    plt.subplots_adjust(bottom=0.001, top=0.9)
    plt.show()