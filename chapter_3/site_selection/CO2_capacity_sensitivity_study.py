"""
Phase 2: Calculate spatial sensitivity of CO2 flux using monthly climatologies
"""

import numpy as np
import xarray as xr
import matplotlib.pylab as plt
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

# ##############################################################
# ##                       USER INPUTS                        ##
# ##############################################################

# # Directory where your monthly climatologies are saved
# data_dir = Ldir['LOo'] / 'chapter_3' / 'data'
# months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# # Assume a constant climatological atmospheric pCO2 for the sensitivity tests (e.g., ~2015-2025 average)
# PCO2air_mean = 412.0 

# ##############################################################
# ##                 Important coefficients                   ##
# ##############################################################

# # Schmidt number transfer coefficients (Wanninkhof, 1992)
# A_CO2, B_CO2, C_CO2, D_CO2, E_CO2 = 2073.1, 125.62, 3.6276, 0.043219, 0.0

# # Surface CO2 solubility coefficients (Weiss, 1974)
# A1, A2, A3 = -60.2409, 93.4517, 23.3585
# B1, B2, B3 = 0.023517, -0.023656, 0.0047036

# ##############################################################
# ##                    HELPER FUNCTIONS                      ##
# ##############################################################

# def calc_CO2_flux(temp, salt, alk_vol, tic_vol, wind2, pres, lons, lats, pco2_air):
#     """ Computes 1D array of CO2 flux. Accepts scalars or 1D arrays """
#     # Density & Units (alk_vol and tic_vol are in mmol/m3, we need umol/kg for PyCO2SYS)
#     SA = gsw.SA_from_SP(salt, pres, lons, lats)
#     CT = gsw.CT_from_pt(SA, temp)
#     rho = gsw.rho(SA, CT, pres)
#     alk_kg = 1000 * alk_vol / rho
#     tic_kg = 1000 * tic_vol / rho
    
#     # PyCO2SYS
#     CO2dict = pyco2.sys(par1=alk_kg, par1_type=1, par2=tic_kg, par2_type=2,
#             salinity=salt, temperature=temp, pressure=pres,
#             total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
#     pCO2_ocean = CO2dict['pCO2']
    
#     # Transfer & Solubility
#     Sc = A_CO2 + B_CO2*temp + C_CO2*temp**2 + D_CO2*temp**3 + E_CO2*temp**4
#     cff3 = 24/100 * 0.31 * wind2 * (Sc/660)**(-1/2) # [m/day]
    
#     tempK0p01 = (temp + 273.15)/100
#     CO2_sol = np.exp(A1 + A2/tempK0p01 + A3*np.log(tempK0p01) + 
#                      salt*(B1 + B2*tempK0p01 + B3*tempK0p01**2))
    
#     # Capacity and Flux calculations
#     capacity_ideal = CO2_sol * (pco2_air - pCO2_ocean)
#     flux = capacity_ideal * cff3
    
#     return capacity_ideal, flux, pCO2_ocean

# def map_to_2d(arr_1d, mask):
#     """ Helper to map 1D valid water cells back to 2D """
#     arr_2d = np.full(mask.shape, np.nan)
#     arr_2d[mask == 1] = arr_1d
#     return arr_2d

# ##############################################################
# ##                      PROCESS DATA                        ##
# ##############################################################

# # Get grid data once
# grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
# mask_rho = grid_ds.mask_rho.values
# lon_rho = grid_ds.lon_rho.values
# lat_rho = grid_ds.lat_rho.values
# water_idx = (mask_rho == 1)

# aLon = lon_rho[water_idx]
# aLat = lat_rho[water_idx]
# aPres = 0.0 # Surface pressure is 0 dbar

# print('Starting monthly sensitivity calculations...\n')

# for month in months:
#     print(f'Processing {month}...')
    
#     # Load the monthly climatology
#     ds_in = xr.open_dataset(data_dir / f'monthly_climatology_freshwaterCO2_{month}.nc')
    
#     # Extract 1D surface mean arrays
#     aTemp = ds_in['surf_temp_mean'].values[water_idx]
#     aSalt = ds_in['surf_salt_mean'].values[water_idx]
#     aALK  = ds_in['surf_alk_mean'].values[water_idx]
#     aTIC  = ds_in['surf_TIC_mean'].values[water_idx]
#     aWind2= ds_in['wind2_mean'].values[water_idx]
    
#     # Calculate True Domain Averages
#     mTemp = aTemp.mean()
#     mSalt = aSalt.mean()
#     mALK  = aALK.mean()
#     mTIC  = aTIC.mean()
#     mWind2= aWind2.mean()
    
#     # ---------------- RUN SCENARIOS ----------------
#     # 1. Fully varying (Actual Climatological Mean Flux)
#     _, flux_actual, pco2_ocean_actual = calc_CO2_flux(aTemp, aSalt, aALK, aTIC, aWind2, aPres, aLon, aLat, PCO2air_mean)
    
#     # 2. Varying ONE variable at a time
#     _, flux_var_temp, _ = calc_CO2_flux(aTemp, mSalt, mALK, mTIC, mWind2, aPres, aLon, aLat, PCO2air_mean)
#     _, flux_var_alk, _  = calc_CO2_flux(mTemp, mSalt, aALK, mTIC, mWind2, aPres, aLon, aLat, PCO2air_mean)
#     _, flux_var_tic, _  = calc_CO2_flux(mTemp, mSalt, mALK, aTIC, mWind2, aPres, aLon, aLat, PCO2air_mean)
#     _, flux_var_salt, _ = calc_CO2_flux(mTemp, aSalt, mALK, mTIC, mWind2, aPres, aLon, aLat, PCO2air_mean)
#     _, flux_var_wind, _ = calc_CO2_flux(mTemp, mSalt, mALK, mTIC, aWind2, aPres, aLon, aLat, PCO2air_mean)
    
#     # 3. Varying ONLY Delta pCO2
#     mTempK0p01 = (mTemp + 273.15) / 100
#     mCO2_sol = np.exp(A1 + A2/mTempK0p01 + A3*np.log(mTempK0p01) + mSalt*(B1 + B2*mTempK0p01 + B3*mTempK0p01**2))
#     mSc = A_CO2 + B_CO2*mTemp + C_CO2*mTemp**2 + D_CO2*mTemp**3 + E_CO2*mTemp**4
#     mcff3 = 24/100 * 0.31 * mWind2 * (mSc/660)**(-1/2)
#     flux_var_dpco2 = mCO2_sol * (PCO2air_mean - pco2_ocean_actual) * mcff3

#     # ---------------- SAVE DATA ----------------
#     ds_out = xr.Dataset(
#         data_vars=dict(
#             flux_actual    = (['eta_rho','xi_rho'], map_to_2d(flux_actual, mask_rho)),
#             flux_var_temp  = (['eta_rho','xi_rho'], map_to_2d(flux_var_temp, mask_rho)),
#             flux_var_alk   = (['eta_rho','xi_rho'], map_to_2d(flux_var_alk, mask_rho)),
#             flux_var_tic   = (['eta_rho','xi_rho'], map_to_2d(flux_var_tic, mask_rho)),
#             flux_var_salt  = (['eta_rho','xi_rho'], map_to_2d(flux_var_salt, mask_rho)),
#             flux_var_wind  = (['eta_rho','xi_rho'], map_to_2d(flux_var_wind, mask_rho)),
#             flux_var_dpco2 = (['eta_rho','xi_rho'], map_to_2d(flux_var_dpco2, mask_rho)),
#         ),
#         coords=dict(eta_rho=grid_ds.eta_rho.values, xi_rho=grid_ds.xi_rho.values)
#     )
    
#     # Save the sensitivity dataset
#     ds_out.to_netcdf(data_dir / f'monthly_sensitivity_CO2_{month}.nc')
#     ds_in.close()

# print('Done processing all months!')

##############################################
# TEST PLOTTING (Example: August)
##############################################
month_to_plot = 'Aug'
ds_plot = xr.open_dataset(data_dir / f'monthly_sensitivity_CO2_{month_to_plot}.nc')

# Calculate the purely spatial effect by subtracting the domain average (baseline)
# The baseline is roughly the mean of the actual flux.
baseline = np.nanmean(ds_plot['flux_actual'].values)

spatial_effect_temp = ds_plot['flux_var_temp'].values - baseline
spatial_effect_tic  = ds_plot['flux_var_tic'].values - baseline

fig, axes = plt.subplots(1, 2, figsize=(14, 7))
px, py = pfun.get_plon_plat(lon_rho, lat_rho)

# Plot Temp Effect
cs1 = axes[0].pcolormesh(px, py, spatial_effect_temp, vmin=-10, vmax=10, cmap=cmc.vik)
axes[0].set_title(f'{month_to_plot}: Spatial Anomaly Driven by Temp', fontweight='bold')
fig.colorbar(cs1, ax=axes[0], orientation='horizontal', pad=0.05).set_label('mmol m$^{-2}$ d$^{-1}$')

# Plot TIC Effect
cs2 = axes[1].pcolormesh(px, py, spatial_effect_tic, vmin=-10, vmax=10, cmap=cmc.vik)
axes[1].set_title(f'{month_to_plot}: Spatial Anomaly Driven by TIC', fontweight='bold')
fig.colorbar(cs2, ax=axes[1], orientation='horizontal', pad=0.05).set_label('mmol m$^{-2}$ d$^{-1}$')

for ax in axes:
    pfun.add_coast(ax, color='silver')
    pfun.dar(ax)
    ax.set_xlim([-126, -122])
    ax.set_ylim([45.5, 50.5])
    ax.tick_params(left=False, bottom=False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

plt.tight_layout()
plt.show()