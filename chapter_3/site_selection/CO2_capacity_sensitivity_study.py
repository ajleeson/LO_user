"""
Full 10-Year CO₂ Flux Sensitivity Plots at Each Site, Using PyCO2SYS On the Fly.
Atmospheric pCO₂ is allowed to vary for all runs; only a single ocean
property is varied in each scenario, to focus the sensitivity on ocean conditions.
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import gsw
import PyCO2SYS as pyco2
from datetime import datetime
from lo_tools import Lfun,zfun

plt.close('all')

# --- USER INPUTS ---

gtagex = 'cas7_t1_x11ab'
years = ['2019']#,'2016','2017','2018','2019','2020','2021','2022','2023','2024']
data_dir = Lfun.Lstart()['LOo'] / 'chapter_3' / 'data'
grid_file = Lfun.Lstart()['data'] / 'grids/cas7/grid.nc'

# locations=['Columbia River Plume',
#         'Saratoga Passage',
#         'Hood Canal',
#         'Van Island Coast',
#         'Quadra Island',
#         'S. SoG (Fraser plume)',
#         'S. of San Juan Islands']
# natural_ys = [464,875,716,1047,1199,1075,898]
# natural_xs = [360,578,504,186, 233, 468, 513]

locations=['Hood Canal']
natural_ys = [716]
natural_xs = [504]

grid_ds = xr.open_dataset(grid_file)
lon_sites = grid_ds['lon_rho'].values[tuple(natural_ys), tuple(natural_xs)]
lat_sites = grid_ds['lat_rho'].values[tuple(natural_ys), tuple(natural_xs)]

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

# Atmospheric secular trend
D0 = 282.6
D1 = 0.125
D2 =-7.18
D3 = 0.862
D4 =-0.99
D5 = 0.28
D6 =-0.80
D7 = 0.06

def get_PCO2air_secular(dates):
    """Vectorized calculation of pCO2_atm for all datetimes."""
    year = np.array([dt.year for dt in dates])
    yearday = np.array([int(dt.strftime("%j")) for dt in dates])
    pmonth = (year - 1951) + yearday/365
    pi2 = 2 * np.pi
    PCO2air = (D0 + D1*pmonth*12 + D2*np.sin(pi2*pmonth+D3)
                + D4*np.sin(pi2*pmonth+D5) + D6*np.sin(pi2*pmonth+D7))
    return PCO2air

def calc_CO2_flux(temp, salt, alk_vol, tic_vol, wind, pres, lon, lat, pco2_air):
    SA = gsw.SA_from_SP(salt, pres, lon, lat)
    CT = gsw.CT_from_pt(SA, temp)
    rho = gsw.rho(SA, CT, pres)
    alk_kg = 1000 * alk_vol / rho
    tic_kg = 1000 * tic_vol / rho
    CO2dict = pyco2.sys(par1=alk_kg, par1_type=1, par2=tic_kg, par2_type=2,
                        salinity=salt, temperature=temp, pressure=pres,
                        total_silicate=50, total_phosphate=2,
                        opt_k_carbonic=10, opt_buffers_mode=0)
    pCO2_ocean = CO2dict['pCO2']
    Sc = A_CO2 + B_CO2*temp + C_CO2*temp**2 + D_CO2*temp**3 + E_CO2*temp**4
    cff3 = 24/100 * 0.31 * wind**2 * (Sc/660)**(-0.5)
    tempK0p01 = (temp + 273.15)/100
    CO2_sol = np.exp(A1 + A2/tempK0p01 + A3*np.log(tempK0p01) +
                     salt*(B1 + B2*tempK0p01 + B3*tempK0p01**2))
    capacity_ideal = CO2_sol * (pco2_air - pCO2_ocean)
    flux = capacity_ideal * cff3
    return flux, pCO2_ocean

site_files = [data_dir / f"{gtagex}_{year}_freshwatercontent_CO2uptake.nc" for year in years]
ds_all = xr.open_mfdataset(site_files, combine='by_coords')
y_idx = xr.DataArray(natural_ys, dims="site", coords={"site": locations})
x_idx = xr.DataArray(natural_xs, dims="site", coords={"site": locations})
ds_sites = ds_all.isel(eta_rho=y_idx, xi_rho=x_idx)  # shape: (time, site)

time = ds_sites["ocean_time"].values  # (time,)
time_dt = pd.to_datetime(time)
T = ds_sites['surf_temp'].values      # (time,site)
S = ds_sites['surf_salt'].values
A = ds_sites['surf_alk'].values
C = ds_sites['surf_TIC'].values
W2 = ds_sites['wind2'].values
W = np.sqrt(W2)
flux_actual = ds_sites['CO2_flux_actual'].values  # (time,site)

pres = 0

for i, site in enumerate(locations):
    lon = lon_sites[i]
    lat = lat_sites[i]
    temp_ts = T[:,i]
    salt_ts = S[:,i]
    alk_ts  = A[:,i]
    tic_ts  = C[:,i]
    wind_ts = W[:,i]

    # Means for sensitivity scenarios (over ALL years at this site)
    mT, mS, mA, mC, mW = map(np.nanmean, [temp_ts, salt_ts, alk_ts, tic_ts, wind_ts])

    # Vectorized pCO₂_air (shape: n_times) for this location/time
    PCO2air_vec = get_PCO2air_secular(time_dt)

    # Sensitivity runs: for each, pCO2_air varies. Only ONE ocean property varies per run.
    flux_vary_temp, _ = calc_CO2_flux(temp_ts, mS, mA, mC, mW, pres, lon, lat, PCO2air_vec)
    flux_vary_salt, _ = calc_CO2_flux(mT, salt_ts, mA, mC, mW, pres, lon, lat, PCO2air_vec)
    flux_vary_alk,  _ = calc_CO2_flux(mT, mS, alk_ts, mC, mW, pres, lon, lat, PCO2air_vec)
    flux_vary_tic,  _ = calc_CO2_flux(mT, mS, mA, tic_ts, mW, pres, lon, lat, PCO2air_vec)
    flux_vary_wind, _ = calc_CO2_flux(mT, mS, mA, mC, wind_ts, pres, lon, lat, PCO2air_vec)
    # Optionally: only ΔpCO₂ varies (here, still allow PCO2air_vec, as in ocean-only sensitivities)
    flux_actual_test, pCO2_ocean = calc_CO2_flux(temp_ts, salt_ts, alk_ts, tic_ts, wind_ts, pres, lon, lat, PCO2air_vec)
    mTempK0p01 = (mT + 273.15)/100
    mCO2_sol = np.exp(A1 + A2/mTempK0p01 + A3*np.log(mTempK0p01) +
            mS*(B1 + B2*mTempK0p01 + B3*mTempK0p01**2))
    mSc = A_CO2 + B_CO2*mT + C_CO2*mT**2 + D_CO2*mT**3 + E_CO2*mT**4
    mcff3 = 24/100 * 0.31 * mW * (mSc/660)**(-0.5)
    flux_var_dpco2 = mCO2_sol * (PCO2air_vec - pCO2_ocean) * mcff3

    # --------- PLOT ---------
    colors = ['#4477AA','#66CCEE','#228833','#CCBB44','#EE6677','#AA3377','#BBBBBB']
    plt.figure(figsize=(8,6))
    nwin=14
    plt.plot(time_dt, zfun.lowpass(flux_actual[:,i],n=nwin), label='Actual flux', color='k', lw=3)
    plt.plot(time_dt, zfun.lowpass(flux_actual_test,n=nwin), label='Actual flux test', color='deeppink', lw=1, linestyle=':')
    plt.plot(time_dt, zfun.lowpass(flux_vary_temp,n=nwin), label='Vary temp', color='#4477AA')
    plt.plot(time_dt, zfun.lowpass(flux_vary_salt,n=nwin), label='Vary salt', color='#66CCEE')
    plt.plot(time_dt, zfun.lowpass(flux_vary_alk,n=nwin),  label='Vary ALK',  color='#228833')
    plt.plot(time_dt, zfun.lowpass(flux_vary_tic,n=nwin),  label='Vary TIC',  color='#CCBB44')
    plt.plot(time_dt, zfun.lowpass(flux_vary_wind,n=nwin), label='Vary wind', color='#EE6677')
    plt.plot(time_dt, zfun.lowpass(flux_var_dpco2,n=nwin), label=r'Vary $\Delta$pCO2', color='#AA3377')
    plt.axhline(0, color='gray', lw=0.7)
    plt.ylabel(r'CO$_2$ Flux [mmol m$^{-2}$ d$^{-1}$]')
    plt.xlabel('Date')
    plt.title(f'CO2 Flux Sensitivity at {site}\n({nwin}-day Hanning Window)', fontsize=13)
    plt.legend(fontsize=11, loc='best', ncol=2)
    plt.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.show()