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
import matplotlib.dates as mdates
import gsw
import PyCO2SYS as pyco2
from datetime import datetime
from lo_tools import Lfun,zfun

plt.close('all')

# --- USER INPUTS ---

gtagex = 'cas7_t1_x11ab'
years = ['2015','2016']#,'2017','2018','2019','2020','2021','2022','2023','2024']
data_dir = Lfun.Lstart()['LOo'] / 'chapter_3' / 'data'
grid_file = Lfun.Lstart()['data'] / 'grids/cas7/grid.nc'

locations=['Columbia River Plume',
        'Saratoga Passage',
        'Hood Canal',
        'Van Island Coast',
        'Quadra Island',
        'S. SoG (Fraser plume)',
        'S. of San Juan Islands']
natural_ys = [464,875,716,1047,1199,1075,898]
natural_xs = [360,578,504,186, 233, 468, 513]

# locations=['Hood Canal']
# natural_ys = [716]
# natural_xs = [504]

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
    mcff3 = 24/100 * 0.31 * mW**2 * (mSc/660)**(-0.5)
    flux_var_dpco2 = mCO2_sol * (PCO2air_vec - pCO2_ocean) * mcff3

    # --------- PLOT ---------
    colors = ['#4477AA','#66CCEE','#228833','#CCBB44','#EE6677','#AA3377','#BBBBBB']
    fig, ax = plt.subplots(1,1,figsize = (13,6))
    nwin=14
    ax.plot(time_dt, zfun.lowpass(flux_actual[:,i],n=nwin), label='Actual flux', color='k', lw=3)
    # plt.plot(time_dt, zfun.lowpass(flux_actual_test,n=nwin), label='Actual flux test', color='deeppink', lw=1, linestyle=':')
    ax.plot(time_dt, zfun.lowpass(flux_vary_temp,n=nwin), label='Vary temp', color='#4477AA')
    ax.plot(time_dt, zfun.lowpass(flux_vary_salt,n=nwin), label='Vary salt', color='#66CCEE')
    ax.plot(time_dt, zfun.lowpass(flux_vary_alk,n=nwin),  label='Vary ALK',  color='#228833')
    ax.plot(time_dt, zfun.lowpass(flux_vary_tic,n=nwin),  label='Vary TIC',  color='#CCBB44')
    ax.plot(time_dt, zfun.lowpass(flux_vary_wind,n=nwin), label='Vary wind', color='#EE6677')
    ax.plot(time_dt, zfun.lowpass(flux_var_dpco2,n=nwin), label=r'Vary $\Delta$pCO2', color='#AA3377')

    # format figure
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.tick_params(axis='x',rotation=30, labelsize=12)
    ax.tick_params(axis='y',labelsize=12)
    ax.axhline(0, color='silver', lw=1, linestyle='--')
    ax.grid(visible=True, axis='x', color='silver', linestyle='--')
    ax.set_ylabel(r'CO$_2$ Flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=12)
    ax.set_xlabel('Year', fontsize=12)
    ax.set_title(f'CO2 Flux Sensitivity at {site}\n({nwin}-day Hanning Window)', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12, loc='best', ncol=2)
    plt.tight_layout()
    plt.show()

# """
# CO2 Flux Sensitivity Analysis on Climatological (1-year) Daily Means, per Site.
# Atmospheric pCO₂ varies according to the algorithm, all variables are annual daily climatologies.
# """
# import numpy as np
# import xarray as xr
# import matplotlib.pyplot as plt
# import pandas as pd
# import gsw
# import PyCO2SYS as pyco2
# from lo_tools import Lfun, zfun

# import sys
# from pathlib import Path
# pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
# if str(pth) not in sys.path:
#     sys.path.append(str(pth))
# import gfun

# Gr = gfun.gstart()

# Ldir = Lfun.Lstart()

# plt.close('all')

# # --- USER INPUTS ---

# data_dir = Lfun.Lstart()['LOo'] / 'chapter_3' / 'data'
# clim_file = data_dir / 'site_daily_climatologies.nc'
# grid_file = Lfun.Lstart()['data'] / 'grids/cas7/grid.nc'

# locations = ['Columbia River Plume', 'Saratoga Passage', 'Hood Canal', 'Tofino',
#     'Quadra Island', 'S. SoG (Fraser plume)', 'S. of San Juan Islands']
# natural_ys = [464, 875, 716, 1047, 1199, 1075, 898]
# natural_xs = [360, 578, 504, 186, 233, 468, 513]

# # Extract correct lon/lat for all sites in order:
# grid_ds = xr.open_dataset(grid_file)
# lon_sites = grid_ds['lon_rho'].values[tuple(natural_ys), tuple(natural_xs)]
# lat_sites = grid_ds['lat_rho'].values[tuple(natural_ys), tuple(natural_xs)]

# # --- Important coefficients ---
# A_CO2, B_CO2, C_CO2, D_CO2, E_CO2 = 2073.1, 125.62, 3.6276, 0.043219, 0.0
# A1, A2, A3 = -60.2409, 93.4517, 23.3585
# B1, B2, B3 = 0.023517, -0.023656, 0.0047036
# D0 = 282.6
# D1 = 0.125
# D2 =-7.18
# D3 = 0.862
# D4 =-0.99
# D5 = 0.28
# D6 =-0.80
# D7 = 0.06

# def get_PCO2air_secular(dates):
#     """Vectorized calculation of pCO2_atm for all datetimes."""
#     year = np.array([dt.year for dt in dates])
#     yearday = np.array([int(dt.strftime("%j")) for dt in dates])
#     pmonth = (year - 1951) + yearday / 365
#     pi2 = 2 * np.pi
#     PCO2air = (D0 + D1*pmonth*12 + D2*np.sin(pi2*pmonth+D3)
#                 + D4*np.sin(pi2*pmonth+D5) + D6*np.sin(pi2*pmonth+D7))
#     return PCO2air

# def calc_CO2_flux(temp, salt, alk_vol, tic_vol, wind, pres, lon, lat, pco2_air):
#     SA = gsw.SA_from_SP(salt, pres, lon, lat)
#     CT = gsw.CT_from_pt(SA, temp)
#     rho = gsw.rho(SA, CT, pres)
#     alk_kg = 1000 * alk_vol / rho
#     tic_kg = 1000 * tic_vol / rho
#     CO2dict = pyco2.sys(par1=alk_kg, par1_type=1, par2=tic_kg, par2_type=2,
#                         salinity=salt, temperature=temp, pressure=pres,
#                         total_silicate=50, total_phosphate=2,
#                         opt_k_carbonic=10, opt_buffers_mode=0)
#     pCO2_ocean = CO2dict['pCO2']
#     Sc = A_CO2 + B_CO2 * temp + C_CO2 * temp**2 + D_CO2 * temp**3 + E_CO2 * temp**4
#     # -- Main bugfix: wind² should be used! --
#     cff3 = 24/100 * 0.31 * (wind**2) * (Sc/660)**(-0.5)
#     tempK0p01 = (temp + 273.15)/100
#     CO2_sol = np.exp(A1 + A2/tempK0p01 + A3*np.log(tempK0p01) +
#                      salt * (B1 + B2*tempK0p01 + B3*tempK0p01**2))
#     capacity_ideal = CO2_sol * (pco2_air - pCO2_ocean)
#     flux = capacity_ideal * cff3
#     return flux, pCO2_ocean

# # --- LOAD CLIMATOLOGY ---
# clim_ds = xr.open_dataset(clim_file)
# days = clim_ds['dayofyear'].values                   # day (1..366)
# dates = pd.to_datetime('2001-01-01') + pd.to_timedelta(days-1, unit='D')  # Dummy year for solar/seasonal cycle

# for i, site in enumerate(locations):
#     lon = lon_sites[i]
#     lat = lat_sites[i]
#     temp = clim_ds['surf_temp_mean'].sel(site=site).values
#     salt = clim_ds['surf_salt_mean'].sel(site=site).values
#     alk = clim_ds['surf_alk_mean'].sel(site=site).values
#     tic = clim_ds['surf_TIC_mean'].sel(site=site).values
#     wind2 = clim_ds['wind2_mean'].sel(site=site).values
#     wind = np.sqrt(wind2)

#     # Means for sensitivity scenarios over climatological year
#     mT, mS, mA, mC, mW = map(np.nanmean, [temp, salt, alk, tic, wind])

#     pres = 0

#     # Atmospheric pCO₂, always seasonally varying
#     PCO2air_vec = get_PCO2air_secular(dates)

#     # Actual flux (recompute, don't use stored, guarantees apples-to-apples for sensitivity)
#     flux_actual_test, _ = calc_CO2_flux(temp, salt, alk, tic, wind, pres, lon, lat, PCO2air_vec)
#     # Sensitivity runs (only one varies)
#     flux_vary_temp, _ = calc_CO2_flux(temp, mS, mA, mC, mW, pres, lon, lat, PCO2air_vec)
#     flux_vary_salt, _ = calc_CO2_flux(mT, salt, mA, mC, mW, pres, lon, lat, PCO2air_vec)
#     flux_vary_alk, _ = calc_CO2_flux(mT, mS, alk, mC, mW, pres, lon, lat, PCO2air_vec)
#     flux_vary_tic, _ = calc_CO2_flux(mT, mS, mA, tic, mW, pres, lon, lat, PCO2air_vec)
#     flux_vary_wind, _ = calc_CO2_flux(mT, mS, mA, mC, wind, pres, lon, lat, PCO2air_vec)
#     # Only delta pCO₂ varies (all other processes mean, but seasonal atm and ocean pCO₂)
#     _, pCO2_ocean = calc_CO2_flux(temp, salt, alk, tic, wind, pres, lon, lat, PCO2air_vec)
#     mTempK0p01 = (mT + 273.15)/100
#     mCO2_sol = np.exp(A1 + A2/mTempK0p01 + A3*np.log(mTempK0p01) +
#             mS*(B1 + B2*mTempK0p01 + B3*mTempK0p01**2))
#     mSc = A_CO2 + B_CO2*mT + C_CO2*mT**2 + D_CO2*mT**3 + E_CO2*mT**4
#     mcff3 = 24/100 * 0.31 * (mW**2) * (mSc/660)**(-0.5)
#     flux_var_dpco2 = mCO2_sol * (PCO2air_vec - pCO2_ocean) * mcff3

#     # get values I calculated
#     out_dir = Ldir['LOo'] / 'chapter_3' / 'data'
#     clim_ds = xr.open_dataset(out_dir / 'site_daily_climatologies.nc')
#     vn = 'CO2_flux_actual'
#     flux_actual = clim_ds[f'{vn}_mean'].sel(site=site).values

#     # --- Plot ---
#     nwin = 14
#     plt.figure(figsize=(10,6))
#     plt.plot(dates, zfun.lowpass(flux_actual, n=nwin), label='Actual', color='k', lw=3)
#     plt.plot(dates, zfun.lowpass(flux_actual_test, n=nwin), label='Actual test', color='deeppink', linestyle=':', lw=2)
#     plt.plot(dates, zfun.lowpass(flux_vary_temp, n=nwin), label='Vary temp', color='#4477AA')
#     plt.plot(dates, zfun.lowpass(flux_vary_salt, n=nwin), label='Vary salt', color='#66CCEE')
#     plt.plot(dates, zfun.lowpass(flux_vary_alk, n=nwin),  label='Vary ALK',  color='#228833')
#     plt.plot(dates, zfun.lowpass(flux_vary_tic, n=nwin),  label='Vary TIC',  color='#CCBB44')
#     plt.plot(dates, zfun.lowpass(flux_vary_wind, n=nwin), label='Vary wind', color='#EE6677')
#     plt.plot(dates, zfun.lowpass(flux_var_dpco2, n=nwin), label=r'Vary $\Delta$pCO2', color='#AA3377')
#     plt.axhline(0, color='gray', lw=0.7)
#     plt.ylabel(r'CO$_2$ Flux [mmol m$^{-2}$ d$^{-1}$]')
#     plt.xlabel('Day of Year')
#     plt.title(f'CO2 Flux Sensitivity at {site} (Climatology, {nwin}-day lowpass)', fontsize=13)
#     plt.legend(fontsize=11, loc='best', ncol=2)
#     plt.grid(True, alpha=0.2)
#     plt.tight_layout()
#     plt.show()