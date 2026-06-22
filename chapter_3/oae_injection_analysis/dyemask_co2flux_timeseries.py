"""
Plot OAE time series
"""

# import things
import numpy as np
import xarray as xr
import csv
import matplotlib.pylab as plt
import matplotlib.colors as colors
from pathlib import Path
from datetime import datetime
import matplotlib.colors as mcolors
import pandas as pd
import matplotlib.dates as mdates
import gsw
from matplotlib.patches import Rectangle
import PyCO2SYS as pyco2
import cmcrameri.cm as cmc
from lo_tools import Lfun, zrfun, zfun
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

# date for map
date = '2020.07.31' 
# convert date to format yyy-mm-dd
date_formatted = date.replace('.','-')

# which  model runs to look at?
basline = 'cas7_t1_x11ab'
perturbation = 'cas7_t1dgeWB_x11abd'

##############################################################
##                      Get map data                        ##
##############################################################

ds_base = xr.open_dataset(Ldir['roms_out'] / basline      / ('f' + date) / 'ocean_avg_0001.nc')
ds_pert = xr.open_dataset(Ldir['roms_out'] / perturbation / ('f' + date) / 'ocean_avg_0001.nc')

surf_alk_base = ds_base['alkalinity'].values[0,-1,:,:]
surf_alk_pert = ds_pert['alkalinity'].values[0,-1,:,:]

# surf_alk_base = ds_base['alkalinity'].values[0,20,:,:]
# surf_alk_pert = ds_pert['alkalinity'].values[0,20,:,:]

# mid-water column DIC
midD_dic_base = ds_base['TIC'].values[0,20,:,:]
midD_dic_pert = ds_pert['TIC'].values[0,20,:,:]

# get surface dye
surf_dye = ds_pert['dye_01'].values[0,-1,:,:]
# surf_dye = ds_pert['dye_01'].values[0,20,:,:]
# convert units from kg/m3 to mmol/m3, to equate to alkaliinty
surf_dye_alk_units = surf_dye / 1.7e-5 

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
        
# GET Z_RHO
# get S for the whole grid
Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
reader = csv.DictReader(open(Sfp))
S_dict = {}
for row in reader:
    S_dict[row['ITEMS']] = row['VALUES']
S = zrfun.get_S(S_dict)
# get cell thickness (can use either ds because hydrodynamics identical)
h = ds_base['h'].values # height of water column
zeta = ds_base['zeta'].values
z_rho, z_w = zrfun.get_z(h, zeta[0, :, :], S)

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
date_obj = datetime.strptime(date, "%Y.%m.%d")
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

##############################################################
##                  Get time series data                    ##
##############################################################

# Using dye mask (threshold = 5e-4 mmol/m3)
# get the first month that has continuous alkalinity and dye addition (cropped subdomain)
fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'dyemask_CO2flux_2020.06.01_2020.06.30.nc'
ds_addition = xr.open_dataset(fp)
# # crop to just first 30 days
# ds_addition = ds_addition.isel(ocean_time=slice(0,30))
# get the remaining time after the first addition
fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'dyemask_allintegral_oae_deltas_2020.07.01_2020.10.31.nc'
ds_later = xr.open_dataset(fp)

# combine the two datasets
ds_combined = xr.concat([ds_addition, ds_later], dim='ocean_time')

time_combined = ds_combined.ocean_time.values
delta_DIC_combined = ds_combined.delta_DIC_full.values # kmol
delta_Alk_combined = ds_combined.delta_Alk_full.values # kmol
total_dye_combined = ds_combined.total_dye_full.values # kmol
CO2_flux_combined  = ds_combined.total_CO2_flux_full.values # kmol

# get index of date in time_combined
target = np.datetime64(date_formatted+'T12:00:00')
t_index = np.where(time_combined == target)[0][0]

###################################################################
##                         Plotting                              ##  
################################################################### 

# get grid info
lons = ds_base.coords['lon_rho'].values
lats = ds_base.coords['lat_rho'].values
px, py = pfun.get_plon_plat(lons,lats)

# injection location
inj_lon = -122.5674
inj_lat = 48.1956

# Initialize figure
plt.close('all')
def generate_axes(fig):
    gridspec = fig.add_gridspec(nrows=6, ncols=15, wspace=0.5, hspace=0.5)
    ax = {}
    ax['alkT'] = fig.add_subplot(gridspec[0:2, 0:3])
    ax['dicT'] = fig.add_subplot(gridspec[2:4, 0:3], sharex=ax['alkT'])
    ax['effT'] = fig.add_subplot(gridspec[4:6, 0:3], sharex=ax['alkT'])
    ax['dyeM'] = fig.add_subplot(gridspec[0:6, 3:7])
    ax['alkM'] = fig.add_subplot(gridspec[0:6, 7:11])
    ax['CO2M'] = fig.add_subplot(gridspec[0:6, 11:15])
    return ax

# Initialize figure
plt.close('all')
fig = plt.figure(figsize=(12,6))
ax = generate_axes(fig)

# color scale for surf dye and alk
vmin = -1e-2
vmax =  1e-2

# colormaps
surf_cmap = cmc.broc_r
surf_cmap.set_bad(color='darkgray')
flux_cmap = cmc.lajolla_r
flux_cmap.set_bad(color='darkgray')

# dye colormaps TO FIGURE OUT THRESHOLD FOR MASKING DYE
threshold = 5e-4
data = surf_dye_alk_units.copy()
data[data < threshold] = 0
cmap = cmc.devon_r.copy()      # choose any Crameri colormap you like
cmap.set_bad('darkgray')    # color for masked values
cs = ax['dyeM'].pcolormesh(px, py, data, cmap=cmap,vmin=0, vmax=0.01)

# add surface dye
# cs = ax['dyeM'].pcolormesh(px,py,surf_dye_alk_units,cmap=surf_cmap,vmin=vmin,vmax=vmax)
cbar = fig.colorbar(cs, ax=ax['dyeM'],orientation='horizontal',
                    location='bottom',pad=0.02,fraction=0.06)
cbar.ax.tick_params(labelsize=12,rotation=30)
cbar.outline.set_visible(False)
# format figure
ax['dyeM'].set_yticklabels([])
ax['dyeM'].set_xticklabels([])
pfun.dar(ax['dyeM'])
# ax[0].axis('off')
ax['dyeM'].set_title(r'Surface dye $>$ '+str(threshold)+'\n'+r'[mmol OH$^-$ m$^{-3}$]', fontsize=14)
# add injection location
ax['dyeM'].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
ax['dyeM'].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)
# for spine in ax['dyeM'].spines.values():
#     spine.set_visible(False)

# add surface alkalinity
# plot difference in surface alkalinity
diff = surf_alk_pert - surf_alk_base
cs = ax['alkM'].pcolormesh(px,py,diff,cmap=surf_cmap,vmin=vmin,vmax=vmax)
cbar = fig.colorbar(cs, ax=ax['alkM'],orientation='horizontal',
                    location='bottom',pad=0.02,fraction=0.06)
cbar.ax.tick_params(labelsize=12,rotation=30)
cbar.outline.set_visible(False)
# format figure
ax['alkM'].set_yticklabels([])
ax['alkM'].set_xticklabels([])
# ax[0].axis('off')
pfun.dar(ax['alkM'])
ax['alkM'].set_title(r'Surface $\Delta$ Alk$_{T}$'+'\n'+r'[meq m$^{-3}$]', fontsize=14)
# add injection location
ax['alkM'].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
ax['alkM'].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)
# for spine in ax['alkM'].spines.values():
#     spine.set_visible(False)

# add difference in surface CO2 flux
diff = CO2_flux_actual_pert - CO2_flux_actual_base
cs = ax['CO2M'].pcolormesh(px,py,diff,cmap=surf_cmap,vmin=-0.005,vmax=0.005)
cbar = fig.colorbar(cs, ax=ax['CO2M'],orientation='horizontal',
                    location='bottom',pad=0.02,fraction=0.06)
cbar.ax.tick_params(labelsize=12,rotation=30)
cbar.outline.set_visible(False)
# format figure
ax['CO2M'].set_yticklabels([])
ax['CO2M'].set_xticklabels([])
# ax[0].axis('off')
pfun.dar(ax['CO2M'])
ax['CO2M'].set_title(r'$\Delta$ CO$_2$ flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
# add injection location
ax['CO2M'].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
ax['CO2M'].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

# draw box around analysis region
# xmin = -124 #-126
# xmax = -122
# ymin = 46.7 #45.5
# ymax = 49 #50.5
xmin = -126
xmax = -122
ymin = 45.5
ymax = 50.5
# # draw box around study domain
# bordercolor = 'black'
# ax['alkM'].add_patch(Rectangle((xmin, ymin), xmax-xmin, ymax-ymin,
#              edgecolor = bordercolor, facecolor='none', lw=1.5))
# ax['dicM'].add_patch(Rectangle((xmin, ymin), xmax-xmin, ymax-ymin,
#              edgecolor = bordercolor, facecolor='none', lw=1.5))


# plot time series ------------------------------------

# delta alkalinity
ax['alkT'].plot(time_combined,delta_Alk_combined, linewidth=3, color='royalblue',alpha=0.3, label='alk')
ax['alkT'].plot(time_combined,total_dye_combined, linewidth=1.5, color='royalblue',linestyle='--', label='dye')
ax['alkT'].scatter(time_combined[t_index],delta_Alk_combined[t_index],color='black',marker='o', s=50)
ax['alkT'].set_ylabel(r'$\Delta$ Alk [kmol]', fontsize=14)
ax['alkT'].set_xticklabels([])
ax['alkT'].tick_params(axis='x', which='both', labelbottom=False) 
ax['alkT'].tick_params(axis='y', labelsize=12, rotation=30) 
ax['alkT'].grid(True, color='silver', linestyle=':')
ax['alkT'].set_ylim([0,10000])
ax['alkT'].legend(fontsize=10,ncol=2, loc='upper left', frameon=False)

# delta DIC
ax['dicT'].plot(time_combined,delta_DIC_combined, linewidth=3, color='royalblue',alpha=0.3, label=r'$\Delta$ DIC')
ax['dicT'].plot(time_combined,np.cumsum(CO2_flux_combined), linewidth=1.5, color='royalblue',linestyle='--', label=r'CO$_2$ flux')
ax['dicT'].scatter(time_combined[t_index],delta_DIC_combined[t_index],color='black',marker='o', s=50)
ax['dicT'].set_ylabel(r'CO$_2$ Uptake' + '\n[kmol]', fontsize=14)
ax['dicT'].set_xticklabels([])
ax['dicT'].tick_params(axis='x', which='both', labelbottom=False) 
ax['dicT'].tick_params(axis='y', labelsize=12, rotation=30)  
ax['dicT'].grid(True, color='silver', linestyle=':')
ax['dicT'].set_ylim([0,2000])
ax['dicT'].legend(fontsize=10,loc='lower right', frameon=False)


alk_hold = delta_Alk_combined.copy()
alk_hold[29:] = delta_Alk_combined[29]
ax['effT'].plot(time_combined,delta_DIC_combined/alk_hold,
                linewidth=3, color='royalblue',alpha=0.3,label=r'$\Delta$ DIC')
ax['effT'].plot(time_combined,np.cumsum(CO2_flux_combined)/alk_hold,
                linewidth=1.5, color='royalblue',linestyle='--',label=r'CO$_2$ flux')
ax['effT'].scatter(time_combined[t_index],
                   delta_DIC_combined[t_index]/alk_hold[t_index],
                   color='black',marker='o', s=50)
ax['effT'].legend(fontsize=10,loc='lower right', frameon=False)

ax['effT'].set_ylabel(r'$\eta =$' + '\n' +  r'$\Delta$ Carbon / $\Delta$ Alk', fontsize=14)
ax['effT'].grid(True, color='silver', linestyle=':')
ax['effT'].tick_params(axis='both', labelsize=12, rotation=30)  
loc = mdates.MonthLocator(interval=1)
ax['effT'].xaxis.set_major_locator(loc)
ax['effT'].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax['effT'].set_xlim([np.datetime64('2020-06-01'), np.datetime64('2020-10-31')])
# ax['effT'].set_ylim([0,1])