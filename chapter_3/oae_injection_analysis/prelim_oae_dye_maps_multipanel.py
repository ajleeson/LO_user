"""
Plot difference in CO2 flux between baseline and perturbation
Takes about 30 seconds to run
"""

# import things
import numpy as np
import xarray as xr
import csv
import matplotlib.pylab as plt
import matplotlib.colors as colors
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
from datetime import datetime
import gsw
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

dates = ['2020.06.07','2020.06.30','2020.07.31','2020.08.31']

# create dict of data
dict_data = {date: {} for date in dates}

# which  model runs to look at?
basline = 'cas7_t1_x11ab'
perturbation = 'cas7_t1dgeWB_x11abd'

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
##        Get CO2 flux in baseline and pertubation          ##
##############################################################

for y,date in enumerate(dates):

    print('Processing {}'.format(date))

    ds_base = xr.open_dataset(Ldir['roms_out'] / basline      / ('f' + date) / 'ocean_avg_0001.nc')
    ds_pert = xr.open_dataset(Ldir['roms_out'] / perturbation / ('f' + date) / 'ocean_avg_0001.nc')
            
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

    # Get surface delta alkalinity and surface dye

    surf_alk_base = ds_base['alkalinity'].values[0,-1,:,:]
    surf_alk_pert = ds_pert['alkalinity'].values[0,-1,:,:]

    # get surface dye
    surf_dye = ds_pert['dye_01'].values[0,-1,:,:]
    # convert units from kg/m3 to mmol/m3, to equate to alkaliinty
    surf_dye_alk_units = surf_dye / 1.7e-5 

    # add all values to df_dict
    dict_data[date]['alk_diff'] = surf_alk_pert - surf_alk_base
    dict_data[date]['CO2_flux_actual_base'] = CO2_flux_actual_base
    dict_data[date]['CO2_flux_diff'] = CO2_flux_actual_pert - CO2_flux_actual_base

# convert dict back into dict
dict_data_final = dict(dict_data)

print(dict_data_final)

###################################################################
##                         Plotting                              ##  
################################################################### 

# # Initialize figure
# plt.close('all')
# fig,ax = plt.subplots(3,4,figsize=(11,8.5), sharex=True, sharey=True) 

# for y,date in enumerate(dates):

#     # get grid info
#     lons = ds_base.coords['lon_rho'].values
#     lats = ds_base.coords['lat_rho'].values
#     px, py = pfun.get_plon_plat(lons,lats)

#     # set colormaps
#     dye_alk_cmap = cmc.devon_r
#     flux_cmap = cmc.vik
#     diff_cmap = cmc.lajolla_r
#     flux_cmap.set_bad(color='gray')
#     diff_cmap.set_bad(color='gray')
#     dye_alk_cmap.set_bad(color='gray')

#     # # Southern Salish Sea
#     # ymin = 46.8
#     # ymax = 48.9
#     # xmin = -125
#     # xmax = -122.0

#     # # Whidbey Basin
#     # ymin = 47.75
#     # ymax = 48.5
#     # xmin = -123.1
#     # xmax = -122.0

#     # Whidbey and SJdF
#     ymin = 47.75
#     ymax = 48.5
#     xmin = -124
#     xmax = -122.0

#     # injection location
#     inj_lon = -122.5674
#     inj_lat = 48.1956

#     # plot difference in surface alkalinity
#     cs = ax[0,y].pcolormesh(px,py,dict_data_final[date]['alk_diff'],cmap=dye_alk_cmap,norm=colors.LogNorm(vmin=1e-3,vmax=1e-1))
#     if y == 3:
#         cax = inset_axes(ax[0,y],width="6%",height="100%",loc="lower left",
#         bbox_to_anchor=(1.02, 0., 1, 1),bbox_transform=ax[0,y].transAxes, borderpad=0)
#         cbar = fig.colorbar(cs, cax=cax)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.outline.set_visible(False)
#     if y == 0:
#         ax[0,y].set_ylabel(r'Surface $\Delta$ Alk$_{T}$ [meq m$^{-3}$]', fontsize=14)
#     # format figure
#     ax[0,y].set_xlim([xmin,xmax])
#     ax[0,y].set_ylim([ymin,ymax])
#     ax[0,y].set_yticklabels([])
#     ax[0,y].set_xticklabels([])
#     # ax[0].axis('off')
#     pfun.dar(ax[0,y])
#     # ax[0,y].set_title(r'Surface $\Delta$ Alk$_{T}$ [meq m$^{-3}$]', fontsize=14)
#     # add injection location
#     ax[0,y].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
#     ax[0,y].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

#     # # plot surface dye
#     # cs = ax[1].pcolormesh(px,py,surf_dye_alk_units,cmap=dye_alk_cmap,norm=colors.LogNorm(vmin=1e-3,vmax=1e-1))
#     # cbar = fig.colorbar(cs)
#     # cbar.ax.tick_params(labelsize=14)
#     # cbar.outline.set_visible(False)
#     # # format figure
#     # ax[1].set_xlim([xmin,xmax])
#     # ax[1].set_ylim([ymin,ymax])
#     # ax[1].set_yticklabels([])
#     # ax[1].set_xticklabels([])
#     # # ax[1].axis('off')
#     # pfun.dar(ax[1,y])
#     # ax[1].set_title(r'Surface dye [mmol OH$^-$ m$^{-3}$]', fontsize=14)
#     # # add injection location
#     # ax[1].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
#     # ax[1].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

#     # plot basline CO2 flux values
#     cs = ax[2,y].pcolormesh(px,py,dict_data_final[date]['CO2_flux_actual_base'],vmin=-40,vmax=40,cmap=flux_cmap)
#     if y == 3:
#         cax = inset_axes(ax[2,y],width="6%",height="100%",loc="lower left",
#         bbox_to_anchor=(1.02, 0., 1, 1),bbox_transform=ax[2,y].transAxes, borderpad=0)
#         cbar = fig.colorbar(cs, cax=cax)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.outline.set_visible(False)
#     if y == 0:
#         ax[2,y].set_ylabel(r'Baseline CO$_2$ flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
#     # format figure
#     ax[2,y].set_xlim([xmin,xmax])
#     ax[2,y].set_ylim([ymin,ymax])
#     ax[2,y].set_yticklabels([])
#     ax[2,y].set_xticklabels([])
#     # ax[2].axis('off')
#     pfun.dar(ax[2,y])
#     # ax[2,y].set_title(r'Baseline CO$_2$ flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
#     # add injection location
#     ax[2,y].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
#     ax[2,y].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

#     # plot difference in CO2 flux values
#     cs = ax[1,y].pcolormesh(px,py,dict_data_final[date]['CO2_flux_diff'],cmap=diff_cmap,vmin=0,vmax=0.008)#,norm=colors.LogNorm(vmin=1e-4,vmax=1e-1))
#     if y == 3:
#         cax = inset_axes(ax[1,y],width="6%",height="100%",loc="lower left",
#         bbox_to_anchor=(1.02, 0., 1, 1),bbox_transform=ax[1,y].transAxes, borderpad=0)
#         cbar = fig.colorbar(cs, cax=cax)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.outline.set_visible(False)
#         cbar.outline.set_visible(False)
#     if y == 0:
#         ax[1,y].set_ylabel(r'$\Delta$ CO$_2$ flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
#     # format figure
#     ax[1,y].set_xlim([xmin,xmax])
#     ax[1,y].set_ylim([ymin,ymax])
#     ax[1,y].set_yticklabels([])
#     ax[1,y].set_xticklabels([])
#     # ax[3].axis('off')
#     pfun.dar(ax[1,y])
#     # ax[1,y].set_title(r'$\Delta$ CO$_2$ flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
#     # add injection location
#     ax[1,y].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
#     ax[1,y].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

# Initialize figure
plt.close('all')
fig,ax = plt.subplots(4,3,figsize=(8.5,8.5), sharex=True, sharey=True) 

for y,date in enumerate(dates):

    # get grid info
    lons = ds_base.coords['lon_rho'].values
    lats = ds_base.coords['lat_rho'].values
    px, py = pfun.get_plon_plat(lons,lats)

    # set colormaps
    dye_alk_cmap = cmc.devon_r
    flux_cmap = cmc.roma_r #cmc.vik
    diff_cmap = cmc.batlowW_r #cmc.tokyo_r #cmc.lajolla_r #
    flux_cmap.set_bad(color   ='gray')
    diff_cmap.set_bad(color   ='gray')
    dye_alk_cmap.set_bad(color='gray')

    # # Southern Salish Sea
    # ymin = 46.8
    # ymax = 48.9
    # xmin = -125
    # xmax = -122.0

    # # Whidbey Basin
    # ymin = 47.75
    # ymax = 48.5
    # xmin = -123.1
    # xmax = -122.0

    # Whidbey Basin 2.0
    ymin = 47.75
    ymax = 48.5
    xmin = -123.5
    xmax = -122.0

    # # Whidbey and SJdF
    # ymin = 47.75
    # ymax = 48.5
    # xmin = -124
    # xmax = -122.0

    # # Whidbey and SJdF 2.0
    # ymin = 47.3
    # ymax = 48.5
    # xmin = -124.8
    # xmax = -122.0

    # injection location
    inj_lon = -122.5674
    inj_lat = 48.1956

    # plot difference in surface alkalinity
    vmin = 1e-3
    clipped = np.where(dict_data_final[date]['alk_diff'] <= 0, vmin, dict_data_final[date]['alk_diff'])
    cs = ax[y,0].pcolormesh(px,py,clipped,cmap=dye_alk_cmap,norm=colors.LogNorm(vmin=vmin,vmax=1e1))
    # cs = ax[y,0].pcolormesh(px,py,dict_data_final[date]['alk_diff'],vmin=0,vmax=1,cmap=dye_alk_cmap)
    if y == 3:
        cax = inset_axes(ax[y,0],width="100%",height="9%",loc="lower left",
        bbox_to_anchor=(0., -0.12, 1, 1),bbox_transform=ax[y,0].transAxes, borderpad=0)
        cbar = fig.colorbar(cs, cax=cax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=14, rotation=30)
        cbar.outline.set_visible(False)
    if y == 0:
        ax[y,0].set_title(r'Surface $\Delta$ Alk$_{T}$'+'\n'+r'[meq m$^{-3}$]', fontsize=14)
    # format figure
    ax[y,0].set_xlim([xmin,xmax])
    ax[y,0].set_ylim([ymin,ymax])
    ax[y,0].set_yticklabels([])
    ax[y,0].set_xticklabels([])
    # ax[0].axis('off')
    pfun.dar(ax[y,0])
    # ax[0,y].set_title(r'Surface $\Delta$ Alk$_{T}$ [meq m$^{-3}$]', fontsize=14)
    # add injection location
    ax[y,0].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o',    s=100,linewidth=3, zorder=5)
    ax[y,0].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2, zorder=5)
    # add date label
    ax[y,0].set_ylabel(date, fontsize=14)
    pfun.add_coast(ax[y,0], color='black')

    # # plot surface dye
    # cs = ax[1].pcolormesh(px,py,surf_dye_alk_units,cmap=dye_alk_cmap,norm=colors.LogNorm(vmin=1e-3,vmax=1e-1))
    # cbar = fig.colorbar(cs)
    # cbar.ax.tick_params(labelsize=14)
    # cbar.outline.set_visible(False)
    # # format figure
    # ax[1].set_xlim([xmin,xmax])
    # ax[1].set_ylim([ymin,ymax])
    # ax[1].set_yticklabels([])
    # ax[1].set_xticklabels([])
    # # ax[1].axis('off')
    # pfun.dar(ax[1,y])
    # ax[1].set_title(r'Surface dye [mmol OH$^-$ m$^{-3}$]', fontsize=14)
    # # add injection location
    # ax[1].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
    # ax[1].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

    # plot basline CO2 flux values
    cs = ax[y,2].pcolormesh(px,py,dict_data_final[date]['CO2_flux_actual_base'],vmin=-50,vmax=50,cmap=flux_cmap)
    if y == 3:
        cax = inset_axes(ax[y,2],width="100%",height="9%",loc="lower left",
        bbox_to_anchor=(0., -0.12, 1, 1),bbox_transform=ax[y,2].transAxes, borderpad=0)
        cbar = fig.colorbar(cs, cax=cax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=14, rotation=30)
        cbar.outline.set_visible(False)
    if y == 0:
        ax[y,2].set_title(r'Baseline CO$_2$ flux'+'\n'+r'[mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
    # format figure
    ax[y,2].set_xlim([xmin,xmax])
    ax[y,2].set_ylim([ymin,ymax])
    ax[y,2].set_yticklabels([])
    ax[y,2].set_xticklabels([])
    # ax[2].axis('off')
    pfun.dar(ax[y,2])
    # ax[2,y].set_title(r'Baseline CO$_2$ flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
    # add injection location
    ax[y,2].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o',    s=100,linewidth=3, zorder=5)
    ax[y,2].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2, zorder=5)
    pfun.add_coast(ax[y,2], color='black')

    # plot difference in CO2 flux values
    vmin = 1e-4
    clipped = np.where(dict_data_final[date]['CO2_flux_diff'] <= 0, vmin, dict_data_final[date]['CO2_flux_diff'])
    cs = ax[y,1].pcolormesh(px,py,clipped,cmap=diff_cmap,norm=colors.LogNorm(vmin=vmin,vmax=1e-2)) # ,vmin=0,vmax=0.008)#
    # cs = ax[y,1].pcolormesh(px,py,dict_data_final[date]['CO2_flux_diff'],cmap=diff_cmap,vmin=0,vmax=0.008)#
    if y == 3:
        cax = inset_axes(ax[y,1],width="100%",height="9%",loc="lower left",
        bbox_to_anchor=(0., -0.12, 1, 1),bbox_transform=ax[y,1].transAxes, borderpad=0)
        cbar = fig.colorbar(cs, cax=cax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=14, rotation=30)
        cbar.outline.set_visible(False)
        cbar.outline.set_visible(False)
    if y == 0:
        ax[y,1].set_title(r'$\Delta$ CO$_2$ flux'+'\n'+r'mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
    # format figure
    ax[y,1].set_xlim([xmin,xmax])
    ax[y,1].set_ylim([ymin,ymax])
    ax[y,1].set_yticklabels([])
    ax[y,1].set_xticklabels([])
    # ax[3].axis('off')
    pfun.dar(ax[y,1])
    # ax[1,y].set_title(r'$\Delta$ CO$_2$ flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
    # add injection location
    ax[y,1].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o',    s=100,linewidth=3, zorder=5)
    ax[y,1].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2, zorder=5)
    pfun.add_coast(ax[y,1], color='black')

# plt.adjust_subplots(wspace=0.05, hspace=0.05)