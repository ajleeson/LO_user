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
##                        Get data                          ##
##############################################################

fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'oae_deltas_2020.06.01_2020.08.31.nc'
ds = xr.open_dataset(fp)

time = ds.ocean_time.values
delta_DIC = ds.delta_DIC.values # koml

###################################################################
##                         Plotting                              ##  
################################################################### 

# Initialize figure
plt.close('all')
fig,ax = plt.subplots(2,2,figsize=(9,6), sharex=True) 
ax = ax.ravel()

# delta DIC
ax[0].plot(delta_DIC)

# # set colormaps
# dye_alk_cmap = cmc.devon_r
# flux_cmap = cmc.vik
# diff_cmap = cmc.lajolla_r
# flux_cmap.set_bad(color='gray')
# diff_cmap.set_bad(color='gray')
# dye_alk_cmap.set_bad(color='gray')

# # # Southern Salish Sea
# # ymin = 46.8
# # ymax = 48.9
# # xmin = -125
# # xmax = -122.0

# # Whidbey Basin
# ymin = 47.75
# ymax = 48.5
# xmin = -123.1
# xmax = -122.0

# # injection location
# inj_lon = -122.5674
# inj_lat = 48.1956

# # plot difference in surface alkalinity
# diff = surf_alk_pert - surf_alk_base
# cs = ax[0].pcolormesh(px,py,diff,cmap=dye_alk_cmap,norm=colors.LogNorm(vmin=1e-2,vmax=1e1))
# cbar = fig.colorbar(cs)
# cbar.ax.tick_params(labelsize=14)
# cbar.outline.set_visible(False)
# # format figure
# ax[0].set_xlim([xmin,xmax])
# ax[0].set_ylim([ymin,ymax])
# ax[0].set_yticklabels([])
# ax[0].set_xticklabels([])
# # ax[0].axis('off')
# pfun.dar(ax[0])
# ax[0].set_title(r'Surface $\Delta$ Alk$_{T}$ [meq m$^{-3}$]', fontsize=14)
# # add injection location
# ax[0].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
# ax[0].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

# # plot surface dye
# cs = ax[1].pcolormesh(px,py,surf_dye_alk_units,cmap=dye_alk_cmap,norm=colors.LogNorm(vmin=1e-2,vmax=1e1))
# cbar = fig.colorbar(cs)
# cbar.ax.tick_params(labelsize=14)
# cbar.outline.set_visible(False)
# # format figure
# ax[1].set_xlim([xmin,xmax])
# ax[1].set_ylim([ymin,ymax])
# ax[1].set_yticklabels([])
# ax[1].set_xticklabels([])
# # ax[1].axis('off')
# pfun.dar(ax[1])
# ax[1].set_title(r'Surface dye [mmol OH$^-$ m$^{-3}$]', fontsize=14)
# # add injection location
# ax[1].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
# ax[1].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

# # plot basline CO2 flux values
# cs = ax[2].pcolormesh(px,py,CO2_flux_actual_base,vmin=-15,vmax=15,cmap=flux_cmap)
# cbar = fig.colorbar(cs)
# cbar.ax.tick_params(labelsize=14)
# cbar.outline.set_visible(False)
# # format figure
# ax[2].set_xlim([xmin,xmax])
# ax[2].set_ylim([ymin,ymax])
# ax[2].set_yticklabels([])
# ax[2].set_xticklabels([])
# # ax[2].axis('off')
# pfun.dar(ax[2])
# ax[2].set_title(r'Baseline CO$_2$ flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
# # add injection location
# ax[2].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
# ax[2].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

# # plot difference in CO2 flux values
# diff = CO2_flux_actual_pert - CO2_flux_actual_base
# cs = ax[3].pcolormesh(px,py,diff,cmap=diff_cmap,vmin=0,vmax=0.03)#,norm=colors.LogNorm(vmin=1e-4,vmax=1e-1))
# cbar = fig.colorbar(cs)
# cbar.ax.tick_params(labelsize=14)
# cbar.outline.set_visible(False)
# # format figure
# ax[3].set_xlim([xmin,xmax])
# ax[3].set_ylim([ymin,ymax])
# ax[3].set_yticklabels([])
# ax[3].set_xticklabels([])
# # ax[3].axis('off')
# pfun.dar(ax[3])
# ax[3].set_title(r'$\Delta$ CO$_2$ flux [mmol m$^{-2}$ d$^{-1}$]', fontsize=14)
# # add injection location
# ax[3].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
# ax[3].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

# # add 10 km bar
# lat0 = 47.8
# lon0 = -122.28
# lat1 = lat0
# lon1 = -122.08
# distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
# x_dist_km = round(distances_m[0]/1000)
# ax[0].plot([lon0,lon1],[lat0,lat1],color='white',linewidth=2)
# ax[0].text((lon0+lon1)/2,lat0+0.03,'{} km'.format(x_dist_km),color='white',
#         horizontalalignment='center', fontsize=12)

# # letter labels
# ax[0].text(0.95, 0.95, '(a)', transform=ax[0].transAxes, color='white',
#            fontsize=15, fontweight='bold', va='top', ha='right')
# ax[1].text(0.95, 0.95, '(b)', transform=ax[1].transAxes, color='white',
#            fontsize=15, fontweight='bold', va='top', ha='right')
# ax[2].text(0.95, 0.95, '(c)', transform=ax[2].transAxes, color='white',
#            fontsize=15, fontweight='bold', va='top', ha='right')
# ax[3].text(0.95, 0.95, '(d)', transform=ax[3].transAxes, color='white',
#            fontsize=15, fontweight='bold', va='top', ha='right')

# # Generate plot
# plt.tight_layout
# plt.suptitle('OAE daily averages on {}\n(after 3 months of daily alkalinity release)'.format(date), fontsize=14, fontweight='bold')