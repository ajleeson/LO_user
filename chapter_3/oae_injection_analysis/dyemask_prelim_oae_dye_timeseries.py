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
date = '2020.09.20' 
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
##                  Get time series data                    ##
##############################################################

# SMALLER SUB-DOMAIN
# # get the first month that has continuous alkalinity and dye addition (cropped subdomain)
# fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'onemonimpulse_oae_deltas_superSUBDOMAIN_2020.06.01_2020.06.30.nc'
# ds_addition = xr.open_dataset(fp)
# # # crop to just first 30 days
# # ds_addition = ds_addition.isel(ocean_time=slice(0,30))
# # get the remaining time after the first addition
# fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'onemonimpulse_oae_deltas_superSUBDOMAIN_2020.07.01_2020.10.31.nc'
# ds_later = xr.open_dataset(fp)

# BIGGER SUB-DOMAIN
# get the first month that has continuous alkalinity and dye addition (cropped subdomain)
fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'oae_deltas_SUBDOMAIN_2020.06.01_2020.08.31.nc'
ds_addition = xr.open_dataset(fp)
# crop to just first 30 days
ds_addition = ds_addition.isel(ocean_time=slice(0,30))
# get the remaining time after the first addition
fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'onemonimpulse_oae_deltas_SUBDOMAIN_2020.07.01_2020.10.31.nc'
ds_later = xr.open_dataset(fp)

# combine the two datasets
ds_combined = xr.concat([ds_addition, ds_later], dim='ocean_time')

time_combined = ds_combined.ocean_time.values
delta_DIC_combined = ds_combined.delta_DIC.values # kmol
delta_Alk_combined = ds_combined.delta_Alk.values # kmol
total_dye_combined = ds_combined.total_dye.values # kmol
surf_dye_combined  = ds_combined.surf_dye.values  # kmol

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
    ax['dicM'] = fig.add_subplot(gridspec[0:6, 11:15])
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

# add mid-column DIC
# plot difference in mid-column alkalinity
diff = midD_dic_pert - midD_dic_base
cs = ax['dicM'].pcolormesh(px,py,diff,cmap=surf_cmap,vmin=vmin,vmax=vmax)
cbar = fig.colorbar(cs, ax=ax['dicM'],orientation='horizontal',
                    location='bottom',pad=0.02,fraction=0.06)
cbar.ax.tick_params(labelsize=12,rotation=30)
cbar.outline.set_visible(False)
# format figure
ax['dicM'].set_yticklabels([])
ax['dicM'].set_xticklabels([])
# ax[0].axis('off')
pfun.dar(ax['dicM'])
ax['dicM'].set_title(r'Mid-water column ($\sigma$=20)'+'\n'+r'$\Delta$ DIC [mmol m$^{-3}$]', fontsize=14)
# add injection location
ax['dicM'].scatter(inj_lon, inj_lat, color='none', edgecolor='pink',marker='o', s=100,linewidth=3)
ax['dicM'].scatter(inj_lon, inj_lat, color='none', edgecolor='crimson',marker='o', s=100,linewidth=2)

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
ax['dicT'].plot(time_combined,delta_DIC_combined, linewidth=3, color='royalblue',alpha=0.3)
ax['dicT'].scatter(time_combined[t_index],delta_DIC_combined[t_index],color='black',marker='o', s=50)
ax['dicT'].set_ylabel(r'$\Delta$ DIC [kmol]', fontsize=14)
ax['dicT'].set_xticklabels([])
ax['dicT'].tick_params(axis='x', which='both', labelbottom=False) 
ax['dicT'].tick_params(axis='y', labelsize=12, rotation=30)  
ax['dicT'].grid(True, color='silver', linestyle=':')
ax['dicT'].set_ylim([0,6000])

# efficiency
# # get cumulative alkalinity up until Jun 30, then hold cumsum constant
# totalalk = np.cumsum(delta_Alk_combined)
# totalalk29 = totalalk[29]                          # sum of indices 0..29 (first 30 points)
# totalalk_hold = np.concatenate([totalalk[:30],     # keep through index 29
#                           np.full(len(totalalk) - 30, totalalk29)])  # hold constant after
# ax['effT'].plot(time_combined,delta_DIC_combined/totalalk_hold,
#                 linewidth=3, color='royalblue',alpha=0.3)
# ax['effT'].scatter(time_combined[t_index],
#                    delta_DIC_combined[t_index]/totalalk_hold[t_index],
#                    color='black',marker='o', s=50)

# # get cumulative alkalinity up until Jun 30, then hold cumsum constant
# totalalk = np.cumsum(delta_Alk_combined)
# totalalk29 = totalalk[29]                          # sum of indices 0..29 (first 30 points)
# totalalk_hold = np.concatenate([totalalk[:30],     # keep through index 29
#                           np.full(len(totalalk) - 30, totalalk29)])  # hold constant after
# ax['effT'].plot(time_combined,np.cumsum(delta_DIC_combined)/totalalk_hold,
#                 linewidth=3, color='royalblue',alpha=0.3)
# ax['effT'].scatter(time_combined[t_index],
#                    np.cumsum(delta_DIC_combined)[t_index]/totalalk_hold[t_index],
#                    color='black',marker='o', s=50)

# ax['effT'].plot(time_combined,np.cumsum(delta_DIC_combined)/np.nansum(delta_Alk_combined[0:30]),
#                 linewidth=3, color='royalblue',alpha=0.3)
# ax['effT'].scatter(time_combined[t_index],
#                    np.cumsum(delta_DIC_combined)[t_index]/np.nansum(delta_Alk_combined[0:30]),
#                    color='black',marker='o', s=50)

alk_hold = delta_Alk_combined.copy()
alk_hold[29:] = delta_Alk_combined[29]
ax['effT'].plot(time_combined,delta_DIC_combined/alk_hold,
                linewidth=3, color='royalblue',alpha=0.3)
ax['effT'].scatter(time_combined[t_index],
                   delta_DIC_combined[t_index]/alk_hold[t_index],
                   color='black',marker='o', s=50)

# ax['effT'].plot(time_combined,delta_DIC_combined/delta_Alk_combined,
#                 linewidth=3, color='royalblue',alpha=0.3)
# ax['effT'].scatter(time_combined[t_index],
#                    delta_DIC_combined[t_index]/delta_Alk_combined[t_index],
#                    color='black',marker='o', s=50)

ax['effT'].set_ylabel(r'$\eta = \Delta$ DIC / $\Delta$ Alk', fontsize=14)
ax['effT'].grid(True, color='silver', linestyle=':')
ax['effT'].tick_params(axis='both', labelsize=12, rotation=30)  
loc = mdates.MonthLocator(interval=1)
ax['effT'].xaxis.set_major_locator(loc)
ax['effT'].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax['effT'].set_xlim([np.datetime64('2020-06-01'), np.datetime64('2020-10-31')])
# ax['effT'].set_ylim([0,1])