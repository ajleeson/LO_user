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

# get the first month that has continuous alkalinity and dye addition (cropped subdomain)
fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'oae_deltas_SUBDOMAIN_2020.06.01_2020.08.31.nc'
ds_addition = xr.open_dataset(fp)
# crop to just first 30 days
ds_addition = ds_addition.isel(ocean_time=slice(0,30))

# get the remaining time after the first addition
fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'onemonimpulse_oae_deltas_SUBDOMAIN_2020.07.01_2020.08.06.nc'
ds_later = xr.open_dataset(fp)

# combine the two datasets
ds_combined = xr.concat([ds_addition, ds_later], dim='ocean_time')

time_combined = ds_combined.ocean_time.values
delta_DIC_combined = ds_combined.delta_DIC.values # kmol
delta_Alk_combined = ds_combined.delta_Alk.values # kmol
total_dye_combined = ds_combined.total_dye.values # kmol
surf_dye_combined  = ds_combined.surf_dye.values  # kmol

# -------------------------

# # get data
# fp = Ldir['LOo'] / 'chapter_3' / 'data' / 'oae_deltas_SUBDOMAIN_2020.06.01_2020.08.31.nc'
# ds_cropped = xr.open_dataset(fp)

# fp = Ldir['LOo'] / 'chapter_3' / 'data' / '3moncont_oae_deltas_2020.06.01_2020.08.31.nc'
# ds_fulldomain = xr.open_dataset(fp)

# time_crop = ds_cropped.ocean_time.values
# delta_DIC_crop = ds_cropped.delta_DIC.values # kmol
# delta_Alk_crop = ds_cropped.delta_Alk.values # kmol
# total_dye_crop = ds_cropped.total_dye.values # kmol
# surf_dye_crop  = ds_cropped.surf_dye.values  # kmol

# time_full = ds_fulldomain.ocean_time.values
# delta_DIC_full = ds_fulldomain.delta_DIC.values # kmol
# delta_Alk_full = ds_fulldomain.delta_Alk.values # kmol
# total_dye_full = ds_fulldomain.total_dye.values # kmol
# surf_dye_full  = ds_fulldomain.surf_dye.values  # kmol

###################################################################
##                         Plotting                              ##  
################################################################### 

# Initialize figure
plt.close('all')
fig,ax = plt.subplots(2,2,figsize=(9,6), sharex=True) 
ax = ax.ravel()

# delta DIC
# ax[0].plot(delta_DIC_full)
# ax[0].plot(delta_DIC_crop)
ax[0].plot(delta_DIC_combined)
ax[0].set_title(r'$\Delta$ DIC [kmol]', fontsize=14)

# delta alkalinity
# ax[1].plot(delta_Alk_full)
# ax[1].plot(delta_Alk_crop)
# ax[1].plot(total_dye_crop, linestyle='--')
ax[1].plot(delta_Alk_combined)
ax[1].plot(total_dye_combined, linestyle='--')
ax[1].set_title(r'$\Delta$ Alk [kmol]', fontsize=14)

# dye
# ax[2].plot(total_dye[0:n])
# ax[2].plot(surf_dye[0:n])
# ax[2].plot(surf_dye_crop/total_dye_crop)
ax[2].plot(surf_dye_combined/total_dye_combined)
ax[2].set_title(r'Surface dye / total dye', fontsize=14)

# efficiency
# ax[3].plot(delta_DIC_full/delta_Alk_full)
# ax[3].plot(delta_DIC_crop/delta_Alk_crop)
# ax[3].plot(np.cumsum(delta_DIC_crop)/np.cumsum(delta_Alk_crop))
ax[3].plot(np.cumsum(delta_DIC_combined)/np.cumsum(delta_Alk_combined))
ax[3].set_title(r'$\Delta$ DIC / $\Delta$ Alk', fontsize=14)