"""

Compare surface alkalinity concentration
with and without the oae module

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import pandas as pd
import cmocean
import matplotlib.pylab as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
from PyCO2SYS import CO2SYS
import seawater as sw
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

vn = 'alkalinity'
his_file = 'ocean_his_0025.nc'

##############################################################
##                         MAKE MAP                         ##
##############################################################

# zoom into columbia river plume
xmin = -124.6
xmax = -123.6
ymin = 45.8
ymax = 46.8

# set axes range
vmin = 500
vmax = 2500
cmap = 'plasma'

# scale variable & get units
# scale =  pinfo.fac_dict[vn]
units = pinfo.units_dict[vn]

# Initialize figure
fig, ax = plt.subplots(1,2,figsize = (10.5,5), sharex = True, sharey = True)

# get data
fp = Ldir['roms_out'] / 'cas7_t1jxoae_x11bjx' / 'f2013.07.01' / his_file
G = zrfun.get_basic_info(fp, only_G=True)
ds = xr.open_dataset(fp)

fp = Ldir['roms_out'] / 'cas7_t1jxoae_x11ecb' / 'f2013.07.01' / his_file
ds_withmod = xr.open_dataset(fp)

# get lat and lon
lons = ds.coords['lon_rho'].values
lats = ds.coords['lat_rho'].values
px, py = pfun.get_plon_plat(lons,lats)

# get surface alkalinity
surf_alk_nomodule = ds['alkalinity'].values[0,-1,:,:] # (time, surface, all y, all x)
surf_alk_withmodule = ds_withmod['alkalinity'].values[0,-1,:,:] # (time, surface, all y, all x)

# # subtract ambient
# fn = Ldir['roms_out'] / 'cas7_t1jxoae_x11b' / 'f2013.07.01' / his_file
# ds_ambient = xr.open_dataset(fn)
# ambient_surf_alk = ds_ambient['alkalinity'].values[0,-1,:,:]

# surf_alk_nomodule = surf_alk_nomodule - ambient_surf_alk
# surf_alk_withmodule = surf_alk_withmodule - ambient_surf_alk

# # plot estimate of aragonite saturation
# v_dict = {}
# vn_in_list = ['temp', 'salt' , 'alkalinity', 'TIC']#['temp', 'salt' , 'rho', 'alkalinity', 'TIC']
# for cvn in vn_in_list:
#     L = ds[cvn][0,-1,:,:].values
#     v_dict[cvn] = L
# # ------------- the CO2SYS steps -------------------------
# # create pressure
# Ld = G['h']
# Lpres = sw.pres(Ld, lats)
# # get in situ temperature from potential temperature
# Ltemp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, Lpres)
# # convert from umol/L to umol/kg using in situ dentity
# Lalkalinity = v_dict['alkalinity']#1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
# Lalkalinity[Lalkalinity < 100] = np.nan
# LTIC = v_dict['TIC']#1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
# LTIC[LTIC < 100] = np.nan
# CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
#     Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
# # PH = CO2dict['pHout']
# # PH = zfun.fillit(PH.reshape((v_dict['salt'].shape)))
# ARAG = CO2dict['OmegaARout']
# ARAG = ARAG.reshape((v_dict['salt'].shape))

# plot without module
cs = ax[0].pcolormesh(px,py,surf_alk_nomodule, vmin=vmin, vmax=vmax, cmap=cmap)#cmocean.cm.balance_r)
# cs = ax.pcolormesh(px,py,ARAG)
cbar = fig.colorbar(cs, location='left')
cbar.ax.tick_params(labelsize=12)#,length=10, width=2)
cbar.outline.set_visible(False)
ax[0].set_title('No Module', fontsize=14)
# add injection site
ax[0].scatter(-124.08447696,46.24604386,
              marker='*',color='cyan',s=130)
# format figure
ax[0].set_xlim([xmin,xmax])
ax[0].set_ylim([ymin,ymax])
cbar.set_label(r'Surface Alkalinity [meq/m$^3$]', fontsize=14)
pfun.dar(ax[0])

# # plot with module
# cs = ax[1].pcolormesh(px,py,surf_alk_withmodule, vmin=vmin, vmax=vmax, cmap=cmap)#cmocean.cm.balance_r)
# # cs = ax.pcolormesh(px,py,ARAG)
# cbar = fig.colorbar(cs, location='left')
# cbar.ax.tick_params(labelsize=12)#,length=10, width=2)
# cbar.outline.set_visible(False)
# ax[1].set_title('With Module', fontsize=14)
# # format figure
# ax[1].set_xlim([xmin,xmax])
# ax[1].set_ylim([ymin,ymax])
# pfun.dar(ax[1])

print(ds['alkalinity'])

# plot difference
difference = surf_alk_nomodule - surf_alk_withmodule
cs = ax[1].pcolormesh(px,py,difference, vmin=-100, vmax=100, cmap=cmocean.cm.balance_r)
# cs = ax.pcolormesh(px,py,ARAG)
cbar = fig.colorbar(cs, location='left')
cbar.ax.tick_params(labelsize=12)#,length=10, width=2)
cbar.outline.set_visible(False)
ax[1].set_title('No Module minus With Module', fontsize=14)
# format figure
ax[1].set_xlim([xmin,xmax])
ax[1].set_ylim([ymin,ymax])
cbar.set_label(r'Difference in Surface Alkalinity [meq/m$^3$]', fontsize=14)
pfun.dar(ax[1])




plt.tight_layout()