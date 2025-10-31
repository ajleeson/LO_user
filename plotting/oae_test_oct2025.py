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
fig, ax = plt.subplots(1,3,figsize = (15,6), sharex = True, sharey = True)

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
# format figure
ax[0].set_xlim([xmin,xmax])
ax[0].set_ylim([ymin,ymax])
pfun.dar(ax[0])

# plot with module
cs = ax[1].pcolormesh(px,py,surf_alk_withmodule, vmin=vmin, vmax=vmax, cmap=cmap)#cmocean.cm.balance_r)
# cs = ax.pcolormesh(px,py,ARAG)
cbar = fig.colorbar(cs, location='left')
cbar.ax.tick_params(labelsize=12)#,length=10, width=2)
cbar.outline.set_visible(False)
ax[1].set_title('With Module', fontsize=14)
# format figure
ax[1].set_xlim([xmin,xmax])
ax[1].set_ylim([ymin,ymax])
pfun.dar(ax[1])

print(ds['alkalinity'])

# plot difference
difference = surf_alk_nomodule - surf_alk_withmodule
cs = ax[2].pcolormesh(px,py,difference, vmin=-5, vmax=5, cmap=cmocean.cm.balance_r)
# cs = ax.pcolormesh(px,py,ARAG)
cbar = fig.colorbar(cs, location='left')
cbar.ax.tick_params(labelsize=12)#,length=10, width=2)
cbar.outline.set_visible(False)
ax[2].set_title('No Module minus With Module', fontsize=14)
# format figure
ax[2].set_xlim([xmin,xmax])
ax[2].set_ylim([ymin,ymax])
pfun.dar(ax[2])

# # loop through and plot both conditions
# for i,year in enumerate(years):
                
#     v  = val_dict[year] 

#     # plot natural condition
#     if i == 0:
#         cs = ax[i].pcolormesh(px,py,v, vmin=vmin, vmax=vmax, cmap=cmap)#cmocean.cm.balance_r)
#         cbar = fig.colorbar(cs, location='left')
#         cbar.ax.tick_params(labelsize=32)#,length=10, width=2)
#         cbar.outline.set_visible(False)
#         ax[i].set_title('Natural', fontsize=38)


#     if i == 1:
#         diff = (val_dict[years[i]] - val_dict[years[i-1]])
#         mindiff = np.nanmin(diff)
#         maxdiff = np.nanmax(diff)
#         # make sure colorbar axis contains zero
#         if mindiff > 0 and maxdiff > 0:
#             mindiff = maxdiff*-1.01
#         if mindiff < 0 and maxdiff < 0:
#             maxdiff = mindiff*-1.01
#         # don't let colorbar axis scale get too large
#         if maxdiff > vmax:
#             maxdiff = vmax
#         # make sure the colorbar is always centered about zero
#         cmap = cmocean.tools.crop(cmocean.cm.balance_r, mindiff, maxdiff, 0)
#         cs = ax[i].pcolormesh(px,py,diff, vmin=mindiff, vmax=maxdiff, cmap=cmap)
#         ax[i].set_title('Anthropogenic - Natural', fontsize=38)
#         cbar = fig.colorbar(cs, location='right')
#         cbar.ax.tick_params(labelsize=32)#,length=10, width=2)
#         cbar.outline.set_visible(False)

#     # format figure
#     ax[i].set_xlim([xmin,xmax])
#     ax[i].set_ylim([ymin,ymax])
#     ax[i].set_yticklabels([])
#     ax[i].set_xticklabels([])
#     ax[i].axis('off')
#     pfun.dar(ax[i])

# # add wwtp locations
# if WWTP_loc == True:
#     ax[1].scatter(lon_wwtps,lat_wwtps,color='none', edgecolors='k', linewidth=3, s=sizes_wwtps, label='WWTPs')
#     leg_szs = [100, 1000, 10000]
#     szs = [0.3*(leg_sz) for leg_sz in leg_szs]
#     l0 = plt.scatter([],[], s=szs[0], color='none', edgecolors='k', linewidth=3)
#     l1 = plt.scatter([],[], s=szs[1], color='none', edgecolors='k', linewidth=3)
#     l2 = plt.scatter([],[], s=szs[2], color='none', edgecolors='k', linewidth=3)
#     labels = ['< 100', '1,000', '10,000']
#     legend = ax[1].legend([l0, l1, l2], labels, fontsize = 18, markerfirst=False,
#         title='WWTP loading \n'+r' (kg N d$^{-1}$)',loc='lower right', labelspacing=1, borderpad=0.8)
#     plt.setp(legend.get_title(),fontsize=20)

# # add 10 km bar
# lat0 = 46.94
# lon0 = -123.05
# lat1 = lat0
# lon1 = -122.91825
# distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
# x_dist_km = round(distances_m[0]/1000)
# ax[0].plot([lon0,lon1],[lat0,lat1],color='k',linewidth=8)
# ax[0].text(lon0-0.04,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=24)

# # format figure
# ax[0].set_xlim([xmin,xmax])
# ax[0].set_ylim([ymin,ymax])
# ax[0].set_yticklabels([])
# ax[0].set_xticklabels([])
# ax[0].axis('off')
                                
# # Add colormap title
# plt.suptitle('2014 ' + start+' to '+end+' average ' + stext + ' ' + vn + ' ' + units,
#             fontsize=44, fontweight='bold', y=0.95)

# # Generate plot
# plt.tight_layout
# plt.subplots_adjust(left=0.05, right=0.95, top=0.85, wspace=0.02)
# plt.savefig(out_dir / (vn+'_'+stext+'_avg_2014A-N_'+ start + 'THRU'+end+'.png'))