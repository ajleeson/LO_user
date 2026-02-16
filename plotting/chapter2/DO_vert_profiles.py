"""
Compare vertical profiles of DO
in Penn Cove and Hood Canal
"""

# import things
import numpy as np
import xarray as xr
import pandas as pd
import cmocean
import matplotlib.pylab as plt

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
print(pth)
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

plt.close('all')

##############################################################
##                       USER INPUTS                        ##
##############################################################

# Show WWTP locations?
WWTP_loc = True

# years = ['2015','2016','2017','2018','2019','2020']
years = ['2017']

# which  model run to look at?
gtagexes = ['cas7_t1noDIN_x11ab','cas7_t1_x11ab']

# hypoxic season
# start = '09.01'
# end = '10.31'

months = ['January','April','July','October']
colors = ['cornflowerblue','hotpink','yellowgreen','orange']
colors_dark = ['navy','deeppink','darkgreen','orangered']

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# open dataset for every year, and add to dictionary, with year as key

# open datasets
# initialize empty dictionary
ds_loading_dict_hc = {}
ds_loading_dict_pc = {}
ds_noloading_dict_hc = {}
ds_noloading_dict_pc = {}
for year in years:
    for month in months:
        if month == 'January':
            start = '01.01'
            end   = '01.31'
        elif month == 'April':
            start = '04.01'
            end   = '04.30'
        elif month == 'July':
            start = '07.01'
            end   = '07.31'
        elif month == 'October':
            start = '10.01'
            end   = '10.31'
        # open datasets
        ds_loading_pc = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / 'moor_extractions' / 'cas7_t1_x11ab' / ('penncove_2017.' + start + '_2017.' + end + '.nc'))
        ds_loading_hc = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / 'moor_extractions' / 'cas7_t1_x11ab' / ('hoodcanal_2017.' + start + '_2017.' + end + '.nc'))
        ds_noloading_pc = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / 'moor_extractions' / 'cas7_t1noDIN_x11ab' / ('penncove_2017.' + start + '_2017.' + end + '.nc'))
        ds_noloading_hc = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / 'moor_extractions' / 'cas7_t1noDIN_x11ab' / ('hoodcanal_2017.' + start + '_2017.' + end + '.nc'))
        # add ds to dictionary
        ds_loading_dict_hc[year+month] = ds_loading_hc
        ds_loading_dict_pc[year+month] = ds_loading_pc
        ds_noloading_dict_hc[year+month] = ds_noloading_hc
        ds_noloading_dict_pc[year+month] = ds_noloading_pc

# # get climatology of bottom DO for all seven years
# botDO_loading_arrays = [ds['DO_bot'].data for ds in ds_loading_dict.values()]  # [365, y, x]
# botDO_noloading_arrays = [ds['DO_bot'].data for ds in ds_noloading_dict.values()]  # [365, y, x]
# # stack the arrays
# botDO_loading_stacked = np.stack(botDO_loading_arrays, axis=0)  # [7, 365, y, x]
# botDO_noloading_stacked = np.stack(botDO_noloading_arrays, axis=0)  # [7, 365, y, x]
# # average over all seven years
# bottDO_clim_loading = botDO_loading_stacked.mean(axis=0)  # [365, y, x]
# bottDO_clim_noloading = botDO_noloading_stacked.mean(axis=0)  # [365, y, x]


##############################################################
##                    AVERAGE BOTTOM DO                     ##
##############################################################

# Initialize figure
fig, axes = plt.subplots(2,2,figsize = (6,7), sharey='row', sharex='col')
ax = axes.ravel()

for year in years:
    for i,month in enumerate(months):
        ds_pc_loading = ds_loading_dict_pc[year+month]
        ds_hc_loading = ds_loading_dict_hc[year+month]
        ds_pc_noloading = ds_noloading_dict_pc[year+month]
        ds_hc_noloading = ds_noloading_dict_hc[year+month]

        # a single day
        # # plot noloading
        # day = 1
        # ax[0].plot(ds_pc_noloading['oxygen'][day,:]*32/1000,ds_pc_noloading['z_rho'][day,:],color='hotpink',linewidth=3)
        # ax[1].plot(ds_hc_noloading['oxygen'][day,:]*32/1000,ds_hc_noloading['z_rho'][day,:],color='hotpink',linewidth=3)
        # # plot loading
        # ax[0].plot(ds_pc_loading['oxygen'][day,:]*32/1000,ds_pc_loading['z_rho'][day,:],color='navy',linewidth=3)
        # ax[1].plot(ds_hc_loading['oxygen'][day,:]*32/1000,ds_hc_loading['z_rho'][day,:],color='navy',linewidth=3)


        # average profile over hypoxic season
        do_pc_noloading = np.nanmean(ds_pc_noloading['oxygen'],axis=0)*32/1000
        zrho_pc_noloading = np.nanmean(ds_pc_noloading['z_rho'],axis=0)
        do_hc_noloading = np.nanmean(ds_hc_noloading['oxygen'],axis=0)*32/1000
        zrho_hc_noloading = np.nanmean(ds_hc_noloading['z_rho'],axis=0)

        do_pc_loading = np.nanmean(ds_pc_loading['oxygen'],axis=0)*32/1000
        zrho_pc_loading = np.nanmean(ds_pc_loading['z_rho'],axis=0)
        do_hc_loading = np.nanmean(ds_hc_loading['oxygen'],axis=0)*32/1000
        zrho_hc_loading = np.nanmean(ds_hc_loading['z_rho'],axis=0)

        # plot noloading
        ax[0].plot(do_pc_noloading,zrho_pc_noloading,color=colors[i],linewidth=3,alpha=0.4, label=month)
        ax[2].plot(do_hc_noloading,zrho_hc_noloading,color=colors[i],linewidth=3,alpha=0.4)
        # plot loading
        if i == 3:
            ax[0].plot(do_pc_loading,zrho_pc_loading,color='black',linewidth=1, label='Loading')
        else:
            ax[0].plot(do_pc_loading,zrho_pc_loading,color='black',linewidth=1)
        ax[2].plot(do_hc_loading,zrho_hc_loading,color='black',linewidth=1)
        # plot difference
        ax[1].plot(do_pc_loading-do_pc_noloading,zrho_pc_loading,color=colors[i],linewidth=3,alpha=0.5)
        ax[3].plot(do_hc_loading-do_hc_noloading,zrho_hc_loading,color=colors[i],linewidth=3,alpha=0.5)

# format difference panels
# ax[2].set_xlim([0,14])
ax[2].xaxis.set_ticks(np.arange(0, 13, 2))
ax[3].set_xlim([-0.2,0.2])
ax[1].axvline(0, color='silver',linestyle=':')
ax[3].axvline(0, color='silver',linestyle=':')
ax[0].axvline(2, color='silver',linestyle=':')
ax[2].axvline(2, color='silver',linestyle=':')

# add titles
ax[0].set_title('Monthly DO Profiles',fontweight='bold',fontsize=14)
ax[1].set_title(r'Loading $-$ No-loading',fontweight='bold',fontsize=14)
ax[2].set_xlabel('DO [mg/L]', fontsize=12)
ax[3].set_xlabel('DO difference [mg/L]', fontsize=12)
for axis in ax:
    axis.tick_params(axis='both', labelsize=12)
ax[0].set_ylabel('z [m]', fontsize=12)
ax[2].set_ylabel('z [m]', fontsize=12)
# format figure
ax[0].legend(loc='upper left')
fig.supylabel('Hood Canal                               Penn Cove',fontweight='bold',fontsize=14)
ax[0].text(0.8, 0.05, '(a)', fontsize = 12, fontweight='bold', transform=ax[0].transAxes)
ax[2].text(0.8, 0.05, '(b)', fontsize = 12, fontweight='bold', transform=ax[2].transAxes)
ax[1].text(0.8, 0.05, '(c)', fontsize = 12, fontweight='bold', transform=ax[1].transAxes)
ax[3].text(0.8, 0.05, '(d)', fontsize = 12, fontweight='bold', transform=ax[3].transAxes)

# Generate plot
plt.tight_layout()
plt.show()