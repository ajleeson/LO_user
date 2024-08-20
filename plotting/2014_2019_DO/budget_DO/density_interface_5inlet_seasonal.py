"""
Plot average density profile at 5 inlets, by season
To determine depth of interface
For 2-layer simplification
"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import gsw
import pickle

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
jobname = 'twentyoneinlets'
startdate = '2014.01.01'
enddate = '2014.12.31'
year = '2014' # for making a date label

# hypoxic season
start = '08-01'
end = '09-30'


# figure settings
fs = 12 # figure font size
ls = 11 # label size
ts = 14 # title size

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures' / 'twentyone' 
Lfun.make_dir(out_dir)

# get observations
in_dir = Ldir['parent'] / 'LO_output' / 'obsmod'
in_fn = in_dir / ('multi_ctd_' + year + '.p')
df_dict = pickle.load(open(in_fn, 'rb'))
# only look at ecology stations
source = 'ecology_nc'
df_obs = df_dict['obs'].loc[df_dict['obs'].source==source,:]
df_model = df_dict[gtagex].loc[df_dict[gtagex].source==source,:]

##########################################################
##                 Plot station locations               ##
##########################################################

plt.close('all')

# get the grid data
ds = xr.open_dataset('../../../../LO_data/grids/cas7/grid.nc')
z = -ds.h.values
mask_rho = np.transpose(ds.mask_rho.values)
lon = ds.lon_rho.values
lat = ds.lat_rho.values
# make a version of z with nans where masked
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan
zm[np.transpose(mask_rho) != 0] = -1

fig, ax = plt.subplots(5,4,figsize = (10,8), sharex=True)
# ax = axes.ravel()

stations = ['lynchcove', 'carr', 'penn', 'case', 'budd']
seasons = ['Win (Jan-Mar)','Spr (Apr-Jun)','Sum (Jul-Sep)','Aut (Oct-Dec)']
season_start = ['01-01','04-01','07-01','10-01']
season_end = ['03-31','06-30','09-30','12-31']
season_colors = ['aliceblue','lavenderblush','honeydew','antiquewhite']
season_colors_dark = ['cornflowerblue','hotpink','yellowgreen','sandybrown']

for i,station in enumerate(stations):
    ax[i,0].text(0.1,0.85,station,transform=ax[i,0].transAxes, fontweight='bold', fontsize=14)
    
    # calculate lat/lon for station
    lon = sta_dict[station][0]
    lat = sta_dict[station][1]

    # download .nc files
    fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds_orig = xr.open_dataset(fn)

    for j,season in enumerate(seasons):
        # crop to season
        ds = ds_orig.sel(ocean_time=slice(np.datetime64(year+'-'+season_start[j]),np.datetime64(year+'-'+season_end[j])))

        # get density
        # calculate absolute salinity from practical salinity
        salt_abs_all = gsw.conversions.SA_from_SP(ds['salt'].values, ds['z_rho'].values, lon, lat)
        # calculate conservative temperature from potential temperature
        cons_temp_all = gsw.conversions.CT_from_pt(salt_abs_all, ds['temp'].values)
        # calculate potential density
        rho_all = gsw.density.sigma0(salt_abs_all,cons_temp_all)
        # average in time
        rho_ann_avg = rho_all.mean(axis=0)

        # get DO
        DO = ds['oxygen'].values * 32/1000 # mg/L
        DO_ann_avg = DO.mean(axis=0)

        # get depth
        z_rho = ds['z_rho'].values
        z_rho_ann_avg = z_rho.mean(axis=0)

        # plot annual average density profile
        ax[i,j].plot(rho_ann_avg,z_rho_ann_avg, color='deepskyblue',label=r'$\rho$')
        ax[i,j].set_xlim([14,26])
        ax[i,j].set_xticks(np.arange(14, 28, 4))
        # plot annual average oxygen profile
        ax_oxy = ax[i,j].twiny()
        ax_oxy.plot(DO_ann_avg,z_rho_ann_avg, color='mediumorchid',
                    linewidth=3,alpha=0.5,label='DO')
        ax_oxy.set_xlim([0,12])
        ax_oxy.set_xticks(np.arange(0, 14, 4))
        if i != 0:
            ax_oxy.set_xticks([])
        

        # decide if one or two layers
        if i == 0:
            ax[i,0].text(-0.3,1.25,'Two-layers if:\n'+r'$\rho_{bot}-\rho_{top}$ > 1 kg/m$^3$',
                    transform=ax[i,0].transAxes)
        #     if i == 0:
        #         ax[0].text(0.1,1.16,'Two-layers if:\n'+r'$DO_{bot}-DO_{top}$ > 2 mg/L',
        #                 transform=ax[i].transAxes)

        if rho_ann_avg[0] - rho_ann_avg[-1] <= 1:
            ax[i,j].text(0.1,0.8,'One-layer',transform=ax[i,j].transAxes)
        # if max(DO_ann_avg) - min(DO_ann_avg) < 2:
        #     ax[i].text(0.1,0.8,'One-layer',transform=ax[i].transAxes)

        else: # two layers
            # calculate interface depth (based on density)
            drho = np.diff(rho_ann_avg) # delta rho
            dz = np.diff(z_rho_ann_avg) # delta z
            drho_dz = drho/dz # drho/dz
            indmax_drho_dz = np.argmin(drho_dz) # largest magnitude drho/dz (min because negative)
            rho_max_slope = np.mean([rho_ann_avg[indmax_drho_dz],rho_ann_avg[indmax_drho_dz+1]])
            z_max_slope = np.mean([z_rho_ann_avg[indmax_drho_dz],z_rho_ann_avg[indmax_drho_dz+1]])
            # plot location of max slope (interface)
            ax[i,j].scatter(rho_max_slope,z_max_slope, s=30, color='deepskyblue', zorder=3)
            # label interface depth on plot
            ax[i,j].text(0.05,0.18,'{} m'.format(round(z_max_slope,1)),transform=ax[i,j].transAxes,
                    color='deepskyblue', fontweight='bold')

            # calculate interface depth (based on oxygen)
            dDO = np.diff(DO_ann_avg) # delta rho
            dz = np.diff(z_rho_ann_avg) # delta z
            dDO_dz = dDO/dz # drho/dz
            indmax_dDO_dz = np.argmax(dDO_dz) # largest magnitude dDO/dz (max because positive)
            DO_max_slope = np.mean([DO_ann_avg[indmax_dDO_dz],DO_ann_avg[indmax_dDO_dz+1]])
            z_max_slope = np.mean([z_rho_ann_avg[indmax_dDO_dz],z_rho_ann_avg[indmax_dDO_dz+1]])
            # plot location of max slope (interface)
            ax_oxy.scatter(DO_max_slope,z_max_slope, s=30, color='mediumorchid', zorder=3)
            # label interface depth on plot
            ax_oxy.text(0.05,0.05,'{} m'.format(round(z_max_slope,1)),transform=ax[i,j].transAxes,
                    color='mediumorchid', fontweight='bold')
            
            # add season labels
            if i == 0:
                ax[i,j].text(0.9,0.05,season,fontweight='bold',transform=ax[i,j].transAxes,
                            fontsize=10,color=season_colors_dark[j],rotation=90)
        

        # format figure
        ax[i,j].set_ylim([z_rho_ann_avg.min(),-1*z_rho_ann_avg.min()*0.15])
        # format grid
        ax[i,j].grid(True,color='lightgray',linewidth=1)
        if j > 0:
            ax[i,j].yaxis.set_tick_params(labelleft=False)
        else:
            ax[i,j].set_ylabel('Depth [m]')
        # format colors
        ax[i,j].set_facecolor(season_colors[j])#'#EEEEEE')
        ax[i,j].tick_params(axis='x', labelcolor='deepskyblue')
        if i < 5:
            ax_oxy.tick_params(axis='x', labelcolor='mediumorchid')
        else:
            ax_oxy.set_xticklabels([])
        for border in ['top','right','bottom','left']:
            ax[i,j].spines[border].set_visible(False)
            ax_oxy.spines[border].set_visible(False)


plt.suptitle(r'   Depth [m] vs. seasonal average DO [mg L$^{-1}$] and $\rho$ [kg m$^{-3}$]',
             fontsize=16)
# highlight labels in title
rectDO = Rectangle((2.26, 1.25), width=0.7, height=0.25, color='mediumorchid',
    transform=ax[0,0].transAxes, clip_on=False, alpha=0.2)
rectRHO = Rectangle((3.24, 1.25), width=0.62, height=0.25, color='deepskyblue',
    transform=ax[0,0].transAxes, clip_on=False, alpha=0.2)
ax[0,0].add_patch(rectDO)
ax[0,0].add_patch(rectRHO)
plt.tight_layout
plt.subplots_adjust(hspace=0.08, wspace=0.12, top=0.9,
                     bottom=0.05, left=0.08, right=0.95)
plt.savefig(out_dir / 'density_5inlet_season.png')
