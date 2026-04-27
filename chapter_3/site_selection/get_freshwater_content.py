"""
Get freshwater content from model history file

.nc files are saved in LO_output/chapter_3/data
"""

# import things
import numpy as np
import xarray as xr
import csv
import matplotlib.pylab as plt
from pathlib import Path
import gsw
import cmcrameri.cm as cmc
from lo_tools import Lfun, zrfun
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

years = ['2015','2016','2017','2018','2019','2020','2021','2022','2023','2024','2025']

# which  model run to look at?
gtagex = 'cas7_t1_x11ab' # 'cas7_t1d_x11ad'

# where to put output files
out_dir = Ldir['LOo'] / 'chapter_3' / 'data'
Lfun.make_dir(out_dir)

##############################################################
##                    HELPER FUNCTIONS                      ##
##############################################################

def start_ds(ocean_time,eta_rho,xi_rho):
    '''
    Initialize dataset to store processed freshwater content
    ocean_time = ocean time vector
    eta_rho = eta rho vector
    xi_rho = xi_rho vector
    '''
    Ndays = len(ocean_time.values)
    Neta = len(eta_rho.values)
    Nxi = len(xi_rho.values)

    ds = xr.Dataset(data_vars=dict(
        # Freshwater content [m]
        Fs_32p5      = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        Fs_31        = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        Fs_30        = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        # Freshwater content, normalized by depth
        Fs_32p5_norm   = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        Fs_31_norm     = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        Fs_30_norm     = (['ocean_time','eta_rho','xi_rho'], np.zeros((Ndays,Neta,Nxi))),
        ),
    coords=dict(ocean_time=ocean_time, eta_rho=eta_rho, xi_rho=xi_rho,),)
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for processed freshwater content data
    '''
    ds['Fs_32p5'].attrs['long_name'] = 'freshwater content (s0 = 32.5)'
    ds['Fs_32p5'].attrs['units'] = 'm'

    ds['Fs_32p5_norm'].attrs['long_name'] = 'freshwater content normalized by depth (s0 = 32.5)'
    ds['Fs_32p5_norm'].attrs['units'] = 'dimensionless'

    ds['Fs_31'].attrs['long_name'] = 'freshwater content (s0 = 31)'
    ds['Fs_31'].attrs['units'] = 'm'

    ds['Fs_31_norm'].attrs['long_name'] = 'freshwater content normalized by depth (s0 = 31)'
    ds['Fs_31_norm'].attrs['units'] = 'dimensionless'

    ds['Fs_30'].attrs['long_name'] = 'freshwater content (s0 = 30)'
    ds['Fs_30'].attrs['units'] = 'm'

    ds['Fs_30_norm'].attrs['long_name'] = 'freshwater content normalized by depth (s0 = 30)'
    ds['Fs_30_norm'].attrs['units'] = 'dimensionless'

    return ds


##############################################################
##                      PROCESS DATA                        ##
##############################################################

print('Processing started...\n')

for year in years:

    print('{}'.format(year))

    # Initialize empty list of datasets (meant for daily datasets)
    ds_list = []

    # get info to find history files
    gridname, tag, ex_name = gtagex.split('_')
    # get the dict Ldir
    Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
    # add more entries to Ldir
    Ldir['roms_out_num'] = 5
    Ldir['ds0'] = year + '.01.01'
    Ldir['ds1'] = year + '.12.31'
    Ldir['list_type'] = 'average'
    # get history files
    fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir, Ldir['ds0'], Ldir['ds1'])

    # loop through history files for the year
    for i,fn in enumerate(fn_list):

        date_str = Path(fn).parent.name
        print('    ' + date_str) 
        # get data
        ds_raw = xr.open_dataset(fn)
        # initalize dataset to store daily data
        daily_ds = start_ds(ds_raw['ocean_time'],
                        ds_raw['eta_rho'],
                        ds_raw['xi_rho'])
        
        # GET THE CELL THICKNESS FOR VERTICAL INTEGRALS
        # get S for the whole grid
        Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
        reader = csv.DictReader(open(Sfp))
        S_dict = {}
        for row in reader:
            S_dict[row['ITEMS']] = row['VALUES']
        S = zrfun.get_S(S_dict)
        # get cell thickness
        h = ds_raw['h'].values # height of water column
        # loop over time to get z_rho and z_w:
        zeta = ds_raw['zeta'].values
        Nt = zeta.shape[0]
        Nz = S['N']
        dzr_all = np.empty((Nt, Nz, *h.shape))
        z_w_all = np.empty((Nt, Nz + 1, *h.shape))
        z_rho_all = np.empty((Nt, Nz, *h.shape))
        # loop over time -- but in this case, there is only one time step, so this loop will only run once
        for t in range(Nt):
            # make sure to use zeta as an input to account for SSH variability!!!
            z_rho, z_w = zrfun.get_z(h, zeta[t, :, :], S)
            dzr_all[t, :, :, :] = np.diff(z_w, axis=0) # sigma layer thickness at one time
            z_w_all[t, :, :, :] = z_w # sigma layer boundaries at one time
            z_rho_all[t, :, :, :] = z_rho # sigma layer centers at one time

        # convert practical salinity to absolute salinity
        SP = ds_raw['salt'].values # practical salinity
        # z = z_rho_all
        # lons = ds_raw['lon_rho'].values
        # lats = ds_raw['lat_rho'].values
        # p = gsw.p_from_z(z, lats)
        # SA = gsw.SA_from_SP(SP, p, lons, lats)

        # CALCULATE SO - S(x,z) / S0 (use values of s0 = 32.5, 31, and 30)
        deltaS_32p5 = (32.5 - SP) / 32.5
        deltaS_31   = (31   - SP) / 31
        deltaS_30   = (30   - SP) / 30
        # set values < 0 to 0, so that we only sum where s < s0
        deltaS_32p5[deltaS_32p5 < 0] = 0
        deltaS_31[deltaS_31 < 0]     = 0
        deltaS_30[deltaS_30 < 0]     = 0

        # vertically integrate to get freshwater content
        freshwater_content_32p5 = np.nansum(deltaS_32p5 * dzr_all,axis=1)
        freshwater_content_31   = np.nansum(deltaS_31   * dzr_all,axis=1)
        freshwater_content_30   = np.nansum(deltaS_30   * dzr_all,axis=1)
        # apply land mask
        land_mask_2d = ds_raw['mask_rho'].values == 0
        land_mask_3d = land_mask_2d[np.newaxis, :, :]
        freshwater_content_32p5 = np.ma.masked_where(land_mask_3d, freshwater_content_32p5)
        freshwater_content_31   = np.ma.masked_where(land_mask_3d, freshwater_content_31)
        freshwater_content_30   = np.ma.masked_where(land_mask_3d, freshwater_content_30)

        # normalize by depth
        freshwater_content_normalized_32p5 = freshwater_content_32p5 / (z_w_all[0,-1, :, :] - z_w_all[0,0, :, :])
        freshwater_content_normalized_31   = freshwater_content_31   / (z_w_all[0,-1, :, :] - z_w_all[0,0, :, :])
        freshwater_content_normalized_30   = freshwater_content_30   / (z_w_all[0,-1, :, :] - z_w_all[0,0, :, :])

        # add data to daily dataset
        daily_ds['Fs_32p5'] = xr.DataArray(freshwater_content_32p5,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['Fs_31'] = xr.DataArray(freshwater_content_31,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['Fs_30'] = xr.DataArray(freshwater_content_30,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['Fs_32p5_norm'] = xr.DataArray(freshwater_content_normalized_32p5, 
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['Fs_31_norm'] = xr.DataArray(freshwater_content_normalized_31,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        daily_ds['Fs_30_norm'] = xr.DataArray(freshwater_content_normalized_30,
                                    coords={'ocean_time': ds_raw['ocean_time'].values,
                                            'eta_rho': ds_raw['eta_rho'].values,
                                            'xi_rho': ds_raw['xi_rho'].values},
                                    dims=['ocean_time','eta_rho', 'xi_rho'])
        # Append to list
        ds_list.append(daily_ds)

    # PREPARE DATA FOR SAVING
    # Combine all of the daily datasets into one dataset
    ds = xr.concat(ds_list, dim='ocean_time')
    # Add your metadata once at the end
    ds = add_metadata(ds)
    print('    Saving dataset')
    ds.to_netcdf(out_dir / (gtagex + '_' + year + '_FreshwaterContent.nc'))


    # print(ds)

print('Done')


# ##############################################
# # TEST PLOTTING

# # get pcolormesh values
# Fs_32p5 = ds['Fs_32p5'][-1,:,:].values

# # make list of vars
# vars = [Fs_32p5]

# # list of vmins and vmax
# vmins = [0,]
# vmaxs = [15]

# # get colormaps
# cmaps = [cmc.batlowW_r]

# # titles
# titles = ['Freshwater content\n'+r'(s$_0$ = 32.5) [m]']

# # Get grid data
# G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
# grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
# lon = grid_ds.lon_rho.values
# lat = grid_ds.lat_rho.values
# lon_u = grid_ds.lon_u.values
# lat_u = grid_ds.lat_u.values
# lon_v = grid_ds.lon_v.values
# lat_v = grid_ds.lat_v.values
# px, py = pfun.get_plon_plat(lon,lat)

# # lon/lat limits (Study Domain)
# xmin = -126
# xmax = -122
# ymin = 45.5
# ymax = 50.5

# for var, vmin, vmax, cmap, title in zip(vars, vmins, vmaxs, cmaps, titles):

#     # Initialize figure
#     fig,ax = plt.subplots(1,1, figsize=(7,9))

#     # plot values
#     cs = ax.pcolormesh(px,py,var,vmin=vmin, vmax=vmax, cmap=cmap)

#     # Add Puget Sound Inset
#     # [x0, y0, width, height]
#     axins = ax.inset_axes([0.71, 0.0, 0.45, 0.6])
#     # plot values in inset
#     axins.pcolormesh(px, py, var, vmin=vmin, vmax=vmax, cmap=cmap)
#     # Puget Sound limits
#     axins.set_xlim(-123.2, -122.1)
#     axins.set_ylim(46.95, 48.4)
#     axins.tick_params(left=False, bottom=False)
#     # format
#     axins.set_xticklabels([])
#     axins.set_yticklabels([])
#     for spine in axins.spines.values():
#         spine.set_edgecolor('grey')
#         spine.set_linewidth(2)

#     # add colorbar
#     cbar = fig.colorbar(cs, ax=ax, location='bottom', shrink=0.7, pad=0.03)
#     cbar.ax.tick_params(labelsize=18, rotation=30)
#     cbar.outline.set_visible(False)

#     # format figure
#     ax.set_xlim([xmin,xmax])
#     ax.set_ylim([ymin,ymax])
#     ax.set_yticklabels([])
#     ax.set_xticklabels([])
#     ax.tick_params(left=False, bottom=False)
#     pfun.add_coast(ax, color='silver')
#     pfun.add_coast(axins, color='silver')
#     pfun.dar(ax)
#     pfun.dar(axins)
#     ax.set_title(title, fontsize=20,
#                 loc='Left', fontweight='bold')

#     # Generate plot
#     plt.tight_layout
#     plt.subplots_adjust(bottom=0.001, top=0.9)
#     plt.show()