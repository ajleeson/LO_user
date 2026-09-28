"""
This script determines the average DO concentration
in the shallow (<= 10 m) and deep layer of the water column and
saves the data in a new .nc file.

This script searches for yearly box extractions in LO_output, for the
region "pugetsoundDO"


.nc files are saved in LO_output/chapter_2/data
"""

# import things
import numpy as np
import xarray as xr
import csv
import pinfo
from lo_tools import Lfun, zrfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()


##############################################################
##                       USER INPUTS                        ##
##############################################################

years = ['2015','2016','2017','2018','2019','2020']

gtagexes = ['cas7_t1_x11ab','cas7_t1noDIN_x11ab']

out_dir = Ldir['LOo'] / 'chapter_2' / 'data'
Lfun.make_dir(out_dir)

##############################################################
##                    HELPER FUNCTIONS                      ##
##############################################################

def start_ds(ocean_time):
    '''
    Initialize dataset to store processed DO data
    ocean_time = ocean time vector
    eta_rho = eta rho vector
    xi_rho = xi_rho vector
    '''
    Ndays = len(ocean_time.values)

    ds = xr.Dataset(data_vars=dict(

        # average shallow DO concentration (z <= 10 m)
        shallow_DO_mgL_ps = (['ocean_time'], np.zeros((Ndays))),
        shallow_DO_mgL_wb = (['ocean_time'], np.zeros((Ndays))),
        shallow_DO_mgL_ss = (['ocean_time'], np.zeros((Ndays))),
        shallow_DO_mgL_mb = (['ocean_time'], np.zeros((Ndays))),
        shallow_DO_mgL_hc = (['ocean_time'], np.zeros((Ndays))),

        # average deep DO concentration (z > 10 m)
        deep_DO_mgL_ps  = (['ocean_time'], np.zeros((Ndays))),
        deep_DO_mgL_wb  = (['ocean_time'], np.zeros((Ndays))),
        deep_DO_mgL_ss  = (['ocean_time'], np.zeros((Ndays))),
        deep_DO_mgL_mb  = (['ocean_time'], np.zeros((Ndays))),
        deep_DO_mgL_hc  = (['ocean_time'], np.zeros((Ndays))),),

    coords=dict(ocean_time=ocean_time),)
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for processed DO data
    '''

    ds['shallow_DO_mgL_ps'].attrs['long_name'] = 'average DO concentration in the shallow (z <= 10 m) layer of Puget Sound'
    ds['shallow_DO_mgL_ps'].attrs['units'] = 'mg/L'
    ds['shallow_DO_mgL_wb'].attrs['long_name'] = 'average DO concentration in the shallow (z <= 10 m) layer of Whidbey Basin'
    ds['shallow_DO_mgL_wb'].attrs['units'] = 'mg/L'
    ds['shallow_DO_mgL_ss'].attrs['long_name'] = 'average DO concentration in the shallow (z <= 10 m) layer of South Sound'
    ds['shallow_DO_mgL_ss'].attrs['units'] = 'mg/L'
    ds['shallow_DO_mgL_mb'].attrs['long_name'] = 'average DO concentration in the shallow (z <= 10 m) layer of Mabin Basin'
    ds['shallow_DO_mgL_mb'].attrs['units'] = 'mg/L'
    ds['shallow_DO_mgL_hc'].attrs['long_name'] = 'average DO concentration in the shallow (z <= 10 m) layer of Hood Canal'
    ds['shallow_DO_mgL_hc'].attrs['units'] = 'mg/L'

    ds['deep_DO_mgL_ps'].attrs['long_name'] = 'average DO concentration in the deep (z > 10 m) layer of Puget Sound'
    ds['deep_DO_mgL_ps'].attrs['units'] = 'mg/L'
    ds['deep_DO_mgL_wb'].attrs['long_name'] = 'average DO concentration in the deep (z <= 10 m) layer of Whidbey Basin'
    ds['deep_DO_mgL_wb'].attrs['units'] = 'mg/L'
    ds['deep_DO_mgL_ss'].attrs['long_name'] = 'average DO concentration in the deep (z <= 10 m) layer of South Sound'
    ds['deep_DO_mgL_ss'].attrs['units'] = 'mg/L'
    ds['deep_DO_mgL_mb'].attrs['long_name'] = 'average DO concentration in the deep (z <= 10 m) layer of Mabin Basin'
    ds['deep_DO_mgL_mb'].attrs['units'] = 'mg/L'
    ds['deep_DO_mgL_hc'].attrs['long_name'] = 'average DO concentration in the deep (z <= 10 m) layer of Hood Canal'
    ds['deep_DO_mgL_hc'].attrs['units'] = 'mg/L'

    return ds


##############################################################
##                      PROCESS DATA                        ##
##############################################################

print('Processing started...\n')

# read in masks
basin_mask_ds = grid_ds = xr.open_dataset('../../../LO_output/chapter_2/data/basin_masks_from_pugetsoundDObox.nc')
mask_rho = basin_mask_ds.mask_rho.values
mask_hc = basin_mask_ds.mask_hoodcanal.values
mask_ss = basin_mask_ds.mask_southsound.values
mask_wb = basin_mask_ds.mask_whidbeybasin.values
mask_mb = basin_mask_ds.mask_mainbasin.values
mask_ps = basin_mask_ds.mask_pugetsound.values

# get horizontal area
fp = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'box' / 'pugetsoundDO_2014.01.01_2014.12.31.nc'
box_ds = xr.open_dataset(fp)
DX = (box_ds.pm.values)**-1
DY = (box_ds.pn.values)**-1
DA = DX*DY # get area in m2

for gtagex in gtagexes:
        for year in years:
            print('{}, {}'.format(gtagex,year))

            # get data
            fp = Ldir['LOo'] / 'extract' / gtagex / 'box' / ('pugetsoundDO_'+year+'.01.01_'+year+'.12.31.nc')
            ds_raw = xr.open_dataset(fp)

            # initialize dataset
            ds = start_ds(ds_raw['ocean_time'])
            # add metadata
            ds = add_metadata(ds)

            # get water column depths
            Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
            reader = csv.DictReader(open(Sfp))
            S_dict = {}
            for row in reader:
                S_dict[row['ITEMS']] = row['VALUES']
            S = zrfun.get_S(S_dict)
            # get cell thickness
            h = ds_raw['h'].values # height of water column

            # loop over time to get z_rho and z_w:
            print('    Get sigma layer thicknesses')
            zeta = ds_raw['zeta'].values
            Nt = zeta.shape[0]
            Nz = S['N']
            dzr_all = np.empty((Nt, Nz, *h.shape))
            zr_all = np.empty((Nt, Nz, *h.shape))
            zw_all = np.empty((Nt, Nz + 1, *h.shape))
            for t in range(Nt):
                # make sure to use zeta as an input to account for SSH variability!!!
                z_rho, z_w = zrfun.get_z(h, zeta[t, :, :], S)
                dzr_all[t, :, :, :] = np.diff(z_w, axis=0) # sigma layer thickness at one time
                zr_all[t, :, :, :] = z_rho # sigma layer center point at one time
                zw_all[t, :, :, :] = z_w # sigma layer boundaries at one time

            # split the water column into shallow and deep
            ix_upper10 = zr_all >= -10 # split into upper 10 m (shallow)
            ix_deeper = ~ix_upper10
            dzr_shallow = dzr_all * ix_upper10
            dzr_deeper  = dzr_all * ix_deeper

            # get total volume over time
            # water_depth = np.nansum(dzr_all, axis=1) # [m], with shape t,y,x
            water_depth_shallow = np.nansum(dzr_shallow, axis=1) # [m], with shape t,y,x
            water_depth_deeper  = np.nansum(dzr_deeper, axis=1) # [m], with shape t,y,x

            # apply basin masks and get volume
            # PS_water_depth = water_depth * mask_ps # [m], with shape t,y,x
            # Puget Sound
            PS_water_depth_shallow = water_depth_shallow * mask_ps # [m], with shape t,y,x
            PS_water_depth_deeper  = water_depth_deeper  * mask_ps # [m], with shape t,y,x
            # South Sound
            SS_water_depth_shallow = water_depth_shallow * mask_ss # [m], with shape t,y,x
            SS_water_depth_deeper  = water_depth_deeper  * mask_ss # [m], with shape t,y,x
            # Whidbey Basin
            WB_water_depth_shallow = water_depth_shallow * mask_wb # [m], with shape t,y,x
            WB_water_depth_deeper  = water_depth_deeper  * mask_wb # [m], with shape t,y,x
            # Main Basin
            MB_water_depth_shallow = water_depth_shallow * mask_mb # [m], with shape t,y,x
            MB_water_depth_deeper  = water_depth_deeper  * mask_mb # [m], with shape t,y,x
            # Hood Canal
            HC_water_depth_shallow = water_depth_shallow * mask_hc # [m], with shape t,y,x
            HC_water_depth_deeper  = water_depth_deeper  * mask_hc # [m], with shape t,y,x

            # multiply by area to get volume of each water column
            # PS_volume_per_column = PS_water_depth * DA # [m3], with shape t,y,x
            # Puget Sound
            PS_volume_per_column_shallow = PS_water_depth_shallow * DA # [m3], with shape t,y,x
            PS_volume_per_column_deeper = PS_water_depth_deeper  * DA # [m3], with shape t,y,x
            # South Sound
            SS_volume_per_column_shallow = SS_water_depth_shallow * DA # [m3], with shape t,y,x
            SS_volume_per_column_deeper = SS_water_depth_deeper  * DA # [m3], with shape t,y,x
            # Whidbey Basin
            WB_volume_per_column_shallow = WB_water_depth_shallow * DA # [m3], with shape t,y,x
            WB_volume_per_column_deeper = WB_water_depth_deeper  * DA # [m3], with shape t,y,x
            # Main Basin
            MB_volume_per_column_shallow = MB_water_depth_shallow * DA # [m3], with shape t,y,x
            MB_volume_per_column_deeper = MB_water_depth_deeper  * DA # [m3], with shape t,y,x
            # Hood Canal
            HC_volume_per_column_shallow = HC_water_depth_shallow * DA # [m3], with shape t,y,x
            HC_volume_per_column_deeper = HC_water_depth_deeper  * DA # [m3], with shape t,y,x

            # sum over y and x to get total Puget Sound volum
            # PS_total_volume = np.sum(PS_volume_per_column, axis=(1, 2)) # [m3], with shape t
            # Puget Sound
            PS_total_volume_shallow = np.sum(PS_volume_per_column_shallow, axis=(1, 2)) # [m3], with shape t
            PS_total_volume_deeper = np.sum(PS_volume_per_column_deeper, axis=(1, 2)) # [m3], with shape t
            # South Sound
            SS_total_volume_shallow = np.sum(SS_volume_per_column_shallow, axis=(1, 2)) # [m3], with shape t
            SS_total_volume_deeper = np.sum(SS_volume_per_column_deeper, axis=(1, 2)) # [m3], with shape t
            # Whidbey Basin
            WB_total_volume_shallow = np.sum(WB_volume_per_column_shallow, axis=(1, 2)) # [m3], with shape t
            WB_total_volume_deeper = np.sum(WB_volume_per_column_deeper, axis=(1, 2)) # [m3], with shape t
            # Main Basin
            MB_total_volume_shallow = np.sum(MB_volume_per_column_shallow, axis=(1, 2)) # [m3], with shape t
            MB_total_volume_deeper = np.sum(MB_volume_per_column_deeper, axis=(1, 2)) # [m3], with shape t
            # Hood Canal
            HC_total_volume_shallow = np.sum(HC_volume_per_column_shallow, axis=(1, 2)) # [m3], with shape t
            HC_total_volume_deeper = np.sum(HC_volume_per_column_deeper, axis=(1, 2)) # [m3], with shape t


            print('    Calculate average DO concentrations')
            # Now get oxygen values at every grid cell and convert to mg/L
            # just for the month of August
            oxy_mgL = pinfo.fac_dict['oxygen'] * ds_raw['oxygen'].values # [mg/L], with shape t,z,y,x

            # get vertical integral within each sigma layer
            # DO_sigma_m_mgL = dzr_all * oxy_mgL # [m * mg/L], with shape t,z,y,x
            DO_sigma_m_mgL_shallow = dzr_shallow * oxy_mgL # [m * mg/L], with shape t,z,y,x
            DO_sigma_m_mgL_deeper   = dzr_deeper * oxy_mgL # [m * mg/L], with shape t,z,y,x

            # Sum along z to get vertical integrals
            # DO_vert_int = np.nansum(DO_sigma_m_mgL,axis=1) # [m * mg/L], with shape t,y,x
            DO_vert_int_shallow = np.nansum(DO_sigma_m_mgL_shallow,axis=1) # [m * mg/L], with shape t,y,x
            DO_vert_int_deeper  = np.nansum(DO_sigma_m_mgL_deeper,axis=1) # [m * mg/L], with shape t,y,x

            # apply basin masks
            eta_rho = ds_raw['eta_rho']
            xi_rho = ds_raw['xi_rho']
            # PS_DO_m_mgL = DO_vert_int * mask_ps # [m * mg/L], with shape t,y,x
            # Puget Sound
            PS_DO_m_mgL_shallow = DO_vert_int_shallow * mask_ps # [m * mg/L], with shape t,y,x
            PS_DO_m_mgL_deeper  = DO_vert_int_deeper  * mask_ps # [m * mg/L], with shape t,y,x
            # South Sound
            SS_DO_m_mgL_shallow = DO_vert_int_shallow * mask_ss # [m * mg/L], with shape t,y,x
            SS_DO_m_mgL_deeper  = DO_vert_int_deeper  * mask_ss # [m * mg/L], with shape t,y,x
            # Whidbey Basin
            WB_DO_m_mgL_shallow = DO_vert_int_shallow * mask_wb # [m * mg/L], with shape t,y,x
            WB_DO_m_mgL_deeper  = DO_vert_int_deeper  * mask_wb # [m * mg/L], with shape t,y,x
            # Main Basin
            MB_DO_m_mgL_shallow = DO_vert_int_shallow * mask_mb # [m * mg/L], with shape t,y,x
            MB_DO_m_mgL_deeper  = DO_vert_int_deeper  * mask_mb # [m * mg/L], with shape t,y,x
            # Hood Canal
            HC_DO_m_mgL_shallow = DO_vert_int_shallow * mask_hc # [m * mg/L], with shape t,y,x
            HC_DO_m_mgL_deeper  = DO_vert_int_deeper  * mask_hc # [m * mg/L], with shape t,y,x

            # multiply by area to get volume integral
            # PS_DO_vol_int = PS_DO_m_mgL * DA # [m3 * mg/L], with shape t,y,x
            # Puget Sound
            PS_DO_vol_int_shallow = PS_DO_m_mgL_shallow * DA # [m3 * mg/L], with shape t,y,x
            PS_DO_vol_int_deeper  = PS_DO_m_mgL_deeper  * DA # [m3 * mg/L], with shape t,y,x
            # South Sound
            SS_DO_vol_int_shallow = SS_DO_m_mgL_shallow * DA # [m3 * mg/L], with shape t,y,x
            SS_DO_vol_int_deeper  = SS_DO_m_mgL_deeper  * DA # [m3 * mg/L], with shape t,y,x
            # Whidbey Basin
            WB_DO_vol_int_shallow = WB_DO_m_mgL_shallow * DA
            WB_DO_vol_int_deeper  = WB_DO_m_mgL_deeper  * DA
            # Main Basin
            MB_DO_vol_int_shallow = MB_DO_m_mgL_shallow * DA
            MB_DO_vol_int_deeper  = MB_DO_m_mgL_deeper  * DA
            # Hood Canal
            HC_DO_vol_int_shallow = HC_DO_m_mgL_shallow * DA
            HC_DO_vol_int_deeper  = HC_DO_m_mgL_deeper  * DA

            # sum over y and x to get total DO in Puget Sound
            # PS_total_DO = np.sum(PS_DO_vol_int, axis=(1, 2)) # [m3 * mg/L], with shape t
            # Puget Sound
            PS_total_DO_shallow = np.sum(PS_DO_vol_int_shallow, axis=(1, 2)) # [m3 * mg/L], with shape t
            PS_total_DO_deeper  = np.sum(PS_DO_vol_int_deeper, axis=(1, 2)) # [m3 * mg/L], with shape t
            # South Sound
            SS_total_DO_shallow = np.sum(SS_DO_vol_int_shallow, axis=(1, 2))
            SS_total_DO_deeper  = np.sum(SS_DO_vol_int_deeper, axis=(1, 2))
            # Whidbey Basin
            WB_total_DO_shallow = np.sum(WB_DO_vol_int_shallow, axis=(1, 2))
            WB_total_DO_deeper  = np.sum(WB_DO_vol_int_deeper, axis=(1, 2))
            # Main Basin
            MB_total_DO_shallow = np.sum(MB_DO_vol_int_shallow, axis=(1, 2))
            MB_total_DO_deeper  = np.sum(MB_DO_vol_int_deeper, axis=(1, 2))
            # Hood Canal
            HC_total_DO_shallow = np.sum(HC_DO_vol_int_shallow, axis=(1, 2))
            HC_total_DO_deeper  = np.sum(HC_DO_vol_int_deeper, axis=(1, 2))

            # get average concentration over time
            # PS_avg_DO = PS_total_DO / PS_total_volume # [mg/L], with shape t
            # Puget Sound
            PS_avg_DO_shallow = PS_total_DO_shallow / PS_total_volume_shallow # [mg/L], with shape t
            PS_avg_DO_deeper  = PS_total_DO_deeper / PS_total_volume_deeper # [mg/L], with shape t
            # South Sound
            SS_avg_DO_shallow = SS_total_DO_shallow / SS_total_volume_shallow
            SS_avg_DO_deeper  = SS_total_DO_deeper / SS_total_volume_deeper
            # Whidbey Basin
            WB_avg_DO_shallow = WB_total_DO_shallow / WB_total_volume_shallow
            WB_avg_DO_deeper  = WB_total_DO_deeper / WB_total_volume_deeper
            # Main Basin
            MB_avg_DO_shallow = MB_total_DO_shallow / MB_total_volume_shallow
            MB_avg_DO_deeper  = MB_total_DO_deeper / MB_total_volume_deeper
            # Hood Canal
            HC_avg_DO_shallow = HC_total_DO_shallow / HC_total_volume_shallow
            HC_avg_DO_deeper  = HC_total_DO_deeper / HC_total_volume_deeper


            # # outflow of skagit river
            # test_lat = 48.30
            # test_lon = -122.42
            # lons = ds_raw['lon_rho'].values[0,:]
            # lats = ds_raw['lat_rho'].values[:,0]
            # test_y = min(range(len(lats)), key=lambda i: abs(lats[i]-test_lat)) + 1
            # test_x = min(range(len(lons)), key=lambda i: abs(lons[i]-test_lon)) + 1
            # print(oxy_mgL[0,0:2,test_y,test_x])
            # print(dzr_all[0,0:2,test_y,test_x])
            # print(DO_sigma_m_mgL[0,0:2,test_y,test_x])
            # print(DO_vert_int[0,test_y,test_x])
            # print(PS_DO_m_mgL[0,test_y,test_x])


            # add data to ds
            print('    Adding data to dataset')

            # Average Puget Sound DO concentration time series in depths 10 m or shallower
            # Puget Sound
            ds['shallow_DO_mgL_ps'] = xr.DataArray(PS_avg_DO_shallow,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])
            # South Sound
            ds['shallow_DO_mgL_ss'] = xr.DataArray(SS_avg_DO_shallow,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])
            # Whidbey Basin
            ds['shallow_DO_mgL_wb'] = xr.DataArray(WB_avg_DO_shallow,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])
            # Main Basin
            ds['shallow_DO_mgL_mb'] = xr.DataArray(MB_avg_DO_shallow,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])
            # Hood Canal
            ds['shallow_DO_mgL_hc'] = xr.DataArray(HC_avg_DO_shallow,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])
            
            # Average Puget Sound DO concentration time series in deeper than 10 m
            # Puget Sound
            ds['deep_DO_mgL_ps'] = xr.DataArray(PS_avg_DO_deeper,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])
            # South Sound
            ds['deep_DO_mgL_ss'] = xr.DataArray(SS_avg_DO_deeper,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])
            # Whidbey Basin
            ds['deep_DO_mgL_wb'] = xr.DataArray(WB_avg_DO_deeper,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])
            # Main Basin
            ds['deep_DO_mgL_mb'] = xr.DataArray(MB_avg_DO_deeper,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])
            # Hood Canal
            ds['deep_DO_mgL_hc'] = xr.DataArray(HC_avg_DO_deeper,
                                        coords={'ocean_time': ds_raw['ocean_time'].values},
                                        dims=['ocean_time'])

            print('    Saving dataset')
            ds.to_netcdf(out_dir / ('basins_shallow10m_deep_DO_' + gtagex + year + '.nc'))

print('Done')