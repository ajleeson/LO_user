"""
Helper script called by make_forcing_main
to generate forcing for point sources 
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

import sys
import os
import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd
import rivfun
import trapsfun

#################################################################################
#                   Initialize function and empty dataset                       #
#################################################################################

def make_forcing(N,NT,NRIV,NTRIV,NWWTP_moh,dt_ind, yd_ind,ot_vec,Ldir,enable,trapsP,trapsD,ctag):

    # Start Dataset
    wwtp_ds = xr.Dataset()
    NWWTP = 0

    # get year list
    years = [fulldate.year for fulldate in dt_ind]

    # adjust date format
    dt_ind = dt_ind.normalize()

#################################################################################
#                                Get data                                       #
#################################################################################

    # only get data if WWTPs are enabled
    if enable == True:

        # get raw data
        was24_wwtp_fn = Ldir['data'] / trapsD / 'processed_data'/ 'wwtp_data_wasielewski_etal_2024.nc'
        was24_wwtp_data_ds = xr.open_dataset(was24_wwtp_fn)


        # first, make sure file exists
        gwi_fn = Ldir['grid'] / 'was24_wwtp_info.csv'
        if not os.path.isfile(gwi_fn):
            print('***Missing was24_wwtp_info.csv file. Please run traps_placement')
            sys.exit()
        # then get the list of WWTPs and indices for this grid
        gwi_df = pd.read_csv(gwi_fn, index_col='rname')
        # if testing, only look at a few sources
        if Ldir['testing']:
            gwi_df = gwi_df.loc[['King County West Point WWTP', 'BELLINGHAM STP'],:]

#################################################################################
#                       Prepare dataset for data                                #
#################################################################################

        # get number of wwtps after consolidating overlapping ones
        NWWTP = len(gwi_df)

        # Add time coordinate
        wwtp_ds['river_time'] = (('river_time',), ot_vec)
        wwtp_ds['river_time'].attrs['units'] = Lfun.roms_time_units
        wwtp_ds['river_time'].attrs['long_name'] = 'river time'

        # Add river coordinate
        wwtp_ds['river'] = (('river',), np.arange(NRIV+NTRIV+NWWTP_moh+1,NRIV+NTRIV+NWWTP_moh+NWWTP+1))
        wwtp_ds['river'].attrs['long_name'] = 'marine point source identification number'

        # Add river names
        wwtp_ds['river_name'] = (('river',), list(gwi_df.index))
        wwtp_ds['river_name'].attrs['long_name'] = 'point source name'

#################################################################################
#  Add vertical distribution of sources. All WWTPs discharge from bottom layer  #
#################################################################################

        # Add Vshape
        vn = 'river_Vshape'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = ('s_rho', 'river')
        # All discharge coming from the bottom layer
        Vshape = np.zeros((N, NWWTP))
        Vshape[0,:] = 1
        wwtp_ds[vn] = (dims, Vshape)
        wwtp_ds['river_Vshape'].attrs['long_name'] = vinfo['long_name']

#################################################################################
#             Add indices of sources. WWTPs located on the rho-grid             #
#################################################################################

        # Add position and direction
        for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            # set point source diretion to enter vertically (Dsrc = 2)
            if vn == 'river_direction':
                wwtp_direction = 2 * np.ones(NWWTP) 
                wwtp_ds[vn] = (('river',), wwtp_direction)
            # Add X-position (column index)
            elif vn == 'river_Xposition':
                X_vec = np.nan * np.ones(NWWTP)
                for ii,wn in enumerate(gwi_df.index):
                    X_vec[ii] = gwi_df.loc[wn, 'col_py']
                wwtp_ds[vn] = (('river',), X_vec)
            # Add E-position (row index)
            elif vn == 'river_Eposition':
                E_vec = np.nan * np.ones(NWWTP)
                for ii,wn in enumerate(gwi_df.index):
                    E_vec[ii] = gwi_df.loc[wn, 'row_py']
                wwtp_ds[vn] = (('river',), E_vec)
            # add metadata
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']

#################################################################################
#                               Add source flowrate                             #
#################################################################################

        # Add transport
        vn = 'river_transport'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('river',)
        Q_mat = np.zeros((NT, NWWTP))
        for rr,rn in enumerate(gwi_df.index):
            flow = was24_wwtp_data_ds.flow.sel(source=
                                               was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == rn].item(),
                                               date=dt_ind).values
            # update flowrate with open/close date information
            Q_mat[:,rr] = flow
        # add metadata
        wwtp_ds[vn] = (dims, Q_mat)
        wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
        wwtp_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                         Add source salinity and temp                          #
#################################################################################

        # Add salinity and temperature
        for vn in ['river_salt', 'river_temp']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            # salinity is always zero
            if vn == 'river_salt':
                TS_mat = np.zeros((NT, N, NWWTP))
            # get temperature from climatology
            elif vn == 'river_temp':
                TS_mat = np.nan * np.zeros((NT, N, NWWTP))
                for rr,rn in enumerate(gwi_df.index):
                    temps = flow = was24_wwtp_data_ds.temp.sel(source=
                                               was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == rn].item(),
                                               date=dt_ind).values
                    print('-- {}: filled from raw Wasielewski et al. (2024) dataset'.format(rn))
                    for nn in range(N):
                        TS_mat[:, nn, rr] = temps
            # check for nans
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in point source river_temp!')
                sys.exit()
            # add metadata
            wwtp_ds[vn] = (dims, TS_mat)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
            wwtp_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                            Add source biology                                 #
#################################################################################

        # Add biologeochemistry parameters
        for var in ['NO3', 'NH4', 'TIC', 'TAlk', 'Oxyg']:
        # for var in ['NO3', 'NH4', 'TIC', 'Talk', 'DO']:
            vn = 'river_' + var
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NWWTP))
            for rr,rn in enumerate(gwi_df.index):
                # adjust names to get data from dataset
                if var == 'TAlk':
                    var = 'Talk'
                if var == 'Oxyg':
                    var = 'DO'
                bvals = was24_wwtp_data_ds[var].sel(
                        source=was24_wwtp_data_ds.source[was24_wwtp_data_ds.name == rn].item(),
                        date=dt_ind)
                for nn in range(N):
                    B_mat[:, nn, rr] = bvals
            # check for nans
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in tiny river bio!')
                sys.exit()
            # add metadata
            wwtp_ds[vn] = (dims, B_mat)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
            wwtp_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#                  All other biology variables are zero                         #
#################################################################################

        # Add remaining biology (see the lineup near the end of fennel_var.h)
        # Right now, this is simply filling everything with zeros
        bvn_list = ['Phyt', 'Zoop', 'LDeN', 'SDeN', 'Chlo', 'LDeC', 'SDeC']
        for bvn in bvn_list:
            vn = 'river_' + bvn
            vinfo = zrfun.get_varinfo(vn)
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NWWTP))
            # loop through all sources and fill with zeros
            for rr,rn in enumerate(gwi_df.index):
                for nn in range(N):
                    B_mat[:, nn, rr] = rivfun.get_bio_vec(bvn, rn, yd_ind)
            # check for nans
            if np.isnan(B_mat).any():
                print('Error from traps: nans in B_mat for tiny river ' + vn)
                sys.exit()
            # add metadata
            wwtp_ds[vn] = (dims, B_mat)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
            wwtp_ds[vn].attrs['units'] = vinfo['units']

#################################################################################
#          Return WWTP forcing dataset in the form that ROMS expects            #
#################################################################################

    return wwtp_ds, NWWTP
