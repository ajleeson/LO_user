
from datetime import datetime, timedelta
from lo_tools import forcing_argfun2 as ffun
import sys
import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd
import rivfun
import trapsfun

def make_forcing(N,NT,NRIV,NTRIV,dt_ind, yd_ind,ot_vec,Ldir,enable):

    # Start Dataset
    wwtp_ds = xr.Dataset()

    NWWTP = 0

    if enable == True:

        # define directory for point_source climatology
        wwtp_dir = Ldir['LOo'] / 'pre' / 'traps' / 'point_sources'
        traps_type = 'wwtp'  

        # climatological data files
        year0 = 1999
        year1 = 2017
        # climatological data
        Ldir['Cflow_wwtp_fn'] = wwtp_dir / 'Data_historical' / ('CLIM_flow_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['Ctemp_wwtp_fn'] = wwtp_dir / 'Data_historical' / ('CLIM_temp_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CDO_wwtp_fn']   = wwtp_dir / 'Data_historical' / ('CLIM_DO_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CNH4_wwtp_fn']  = wwtp_dir / 'Data_historical' / ('CLIM_NH4_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CNO3_wwtp_fn']  = wwtp_dir / 'Data_historical' / ('CLIM_NO3_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CTalk_wwtp_fn'] = wwtp_dir / 'Data_historical' / ('CLIM_Talk_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CTIC_wwtp_fn']  = wwtp_dir / 'Data_historical' / ('CLIM_TIC_' + str(year0) + '_' + str(year1) + '.p')

        # get the list of point sources and indices for this grid
        gri_fn = Ldir['grid'] / 'wwtp_info.csv'
        gri_df = pd.read_csv(gri_fn, index_col='rname')
        if Ldir['testing']:
            gri_df = gri_df.loc[['West Point', 'Birch Bay', 'Tacoma Central', 'US Oil & Refining'],:]
        # gri_df = gri_df.drop('Birch Bay') # Remove the Birch Bay treatment plant
        NWWTP = len(gri_df)
        
        # get list of overlapping point sources
        overlapping_wwtps = gri_df[gri_df.duplicated(['row_py','col_py'], keep=False) == True].index.values
        # consolidate overlapping point sources
        combined_names = trapsfun.combine_adjacent(overlapping_wwtps)
        gri_df_no_ovrlp = pd.DataFrame(columns=gri_df.columns) 
        gri_df_no_ovrlp.index.name='rname'
        for psname in gri_df.index: # loop through original dataframe
            if psname in overlapping_wwtps: # look for point sources that are in the list of duplicates
                name_index = np.where(overlapping_wwtps == psname)[0][0] # get index in the list of duplicates
                if name_index%2 == 0: # even index means first occurence of duplicate
                    newname = combined_names[int(name_index/2)] # combine names of duplicates
                    gri_df_no_ovrlp.loc[newname] = gri_df.loc[psname] # add combined source to dataframe
                # Note: second duplicate will be dropped
            else:
                gri_df_no_ovrlp.loc[psname] = gri_df.loc[psname] # if not a duplicate, then just copy over original info

        NWWTP = len(gri_df_no_ovrlp)
        # NWWTP = len(gri_df)

        # get the flow, temperature, and nutrient data for these days
        qtbio_wwtp_df_dict = trapsfun.get_qtbio(gri_df, dt_ind, yd_ind, Ldir, traps_type)

        # Add time coordinate
        wwtp_ds['river_time'] = (('river_time',), ot_vec)
        wwtp_ds['river_time'].attrs['units'] = Lfun.roms_time_units
        wwtp_ds['river_time'].attrs['long_name'] = 'river time'

        # Add river coordinate
        wwtp_ds['river'] = (('river',), np.arange(NRIV+NTRIV+1,NRIV+NTRIV+NWWTP+1))
        wwtp_ds['river'].attrs['long_name'] = 'marine point source identification number'

        # Add river names
        wwtp_ds['river_name'] = (('river',), list(gri_df_no_ovrlp.index))
        wwtp_ds['river_name'].attrs['long_name'] = 'point source name'

        # Add Vshape
        vn = 'river_Vshape'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = ('s_rho', 'river')
        # All discharge coming from the bottom layer
        Vshape = np.zeros((N, NWWTP))
        Vshape[0,:] = 1
        wwtp_ds[vn] = (dims, Vshape)
        wwtp_ds['river_Vshape'].attrs['long_name'] = vinfo['long_name']

        # Add position and direction
        for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            if vn == 'river_direction':
                # set point source diretion to enter vertically (2)
                wwtp_direction = 2 * np.ones(NWWTP) 
                wwtp_ds[vn] = (('river',), wwtp_direction)
            elif vn == 'river_Xposition':
                X_vec = np.nan * np.ones(NWWTP)
                for ii,wn in enumerate(gri_df_no_ovrlp.index):
                    X_vec[ii] = gri_df_no_ovrlp.loc[wn, 'col_py']
                wwtp_ds[vn] = (('river',), X_vec)
            elif vn == 'river_Eposition':
                E_vec = np.nan * np.ones(NWWTP)
                # ii = 0
                for ii,wn in enumerate(gri_df_no_ovrlp.index):
                    E_vec[ii] = gri_df_no_ovrlp.loc[wn, 'row_py']
                    # ii += 1
                wwtp_ds[vn] = (('river',), E_vec)
                    # ii += 1
                wwtp_ds[vn] = (('river',), E_vec)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']

        # Add transport
        vn = 'river_transport'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('river',)
        Q_mat = np.zeros((NT, NWWTP))
        rr = 0
        for rn in gri_df_no_ovrlp.index:
            # sum flowrates together if duplicate point sources
            if '+' in rn:
                # split into individual point sources
                [wwtp1,wwtp2] = rn.split('+')
                # get individual point source flowrates
                qtbio_wwtp_df_1 = qtbio_wwtp_df_dict[wwtp1]
                qtbio_wwtp_df_2 = qtbio_wwtp_df_dict[wwtp2]
                flow1 = qtbio_wwtp_df_1['flow'].values
                flow2 = qtbio_wwtp_df_2['flow'].values
                # combine point source flow
                flow = flow1 + flow2
            else:
                qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
                flow = qtbio_wwtp_df['flow'].values
            Q_mat[:,rr] = flow
            rr += 1
        wwtp_ds[vn] = (dims, Q_mat)
        wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
        wwtp_ds[vn].attrs['units'] = vinfo['units']

        # Add salinity and temperature
        for vn in ['river_salt', 'river_temp']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            if vn == 'river_salt':
                TS_mat = np.zeros((NT, N, NWWTP))
            elif vn == 'river_temp':
                TS_mat = np.nan * np.zeros((NT, N, NWWTP))
                rr = 0
                for rn in gri_df_no_ovrlp.index:
                    if '+' in rn:
                        # split into individual point sources
                        [wwtp1,wwtp2] = rn.split('+')
                        # get individual point source dataframe
                        qtbio_wwtp_df_1 = qtbio_wwtp_df_dict[wwtp1]
                        qtbio_wwtp_df_2 = qtbio_wwtp_df_dict[wwtp2]
                        # calculate weighted average
                        temps = trapsfun.weighted_average('temp',qtbio_wwtp_df_1, qtbio_wwtp_df_2)
                    else:
                        qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
                        temps = qtbio_wwtp_df['temp'].values
                    for nn in range(N):
                        TS_mat[:, nn, rr] = temps
                    rr += 1
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in point source river_temp!')
                sys.exit()
            wwtp_ds[vn] = (dims, TS_mat)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
            wwtp_ds[vn].attrs['units'] = vinfo['units']

        # Add biology that have existing climatology
        for var in ['NO3', 'NH4', 'TIC', 'TAlk', 'Oxyg']:
            vn = 'river_' + var
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NWWTP))
            rr = 0
            for rn in gri_df_no_ovrlp.index:
                if '+' in rn:
                    # split into individual point sources
                    [wwtp1,wwtp2] = rn.split('+')
                    # get individual point source dataframe
                    qtbio_wwtp_df_1 = qtbio_wwtp_df_dict[wwtp1]
                    qtbio_wwtp_df_2 = qtbio_wwtp_df_dict[wwtp2]
                    # calculate weighted average
                    bvals = trapsfun.weighted_average(var,qtbio_wwtp_df_1, qtbio_wwtp_df_2)
                else:
                    qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
                    bvals = qtbio_wwtp_df[var].values
                for nn in range(N):
                    B_mat[:, nn, rr] = bvals
                rr += 1
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in tiny river bio!')
                sys.exit()
            wwtp_ds[vn] = (dims, B_mat)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
            wwtp_ds[vn].attrs['units'] = vinfo['units']

        # Add remaining biology (see the lineup near the end of fennel_var.h)
        # I'm pretty sure this is simply filling everything with zeros
        bvn_list = ['Phyt', 'Zoop', 'LDeN', 'SDeN', 'Chlo', 'LDeC', 'SDeC']
        for bvn in bvn_list:
            vn = 'river_' + bvn
            vinfo = zrfun.get_varinfo(vn)
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NWWTP))
            rr = 0
            for rn in gri_df_no_ovrlp.index:
                # qtbio_wwtp_df = qtbio_wwtp_df_dict[rn]
                for nn in range(N):
                    B_mat[:, nn, rr] = rivfun.get_bio_vec(bvn, rn, yd_ind)
                rr += 1
            if np.isnan(B_mat).any():
                print('Error from traps: nans in B_mat for tiny river ' + vn)
                sys.exit()
            wwtp_ds[vn] = (dims, B_mat)
            wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
            wwtp_ds[vn].attrs['units'] = vinfo['units']

    return wwtp_ds, NWWTP
