
from datetime import datetime, timedelta
from lo_tools import forcing_argfun2 as ffun
import sys
import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd
import rivfun
import trapsfun

def make_forcing(N,NT,NRIV,dt_ind, yd_ind,ot_vec,Ldir,enable):
    # Start Dataset
    triv_ds = xr.Dataset()

    # get the list of rivers and indices for this grid
    gri_fn = Ldir['grid'] / 'triv_info.csv'
    gri_df = pd.read_csv(gri_fn, index_col='rname')

    NTRIV = 0

    if enable == True:

        # define directory for tiny river climatology
        tri_dir = Ldir['LOo'] / 'pre' / 'traps' / 'tiny_rivers'
        traps_type = 'triv'

        # climatological data files
        year0 = 1999
        year1 = 2017
        # climatological data
        Ldir['Cflow_triv_fn'] = tri_dir / 'Data_historical' / ('CLIM_flow_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['Ctemp_triv_fn'] = tri_dir / 'Data_historical' / ('CLIM_temp_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CDO_triv_fn']   = tri_dir / 'Data_historical' / ('CLIM_DO_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CNH4_triv_fn']  = tri_dir / 'Data_historical' / ('CLIM_NH4_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CNO3_triv_fn']  = tri_dir / 'Data_historical' / ('CLIM_NO3_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CTalk_triv_fn'] = tri_dir / 'Data_historical' / ('CLIM_Talk_' + str(year0) + '_' + str(year1) + '.p')
        Ldir['CTIC_triv_fn']  = tri_dir / 'Data_historical' / ('CLIM_TIC_' + str(year0) + '_' + str(year1) + '.p')

        # get the list of rivers and indices for this grid
        gri_fn = Ldir['grid'] / 'triv_info.csv'
        gri_df = pd.read_csv(gri_fn, index_col='rname')
        if Ldir['testing']:
            gri_df = gri_df.loc[['Birch Bay', 'Purdy Cr', 'Burley Cr', 'Perry Cr','McLane Cr'],:]

        # get list of overlapping rivers
        overlapping_trivs = gri_df[gri_df.duplicated(['row_py','col_py'], keep=False) == True].index.values
        # consolidate overlapping rivers
        combined_names = trapsfun.combine_adjacent(overlapping_trivs)
        gri_df_no_ovrlp = pd.DataFrame(columns=gri_df.columns) 
        gri_df_no_ovrlp.index.name='rname'
        for trname in gri_df.index: # loop through original dataframe
            if trname in overlapping_trivs: # look for rivers that are in the list of duplicates
                name_index = np.where(overlapping_trivs == trname)[0][0] # get index in the list of duplicates
                if name_index%2 == 0: # even index means first occurence of duplicate
                    newname = combined_names[int(name_index/2)] # combine names of duplicates
                    gri_df_no_ovrlp.loc[newname] = gri_df.loc[trname] # add combined source to dataframe
                # Note: second duplicate will be dropped (so idir, isign, and uv will come from the first duplicate)
            else:
                gri_df_no_ovrlp.loc[trname] = gri_df.loc[trname] # if not a duplicate, then just copy over original info

        NTRIV = len(gri_df_no_ovrlp)
        # NTRIV = len(gri_df)

        # get the flow, temperature, and nutrient data for these days
        qtbio_triv_df_dict = trapsfun.get_qtbio(gri_df, dt_ind, yd_ind, Ldir, traps_type)

        # Add time coordinate
        triv_ds['river_time'] = (('river_time',), ot_vec)
        triv_ds['river_time'].attrs['units'] = Lfun.roms_time_units
        triv_ds['river_time'].attrs['long_name'] = 'river time'

        # Add river coordinate
        triv_ds['river'] = (('river',), np.arange(NRIV+1,NRIV+NTRIV+1))
        triv_ds['river'].attrs['long_name'] = 'tiny river runoff identification number'

        # Add Vshape
        vn = 'river_Vshape'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = ('s_rho', 'river')
        # For Vtransform = 2, even spacing is a good approximation, and
        # we implement this by using 1/N as the fraction in each vertical cell.
        Vshape = (1/N) * np.ones((N, NTRIV))
        triv_ds[vn] = (dims, Vshape)
        triv_ds[vn].attrs['long_name'] = vinfo['long_name']

        # Add position and direction
        for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            if vn == 'river_direction':
                triv_ds[vn] = (('river',), gri_df_no_ovrlp.idir.to_numpy())
            elif vn == 'river_Xposition':
                X_vec = np.nan * np.ones(NTRIV)
                ii = 0
                for rn in gri_df_no_ovrlp.index:
                    if gri_df_no_ovrlp.loc[rn, 'idir'] == 0:
                        X_vec[ii] = gri_df_no_ovrlp.loc[rn, 'col_py'] + 1
                    elif gri_df_no_ovrlp.loc[rn, 'idir'] == 1:
                        X_vec[ii] = gri_df_no_ovrlp.loc[rn, 'col_py']
                    ii += 1
                triv_ds[vn] = (('river',), X_vec)
            elif vn == 'river_Eposition':
                E_vec = np.nan * np.ones(NTRIV)
                ii = 0
                for rn in gri_df_no_ovrlp.index:
                    if gri_df_no_ovrlp.loc[rn, 'idir'] == 0:
                        E_vec[ii] = gri_df_no_ovrlp.loc[rn, 'row_py']
                    elif gri_df_no_ovrlp.loc[rn, 'idir'] == 1:
                        E_vec[ii] = gri_df_no_ovrlp.loc[rn, 'row_py'] + 1
                    ii += 1
                triv_ds[vn] = (('river',), E_vec)
            triv_ds[vn].attrs['long_name'] = vinfo['long_name']

        # Add transport
        vn = 'river_transport'
        vinfo = zrfun.get_varinfo(vn, vartype='climatology')
        dims = (vinfo['time'],) + ('river',)
        Q_mat = np.zeros((NT, NTRIV))
        rr = 0
        for rn in gri_df_no_ovrlp.index:
            # sum flowrates together if duplicate river
            if '+' in rn:
                # split into individual rivers
                [riv1,riv2] = rn.split('+')
                # get individual river flowrates
                qtbio_triv_df_1 = qtbio_triv_df_dict[riv1]
                qtbio_triv_df_2 = qtbio_triv_df_dict[riv2]
                flow1 = qtbio_triv_df_1['flow'].values
                flow2 = qtbio_triv_df_2['flow'].values
                # combine river flow
                flow = flow1 + flow2
            else:
                qtbio_triv_df = qtbio_triv_df_dict[rn]
                flow = qtbio_triv_df['flow'].values
            Q_mat[:,rr] = flow * gri_df_no_ovrlp.loc[rn, 'isign']
            rr += 1
        triv_ds[vn] = (dims, Q_mat)
        triv_ds[vn].attrs['long_name'] = vinfo['long_name']
        triv_ds[vn].attrs['units'] = vinfo['units']

        # Add salinity and temperature
        for vn in ['river_salt', 'river_temp']:
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            if vn == 'river_salt':
                TS_mat = np.zeros((NT, N, NTRIV))
            elif vn == 'river_temp':
                TS_mat = np.nan * np.zeros((NT, N, NTRIV))
                rr = 0
                for rn in gri_df_no_ovrlp.index:
                    if '+' in rn:
                        # split into individual rivers
                        [riv1,riv2] = rn.split('+')
                        # get individual river dataframe
                        qtbio_triv_df_1 = qtbio_triv_df_dict[riv1]
                        qtbio_triv_df_2 = qtbio_triv_df_dict[riv2]
                        # calculate weighted average
                        temps = trapsfun.weighted_average('temp',qtbio_triv_df_1, qtbio_triv_df_2)
                    else:
                        qtbio_triv_df = qtbio_triv_df_dict[rn]
                        temps = qtbio_triv_df['temp'].values
                    for nn in range(N):
                        TS_mat[:, nn, rr] = temps
                    rr += 1
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in tiny river river_temp!')
                sys.exit()
            triv_ds[vn] = (dims, TS_mat)
            triv_ds[vn].attrs['long_name'] = vinfo['long_name']
            triv_ds[vn].attrs['units'] = vinfo['units']

        # Add biology that have existing climatology
        for var in ['NO3', 'NH4', 'TIC', 'TAlk', 'Oxyg']:
            vn = 'river_' + var
            vinfo = zrfun.get_varinfo(vn, vartype='climatology')
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NTRIV))
            rr = 0
            for rn in gri_df_no_ovrlp.index:
                if '+' in rn:
                    # split into individual rivers
                    [riv1,riv2] = rn.split('+')
                    # get individual river dataframe
                    qtbio_triv_df_1 = qtbio_triv_df_dict[riv1]
                    qtbio_triv_df_2 = qtbio_triv_df_dict[riv2]
                    # calculate weighted average
                    bvals = trapsfun.weighted_average(var,qtbio_triv_df_1, qtbio_triv_df_2)
                else:
                    qtbio_triv_df = qtbio_triv_df_dict[rn]
                    bvals = qtbio_triv_df[var].values
                for nn in range(N):
                    B_mat[:, nn, rr] = bvals
                rr += 1
            if np.isnan(TS_mat).any():
                print('Error from traps: nans in tiny river bio!')
                sys.exit()
            triv_ds[vn] = (dims, B_mat)
            triv_ds[vn].attrs['long_name'] = vinfo['long_name']
            triv_ds[vn].attrs['units'] = vinfo['units']

        # Add remaining biology (see the lineup near the end of fennel_var.h)
        # I'm pretty sure this is simply filling everything with zeros
        bvn_list = ['Phyt', 'Zoop', 'LDeN', 'SDeN', 'Chlo', 'LDeC', 'SDeC']
        for bvn in bvn_list:
            vn = 'river_' + bvn
            vinfo = zrfun.get_varinfo(vn)
            dims = (vinfo['time'],) + ('s_rho', 'river')
            B_mat = np.nan * np.zeros((NT, N, NTRIV))
            rr = 0
            for rn in gri_df_no_ovrlp.index:
                # qtbio_triv_df = qtbio_triv_df_dict[rn]
                for nn in range(N):
                    B_mat[:, nn, rr] = rivfun.get_bio_vec(bvn, rn, yd_ind)
                rr += 1
            if np.isnan(B_mat).any():
                print('Error from traps: nans in B_mat for tiny river ' + vn)
                sys.exit()
            triv_ds[vn] = (dims, B_mat)
            triv_ds[vn].attrs['long_name'] = vinfo['long_name']
            triv_ds[vn].attrs['units'] = vinfo['units']

        # Rename rivers that share name with WWTP. This code appends ' R' at the end of the river name
        duplicates = ['Port Angeles', 'Port Townsend', 'Birch Bay', 'Port Gamble', 'Gig Harbor']
        print(gri_df_no_ovrlp.index)
        gri_df_no_ovrlp.index = np.where(gri_df_no_ovrlp.index.isin(duplicates), gri_df_no_ovrlp.index + ' R', gri_df_no_ovrlp.index)
        print(gri_df_no_ovrlp.index) 

        # Add river names
        triv_ds['river_name'] = (('river',), list(gri_df_no_ovrlp.index))
        triv_ds['river_name'].attrs['long_name'] = 'tiny river name'

    return triv_ds, NTRIV