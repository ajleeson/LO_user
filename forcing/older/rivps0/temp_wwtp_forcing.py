# get the list of wwtps and indices for this grid
gwwtp_fn = Ldir['grid'] / 'wwtp_info.csv'
gwwtp_df = pd.read_csv(gwwtp_fn, index_col='wname')
NWWTP = len(gwwtp_df)

# Start Dataset
wwtp_ds = xr.Dataset()

# Add time coordinate
wwtp_ds['river_time'] = (('river_time',), ot_vec)
wwtp_ds['river_time'].attrs['units'] = Lfun.roms_time_units
wwtp_ds['river_time'].attrs['long_name'] = 'river time'

# Add wwtp coordinate
wwtp_ds['river'] = (('river',), np.arange(NRIV,NRIV+NWWTP+1))
wwtp_ds['river'].attrs['long_name'] = 'point source identification number'

# Add river names
wwtp_ds['river_name'] = (('river',), list(gwwtp_fn.index))
wwtp_ds['river_name'].attrs['long_name'] = 'river name'

# Add Vshape
vn = 'river_Vshape'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = ('s_rho', 'river')
# For Vtransform = 2, even spacing is a good approximation, and
# we implement this by using 1/N as the fraction in each vertical cell.
Vshape = (1/N) * np.ones((N, NWWTP))
wwtp_ds[vn] = (dims, Vshape)
wwtp_ds['river_Vshape'].attrs['long_name'] = vinfo['long_name']

# Add position and direction
for vn in ['river_Xposition', 'river_Eposition', 'river_direction']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    if vn == 'river_direction':
        # set point source diretion to enter vertically (2)
        wwtp_dir = 2 * np.ones(NWWTP) 
        wwtp_ds[vn] = (('river',), wwtp_dir)
    elif vn == 'river_Xposition':
        X_vec = np.nan * np.ones(NWWTP)
        for ii,wn in enumerate(gwwtp_df.index):
            X_vec[ii] = gwwtp_df.loc[wn, 'col_py']
        wwtp_ds[vn] = (('river',), X_vec)
    elif vn == 'river_Eposition':
        E_vec = np.nan * np.ones(NWWTP)
        # ii = 0
        for ii,wn in enumerate(gri_df.index):
            E_vec[ii] = gwwtp_df.loc[wn, 'row_py']
            ii += 1
        wwtp_ds[vn] = (('river',), E_vec)
            # ii += 1
        wwtp_ds[vn] = (('river',), E_vec)
    wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
        
# Add transport
vn = 'river_transport'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = (vinfo['time'],) + ('river',)
Q_mat = np.zeros((NT, NWWTP))
# ii = 0
for ii,rn in enumerate(gwwtp_fn.index):
    if rn == 'user_specified_wwtp':
        Q_mat[:,ii] = 3000 * np.ones(NT)
    else:
        # wwtps all have the same flowrate
        Q_mat[:,ii] = 1000 * np.ones(NT)
    # ii += 1
wwtp_ds[vn] = (dims, Q_mat)
wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
wwtp_ds[vn].attrs['units'] = vinfo['units']

# Add salinity and temperature
for vn in ['river_salt', 'river_temp']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = (vinfo['time'],) + ('s_rho', 'river')
    if vn == 'river_salt':
        TR_mat = np.zeros((NT, N, NWWTP))
    elif vn == 'river_temp':
        TR_mat = 10 * np.ones((NT, N, NWWTP))
    wwtp_ds[vn] = (dims, TR_mat)
    wwtp_ds[vn].attrs['long_name'] = vinfo['long_name']
    wwtp_ds[vn].attrs['units'] = vinfo['units']