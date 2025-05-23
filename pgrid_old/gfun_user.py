"""
User-specific code for pgrid.

You would edit the information to reflect whatever grid you are working on.
"""
import numpy as np
import pandas as pd
from lo_tools import zfun, Lfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

# This is the name of the grid that you are working on.
gridname = 'alpe2'

# default s-coordinate info (could override below)
s_dict = {'THETA_S': 4, 'THETA_B': 2, 'TCLINE': 10, 'N': 30,
        'VTRANSFORM': 2, 'VSTRETCHING': 4}

# Set the gridname and tag to use when creating the Ldir paths.
# They are used for accessing the river tracks, which may be developed for one
# grid but reused in others.
if gridname in ['ai0','hc0', 'sal0', 'so0']:
    # these cases reuse (all or some of) the LiveOcean cas6 model rivers
    base_gridname = 'cas6'
    base_tag = 'v3'
elif gridname in ['ae0']:
    # for analytical cases we create the river info and track in
    # make_initial_info() below, but we still assign gridname and tag so
    # that they get saved in the right places
    base_gridname = 'ae0'
    base_tag = 'v0'
elif gridname in ['al0']:
    # AL's first analytical model estuary
    # for analytical cases we create the river info and track in
    # make_initial_info() below, but we still assign gridname and tag so
    # that they get saved in the right places
    base_gridname = 'al0'
    base_tag = 'v0'
elif gridname in ['alpe']:
    # parabolic cross section estuary, but otherwise straight
    # for analytical cases we create the river info and track in
    # make_initial_info() below, but we still assign gridname and tag so
    # that they get saved in the right places
    base_gridname = 'alpe'
    base_tag = 'v0'
elif gridname in ['alpe2']:
    # a shortened version of alpe
    base_gridname = 'alpe2'
    base_tag = 'v0'
elif gridname in ['fsg']:
    # flat shelf grid
    base_gridname = 'fsg'
    base_tag = 'v0'
elif gridname in ['awwtp1']:
    # flat shelf grid
    base_gridname = 'awwtp1'
    base_tag = 'v0'

def make_initial_info(gridname=gridname):
    # Add an elif section for your grid.

    if gridname == 'sal0':
        # A Salish Sea grid, used as an example.
        dch = gfun.default_choices()
        aa = [-124, -122, 47, 49]
        res = 600 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north', 'west']
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'hc0':
        dch = gfun.default_choices()
        aa = [-123.2, -122.537, 47.3, 47.9]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['t_list'] = [dch['t_dir'] / 'psdem' / 'PS_27m.nc']
        dch['nudging_edges'] = ['north']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'ai0':
        dch = gfun.default_choices()
        aa = [-122.82, -122.36, 47.758, 48.18]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['t_list'] = [dch['t_dir'] / 'psdem_10m' / 'PS_30m.nc']
        dch['nudging_edges'] = ['north', 'south', 'east', 'west']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'so0':
        # South Sound
        dch = gfun.default_choices()
        dch['z_offset'] = -1.3 # NAVD88 is 1.3 m below MSL at Seattle
        dch['excluded_rivers'] = ['skokomish']
        aa = [-123.13, -122.76, 47, 47.42]
        res = 50 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['t_list'] = [dch['t_dir'] / 'srtm15' / 'topo15.nc',
                    dch['t_dir'] / 'psdem_10m' / 'PS_30m.nc']
        dch['nudging_edges'] = ['east']
        dch['nudging_days'] = (0.1, 1.0)
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'ae0':
        # analytical model estuary
        dch = gfun.default_choices()
        lon_list = [-2, 0, 1, 2]
        x_res_list = [2500, 500, 500, 2500]
        lat_list = [43, 44.9, 45.1, 47]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                            lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        dch['analytical'] = True
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['use_z_offset'] = False
        # tidy up dch
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        # make bathymetry by hand
        z = np.zeros(lon.shape)
        x, y = zfun.ll2xy(lon, lat, 0, 45)
        zshelf = x * 1e-3
        zestuary = -20 + 20*x/1e5 + 20/(1e4)*np.abs(y)
        z = zshelf
        mask = zestuary < z
        z[mask] = zestuary[mask]
        
        # create a river file
        Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)
        ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
        Lfun.make_dir(ri_dir)
        gri_fn = ri_dir / 'river_info.csv'
        with open(gri_fn, 'w') as rf:
            rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
            rf.write('creek0,,,,1.0,5.0,m3/s,degC\n')
        # and make a track for the river
        track_dir = ri_dir / 'tracks'
        Lfun.make_dir(track_dir)
        track_fn = track_dir / 'creek0.p'
        track_df = pd.DataFrame()
        NTR = 100
        track_df['lon'] = np.linspace(0,4,NTR) # OK to go past edge of domain
        track_df['lat'] = 45*np.ones(NTR)
        track_df.to_pickle(track_fn)
        # NOTE: tracks go from ocean to land

    elif gridname == 'al0':
        # get list of default choices
        dch = gfun.default_choices()

        # So far I have left this unchanged because I don't understand it
        lon_list = [-2, 0, 1, 2]
        x_res_list = [2500, 500, 500, 2500]
        lat_list = [43, 44.9, 45.1, 47]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                            lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # defining bathymetry analytically
        dch['analytical'] = True
        # I want my river to come from the north and ocean to come from the south
        dch['nudging_edges'] = ['south', 'east', 'west']
        # Let mean sea level equal NAVD88
        dch['use_z_offset'] = False
        # tidy up dch (make sure there is no carryover from prior runs?)
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        # make bathymetry by hand
        z = np.zeros(lon.shape)
        x, y = zfun.ll2xy(lon, lat, 0, 45)

        # sigmoid-shaped shelf
        zshelf = 240*(1/(1+np.exp(-y/5e4))) - 130

        # sinusoid-shaped estuary (amplitude scaled by chi-squared distribution)
        # with parabolic cross section, and slope that drops linearly
        k = 4
        gamma = 0.2
        yscale = 1e4
        zestuary = -20 + 3*(x/5e3 -\
        (((y/yscale)**(k/2-1)*np.exp(-(y/yscale)/2)/(2**(k/2)*gamma)))*np.sin(y/2e3))**2 + 25*y/1e5

        z = zshelf
        mask = zestuary < z
        z[mask] = zestuary[mask]
        
        # create a river file
        Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)
        ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
        Lfun.make_dir(ri_dir)
        gri_fn = ri_dir / 'river_info.csv'
        with open(gri_fn, 'w') as rf:
            rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
            rf.write('creek0,,,,1.0,5.0,m3/s,degC\n')
        # and make a track for the river
        track_dir = ri_dir / 'tracks'
        Lfun.make_dir(track_dir)
        track_fn = track_dir / 'creek0.p'
        track_df = pd.DataFrame()
        NTR = 100
        track_df['lat'] = np.linspace(45,50,NTR) # OK to go past edge of domain
        track_df['lon'] = 0*np.ones(NTR)
        track_df.to_pickle(track_fn)

    elif gridname == 'alpe':
        # get list of default choices
        dch = gfun.default_choices()

        # So far I have left this unchanged because I don't understand it
        lon_list = [-2, 0, 1, 2]
        x_res_list = [2500, 500, 500, 2500]
        lat_list = [43, 44.9, 45.1, 47]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                            lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # defining bathymetry analytically
        dch['analytical'] = True
        # I want my river to come from the north and ocean to come from the south
        dch['nudging_edges'] = ['south', 'east', 'west']
        # Let mean sea level equal NAVD88
        dch['use_z_offset'] = False
        # tidy up dch (make sure there is no carryover from prior runs?)
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        # make bathymetry by hand
        z = np.zeros(lon.shape)
        x, y = zfun.ll2xy(lon, lat, 0, 45)

        # sigmoid-shaped shelf
        zshelf = 240*(1/(1+np.exp(-y/5e4))) - 130

        # parabolic cross section, and slope that drops linearly
        zestuary = -20 + 3*((x)/5e3)**2 + 25*y/1e5

        z = zshelf
        mask = zestuary < z
        z[mask] = zestuary[mask]
        
        # create a river file
        Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)
        ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
        Lfun.make_dir(ri_dir)
        gri_fn = ri_dir / 'river_info.csv'
        with open(gri_fn, 'w') as rf:
            rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
            rf.write('creek0,,,,1.0,5.0,m3/s,degC\n')
        # and make a track for the river
        track_dir = ri_dir / 'tracks'
        Lfun.make_dir(track_dir)
        track_fn = track_dir / 'creek0.p'
        track_df = pd.DataFrame()
        NTR = 100
        track_df['lat'] = np.linspace(45,50,NTR) # OK to go past edge of domain
        track_df['lon'] = 0*np.ones(NTR)
        track_df.to_pickle(track_fn)

    elif gridname == 'alpe2':
        # get list of default choices
        dch = gfun.default_choices()

        # So far I have left this unchanged because I don't understand it
        lon_list = [-1.5, -0.5, 0.5, 1.5]
        x_res_list = [2500, 500, 500, 2500]
        lat_list = [44, 44.9, 45.1, 46]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                            lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # defining bathymetry analytically
        dch['analytical'] = True
        # I want my river to come from the north and ocean to come from the south
        dch['nudging_edges'] = ['south', 'east', 'west']
        # Let mean sea level equal NAVD88
        dch['use_z_offset'] = False
        # tidy up dch (make sure there is no carryover from prior runs?)
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        # make bathymetry by hand
        z = np.zeros(lon.shape)
        x, y = zfun.ll2xy(lon, lat, 0, 45)

        # sigmoid-shaped shelf
        zshelf = 240*(1/(1+np.exp(-y/5e4))) - 130

        # parabolic cross section, and slope that drops linearly
        zestuary = -20 + 3*((x)/5e3)**2 + 25*y/1e5

        z = zshelf
        mask = zestuary < z
        z[mask] = zestuary[mask]
        
        # create a river file
        Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)
        ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
        Lfun.make_dir(ri_dir)
        gri_fn = ri_dir / 'river_info.csv'
        with open(gri_fn, 'w') as rf:
            rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
            rf.write('creek0,,,,1.0,5.0,m3/s,degC\n')
        # and make a track for the river
        track_dir = ri_dir / 'tracks'
        Lfun.make_dir(track_dir)
        track_fn = track_dir / 'creek0.p'
        track_df = pd.DataFrame()
        NTR = 100
        track_df['lat'] = np.linspace(45,48,NTR) # OK to go past edge of domain
        track_df['lon'] = 0*np.ones(NTR)
        track_df.to_pickle(track_fn)

    elif gridname == 'fsg':
        # get list of default choices
        dch = gfun.default_choices()

        # So far I have left this unchanged because I don't understand it
        lon_list = [-1.5, -0.5, 0.5, 1.5]
        x_res_list = [2500, 500, 500, 2500]
        lat_list = [44, 44.4, 44.6, 45]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                            lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # defining bathymetry analytically
        dch['analytical'] = True
        # Land is toward the north
        dch['nudging_edges'] = ['south', 'east', 'west']
        # Let mean sea level equal NAVD88
        dch['use_z_offset'] = False
        # tidy up dch (make sure there is no carryover from prior runs?)
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        # make bathymetry by hand
        z = np.zeros(lon.shape)

        # flat bottomed shelf
        z[lat < 44.75] = -100

    elif gridname == 'awwtp1':
        # get list of default choices
        dch = gfun.default_choices()

        # So far I have left this unchanged because I don't understand it
        lon_list = [-1.5, -0.5, 0.5, 1.5]
        x_res_list = [2500, 500, 500, 2500]
        lat_list = [44, 44.9, 45.1, 46]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                            lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # defining bathymetry analytically
        dch['analytical'] = True
        # I want my river to come from the north and ocean to come from the south
        dch['nudging_edges'] = ['south', 'east', 'west']
        # Let mean sea level equal NAVD88
        dch['use_z_offset'] = False
        # tidy up dch (make sure there is no carryover from prior runs?)
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        # make bathymetry by hand
        z = np.zeros(lon.shape)
        x, y = zfun.ll2xy(lon, lat, 0, 45)

        # sigmoid-shaped shelf
        zshelf = 240*(1/(1+np.exp(-y/5e4))) - 130

        # parabolic cross section, and slope that drops linearly
        zestuary = -20 + 3*((x)/5e3)**2 + 25*y/1e5

        z = zshelf
        mask = zestuary < z
        z[mask] = zestuary[mask]
      
    else:
        print('Error from make_initial_info: unsupported gridname')
        return
        
    # check for odd size of grid and trim if needed
    NR, NC = lon.shape
    if np.mod(NR,2) != 0:
        print('- trimming row from grid')
        lon = lon[:-1,:]
        lat = lat[:-1,:]
        z = z[:-1,:]
    if np.mod(NC,2) != 0:
        print('- trimming column from grid')
        lon = lon[:,:-1]
        lat = lat[:,:-1]
        z = z[:,:-1]

    return lon, lat, z, dch
    

