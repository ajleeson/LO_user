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
gridname = 'hc_al'

# default s-coordinate info (could override below)
s_dict = {'THETA_S': 4, 'THETA_B': 2, 'TCLINE': 10, 'N': 30,
        'VTRANSFORM': 2, 'VSTRETCHING': 4}

def make_initial_info(gridname=gridname):
    # Add an elif section for your grid.

    if gridname == 'test0':
        # A large grid, used as a test.
        dch = gfun.default_choices()
        aa = [-130, -122, 42, 52]
        res = 1000 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north','south','west']
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['srtm15plus','cascadia','nw_pacific','psdem',
               'ttp_patch','grays_harbor','willapa_bay']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'hc0':
        dch = gfun.default_choices()
        aa = [-123.2, -122.537, 47.3, 47.9]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['psdem']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'ai0':
        dch = gfun.default_choices()
        aa = [-122.82, -122.36, 47.758, 48.18]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north', 'south', 'east', 'west']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['psdem']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'so1':
        # South Sound, new version, 2023.04.12
        dch = gfun.default_choices()
        dch['z_offset'] = -1.3 # NAVD88 is 1.3 m below MSL at Seattle
        dch['excluded_rivers'] = ['skokomish']
        aa = [-123.13, -122.76, 47, 47.42]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['east']
        dch['nudging_days'] = (0.1, 1.0)
        
        # by setting a small min_depth were are planning to use
        # wetting and drying in ROMS, but maintaining positive depth
        # for all water cells
        dch['min_depth'] = 0.2 # meters (positive down)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['psdem']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
                    
    elif gridname == 'wgh2':
        # Willapa Bay and Grays Harbor nest
        dch = gfun.default_choices()
        aa = [-124.4,-123.7,46.35,47.1]
        res = 200 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        
        dch['z_offset'] = -2
        # The docs for nw_pacific say the vertical datum is "sea level" and for Willapa
        # Bay and Grays Harbor it is MLLW so to match
        # this we would use z_offset = 0 or -1, but the intention here is to make the z=0
        # level be higher up, so that we catch more of the intertidal when using
        # WET_DRY. This should be matched by a similar intervention to zeta in ocnN.
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['nudging_days'] = (0.1, 1.0)
        
        # by setting a small min_depth were are planning to use
        # WET_DRY in ROMS, but maintaining positive depth
        # for all water cells
        dch['min_depth'] = 0.2 # meters (positive down)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['nw_pacific','grays_harbor','willapa_bay']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'cas2k':
        # cas6 domain but with 2 km resolution
        dch = gfun.default_choices()
        aa = [-130, -122, 42, 52]
        res = 2000 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['nudging_days'] = (3.0, 60.0)
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['srtm15plus','cascadia','nw_pacific','psdem']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'cas7':
        # based completely on cas6 except we carve out Agate Pass and
        # Swinomish Channel by hand. This is an example of working from an
        # existing grid.
        dch = gfun.default_choices()
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['nudging_days'] = (3.0, 60.0)
        Ldir = Lfun.Lstart()
        fn = Ldir['parent'] / 'LO_output' / 'pgrid' / 'cas6' / 'grid.nc'
        dch['maskfile_to_copy'] = fn
        dch['remove_islands'] = False
        dch['trim_grid'] = False

        import xarray as xr
        ds = xr.open_dataset(fn)
        z = -ds.h.values
        lon = ds.lon_rho.values
        lat = ds.lat_rho.values
        
        # The plan is to only run:
        # start_grid
        # make_mask
        # edit_mask
        # (don't run carve_rivers - just copy the file from cas6)
        # smooth_grid
        # make_extras
        # grid_to_LO
            
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
        
        # create a river file by hand
        Ldir = Lfun.Lstart()
        dch['ctag'] = 'ae0_v0'
        ri_dir = Ldir['LOo'] / 'pre' / 'river1' / dch['ctag']
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
        if True:
            track_df['lon'] = np.linspace(0,4,NTR) # OK to go past edge of domain
            track_df['lat'] = 45*np.ones(NTR)
        else: # Debugging with N/S river channel
            track_df['lon'] = 0.25*np.ones(NTR)
            track_df['lat'] = np.linspace(45,44,NTR) # South to North river
        track_df.to_pickle(track_fn)
        # *** NOTE: TRACKS MUST GO FROM OCEAN TO LAND ***

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
        

    # my idalized hood canal estuary
    elif gridname == 'hc_al':
        # get list of default choices
        dch = gfun.default_choices()

        # 200 m resolution inside of Hood Canal, and 2500 far away
        # lon_list = [-124.1, -122.7, -122.5, -121.1]
        lon_list = [-123.5, -122.7, -122.5, -121.7]
        # x_res_list = [2500, 200, 200, 2500]
        x_res_list = [3000, 200, 200, 3000]
        lat_list = [47, 48.1, 48.5]
        # lat_list = [47, 48.2, 49]
        y_res_list = [500, 500, 3000]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                            lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # defining bathymetry analytically
        dch['analytical'] = True
        # I want my river to come from the south and ocean to come from the north
        dch['nudging_edges'] = ['north', 'east', 'west']
        # Let mean sea level equal NAVD88
        dch['use_z_offset'] = False
        # tidy up dch (make sure there is no carryover from prior runs?)
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        # make bathymetry by hand
        z = np.zeros(lon.shape)
        x, y = zfun.ll2xy(lon, lat, -122.6, 48)

        # parabolic shelf
        y0 = 0
        z0 = -50
        y1 = np.max(y)
        z1 = -200
        a = (z0-z1)/((y0**2-y1**2)-2*y1*(y0-y1))
        b = -2*a*y1
        c = z0 - a*y0**2 - b*y0
        zshelf = a*y**2 + b*y + c

        print('domain size')
        print('x = {} km'.format((np.max(x)-np.min(x))/1000))
        print('y = {} km'.format((np.max(y)-np.min(y))/1000))

        # rectangular cross section, and slope that is a 6th order polynomial fit to hood canal-ish bathy
        # coefficients
        A = 1.35760755e-26
        B = 3.58266838e-21
        C = 3.12590423e-16
        D = 8.74304619e-12
        E = -2.08645813e-08
        F = 4.33905852e-05
        G = -50
        # set x and to arbitrarily large value outside of estuary so it will get masked out
        x_temp = x
        y_temp = y
        # hood canal is 3 km wide
        x_temp[np.abs(x_temp)>1500] = 1e7
        y_temp[y>0] = -1e3
        zestuary = (A*y**6 + B*y**5 + C*y**4 + D*y**3 + E*y**2 + F*y + G) +\
                    80*((x_temp)/5E3)**2


        z = zshelf
        mask = zestuary < z
        z[mask] = zestuary[mask]

        # create a river file by hand
        Ldir = Lfun.Lstart()
        # update ctag for river
        dch['ctag'] = 'hc_al_rivs'
        ri_dir = Ldir['LOo'] / 'pre' / 'river1' / dch['ctag']
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
        if True:
            track_df['lon'] = -122.6*np.ones(NTR)
            track_df['lat'] = np.linspace(47.12,46.9,NTR) # OK to go past edge of domain
        track_df.to_pickle(track_fn)
        # *** NOTE: TRACKS MUST GO FROM OCEAN TO LAND ***

    else:
        print('Error from make_initial_info: unsupported gridname')
        return
        
    if dch['trim_grid']:
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
    

