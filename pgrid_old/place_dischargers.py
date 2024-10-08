"""
This script reads in lat/lon wwtp coordinates and places a point source at the nearest coastal grid cell.
"""

import gfun
Gr =gfun.gstart()

from lo_tools import zfun
from lo_tools import plotting_functions as pfun

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seawater as sw
import pickle

# load the default choices
# dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))

# get discharger info
wwtp_fn = Gr['wwtp_dir'] / 'wwtp_loc_info.csv'
wwtp_df = pd.read_csv(wwtp_fn, index_col='dname')

# select and increment grid file (the new tag is called _d for dischargers)
in_fn = gfun.select_file(Gr)
out_fn = gfun.increment_filename(in_fn, '_d')

# get the grid data
ds = xr.open_dataset(in_fn)
z = -ds.h.values
mask_rho = np.transpose(ds.mask_rho.values)
lon = ds.lon_rho.values
lat = ds.lat_rho.values

X = lon[0,:] # grid cell X values
Y = lat[:,0] # grid cell Y values

def in_domain(x, y, X, Y):
    # Utility function to make sure that a point (x, y) is
    # in a domain specified by vectors X and Y.
    # We actually require the point to be 'pad' in from the edge.
    pad = 1
    if x>=X[0+pad] and x<=X[-1-pad] and y>=Y[0+pad] and y<=Y[-1-pad]:
        return True
    else:
        return False

def cell_in_domain(ival,jval,II,JJ):
    # Utility function to make sure that a grid cell with index location (ival,jval)
    # in the domain specified by domain sizes II and JJ.
    if ival>=0 and ival<=II and jval>=0 and jval<=JJ:
        return True
    else:
        return False

def get_cell_info(I_ind,J_ind,X,Y,x,y):
    # function that checks whether a grid cell is in the domain and if it is on water.
    # If both, then records the index of the grid cell, and the distance from
    # the cell to the wwtp
    ii_list = []
    jj_list = []
    distance_list = []
    # if cell is in domain and it is a water cell (mask = 0),
    # then record information
    if cell_in_domain(I_ind,J_ind,len(X)-1,len(Y)-1):
        if mask_rho[I_ind,J_ind]:
            # record cell indices
            II_list.append(I_ind)
            JJ_list.append(J_ind)
            # record distance from center of cell to wwtp location
            xmeters, ymeters = zfun.ll2xy(X[I_ind], Y[J_ind], x, y)
            distance = np.sqrt(xmeters**2 + ymeters**2)
            distance_list.append(distance)
    return ii_list, jj_list, distance_list


# prepare things for looping over the dischargers
ilon_dict = dict()
ilat_dict = dict()
dir_dict = dict()
good_wwtp = []
roms_info_df = pd.DataFrame()

# loop over dischargers
for wn in wwtp_df.index:

    try:
        # get lon and lat (x and y) coordinate of wwtp
        x = wwtp_df._get_value(wn,'lon')
        y = wwtp_df._get_value(wn,'lat')
    
        if in_domain(x, y, X, Y):
            # we only consider dischargers in the domain
            good_wwtp.append(wn)
            print('including ' + wn.title())

            # initialize a boolean to track whether a grid cell has been selected for point source
            ps_located = False

            # figure out in which grid cell the wwtp is located
            ix = zfun.find_nearest_ind(X,x)
            iy = zfun.find_nearest_ind(Y,y)
            wwtpII = np.array(ix)
            wwtpJJ = np.array(iy)

            # search for best grid cell to place point source -----------------------
        
            # First, check if the wwtp is located in the ocean
            if mask_rho[wwtpII,wwtpJJ]: # mask of 1 means water
                # point source located at same location as wwtp
                psII = wwtpII
                psJJ = wwtpJJ
                ps_located = True
                print('{}:in water, PS lat/lon:({},{})'.format(wn,round(X[psII],2),round(Y[psJJ],2)))

            # Otherwise, search in square rings around wwtp for the nearest coastal cell
            else:
                # initialize counters and arrays
                ringnum = 1
                II_list = []
                JJ_list = []
                dist_list = []
                min_dist = 1e9 # arbitratily large number so we can update with min value later

                while not ps_located:

                    # record all water cells that are in top row of ring
                    for iter in range(-ringnum,ringnum+1):
                        I_ind = wwtpII + iter
                        J_ind = wwtpJJ + ringnum
                        # get cell info if it is in water
                        ii_list, jj_list, distance_list = get_cell_info(I_ind,J_ind,X,Y,x,y)
                        # add the cell info to arrays if in water
                        II_list = II_list + ii_list
                        JJ_list = JJ_list + jj_list
                        dist_list = dist_list + distance_list

                    # record all water cells that are in bottom row of ring
                    for iter in range(-ringnum,ringnum+1):
                        I_ind = wwtpII + iter
                        J_ind = wwtpJJ - ringnum
                        # get cell info if it is in water
                        ii_list, jj_list, distance_list = get_cell_info(I_ind,J_ind,X,Y,x,y)
                        # add the cell info to arrays if in water
                        II_list = II_list + ii_list
                        JJ_list = JJ_list + jj_list
                        dist_list = dist_list + distance_list

                    # record all water cells that are in left column of ring (exclude diagonals)
                    for iter in range(-(ringnum-1),ringnum):
                        J_ind = wwtpJJ + iter
                        I_ind = wwtpII - ringnum
                        # get cell info if it is in water
                        ii_list, jj_list, distance_list = get_cell_info(I_ind,J_ind,X,Y,x,y)
                        # add the cell info to arrays if in water
                        II_list = II_list + ii_list
                        JJ_list = JJ_list + jj_list
                        dist_list = dist_list + distance_list

                    # record all water cells that are in right column of ring (exclude diagonals)
                    for iter in range(-(ringnum-1),ringnum):
                        J_ind = wwtpJJ + iter
                        I_ind = wwtpII + ringnum
                        # get cell info if it is in water
                        ii_list, jj_list, distance_list = get_cell_info(I_ind,J_ind,X,Y,x,y)
                        # add the cell info to arrays if in water
                        II_list = II_list + ii_list
                        JJ_list = JJ_list + jj_list
                        dist_list = dist_list + distance_list

                    # calculate minimum distance coastal cell (if there were any cells in water)
                    if len(dist_list) > 0:
                        # minimum distance
                        min_dist_new = np.min(dist_list)
                        # index of minimum distance
                        min_dist_ind = dist_list.index(min_dist_new)
                        # corresponding i an j indices
                        min_II = II_list[min_dist_ind]
                        min_JJ = JJ_list[min_dist_ind]

                        # check if new min_dist is smaller than previous min_dist
                        if min_dist_new < min_dist:
                            # save new min dist and grid cell index
                            min_dist = min_dist_new
                            psII = min_II
                            psJJ = min_JJ

                            # check next ring if there is a closer coastal grid cell
                            ringnum +=1

                        # if new min_dist is larger than before, then we have found the nearest coastal grid cell
                        else:
                            ps_located = True
                            print('{}:on land, PS lat/lon:({},{})'.format(wn,round(X[psII],2),round(Y[psJJ],2)))

                    # if all cells in ring are on land, then go into next ring
                    else:
                        # iterate into next level ring
                        ringnum += 1

            # NOTE: the row_py and col_py that we write to river_info.csv use
            # python indexing conventions, meaning that the rho, u, and v grids
            # start at [0,0].  Hence when we plot things, as in plot_grid.py, we
            # should be able to use [row_py, col_py] directly with the u and v grids.
            #
            # However, when we later create the source positions for running ROMS,
            # e.g. in forcing/riv0/make_forcing_main.py => rivers.nc, we convert these
            # values to ROMS indexing conventions, meaning:
            # - add 1 to row of N/S sources
            # - add 1 to column of E/W sources.

            # save location of point source
            roms_info_df.loc[wn, 'row_py'] = psJJ 
            roms_info_df.loc[wn, 'col_py'] = psII 
        
        else:
            print(' >> excluding ' + wn.title())
        
    except FileNotFoundError:
        pass


# In reality, we don't need to save a new grid file
# but it helps to keep track of what steps have already been run
print('Creating ' + str(out_fn))
ds.to_netcdf(out_fn)
ds.close()

# save the roms point source info
out_wfn = Gr['gdir'] / 'roms_wwtp_info.csv'
print('\nCreating ' + str(out_wfn))
roms_info_df.index.name = 'wname'
roms_info_df.to_csv(out_wfn)

