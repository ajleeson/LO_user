"""
tiny river and point source (traps) functions
"""

import xarray as xr
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt
import cmocean
import sys
from lo_tools import forcing_argfun2 as ffun
from lo_tools import zfun
from lo_tools import plotting_functions as pfun
from matplotlib.patches import Rectangle

Ldir = ffun.intro() # this handles all the argument passing

def in_domain(x, y, X, Y):
    '''
    Note: code borrowed from pgrid/carve_rivers.py
    Utility function to make sure that a point (x, y) is
    in a domain specified by vectors X and Y.
    We actually require the point to be 'pad' in from the edge.

    Inputs:
        x:  point of interest lon coordinates
        y:  point of interest lat coordinates
        X:  domain lon coordinates
        Y:  domain lat coordinates
    Outputs:
        boolean (T if in domain, F if not in domain)
    '''
    pad = 1
    if x>=X[0+pad] and x<=X[-1-pad] and y>=Y[0+pad] and y<=Y[-1-pad]:
        return True
    else:
        return False

def cell_in_domain(ival,jval,II,JJ):
    '''
    Utility function to make sure that a grid cell with index location (ival,jval)
    in the domain specified by domain sizes II and JJ.

    Inputs:
        ival:   grid cell of interest i-index
        jval:   grid cell of interest j-index
        II:     domain max i index
        JJ:     domain max j index
    Outputs:
        boolean (T if in domain, F if not in domain)
    '''
    if ival>=0 and ival<=II and jval>=0 and jval<=JJ:
        return True
    else:
        return False

def get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho):
    '''
    Function that checks whether a grid cell is in the model domain and if it is on water.
    If both, then records the index of the grid cell, and the distance from
    the grid cell to the original WWTP lat/lon coordinates.

    Inputs:
        I_ind:      i-index of WWTP
        J_ind:      j-index of WWTP
        X:          domain lon coordinates
        Y:          domain lat coordinates
        x:          WWTP lon coordinates
        y:          WWTP lat coordinates
        mask_rho:   used to check if grid cell is water/land

    Outputs:
        ii_list:        i-index of WWTP (only saved if cell is in model domain)
        jj_list:        j-index of WWTP (only saved if cell is in model domain)
        distance_list:  distance from grid cell to original WWTP location [m]  
    '''
    ii_list = []
    jj_list = []
    distance_list = []
    # if cell is in domain and it is a water cell (mask = 0),
    # then record information
    if cell_in_domain(I_ind,J_ind,len(X)-1,len(Y)-1):
        if mask_rho[I_ind,J_ind]:
            # record cell indices
            ii_list.append(I_ind)
            jj_list.append(J_ind)
            # record distance from center of cell to wwtp location
            xmeters, ymeters = zfun.ll2xy(X[I_ind], Y[J_ind], x, y)
            distance = np.sqrt(xmeters**2 + ymeters**2)
            distance_list.append(distance)
    return ii_list, jj_list, distance_list

def get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho):
    
    '''
    Function that check whether a grid cell is a coastal water cell.
    If it is coastal, returns updated row, col, idir, isign, and uv, 
    and the distance from the coastal cell to the river mouth. 

    River direction is decided based on the position of the closest 
    coastal land cell to the original river mouth lat/lon coordinates.

    Inputs:
        I_ind:      i-index of river mouth
        J_ind:      j-index of river mouth
        X:          domain lon coordinates
        Y:          domain lat coordinates
        x:          river mouth lon coordinates
        y:          river mouth lat coordinates
        mask_rho:   used to check if grid cell is water/land

    Outputs:
        ii_list:        i-index of river (only saved if cell is a costal cell)
        jj_list:        j-index of river (only saved if cell is a costal cell)
        distance_list:  distance from grid cell to original river mouth location [m]  
        idir_list:      0 = E-W river, 1 = N-S river
        isign_list:     1 = flowing N/E, -1 = flowing S,W
        uv_list:        'u' = E-W river, 'v' = N-S river
    '''
    ii_list = []
    jj_list = []
    distance_list = []
    idir_list = []
    isign_list = []
    uv_list = []
    # if cell is in domain and it is a water cell (mask = 1),
    # then check if it is coastal
    if cell_in_domain(I_ind,J_ind,len(X)-1,len(Y)-1):
        if mask_rho[I_ind,J_ind]: # is water cell

            # Get coordinates of adjacent cells
            # north
            NI = I_ind
            NJ = J_ind + 1
            # south
            SI = I_ind
            SJ = J_ind - 1
            # west
            WI = I_ind - 1
            WJ = J_ind
            # east
            EI = I_ind + 1
            EJ = J_ind

            # tracker if there's a coastal cell or not
            is_coastal = False
            coastal_type = []

            # check if north cell is land
            if cell_in_domain(NI,NJ,len(X)-1,len(Y)-1):
                # check if north cell is land (meaning original cell is coastal)
                if not mask_rho[NI,NJ]:
                    # original cell is coastal
                    is_coastal = True
                    coastal_type = coastal_type + ['N']
            # check if south cell is land
            if cell_in_domain(SI,SJ,len(X)-1,len(Y)-1):
                # check if south cell is land (meaning original cell is coastal)
                if not mask_rho[SI,SJ]:
                    # original cell is coastal
                    is_coastal = True
                    coastal_type = coastal_type + ['S']
            # check if west cell is land
            if cell_in_domain(WI,WJ,len(X)-1,len(Y)-1):
                # check if west cell is land (meaning original cell is coastal)
                if not mask_rho[WI,WJ]:
                    # original cell is coastal
                    is_coastal = True
                    coastal_type = coastal_type + ['W']
            # check if east cell is land
            if cell_in_domain(EI,EJ,len(X)-1,len(Y)-1):
                # check if east cell is land (meaning original cell is coastal)
                if not mask_rho[EI,EJ]:
                    # original cell is coastal
                    is_coastal = True
                    coastal_type = coastal_type + ['E']

            # write data if the cell is coastal
            if is_coastal:
                # record distance from center of cell to source location
                xmeters, ymeters = zfun.ll2xy(X[I_ind], Y[J_ind], x, y)
                distance = np.sqrt(xmeters**2 + ymeters**2)

                # Figure out which adjacent land cell to use
                if len(coastal_type) == 1: # easy if there is only one land cell
                    source_side = coastal_type[0]
                # if there's more than one land cell, pick the one that is closest to the original source
                else:
                    dist2land = []
                    if 'N' in coastal_type:
                        xdist,ydist = zfun.ll2xy(X[NI], Y[NJ], x, y)
                        dist2land = dist2land + [np.sqrt(xdist**2 + ydist**2)]
                    if 'S' in coastal_type:
                        xdist,ydist = zfun.ll2xy(X[SI], Y[SJ], x, y)
                        dist2land = dist2land + [np.sqrt(xdist**2 + ydist**2)]
                    if 'W' in coastal_type:
                        xdist,ydist = zfun.ll2xy(X[WI], Y[WJ], x, y)
                        dist2land = dist2land + [np.sqrt(xdist**2 + ydist**2)]
                    if 'E' in coastal_type:
                        xdist,ydist = zfun.ll2xy(X[EI], Y[EJ], x, y)
                        dist2land = dist2land + [np.sqrt(xdist**2 + ydist**2)]
                    # minimum distance
                    min_dist2land = np.min(dist2land)
                    # index of minimum distance
                    ind_dist2land = dist2land.index(min_dist2land)
                    # corresponding river direction
                    source_side = coastal_type[ind_dist2land]

                # get river direction values
                xoff = 0
                yoff = 0
                # offsets follow same convention as in pgrid code
                if source_side == 'N':
                        idir = 1
                        isign = -1
                        uv = 'v'
                elif source_side == 'S':
                        idir = 1
                        isign = 1
                        uv = 'v'
                        yoff = -1
                elif source_side == 'W':
                        idir = 0
                        isign = 1
                        uv = 'u'
                        xoff = -1
                elif source_side == 'E':
                        idir = 0
                        isign = -1
                        uv = 'u'

                # record distance
                distance_list.append(distance)
                # record cell indices
                ii_list.append(I_ind + xoff)
                jj_list.append(J_ind + yoff)
                # record cell flow directions
                idir_list.append(idir)
                isign_list.append(isign)
                uv_list.append(uv)

    return ii_list, jj_list, distance_list, idir_list, isign_list, uv_list

def get_nearest_coastal_cell_wwtp(sname,x,y,X,Y,mask_rho):
    """
    Algorithm that calculates the nearest coastal grid cell (row/col)
    to a point source,
    given its lat/lon coordinates, the lat/lon coordinates of the grid,
    and the mask of the rho-cells (to identify water vs land)

    Inputs:
        sname:      source name
        x:          WWTP lon coordinates
        y:          WWTP lat coordinates
        X:          domain lon coordinates
        Y:          domain lat coordinates
        mask_rho:   used to check if grid cell is water/land

    Outputs:
        row:        row index of nearest coastal cell    
        col:        col index of nearest coastal cell
        ringnum:    ring that contains nearest coastal cell  
    """

    if in_domain(x, y, X, Y):

        # initialize a boolean to track whether a grid cell has been selected for point source
        ps_located = False

        # figure out in which grid cell the wwtp is located
        ix = zfun.find_nearest_ind(X,x)
        iy = zfun.find_nearest_ind(Y,y)
        wwtpII = np.array(ix)
        wwtpJJ = np.array(iy)

        # search for best grid cell to place wwtp-----------------------
    
        # First, check if the wwtp is located in the ocean
        if mask_rho[wwtpII,wwtpJJ]: # mask of 1 means water
            # point source located at same location as wwtp
            ringnum = 0
            psII = wwtpII
            psJJ = wwtpJJ
            ps_located = True

        # If on land, search in square rings around the source for the nearest coastal cell
        elif not mask_rho[wwtpII,wwtpJJ]:
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
                    ii_list, jj_list, distance_list = get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho)
                    # add the cell info to arrays if in water
                    II_list = II_list + ii_list
                    JJ_list = JJ_list + jj_list
                    dist_list = dist_list + distance_list

                # record all water cells that are in bottom row of ring
                for iter in range(-ringnum,ringnum+1):
                    I_ind = wwtpII + iter
                    J_ind = wwtpJJ - ringnum
                    # get cell info if it is in water
                    ii_list, jj_list, distance_list = get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho)
                    # add the cell info to arrays if in water
                    II_list = II_list + ii_list
                    JJ_list = JJ_list + jj_list
                    dist_list = dist_list + distance_list

                # record all water cells that are in left column of ring (exclude diagonals)
                for iter in range(-(ringnum-1),ringnum):
                    J_ind = wwtpJJ + iter
                    I_ind = wwtpII - ringnum
                    # get cell info if it is in water
                    ii_list, jj_list, distance_list = get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho)
                    # add the cell info to arrays if in water
                    II_list = II_list + ii_list
                    JJ_list = JJ_list + jj_list
                    dist_list = dist_list + distance_list

                # record all water cells that are in right column of ring (exclude diagonals)
                for iter in range(-(ringnum-1),ringnum):
                    J_ind = wwtpJJ + iter
                    I_ind = wwtpII + ringnum
                    # get cell info if it is in water
                    ii_list, jj_list, distance_list = get_cell_info_wwtp(I_ind,J_ind,X,Y,x,y,mask_rho)
                    # add the cell info to arrays if in water
                    II_list = II_list + ii_list
                    JJ_list = JJ_list + jj_list
                    dist_list = dist_list + distance_list

                # calculate minimum distance coastal cell (if there were any coastal cells)
                if len(dist_list) > 0:
                    # minimum distance
                    min_dist_new = np.min(dist_list)
                    # index of minimum distance
                    min_dist_ind = dist_list.index(min_dist_new)
                    # corresponding i and j indices
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

                # if all cells in ring are on land, then go into next ring
                else:
                    # iterate into next level ring
                    ringnum += 1

        row = psJJ
        col = psII
        return [row,col,ringnum]
    
    else:
        print(' >> excluding ' + sname)
        row = np.nan
        col = np.nan
        ringnum = np.nan
        return [row,col,ringnum]
        
def get_nearest_coastal_cell_riv(sname,x,y,X,Y,mask_rho):
    """
    Algorithm that calculates the nearest coastal grid cell (row/col)
    to a rivermouth,
    given its lat/lon coordinates, the lat/lon coordinates of the grid,
    and the mask of the rho-cells (to identify water vs land)
    
    Inputs:
        sname:      source name
        x:          river mouth lon coordinates
        y:          river mouth lat coordinates
        X:          domain lon coordinates
        Y:          domain lat coordinates
        mask_rho:   used to check if grid cell is water/land

    Outputs:
        row:        row index of nearest coastal cell    
        col:        col index of nearest coastal cell
        idir:       0 = E-W river, 1 = N-S river
        isign:      1 = flowing N/E, -1 = flowing S,W
        uv:         'u' = E-W river, 'v' = N-S river
        ringnum:    ring that contains nearest coastal cell  
    """

    if in_domain(x, y, X, Y):

        # initialize a boolean to track whether a grid cell has been selected for point source
        rm_located = False

        # figure out in which grid cell the source is located
        ix = zfun.find_nearest_ind(X,x)
        iy = zfun.find_nearest_ind(Y,y)
        rivII = np.array(ix)
        rivJJ = np.array(iy)
    
        # search for best grid cell to place river mouth -----------------------
    
        # First, check if the river is already located in a coastal cell
        if mask_rho[rivII,rivJJ]: # mask of 1 means water
            ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(rivII,rivJJ,X,Y,x,y,mask_rho)
            # if list is not empty, then the cell is coastal
            if len(ii_list) > 0:
                # save the coastal information
                ringnum = 0
                rmII = ii_list[0]
                rmJJ = jj_list[0]
                idir = idir_list[0]
                isign = isign_list[0]
                uv = uv_list[0]
                rm_located = True

        # If original cell wasn't coastal,
        # Search in square rings around the source for the nearest coastal cell
        # initialize counters and arrays
        if not rm_located:
            ringnum = 1
            II_list = []
            JJ_list = []
            dist_list = []
            IDIR_list = []
            ISIGN_list = []
            UV_list = []
            min_dist = 1e9 # arbitratily large number so we can update with min value later

        while not rm_located:

            # record all water cells that are in top row of ring
            for iter in range(-ringnum,ringnum+1):
                I_ind = rivII + iter
                J_ind = rivJJ + ringnum
                # get cell info if it is in water
                ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho)
                # add the cell info to arrays if in water
                II_list = II_list + ii_list
                JJ_list = JJ_list + jj_list
                dist_list = dist_list + distance_list
                IDIR_list = IDIR_list + idir_list
                ISIGN_list = ISIGN_list + isign_list
                UV_list = UV_list + uv_list

            # record all water cells that are in bottom row of ring
            for iter in range(-ringnum,ringnum+1):
                I_ind = rivII + iter
                J_ind = rivJJ - ringnum
                # get cell info if it is in water
                ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho)
                # add the cell info to arrays if in water
                II_list = II_list + ii_list
                JJ_list = JJ_list + jj_list
                dist_list = dist_list + distance_list
                IDIR_list = IDIR_list + idir_list
                ISIGN_list = ISIGN_list + isign_list
                UV_list = UV_list + uv_list

            # record all water cells that are in left column of ring (exclude diagonals)
            for iter in range(-(ringnum-1),ringnum):
                J_ind = rivJJ + iter
                I_ind = rivII - ringnum
                # get cell info if it is in water
                ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho)
                # add the cell info to arrays if in water
                II_list = II_list + ii_list
                JJ_list = JJ_list + jj_list
                dist_list = dist_list + distance_list
                IDIR_list = IDIR_list + idir_list
                ISIGN_list = ISIGN_list + isign_list
                UV_list = UV_list + uv_list

            # record all water cells that are in right column of ring (exclude diagonals)
            for iter in range(-(ringnum-1),ringnum):
                J_ind = rivJJ + iter
                I_ind = rivII + ringnum
                # get cell info if it is in water
                ii_list, jj_list, distance_list, idir_list, isign_list, uv_list = get_cell_info_riv(I_ind,J_ind,X,Y,x,y,mask_rho)
                # add the cell info to arrays if in water
                II_list = II_list + ii_list
                JJ_list = JJ_list + jj_list
                dist_list = dist_list + distance_list
                IDIR_list = IDIR_list + idir_list
                ISIGN_list = ISIGN_list + isign_list
                UV_list = UV_list + uv_list

            # calculate minimum distance coastal cell (if there were any cells in water)
            if len(dist_list) > 0:
                # minimum distance
                min_dist_new = np.min(dist_list)
                # index of minimum distance
                min_dist_ind = dist_list.index(min_dist_new)
                # corresponding i and j indices
                min_II = II_list[min_dist_ind]
                min_JJ = JJ_list[min_dist_ind]
                # corresponding river directions
                min_IDIR = IDIR_list[min_dist_ind]
                min_ISIGN = ISIGN_list[min_dist_ind]
                min_UV = UV_list[min_dist_ind]

                # check if new min_dist is smaller than previous min_dist
                if min_dist_new < min_dist:
                    # save new min dist and grid cell index
                    min_dist = min_dist_new
                    rmII = min_II
                    rmJJ = min_JJ
                    idir = min_IDIR
                    isign = min_ISIGN
                    uv = min_UV
                    # check next ring if there is a closer coastal grid cell
                    ringnum +=1

                # if new min_dist is larger than before, then we have found the nearest coastal grid cell
                else:
                    rm_located = True
                    # print('Originally on land')
                    # # print('New lat/lon:({},{})'.format(round(Y[psJJ],2),round(X[psII],2)))

            # if all cells in ring are on land, then go into next ring
            else:
                # iterate into next level ring
                ringnum += 1

        # save location of river mouth
        # print('New lat/lon:({},{})'.format(round(Y[rmJJ],2),round(X[rmII],2)))
        row = rmJJ
        col = rmII
        # print('Corresponding grid row/col:({},{})'.format(row,col))
        return [row,col,idir,isign,uv,ringnum]
    
    else:
        print(' >> excluding ' + sname)
        row = np.nan
        col = np.nan
        ringnum = np.nan
        idir = np.nan
        isign = np.nan
        uv = np.nan
        return [row,col,idir,isign,uv,ringnum]

def traps_placement(source_type):
    '''
    Function that looks at all of the rivers and marine point sources 
    in SSM, and identified where to place these tiny rivers and point sources (TRAPS)
    in the LiveOcean grid. 

    Input:
        source_type:    'riv' or 'wwtp'

    Output:
        wwtp_info.csv or triv_info.csv saved in LO_data/grids/[gridname]
        Figures with source locations saved in LO_data/grids/[gridname]
    '''

    if source_type == 'wwtp':
        output_fn = 'wwtp_info.csv'
        fig_fn = 'wwtp_loc.png'
        inflow_type = 'Point Source'
    elif source_type == 'riv':
        output_fn = 'triv_info.csv'
        fig_fn = 'triv_loc.png'
        inflow_type = 'River'

    # get the grid data
    grid_fn = Ldir['grid'] / 'grid.nc'
    ds = xr.open_dataset(grid_fn)
    z = -ds.h.values
    mask_rho = np.transpose(ds.mask_rho.values)
    lon = ds.lon_rho.values
    lat = ds.lat_rho.values
    X = lon[0,:] # grid cell X values
    Y = lat[:,0] # grid cell Y values

    # read overlapping rivers
    repeatrivs_fn = Ldir['data'] / 'traps' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)

    # read lat/lon coordinates for the traps
    trapsll_fn = Ldir['data'] / 'traps' / 'SSM_source_info.xlsx'
    latlon_all_df = pd.read_excel(trapsll_fn,usecols='D,E,F,G,N,O')

    # initialize dataframe to save results
    rowcol_df = pd.DataFrame()
    SSMrivll_df = pd.DataFrame()
    snames = []
    ringnums = []

    # look at either only rivers or point sources
    latlon_df = latlon_all_df.loc[latlon_all_df['Inflow_Typ'] == inflow_type]

    # loop through all of the traps
    for source in latlon_df.index:
        sname = latlon_df._get_value(source,'Name') # source name
        
        if inflow_type == 'River':

            # check if river already in LiveOcean
            SSM_repeats = repeatrivs_df['SSM_rname'] # get names of repeat rivers
            checkname = sname
            # remove 1 and 2 from river name
            if '- 1' in sname:
                checkname = sname.replace(' - 1','')
            elif '- 2' in sname:
                checkname = sname.replace(' - 2','')
            # don't add rivers that already exist
            if SSM_repeats.str.contains(checkname).any():
                continue 

            # add river to LiveOcean if not already present
            if '1' in sname: # some large river mouths have two lat/lon coords
                x1 = latlon_df._get_value(source,'Lon')
                y1 = latlon_df._get_value(source,'Lat')
                # search for corresponding source containing '2' in its name
                sname2 = sname.replace('1','2')
                x2 = latlon_df['Lon'].loc[latlon_df['Name'] == sname2].values[0]
                y2 = latlon_df['Lat'].loc[latlon_df['Name'] == sname2].values[0]
                # average lat/lon
                x = np.mean([x1, x2]) # lon = average of both lons
                y = np.mean([y1, y2]) # lat = average of both lats
                # remove the number from the river name
                sname = sname.replace(' - 1','')
                # get nearest coastal grid cell
                [row,col,idir,isign,uv,ringnum] = get_nearest_coastal_cell_riv(sname,x,y,X,Y,mask_rho)
                # save coordinates to rowcol_df
                rowcol_df.loc[sname, 'row_py'] = row 
                rowcol_df.loc[sname, 'col_py'] = col
                rowcol_df.loc[sname, 'idir']   = idir
                rowcol_df.loc[sname, 'isign']  = isign
                rowcol_df.loc[sname, 'uv']     = uv
                SSMrivll_df.loc[sname,'Lon'] = x
                SSMrivll_df.loc[sname,'Lat'] = y
                snames = snames + [sname]
                # save ringnum for tracking. ringnum = zero means the source didn't move
                ringnums = ringnums + [ringnum]
           
            elif '- 2' in sname:
                pass # already accounted in the above if case

            else:
                # get lat/lon
                x = latlon_df._get_value(source,'Lon')
                y = latlon_df._get_value(source,'Lat')
                # get nearest coastal gridcell
                [row,col,idir,isign,uv,ringnum] = get_nearest_coastal_cell_riv(sname,x,y,X,Y,mask_rho)
                # save coordinates to rowcol_df
                rowcol_df.loc[sname, 'row_py'] = row 
                rowcol_df.loc[sname, 'col_py'] = col
                rowcol_df.loc[sname, 'idir']   = idir
                rowcol_df.loc[sname, 'isign']  = isign
                rowcol_df.loc[sname, 'uv']     = uv
                SSMrivll_df.loc[sname,'Lon'] = x
                SSMrivll_df.loc[sname,'Lat'] = y
                snames = snames + [sname]
                # save ringnum for tracking. ringnum = zero means the source didn't move
                ringnums = ringnums + [ringnum]
        
        else: # wwtps are easier
            # get lat/lon
            x = latlon_df._get_value(source,'Lon')
            y = latlon_df._get_value(source,'Lat')
            # get nearest coastal gridcell
            [row,col,ringnum] = get_nearest_coastal_cell_wwtp(sname,x,y,X,Y,mask_rho)
            # save coordinates to rowcol_df
            rowcol_df.loc[sname, 'row_py'] = row 
            rowcol_df.loc[sname, 'col_py'] = col
            snames = snames + [sname]
            # save ringnum for tracking. ringnum = zero means the source didn't move
            ringnums = ringnums + [ringnum]

        # update name of index column
        rowcol_df.index.name = 'rname'

    # save the source info
    out_rfn = Ldir['grid'] / output_fn
    print('\nCreating ' + str(out_rfn))
    # drop any sources that were outside of the domain, and had values padded with nans
    rowcol_df = rowcol_df.dropna()
    rowcol_df.to_csv(out_rfn)

    plotting = True

    if plotting == True:
        # PLOTTING FOR TESTING ------------------------------------------------------------------------
        plon, plat = pfun.get_plon_plat(lon,lat)

        # make a version of z with nans where masked
        zm = z.copy()
        zm[np.transpose(mask_rho) == 0] = np.nan
        zm[np.transpose(mask_rho) != 0] = -1

        # bathymetry
        fig = plt.figure(figsize=(17,7))
        plt.tight_layout()
        ax0 = plt.subplot2grid((2,5),(0,0), rowspan=2, colspan=2) #fig.add_subplot(141)
        # pfun.add_coast(ax0,color='black')
        ax0.pcolormesh(plon, plat, zm, edgecolor='none', linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

        ax1 = plt.subplot2grid((2,5),(0,2),rowspan=2, colspan=1) #fig.add_subplot(246)
        ax1.pcolormesh(plon, plat, zm, edgecolor='none', linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

        ax2 = plt.subplot2grid((2,5),(0,3),rowspan=2, colspan=2) #fig.add_subplot(122)
        ax2.pcolormesh(plon, plat, zm, edgecolor='none', linewidth=0.5, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

        # All new sources
        x0_1 = -128.5
        x0_2 = -121.7
        y0_1 = 46.75
        y0_2 = 51.25
        # North Sound
        x1_1 = -123.1
        x1_2 = -122.2
        y1_1 = 47.5
        y1_2 = 48.8
        # South Sound
        x2_1 = -123.2
        x2_2 = -122.3
        y2_1 = 46.9
        y2_2 = y1_1

        ax0.set_xlim(x0_1,x0_2) # All new sources
        ax0.set_ylim(y0_1,y0_2) # All new sources
        ax1.set_xlim(x1_1,x1_2) # North Sound
        ax1.set_ylim(y1_1,y1_2) # North Sound
        ax2.set_xlim(x2_1,x2_2) # South Sound
        ax2.set_ylim(y2_1,y2_2) # South Sound

        # Add river track directions -------------------------------------------------------
        if inflow_type == 'River':
            rri_df = rowcol_df
            lon = ds.lon_rho.values
            lat = ds.lat_rho.values
            lon_u = ds.lon_u.values
            lat_u = ds.lat_u.values
            lon_v = ds.lon_v.values
            lat_v = ds.lat_v.values

            # label counters
            first_label = True
            first_label_LO = True

            for i,rn in enumerate(rri_df.index):
                # These are indices (python, zero-based) into either the
                # u or v grids.
                ii = int(rri_df.loc[rn,'col_py'])
                jj = int(rri_df.loc[rn,'row_py'])
                
                uv = rri_df.loc[rn,'uv']
                isign = rri_df.loc[rn,'isign']
                idir = rri_df.loc[rn,'idir']
                
                if uv == 'u' and isign == 1:
                    # River source on W side of rho cell
                    if first_label:
                        ax0.scatter(lon[jj,ii+1], lat[jj,ii+1], color='k', marker='o',s = 10, label='Tiny River')
                        first_label = False
                    else:
                        ax0.scatter(lon[jj,ii+1], lat[jj,ii+1], color='k', marker='o',s = 10)
                    for axis in [ax1,ax2]:
                        ax2.plot(lon_u[jj,ii], lat_u[jj,ii],marker='>',color='firebrick')
                if uv == 'u' and isign == -1:
                    # River source on E side of rho cell
                    ax0.scatter(lon[jj,ii], lat[jj,ii+1], color='k', marker='o',s = 10)
                    for axis in [ax1,ax2]:
                            axis.plot(lon_u[jj,ii], lat_u[jj,ii],marker='<',color='firebrick')
                if uv == 'v' and isign == 1:
                    # River source on S side of rho cell
                    ax0.scatter(lon[jj+1,ii], lat[jj+1,ii], color='k', marker='o',s = 10)
                    for axis in [ax1,ax2]:
                          axis.plot(lon_v[jj,ii], lat_v[jj,ii],marker='^',color='royalblue')  
                if uv == 'v' and isign == -1:
                    # River source on N side of rho cell
                    ax0.scatter(lon[jj,ii], lat[jj,ii], color='k', marker='o',s = 10)
                    for axis in [ax1,ax2]:
                            axis.plot(lon_v[jj,ii], lat_v[jj,ii],marker='v',color='royalblue')

        # -----------------------------------------------------

        # plot point sources
        ps_lon = []
        ps_lat = []
        for i,ind in enumerate(rowcol_df['col_py']):
            if math.isnan(ind):
                ps_lon = ps_lon + [np.nan]
                ps_lat = ps_lat + [np.nan]
            else:
                ps_lon = ps_lon + [X[int(ind)]]
                indy = rowcol_df['row_py'][i]
                ps_lat = ps_lat + [Y[int(indy)]]
        
        if inflow_type == 'Point Source':
            for axis in [ax0,ax1,ax2]:
                axis.scatter(ps_lon,ps_lat, color='k', marker='o', s=20, label='Point Source')

        if inflow_type == 'River':
            # plot pre-existing river locations
            LOrivs_fn = Ldir['grid'] / 'river_info.csv'
            LOrivs_df = pd.read_csv(LOrivs_fn)
            for rn in LOrivs_df.index:
                # These are indices (python, zero-based) into either the
                # u or v grids.
                ii = int(LOrivs_df.loc[rn,'col_py'])
                jj = int(LOrivs_df.loc[rn,'row_py'])
                
                uv = LOrivs_df.loc[rn,'uv']
                isign = LOrivs_df.loc[rn,'isign']
                idir = LOrivs_df.loc[rn,'idir']

                if uv == 'u' and isign == 1:
                    # River source on W side of rho cell
                    if first_label_LO:
                       ax0.scatter(lon[jj,ii+1], lat[jj,ii+1], color='darkorange', marker='*',
                                   edgecolor = 'darkorange', s = 20, label='Pre-existing LO River')
                       first_label_LO = False
                    else:
                        ax0.scatter(lon[jj,ii+1], lat[jj,ii+1], color='darkorange', marker='*',
                                   edgecolor = 'darkorange', s = 20)
                    for axis in [ax1,ax2]:
                        axis.plot(lon_u[jj,ii], lat_u[jj,ii],marker='>',color='firebrick')

                if uv == 'u' and isign == -1:
                    # River source on E side of rho cell
                    ax0.scatter(lon[jj,ii], lat[jj,ii], color='darkorange',
                    marker='*', edgecolor = 'darkorange', s = 20)
                    for axis in [ax1,ax2]:
                        axis.plot(lon_u[jj,ii], lat_u[jj,ii],marker='<',color='firebrick')

                if uv == 'v' and isign == 1:
                    # River source on S side of rho cell
                    ax0.scatter(lon[jj+1,ii], lat[jj+1,ii], color='darkorange',
                    marker='*', edgecolor = 'darkorange', s = 20)
                    for axis in [ax1,ax2]:
                        axis.plot(lon_v[jj,ii], lat_v[jj,ii],marker='^',color='royalblue')
                            
                if uv == 'v' and isign == -1:
                    # River source on N side of rho cell
                    ax0.scatter(lon[jj,ii], lat[jj,ii], color='darkorange',
                    marker='*', edgecolor = 'darkorange', s = 20)
                    for axis in [ax1,ax2]:
                        axis.plot(lon_v[jj,ii], lat_v[jj,ii],marker='v',color='royalblue')
                            

        ax0.add_patch(Rectangle((x1_1,y1_1),x1_2-x1_1,y1_2-y1_1,edgecolor='deeppink',facecolor='none',linewidth=2))
        for spine in ax1.spines.values():
            spine.set_edgecolor('deeppink')
            spine.set_linewidth(2)

        ax0.add_patch(Rectangle((x2_1,y2_1),x2_2-x2_1,y2_2-y2_1,edgecolor='lightslategrey',facecolor='none', linewidth=2))
        for spine in ax2.spines.values():
            spine.set_edgecolor('lightslategrey')
            spine.set_linewidth(2)

        # finalize plot
        ax0.set_ylabel('Lat', fontsize = 12)
        ax0.set_xlabel('Lon', fontsize = 12)
        ax0.legend(loc='lower left', fontsize = 12)
        ax0.tick_params(axis='both', which='major', labelsize=10)
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax2.set_xticks([])
        ax2.set_yticks([])
        plt.suptitle('Algorithm Placement of {}s'.format(inflow_type), fontsize = 16)
        plt.locator_params(nbins=4)
        pfun.dar(ax0) 
        pfun.dar(ax1)
        pfun.dar(ax2)
        # plt.show()
        # save figure
        out_fig = Ldir['grid'] / fig_fn
        plt.savefig(out_fig)

    return

def get_qtbio(gri_df, dt_ind, yd_ind, Ldir, traps_type):

    # Only add biology to pre-existing LO river if Ecology has data
    if traps_type == 'LOriv':
        # get names of duplicate rivers
        repeatrivs_fn = Ldir['data'] / 'traps' / 'LiveOcean_SSM_rivers.xlsx'
        repeatrivs_df = pd.read_excel(repeatrivs_fn)
        LObio_names_all = list(repeatrivs_df.loc[repeatrivs_df['in_both'] == 1, 'LO_rname'])
        # remove the weird rivers
        weird_duplicate_rivers = ['Alberni Inlet', 'Chehalis R', 'Gold River', 'Willapa R', 'Columbia R', 'Comox']
        # Note that these are the names that LO calls the rivers
        LObio_names = [rname for rname in LObio_names_all if LO2SSM_name(rname) not in weird_duplicate_rivers]

    # load climatological data
    if traps_type != 'LOriv':
        # don't need flow and temp for pre-existing LO rivers
        Cflow_df = pd.read_pickle(Ldir['Cflow_'+traps_type+'_fn'])
        Ctemp_df = pd.read_pickle(Ldir['Ctemp_'+traps_type+'_fn'])
    CDO_df   = pd.read_pickle(Ldir['CDO_'+traps_type+'_fn'])
    CNH4_df  = pd.read_pickle(Ldir['CNH4_'+traps_type+'_fn'])
    CNO3_df  = pd.read_pickle(Ldir['CNO3_'+traps_type+'_fn'])
    CTalk_df = pd.read_pickle(Ldir['CTalk_'+traps_type+'_fn'])
    CTIC_df  = pd.read_pickle(Ldir['CTIC_'+traps_type+'_fn'])

    # year day index starts from 1. Convert to start from 0 to work with python
    yd_ind = yd_ind - 1

    # initialize output dict
    qtbio_df_dict = dict()
    for rn in gri_df.index:

        # convert LO river name to SSM river name
        if traps_type == 'LOriv':
            if rn in LObio_names:
                rn = LO2SSM_name(rn)
            else:
                # skips rivers for which Ecology does not have data
                continue    
        
        # initialize screen output
        sout = '-- ' + str(rn) + ': '
        
        # initialize a qtbio (flow and temperature vs. time) DataFrame for this river
        qtbio_df = pd.DataFrame(index=dt_ind, columns=['flow','temp','Oxyg','NH4','NO3','TAlk','TIC'])
        # fill with climatological fields
        if traps_type != 'LOriv':
            qtbio_df.loc[:, 'flow'] = Cflow_df.loc[yd_ind,rn].values
            qtbio_df.loc[:, 'temp'] = Ctemp_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'Oxyg']   = CDO_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'NH4']  = CNH4_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'NO3']  = CNO3_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'TAlk'] = CTalk_df.loc[yd_ind,rn].values
        qtbio_df.loc[:, 'TIC']  = CTIC_df.loc[yd_ind,rn].values
        sout += 'filled from climatology (Ecology dataset)'
          
        # save in the dict
        qtbio_df_dict[rn] = qtbio_df
        
        # screen output
        print(sout)
        sys.stdout.flush()
        
    return qtbio_df_dict

def combine_adjacent(lst):
    """
    Given a list, e.g. ['a','b','c','d']
    returns: ['a+b', 'c+d']
    """
    combined = [x + '+' + y for x, y in zip(lst[::2],lst[1::2])]
    return combined

def weighted_average(vn,qtbio_df_1, qtbio_df_2):
    '''
    Calculate the weighted average properties based on flowrate of two overlapping sources
    '''
    # get flowrates
    flow1 = qtbio_df_1['flow'].values
    flow2 = qtbio_df_2['flow'].values
    # get variable
    var1 = qtbio_df_1[vn].values
    var2 = qtbio_df_2[vn].values
    # calculate weighted average based on flowrate
    waverage = [np.average([var1[i], var2[i]], weights = [flow1[i], flow2[i]]) for i in range(len(flow1))]
    return waverage

def LO2SSM_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = Ldir['data'] / 'traps' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_SSM = repeatrivs_df.loc[repeatrivs_df['LO_rname'] == rname, 'SSM_rname'].values[0]

    return rname_SSM