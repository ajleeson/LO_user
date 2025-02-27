"""
Plots a comparison of bottom DO for two identical grids with different oxygen value.

This is a custom function for a particular experiment, but code can be adapted for other use cases in the future.

From the terminal: python bottom_DO_nutrient_comparison.py

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

# ----------------------------------------------------------------------

# Provide information about models to compare
# basecase (has nutrients)
bc_gtagex = 'cas6_traps2_x2b'
# test condition (no nutrients)
c1_gtagex = 'cas6_traps3_x2b'
gtagexes = [bc_gtagex, c1_gtagex]

hr = '0025'
date = '2017.03.07'

#'Puget Sound','Whidbey Basin','North King County','Lynch Cove'
region = 'Whidbey Basin'

# Variables to compare
vn_list = ['oxygen']#['NO3','phytoplankton','oxygen']

# Show WWTP locations?
WWTP_loc = True

##########################################################
# helper funcions

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = '../../LO_data/traps/LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_LO = repeatrivs_df.loc[repeatrivs_df['SSM_rname'] == rname, 'LO_rname'].values[0]

    return rname_LO

def LO2SSM_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = Ldir['data'] / 'traps' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_SSM = repeatrivs_df.loc[repeatrivs_df['LO_rname'] == rname, 'SSM_rname'].values[0]

    return rname_SSM

# riv fun function to get biology
def get_bio_vec(vn, rn, yd_ind):
    ndt = len(yd_ind)
    yd = yd_ind.values
    ovec = np.ones(ndt)
    if vn == 'NO3':
        if rn == 'fraser':
            vv = 2 + (13/2) + (13/2)*np.cos(2*np.pi*((yd-30)/366))
        elif rn == 'columbia':
            vv = 5 + (35/2) + (35/2)*np.cos(2*np.pi*((yd)/366))
        else:
            vv = 5 * ovec
    elif vn == 'Oxyg':
        vv = 350 * ovec
    elif vn in ['TAlk', 'TIC']:
        if rn in ['columbia', 'deschutes', 'duwamish']:
            vv = 1000 * ovec
        else:
            vv = 300 * ovec
    else:
        vv = 0 * ovec # all others filled with zeros
    return vv

##########################################################
if WWTP_loc == True:
    # set up the time index for the record
    Ldir = Lfun.Lstart()
    dsf = Ldir['ds_fmt']
    dt0 = datetime.strptime('2020.01.01',dsf)
    dt1 = datetime.strptime('2020.12.31',dsf)
    days = (dt0, dt1)
        
    # pandas Index objects
    dt_ind = pd.date_range(start=dt0, end=dt1)
    yd_ind = pd.Index(dt_ind.dayofyear)

    # Get LiveOcean grid info --------------------------------------------------

    # get the grid data
    ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc')
    z = -ds.h.values
    mask_rho = np.transpose(ds.mask_rho.values)
    lon = ds.lon_rho.values
    lat = ds.lat_rho.values
    X = lon[0,:] # grid cell X values
    Y = lat[:,0] # grid cell Y values
    plon, plat = pfun.get_plon_plat(lon,lat)
    # make a version of z with nans where masked
    zm = z.copy()
    zm[np.transpose(mask_rho) == 0] = np.nan
    zm[np.transpose(mask_rho) != 0] = -1

    # Point Sources -------------------------------------------------------------
    # Prepare data for spatial summary plots

    # get flow, nitrate, and ammonium values
    fp_wwtps = '../../LO_output/pre/traps/point_sources/Data_historical/'
    flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow_1999_2017.p')    # m3/s
    no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3_1999_2017.p')      # mmol/m3
    nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4_1999_2017.p')      # mmol/m3

    # calculate total DIN concentration in mg/L
    dindf_wwtps = (no3df_wwtps + nh4df_wwtps)/71.4    # mg/L

    # calculate daily loading timeseries in kg/d
    dailyloaddf_wwtps = 86.4*dindf_wwtps*flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s

    # calculate average daily load over the year (kg/d)
    avgload_wwtps = dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

    # add row and col index for plotting on LiveOcean grid
    griddf0_wwtps = pd.read_csv('../../LO_data/grids/cas6/wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    avgload_wwtps = avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    avgload_wwtps = avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    # get point source lat and lon
    lon_wwtps = [X[int(col)] for col in avgload_wwtps['col_py']]
    lat_wwtps = [Y[int(row)] for row in avgload_wwtps['row_py']]

    # Rivers -------------------------------------------------------------
    # Prepare data for spatial summary plots

    # get flow, nitrate, and ammonium values
    # fp_trivs = '../../LO_output/pre/traps/all_rivers/Data_historical/'
    fp_trivs = '../../LO_output/pre/traps/tiny_rivers/Data_historical/'
    fp_LOrivs = '../../LO_output/pre/river/cas6/Data_historical/'
    fp_LObio = '../../LO_output/pre/traps/LO_rivbio/Data_historical/'
    flowdf_rivs = pd.read_pickle(fp_trivs+'CLIM_flow_1999_2017.p')    # m3/s
    no3df_rivs = pd.read_pickle(fp_trivs+'CLIM_NO3_1999_2017.p')      # mmol/m3
    nh4df_rivs = pd.read_pickle(fp_trivs+'CLIM_NH4_1999_2017.p')      # mmol/m3
    # pre-existing LO river flowrates
    flowdf_LOrivs = pd.read_pickle(fp_LOrivs+'CLIM_flow_1980_2020.p')    # m3/s
    flowdf_LOrivs = flowdf_LOrivs.reset_index(drop=True)
    # Ecology data for pre-existing LO rivers
    no3df_LOrivs_ecol = pd.read_pickle(fp_LObio+'CLIM_NO3_1999_2017.p')      # mmol/m3
    nh4df_LOrivs_ecol = pd.read_pickle(fp_LObio+'CLIM_NH4_1999_2017.p')      # mmol/m3
    # get names of all pre-existing rivers for which Ecology has data (use LO naming convention)
    LObio_names = [SSM2LO_name(riv) for riv in no3df_LOrivs_ecol.columns]

    # get biology data for pre-existing rivers
    no3df_LOrivs = pd.DataFrame()
    nh4df_LOrivs = pd.DataFrame()
    for rn in flowdf_LOrivs:
        # Use ecology data if there exists any
        if rn in LObio_names:
            no3df_LOrivs[rn] = no3df_LOrivs_ecol[LO2SSM_name(rn)]
            nh4df_LOrivs[rn] = nh4df_LOrivs_ecol[LO2SSM_name(rn)]
        # Otherwise use the rivfun method to guess
        else:
            no3df_LOrivs[rn] = get_bio_vec('NO3', rn, yd_ind)
            nh4df_LOrivs[rn] = get_bio_vec('NH4', rn, yd_ind)

    # calculate total DIN concentration in mg/L
    dindf_rivs = (no3df_rivs + nh4df_rivs)/71.4          # mg/L
    dindf_LOrivs = (no3df_LOrivs + nh4df_LOrivs)/71.4    # mg/L

    # calculate daily loading timeseries in kg/d
    dailyloaddf_rivs = 86.4*dindf_rivs*flowdf_rivs       # kg/d = 86.4 * mg/L * m3/s
    dailyloaddf_LOrivs = 86.4*dindf_LOrivs*flowdf_LOrivs # kg/d = 86.4 * mg/L * m3/s

    # calculate average daily load over the year (kg/d)
    avgload_trivs = dailyloaddf_rivs.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')
    avgload_LOrivs = dailyloaddf_LOrivs.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

    # add row and col index for plotting on LiveOcean grid (tiny rivers)
    griddf0_trivs = pd.read_csv('../../LO_data/grids/cas6/triv_info.csv')
    griddf_trivs = griddf0_trivs.set_index('rname') # use river name as index
    avgload_trivs = avgload_trivs.join(griddf_trivs['row_py']) # add row to avg load df (uses rname to index)
    avgload_trivs = avgload_trivs.join(griddf_trivs['col_py']) # do the same for cols

    # add row and col index for plotting on LiveOcean grid (pre-existing rivers)
    griddf0_LOrivs = pd.read_csv('../../LO_data/grids/cas6/river_info.csv')
    griddf_LOrivs = griddf0_LOrivs.set_index('rname') # use river name as index
    avgload_LOrivs = avgload_LOrivs.join(griddf_LOrivs['row_py']) # add row to avg load df (uses rname to index)
    avgload_LOrivs = avgload_LOrivs.join(griddf_LOrivs['col_py']) # do the same for cols

    # get average load of all rivers
    avgload_allrivs = pd.concat([avgload_trivs, avgload_LOrivs])
    # drop nans
    avgload_allrivs = avgload_allrivs.dropna()

    # get trivs lat and lon
    lon_riv = [X[int(col)] for col in avgload_allrivs['col_py']]
    lat_riv = [Y[int(row)] for row in avgload_allrivs['row_py']]

    # PLOTTING ----------------------------------------------------

    # define marker sizes (minimum size is 10 so dots don't get too small)
    sizes_wwtps = [max(0.3*load,30) for load in avgload_wwtps['avg-daily-load(kg/d)']]
    sizes_rivs = [max(0.3*load,30) for load in avgload_allrivs['avg-daily-load(kg/d)']]

##########################################################

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# get plotting limits based on region
if region == 'North King County':
    xmin = -122.6
    xmax = -122.3
    ymin = 47.55
    ymax = 47.75
elif region == 'Lynch Cove':
    xmin = -123.2
    xmax = -122.8
    ymin = 47.3
    ymax = 47.5
elif region == 'Whidbey Basin':
    xmin = -122.8
    xmax = -122.1
    ymin = 48
    ymax = 48.4
elif region == 'Puget Sound':
    xmin = -123.5
    xmax = -122
    ymin = 47
    ymax = 48.5

# Get grid data
G = zrfun.get_basic_info('/home/aleeson/LO_data/grids/cas6/grid.nc', only_G=True)
grid_ds = xr.open_dataset('/home/aleeson/LO_data/grids/cas6/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values
px, py = pfun.get_plon_plat(G['lon_rho'],G['lat_rho'])
# get indices of min and max
imin_x = find_nearest(px[0,:],xmin)
imax_x = find_nearest(px[0,:],xmax)
imin_y = find_nearest(py[:,0],ymin)
imax_y = find_nearest(py[:,0],ymax)

# Loop through variable to compare
for vn in vn_list:
    if vn == 'NO3':
        slev = -1
        stext = 'Surface'
        vmin = 0
        vmax = 44
        cmap = cmocean.cm.matter
    elif vn == 'NH4':
        slev = -1
        stext = 'Surface'
        vmin = 0
        vmax = 5
        cmap = cmocean.cm.matter
    elif vn == 'phytoplankton':
        slev = -1
        stext = 'Surface'
        vmin = 0
        vmax = 10
        cmap = cmocean.cm.algae
    elif vn == 'oxygen':
        slev = 0
        stext = 'Bottom'
        vmin = 0
        vmax = 10
        cmap = cmocean.cm.oxy

    # scale variable
    scale =  pinfo.fac_dict[vn]

    # Initialize figure
    fs = 10
    pfun.start_plot(fs=fs, figsize=(8,9))
    fig = plt.figure()

    # loop through and plot both conditions
    for i,gtagex in enumerate(gtagexes):

        # get data
        fp = Ldir['roms_out'] / gtagex / ('f' + date) / ('ocean_his_'+ hr +'.nc')
        ds = xr.open_dataset(fp)

        # Plot individual map field
        ax = fig.add_subplot(3,2, i+1)
        # cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
        #             cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
        #             vlims_fac=pinfo.range_dict[vn], do_mask_edges=True)
        v = ds[vn][0, slev,:,:].values * scale
        cs = ax.pcolormesh(px,py,v, vmin=vmin, vmax=vmax, cmap=cmap)
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        pfun.add_coast(ax)
        # ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        plt.suptitle('%s \n %s %s %s' % (region,stext, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.5*fs)
        ax.set_xlabel('Longitude')
        if i == 0:
            ax.set_title('With WWTP Nutrients')
            ax.set_ylabel('Latitude')
            # save dataset for later use
            bc_ds = ds
        elif i ==1:
            ax.set_title('No WWTP Nutrients')
            # fig.colorbar(cs)
            # add colorbar
            cb_ax = fig.add_axes([.91,.67,.03,.23])
            fig.colorbar(cs,orientation='vertical',cax=cb_ax)
            ax.set_yticklabels([])
            # save dataset for later use
            c1_ds = ds

    # plot the difference
    plt.subplots_adjust(hspace=0.5)
    ax = fig.add_subplot(3,2,(3,6))
    diff = (bc_ds[vn] - c1_ds[vn]) * scale
    vmin = np.min(diff[0,slev,imin_y:imax_y,imin_x:imax_x])
    vmax = np.max(diff[0,slev,imin_y:imax_y,imin_x:imax_x])
    # make sure the colorbar is always centered about zero
    cmap = cmocean.tools.crop(cmocean.cm.balance, vmin, vmax, 0)
    cs = ax.pcolormesh(px,py,diff[0,slev,:,:], vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.set_title('(With Nutrients Condition) - (No Nutrients Condition)')
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')
    cb_ax = fig.add_axes([.91,.1,.03,.5])
    fig.colorbar(cs,orientation='vertical',cax=cb_ax)

    # add wwtp locations
    if WWTP_loc == True:
        ax.scatter(lon_wwtps,lat_wwtps,color='lawngreen', edgecolors='green', alpha=0.7, s=sizes_wwtps)
        # loc_df = pd.read_csv('../../LO_data/grids/cas6/wwtp_as_triv_info.csv')
        # for rn in loc_df.index:
        # # These are indices (python, zero-based) into either the
        # # u or v grids.
        #     ii = int(loc_df.loc[rn,'col_py'])
        #     jj = int(loc_df.loc[rn,'row_py'])
            
        #     uv = loc_df.loc[rn,'uv']
        #     isign = loc_df.loc[rn,'isign']
        #     idir = loc_df.loc[rn,'idir']

        #     if uv == 'u' and isign == 1:
        #         ax.plot(lon_u[jj,ii], lat_u[jj,ii],marker='>',color='darkorange')

        #     if uv == 'u' and isign == -1:
        #         ax.plot(lon_u[jj,ii], lat_u[jj,ii],marker='<',color='darkorange')

        #     if uv == 'v' and isign == 1:
        #         ax.plot(lon_v[jj,ii], lat_v[jj,ii],marker='^',color='darkorange')
                        
        #     if uv == 'v' and isign == -1:
        #         ax.plot(lon_v[jj,ii], lat_v[jj,ii],marker='v',color='darkorange')
                           

    # Generate plot
    plt.show()

# Calculate total amount of nutrients
S = {'N':30, 'Cs_w':bc_ds['Cs_w'].values, 'hc':bc_ds['hc'].values, 's_w':bc_ds['s_w'].values,
     'Vtransform':bc_ds['Vtransform'].values}
i = 20
j = 20
z_w = zrfun.get_z(bc_ds['h'][i,j], bc_ds['zeta'][:,i,j], S, only_w=True)
dz = np.diff(z_w)