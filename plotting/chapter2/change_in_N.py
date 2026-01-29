"""
Compare average bottom DO between multiple years
(Set up to run for 6 years)

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator
import matplotlib.patches as patches
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import pandas as pd
import cmocean
import matplotlib.pylab as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

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

remove_straits = True

WWTP_loc = True

# Hanning window length
nwin = 20

# years =  ['2014','2015']
years =  ['2014','2015','2016','2017','2018','2019','2020']

# which  model run to look at?
gtagexes = ['cas7_t1_x11ab','cas7_t1noDIN_x11ab'] 

# where to put output figures
out_dir = Ldir['LOo'] / 'chapter_2' / 'figures'
Lfun.make_dir(out_dir)

regions = ['All Puget Sound']
# regions = ['Hood Canal', 'South Sound', 'Whidbey Basin', 'Main Basin', 'All Puget Sound']
# colors = ['hotpink','mediumpurple','dodgerblue','yellowgreen','black']

plt.close('all')

##############################################################
##                    HELPER FUNCTIONS                      ##
##############################################################

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in SSM, find corresponding river name in LO
    """
    repeatrivs_fn = '../../../LO_data/trapsD01/LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_LO = repeatrivs_df.loc[repeatrivs_df['SSM_rname'] == rname, 'LO_rname'].values[0]
    return rname_LO

def LO2SSM_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = Ldir['data'] / 'trapsD01' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_SSM = repeatrivs_df.loc[repeatrivs_df['LO_rname'] == rname, 'SSM_rname'].values[0]
    return rname_SSM


if WWTP_loc == True:

    # get the grid data
    ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
    z = -ds.h.values
    mask_rho = np.transpose(ds.mask_rho.values)
    lon = ds.lon_rho.values
    lat = ds.lat_rho.values
    X = lon[0,:] # grid cell X values
    Y = lat[:,0] # grid cell Y values

    # get flow, nitrate, and ammonium values
    fp_wwtps = '../../../LO_output/pre/trapsP01/moh20_wwtps/lo_base/Data_historical/'
    moh20_flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow.p')    # m3/s
    moh20_no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3.p')      # mmol/m3
    moh20_nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4.p')      # mmol/m3

    fp_wwtps = '../../../LO_output/pre/trapsP01/was24_wwtps/lo_base/Data_historical/'
    was24_flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow.p')    # m3/s
    was24_no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3.p')      # mmol/m3
    was24_nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4.p')      # mmol/m3

    fp_trivs = '../../../LO_output/pre/trapsP01/moh20_tinyrivers/lo_base/Data_historical/'
    flowdf_trivs = pd.read_pickle(fp_trivs+'CLIM_flow.p')    # m3/s
    no3df_trivs = pd.read_pickle(fp_trivs+'CLIM_NO3.p')      # mmol/m3
    nh4df_trivs = pd.read_pickle(fp_trivs+'CLIM_NH4.p')      # mmol/m3

    fp_LOrivbio = '../../../LO_output/pre/trapsP01/moh20_LOrivbio/lo_base/Data_historical/'
    fp_LOrivflo = '../../../LO_output/pre/river1/lo_base/Data_historical/'
    flowdf_LOrivflo = pd.read_pickle(fp_LOrivflo+'CLIM_flow.p')    # m3/s
    no3df_LOrivbio = pd.read_pickle(fp_LOrivbio+'CLIM_NO3.p')      # mmol/m3
    nh4df_LOrivbio = pd.read_pickle(fp_LOrivbio+'CLIM_NH4.p')      # mmol/m3
    # read overlapping rivers
    repeatrivs_fn = Ldir['data'] / 'trapsD01' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    # convert LO name to SSM name
    SSM_repeats = repeatrivs_df['SSM_rname'].values
    # remove nans
    SSM_repeats = [x for x in SSM_repeats if str(x) != 'nan']
    # remove weird rivers
    weird_rivers = ['Alberni Inlet', 'Chehalis R', 'Gold River','Willapa R', 'Columbia R', 'Comox']
    SSM_repeats = [rname for rname in SSM_repeats if rname not in weird_rivers]
    LO_repeats = [SSM2LO_name(rname) for rname in SSM_repeats]
    # Keep only SSM repeat river names from list of river names
    flowdf_LOrivflo.columns = [LO2SSM_name(rname) if rname in LO_repeats else rname for rname in flowdf_LOrivflo.columns]
    # sort the nutrient dfs
    no3df_LOrivbio = no3df_LOrivbio.reindex(columns=flowdf_LOrivflo.columns)
    nh4df_LOrivbio = nh4df_LOrivbio.reindex(columns=flowdf_LOrivflo.columns)
    # add nutrient data
    yd = np.linspace(0,365,366)
    no3df_LOrivbio['columbia'] = 5 + (35/2) + (35/2)*np.cos(2*np.pi*((yd)/366))
    no3df_LOrivbio = no3df_LOrivbio.fillna(5)
    nh4df_LOrivbio = nh4df_LOrivbio.fillna(0)

    # calculate total DIN concentration in mg/L
    moh20_dindf_wwtps = (moh20_no3df_wwtps + moh20_nh4df_wwtps)/71.4    # mg/L
    was24_dindf_wwtps = (was24_no3df_wwtps + was24_nh4df_wwtps)/71.4    # mg/L

    # calculate total NO3 and NH4 concentrations in mg/L
    moh20_no3df_wwtps = (moh20_no3df_wwtps)/71.4    # mg/L
    was24_no3df_wwtps = (was24_no3df_wwtps)/71.4    # mg/L
    moh20_nh4df_wwtps = (moh20_nh4df_wwtps)/71.4    # mg/L
    was24_nh4df_wwtps = (was24_nh4df_wwtps)/71.4    # mg/L

    no3df_trivs = (no3df_trivs)/71.4    # mg/L
    nh4df_trivs = (nh4df_trivs)/71.4    # mg/L
    no3df_LOrivbio = (no3df_LOrivbio)/71.4    # mg/L
    nh4df_LOrivbio = (nh4df_LOrivbio)/71.4    # mg/L

    # calculate daily loading timeseries in kg/d
    moh20_dailyloaddf_wwtps = 86.4*moh20_dindf_wwtps*moh20_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s
    was24_dailyloaddf_wwtps = 86.4*was24_dindf_wwtps*was24_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s

    moh20_no3_dailyloaddf_wwtps = 86.4*moh20_no3df_wwtps*moh20_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s
    was24_no3_dailyloaddf_wwtps = 86.4*was24_no3df_wwtps*was24_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s
    moh20_nh4_dailyloaddf_wwtps = 86.4*moh20_nh4df_wwtps*moh20_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s
    was24_nh4_dailyloaddf_wwtps = 86.4*was24_nh4df_wwtps*was24_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s

    no3_dailyloaddf_triv = 86.4*no3df_trivs*flowdf_trivs # kg/d = 86.4 * mg/L * m3/s
    nh4_dailyloaddf_triv = 86.4*nh4df_trivs*nh4df_trivs # kg/d = 86.4 * mg/L * m3/s
    no3_dailyloaddf_LOrivbio = 86.4*no3df_LOrivbio*flowdf_LOrivflo # kg/d = 86.4 * mg/L * m3/s
    nh4_dailyloaddf_LOrivbio = 86.4*nh4df_LOrivbio*flowdf_LOrivflo # kg/d = 86.4 * mg/L * m3/s

    # calculate average daily load over the year (kg/d)
    moh20_avgload_wwtps = moh20_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')
    was24_avgload_wwtps = was24_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)') 

    moh20_no3_avgload_wwtps = moh20_no3_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-no3-load(kg/d)')
    was24_no3_avgload_wwtps = was24_no3_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-no3-load(kg/d)') 
    moh20_nh4_avgload_wwtps = moh20_nh4_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-nh4-load(kg/d)')
    was24_nh4_avgload_wwtps = was24_nh4_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-nh4-load(kg/d)')

    moh20_avgload_wwtps = moh20_avgload_wwtps.join(moh20_no3_avgload_wwtps['avg-daily-no3-load(kg/d)'])
    moh20_avgload_wwtps = moh20_avgload_wwtps.join(moh20_nh4_avgload_wwtps['avg-daily-nh4-load(kg/d)'])
    was24_avgload_wwtps = was24_avgload_wwtps.join(was24_no3_avgload_wwtps['avg-daily-no3-load(kg/d)'])
    was24_avgload_wwtps = was24_avgload_wwtps.join(was24_nh4_avgload_wwtps['avg-daily-nh4-load(kg/d)'])

    avgload_trivs = no3_dailyloaddf_triv.mean(axis=0).to_frame(name='avg-daily-no3-load(kg/d)')
    nh4_avgload_triv = nh4_dailyloaddf_triv.mean(axis=0).to_frame(name='avg-daily-nh4-load(kg/d)')
    avgload_trivs = avgload_trivs.join(nh4_avgload_triv['avg-daily-nh4-load(kg/d)'])

    avgload_LOriv = no3_dailyloaddf_LOrivbio.mean(axis=0).to_frame(name='avg-daily-no3-load(kg/d)')
    nh4_avgload_LOrivbio = nh4_dailyloaddf_LOrivbio.mean(axis=0).to_frame(name='avg-daily-nh4-load(kg/d)')
    avgload_LOriv = avgload_LOriv.join(nh4_avgload_LOrivbio['avg-daily-nh4-load(kg/d)'])

    # add row and col index for plotting on LiveOcean grid
    griddf0_wwtps = pd.read_csv('../../../LO_data/grids/cas7/moh20_wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    moh20_avgload_wwtps = moh20_avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    moh20_avgload_wwtps = moh20_avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    griddf0_wwtps = pd.read_csv('../../../LO_data/grids/cas7/was24_wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    was24_avgload_wwtps = was24_avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    was24_avgload_wwtps = was24_avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    griddf0_triv = pd.read_csv('../../../LO_data/grids/cas7/triv_info.csv')
    griddf0_triv = griddf0_triv.set_index('rname') # use river name as index
    avgload_trivs = avgload_trivs.join(griddf0_triv['row_py']) # add row to avg load df (uses rname to index)
    avgload_trivs = avgload_trivs.join(griddf0_triv['col_py']) # do the same for cols
    avgload_trivs = avgload_trivs.dropna() # drop Willamette R

    griddf0_LOriv = pd.read_csv('../../../LO_data/grids/cas7/river_info.csv')
    griddf0_LOriv = griddf0_LOriv.set_index('rname') # use river name as index
    griddf0_LOriv.index = [LO2SSM_name(rname) if rname in LO_repeats else rname for rname in griddf0_LOriv.index]
    avgload_LOriv = avgload_LOriv.join(griddf0_LOriv['row_py']) # add row to avg load df (uses rname to index)
    avgload_LOriv = avgload_LOriv.join(griddf0_LOriv['col_py']) # do the same for cols
    avgload_LOriv['row_py'] = avgload_LOriv['row_py'].fillna(678)
    avgload_LOriv['col_py'] = avgload_LOriv['col_py'].fillna(492)

    # # get point source lat and lon
    # moh20_lon_wwtps = [X[int(col)] for col in moh20_avgload_wwtps['col_py']]
    # moh20_lat_wwtps = [Y[int(row)] for row in moh20_avgload_wwtps['row_py']]
    # was24_lon_wwtps = [X[int(col)] for col in was24_avgload_wwtps['col_py']]
    # was24_lat_wwtps = [Y[int(row)] for row in was24_avgload_wwtps['row_py']]

    # add lat and lon to df
    moh20_avgload_wwtps['lon'] = [X[int(col)] for col in moh20_avgload_wwtps['col_py']]
    moh20_avgload_wwtps['lat'] = [Y[int(row)] for row in moh20_avgload_wwtps['row_py']]
    was24_avgload_wwtps['lon'] = [X[int(col)] for col in was24_avgload_wwtps['col_py']]
    was24_avgload_wwtps['lat'] = [Y[int(row)] for row in was24_avgload_wwtps['row_py']]

    avgload_trivs['lon'] = [X[int(col)] for col in avgload_trivs['col_py']]
    avgload_trivs['lat'] = [Y[int(row)] for row in avgload_trivs['row_py']]

    avgload_LOriv['lon'] = [X[int(col)] for col in avgload_LOriv['col_py']]
    avgload_LOriv['lat'] = [Y[int(row)] for row in avgload_LOriv['row_py']]

    # remove WWTPs outside of Puget Sound
    lat_max = 48.49
    lon_min = -123.15
    moh20_avgload_wwtps = moh20_avgload_wwtps[
        (moh20_avgload_wwtps['lat'] <= lat_max) &
        (moh20_avgload_wwtps['lon'] >= lon_min)]
    was24_avgload_wwtps = was24_avgload_wwtps[
        (was24_avgload_wwtps['lat'] <= lat_max) &
        (was24_avgload_wwtps['lon'] >= lon_min)]
    # remove rivers outside of Puget Sound
    avgload_trivs = avgload_trivs[
        (avgload_trivs['lat'] <= lat_max) &
        (avgload_trivs['lon'] >= lon_min)]
    avgload_LOriv = avgload_LOriv[
        (avgload_LOriv['lat'] <= lat_max) &
        (avgload_LOriv['lon'] >= lon_min)]
    # remove WWTPs outside of Puget Sound
    lat_max = 48.5
    lat_min = 47.97
    lon_max = -122.84
    lon_min = -123.3
    moh20_avgload_wwtps = moh20_avgload_wwtps[
        ~((moh20_avgload_wwtps['lat'] <= lat_max) &
        (moh20_avgload_wwtps['lat'] >= lat_min) &
        (moh20_avgload_wwtps['lon'] <= lon_max) &
        (moh20_avgload_wwtps['lon'] >= lon_min))]
    was24_avgload_wwtps = was24_avgload_wwtps[
        ~((was24_avgload_wwtps['lat'] <= lat_max) &
        (was24_avgload_wwtps['lat'] >= lat_min) &
        (was24_avgload_wwtps['lon'] <= lon_max)) &
        (was24_avgload_wwtps['lon'] >= lon_min)]
    # remove rivers outside of Puget Sound
    avgload_trivs = avgload_trivs[
        ~((avgload_trivs['lat'] <= lat_max) &
        (avgload_trivs['lat'] >= lat_min) &
        (avgload_trivs['lon'] <= lon_max) &
        (avgload_trivs['lon'] >= lon_min))]
    avgload_LOriv = avgload_LOriv[
        ~((avgload_LOriv['lat'] <= lat_max) &
        (avgload_LOriv['lat'] >= lat_min) &
        (avgload_LOriv['lon'] <= lon_max) &
        (avgload_LOriv['lon'] >= lon_min))]
    # remove WWTPs outside of Puget Sound
    lat_max = 48.15
    lat_min = 48.14
    lon_max = -122.78
    lon_min = -122.8
    moh20_avgload_wwtps = moh20_avgload_wwtps[
        ~((moh20_avgload_wwtps['lat'] <= lat_max) &
        (moh20_avgload_wwtps['lat'] >= lat_min) &
        (moh20_avgload_wwtps['lon'] <= lon_max) &
        (moh20_avgload_wwtps['lon'] >= lon_min))]
    was24_avgload_wwtps = was24_avgload_wwtps[
        ~((was24_avgload_wwtps['lat'] <= lat_max) &
        (was24_avgload_wwtps['lat'] >= lat_min) &
        (was24_avgload_wwtps['lon'] <= lon_max) &
        (was24_avgload_wwtps['lon'] >= lon_min))]
    # remove WWTPs outside of Puget Sound
    lat_max = 48.40
    lat_min = 48.35
    lon_max = -122.66
    lon_min = -122.67
    moh20_avgload_wwtps = moh20_avgload_wwtps[
        ~((moh20_avgload_wwtps['lat'] <= lat_max) &
        (moh20_avgload_wwtps['lat'] >= lat_min) &
        (moh20_avgload_wwtps['lon'] <= lon_max) &
        (moh20_avgload_wwtps['lon'] >= lon_min))]
    was24_avgload_wwtps = was24_avgload_wwtps[
        ~((was24_avgload_wwtps['lat'] <= lat_max) &
        (was24_avgload_wwtps['lat'] >= lat_min) &
        (was24_avgload_wwtps['lon'] <= lon_max) &
        (was24_avgload_wwtps['lon'] >= lon_min))]


    # get point source lat and lon
    moh20_lon_wwtps = moh20_avgload_wwtps['lon']
    moh20_lat_wwtps = moh20_avgload_wwtps['lat']
    was24_lon_wwtps = was24_avgload_wwtps['lon']
    was24_lat_wwtps = was24_avgload_wwtps['lat']

    triv_lon = avgload_trivs['lon']
    triv_lat = avgload_trivs['lat']

    LOriv_lon = avgload_LOriv['lon']
    LOriv_lat = avgload_LOriv['lat']
    
    # define marker sizes (minimum size is 10 so dots don't get too small)
    moh20_sizes_wwtps = [max(0.05*load,5) for load in moh20_avgload_wwtps['avg-daily-load(kg/d)']]
    was24_sizes_wwtps = [max(0.05*load,5) for load in was24_avgload_wwtps['avg-daily-load(kg/d)']]

    avgload_trivs['avg-daily-load(kg/d)'] = avgload_trivs['avg-daily-no3-load(kg/d)'] + avgload_trivs['avg-daily-nh4-load(kg/d)']
    avgload_LOriv['avg-daily-load(kg/d)'] = avgload_LOriv['avg-daily-no3-load(kg/d)'] + avgload_LOriv['avg-daily-nh4-load(kg/d)']
    triv_sizes = [max(0.05*load,5) for load in avgload_trivs['avg-daily-load(kg/d)']]
    LOriv_sizes = [max(0.05*load,5) for load in avgload_LOriv['avg-daily-load(kg/d)']]

    # print average daily WWTP load to Puget Sound
    total_wwtp_DIN_load = np.sum(moh20_avgload_wwtps['avg-daily-load(kg/d)']) + np.sum(was24_avgload_wwtps['avg-daily-load(kg/d)'])
    total_wwtp_NO3_load = np.sum(moh20_avgload_wwtps['avg-daily-no3-load(kg/d)']) + np.sum(was24_avgload_wwtps['avg-daily-no3-load(kg/d)'])
    total_wwtp_NH4_load = np.sum(moh20_avgload_wwtps['avg-daily-nh4-load(kg/d)']) + np.sum(was24_avgload_wwtps['avg-daily-nh4-load(kg/d)'])

    total_river_NO3_load = np.sum(avgload_trivs['avg-daily-no3-load(kg/d)']) + np.sum(avgload_LOriv['avg-daily-no3-load(kg/d)'])
    total_river_NH4_load = np.sum(avgload_trivs['avg-daily-nh4-load(kg/d)']) + np.sum(avgload_LOriv['avg-daily-nh4-load(kg/d)'])
    total_river_DIN_load = total_river_NO3_load + total_river_NH4_load

    print('Note: TN from WWTPs & rivers is just NO3 + NH4; but TN for the whole Puget Sound is NPZD.')
    print('----------')
    print('Average daily WWTP TN load to Puget Sound: {} kg/day'.format(round(total_wwtp_DIN_load,2)))
    print('Average daily WWTP NO3 load to Puget Sound: {} kg/day'.format(round(total_wwtp_NO3_load,2)))
    print('Average daily WWTP NH4 load to Puget Sound: {} kg/day'.format(round(total_wwtp_NH4_load,2)))
    print('---------')
    print('Average daily river TN load to Puget Sound: {} kg/day'.format(round(total_river_DIN_load,2)))
    print('Average daily river NO3 load to Puget Sound: {} kg/day'.format(round(total_river_NO3_load,2)))
    print('Average daily river NH4 load to Puget Sound: {} kg/day'.format(round(total_river_NH4_load,2)))
    print('---------')

    loading_din = total_river_DIN_load+total_wwtp_DIN_load
    noloading_din = total_river_DIN_load

    loading_no3 = total_river_NO3_load+total_wwtp_NO3_load
    noloading_no3 = total_river_NO3_load

    loading_nh4 = total_river_NH4_load+total_wwtp_NH4_load
    noloading_nh4 = total_river_NH4_load

    print('Land-based TN loads decreased by: {} perc'.format( round( (loading_din-noloading_din)/loading_din * 100 ,2)))
    print('Land-based NO3 loads decreased by: {} perc'.format( round( (loading_no3-noloading_no3)/loading_no3 * 100 ,2)))
    print('Land-based NH4 loads decreased by: {} perc'.format( round( (loading_nh4-noloading_nh4)/loading_nh4 * 100 ,2)))


##############################################################
##                      PROCESS DATA                        ##
##############################################################

# read in masks
basin_mask_ds = grid_ds = xr.open_dataset('../../../LO_output/chapter_2/data/basin_masks_from_pugetsoundDObox.nc')
mask_rho = basin_mask_ds.mask_rho.values
mask_hc = basin_mask_ds.mask_hoodcanal.values
mask_ss = basin_mask_ds.mask_southsound.values
mask_wb = basin_mask_ds.mask_whidbeybasin.values
mask_mb = basin_mask_ds.mask_mainbasin.values
lon = basin_mask_ds['lon_rho'].values
lat = basin_mask_ds['lat_rho'].values
h = basin_mask_ds['h'].values
plon, plat = pfun.get_plon_plat(lon,lat)

##############################################################
# get average concentration per basin

# open datasets
if remove_straits:
    straits = 'noStraits'
else:
    straits = 'withStraits'

# initialize empty dictionaries and fill with vertical integrals
NO3_vert_dict = {}
phyto_vert_dict = {}
zoop_vert_dict = {}
NH4_vert_dict = {}
Ldet_vert_dict = {}
Sdet_vert_dict = {}
DO_vert_dict = {}

for year in years:
    for gtagex in gtagexes:
        ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / (gtagex + '_pugetsoundDO_' + year + '_NPZD_vert_ints_' + straits + '.nc'))
        NO3_vert_int = ds['NO3_vert_int'].values
        phyto_vert_int = ds['phyto_vert_int'].values
        zoop_vert_int = ds['zoop_vert_int'].values
        NH4_vert_int = ds['NH4_vert_int'].values
        LdetritusN_vert_int = ds['LdetritusN_vert_int'].values
        SdetritusN_vert_int = ds['SdetritusN_vert_int'].values
        DO_vert_int = ds['DO_vert_int'].values
        # add data to dictionaries
        NO3_vert_dict[gtagex+year] = NO3_vert_int
        phyto_vert_dict[gtagex+year] = phyto_vert_int
        zoop_vert_dict[gtagex+year] = zoop_vert_int
        NH4_vert_dict[gtagex+year] = NH4_vert_int
        Ldet_vert_dict[gtagex+year] = LdetritusN_vert_int
        Sdet_vert_dict[gtagex+year] = SdetritusN_vert_int
        DO_vert_dict[gtagex+year] = DO_vert_int

# grid cell areas
fp = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
box_ds = xr.open_dataset(fp)
DX = (box_ds.pm.values)**-1
DY = (box_ds.pn.values)**-1
DA = DX*DY # get area in m2


# initialize dictionary for average concentration (volume integrals [mol], normalized by volume)
NO3_vol_norm = {}
phyto_vol_norm = {}
zoop_vol_norm = {}
NH4_vol_norm = {}
Ldet_vol_norm = {}
Sdet_vol_norm = {}
DO_vol_norm = {}

for year in years:
    for region in regions:

        # get mask for the region
        if region == 'Hood Canal':
            mask = mask_hc
        elif region == 'South Sound':
            mask = mask_ss
        elif region == 'Whidbey Basin':
            mask = mask_wb
        elif region == 'Main Basin':
            mask = mask_mb
        elif region == 'All Puget Sound':
            mask = mask_hc + mask_ss + mask_wb + mask_mb

        # basin volume
        h_masked = h * mask
        basin_vol = np.sum(h_masked * DA) # [m3]

        for gtagex in gtagexes:

            NO3_vert_int = NO3_vert_dict[gtagex+year]
            NO3_vert_int_masked = NO3_vert_int * mask
            NO3_vol_timeseries = np.sum(NO3_vert_int_masked * DA, axis=(1, 2)) # [mol]
            NO3_vol_norm[gtagex+region+year] = NO3_vol_timeseries

            NH4_vert_int = NH4_vert_dict[gtagex+year]
            NH4_vert_int_masked = NH4_vert_int * mask
            NH4_vol_timeseries = np.sum(NH4_vert_int_masked * DA, axis=(1, 2)) # [mol]
            NH4_vol_norm[gtagex+region+year] = NH4_vol_timeseries

            phyto_vert_int = phyto_vert_dict[gtagex+year]
            phyto_vert_int_masked = phyto_vert_int * mask
            phyto_vol_timeseries = np.sum(phyto_vert_int_masked * DA, axis=(1, 2)) # [mol]
            phyto_vol_norm[gtagex+region+year] = phyto_vol_timeseries

            zoop_vert_int = zoop_vert_dict[gtagex+year]
            zoop_vert_int_masked = zoop_vert_int * mask
            zoop_vol_timeseries = np.sum(zoop_vert_int_masked * DA, axis=(1, 2)) # [mol]
            zoop_vol_norm[gtagex+region+year] = zoop_vol_timeseries


            Ldet_vert_int = Ldet_vert_dict[gtagex+year]
            Ldet_vert_int_masked = Ldet_vert_int * mask
            Ldet_vol_timeseries = np.sum(Ldet_vert_int_masked * DA, axis=(1, 2)) # [mol]
            Ldet_vol_norm[gtagex+region+year] = Ldet_vol_timeseries

            Sdet_vert_int = Sdet_vert_dict[gtagex+year]
            Sdet_vert_int_masked = Sdet_vert_int * mask
            Sdet_vol_timeseries = np.sum(Sdet_vert_int_masked * DA, axis=(1, 2)) # [mol]
            Sdet_vol_norm[gtagex+region+year] = Sdet_vol_timeseries

            DO_vert_int = DO_vert_dict[gtagex+year]
            DO_vert_int_masked = DO_vert_int * mask
            DO_vol_timeseries = np.sum(DO_vert_int_masked * DA, axis=(1, 2)) # [mol]
            DO_vol_norm[gtagex+region+year] = DO_vol_timeseries


NO3_timeseries_noloading = []
NH4_timeseries_noloading = []
phyto_timeseries_noloading = []
zoop_timeseries_noloading = []
Ldet_timeseries_noloading = []
Sdet_timeseries_noloading = []
DO_timeseries_noloading = []

NO3_timeseries_loading = []
NH4_timeseries_loading = []
phyto_timeseries_loading = []
zoop_timeseries_loading = []
Ldet_timeseries_loading = []
Sdet_timeseries_loading = []
DO_timeseries_loading = []

# get full time series
for year in years:
    region = 'All Puget Sound'

    NO3_timeseries_noloading.extend(NO3_vol_norm['cas7_t1noDIN_x11ab'+region+year])
    NH4_timeseries_noloading.extend(NH4_vol_norm['cas7_t1noDIN_x11ab'+region+year])
    phyto_timeseries_noloading.extend(phyto_vol_norm['cas7_t1noDIN_x11ab'+region+year])
    zoop_timeseries_noloading.extend(zoop_vol_norm['cas7_t1noDIN_x11ab'+region+year])
    Sdet_timeseries_noloading.extend(Sdet_vol_norm['cas7_t1noDIN_x11ab'+region+year])
    Ldet_timeseries_noloading.extend(Ldet_vol_norm['cas7_t1noDIN_x11ab'+region+year])
    DO_timeseries_noloading.extend(DO_vol_norm['cas7_t1noDIN_x11ab'+region+year])

    NO3_timeseries_loading.extend(NO3_vol_norm['cas7_t1_x11ab'+region+year])
    NH4_timeseries_loading.extend(NH4_vol_norm['cas7_t1_x11ab'+region+year])
    phyto_timeseries_loading.extend(phyto_vol_norm['cas7_t1_x11ab'+region+year])
    zoop_timeseries_loading.extend(zoop_vol_norm['cas7_t1_x11ab'+region+year])
    Sdet_timeseries_loading.extend(Sdet_vol_norm['cas7_t1_x11ab'+region+year])
    Ldet_timeseries_loading.extend(Ldet_vol_norm['cas7_t1_x11ab'+region+year])
    DO_timeseries_loading.extend(DO_vol_norm['cas7_t1_x11ab'+region+year])

# get average molar quantity
NO3_mols_daily_avg_noloading = np.nanmean(NO3_timeseries_noloading)
NH4_mols_daily_avg_noloading = np.nanmean(NH4_timeseries_noloading)
phyto_mols_daily_avg_noloading = np.nanmean(phyto_timeseries_noloading)
zoop_mols_daily_avg_noloading = np.nanmean(zoop_timeseries_noloading)
Sdet_mols_daily_avg_noloading = np.nanmean(Sdet_timeseries_noloading)
Ldet_mols_daily_avg_noloading = np.nanmean(Ldet_timeseries_noloading)
DO_mols_daily_avg_noloading = np.nanmean(DO_timeseries_noloading)

NO3_mols_daily_avg_loading = np.nanmean(NO3_timeseries_loading)
NH4_mols_daily_avg_loading = np.nanmean(NH4_timeseries_loading)
phyto_mols_daily_avg_loading = np.nanmean(phyto_timeseries_loading)
zoop_mols_daily_avg_loading = np.nanmean(zoop_timeseries_loading)
Sdet_mols_daily_avg_loading = np.nanmean(Sdet_timeseries_loading)
Ldet_mols_daily_avg_loading = np.nanmean(Ldet_timeseries_loading)
DO_mols_daily_avg_loading = np.nanmean(DO_timeseries_loading)

# get total TN
TN_noloading = (NO3_mols_daily_avg_noloading +
                NH4_mols_daily_avg_noloading +
                phyto_mols_daily_avg_noloading +
                zoop_mols_daily_avg_noloading +
                Sdet_mols_daily_avg_noloading +
                Ldet_mols_daily_avg_noloading)
TN_loading = (NO3_mols_daily_avg_loading +
                NH4_mols_daily_avg_loading +
                phyto_mols_daily_avg_loading +
                zoop_mols_daily_avg_loading +
                Sdet_mols_daily_avg_loading +
                Ldet_mols_daily_avg_loading)

print('----------')
print('Vol-integrated TN in Puget Sound decreased by: {} perc'.format( round(
    (TN_loading-TN_noloading)/TN_noloading * 100 ,2)))
print('Vol-integrated NO3 in Puget Sound decreased by: {} perc'.format( round(
    (NO3_mols_daily_avg_loading-NO3_mols_daily_avg_noloading)/NO3_mols_daily_avg_noloading * 100 ,2)))
print('Vol-integrated NH4 in Puget Sound decreased by: {} perc'.format( round(
    (NH4_mols_daily_avg_loading-NH4_mols_daily_avg_noloading)/NH4_mols_daily_avg_noloading * 100 ,2)))
print('Vol-integrated DO in Puget Sound decreased by: {} perc'.format( round(
    (DO_mols_daily_avg_noloading-DO_mols_daily_avg_loading)/DO_mols_daily_avg_noloading * 100 ,2)))




##############################################################
##                    Plot basin map                        ##
##############################################################

# Puget Sound bounds
xmin = -123.29 #-125
xmax = -122.1
ymin = 46.95
ymax = 48.50 #49.5


# initialize figure
fig, ax0 = plt.subplots(1,1,figsize = (8,8))

# # Hood Canal
# ax0.pcolormesh(plon, plat, np.where(mask_hc == 0, np.nan, mask_hc),
#             vmin=0, vmax=2.5, cmap='RdPu', alpha=0.4 )
# # South Sound
# ax0.pcolormesh(plon, plat, np.where(mask_ss == 0, np.nan, mask_ss),
#             vmin=0, vmax=2, cmap='Purples', alpha=0.4  )
# # Whidbey Basin
# ax0.pcolormesh(plon, plat, np.where(mask_wb == 0, np.nan, mask_wb),
#             vmin=0, vmax=3, cmap='cool', alpha=0.4  )
# # Main Basin
# ax0.pcolormesh(plon, plat, np.where(mask_mb == 0, np.nan, mask_mb),
#             vmin=0, vmax=1.5, cmap='summer', alpha=0.4  )
# Full region
ax0.pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho),
            vmin=0, vmax=3, cmap='Greys', alpha=0.4 )
# Hood Canal
ax0.pcolormesh(plon, plat, np.where(mask_hc == 0, np.nan, mask_hc),
            vmin=0, vmax=1.5, cmap='Greys', alpha=0.4 )
# South Sound
ax0.pcolormesh(plon, plat, np.where(mask_ss == 0, np.nan, mask_ss),
            vmin=0, vmax=1.5, cmap='Greys', alpha=0.4  )
# Whidbey Basin
ax0.pcolormesh(plon, plat, np.where(mask_wb == 0, np.nan, mask_wb),
            vmin=0, vmax=1.5, cmap='Greys', alpha=0.4  )
# Main Basin
ax0.pcolormesh(plon, plat, np.where(mask_mb == 0, np.nan, mask_mb),
            vmin=0, vmax=1.5, cmap='Greys', alpha=0.4  )


# format figure
ax0.set_xlim([xmin,xmax])
ax0.set_ylim([ymin,ymax])
ax0.set_ylabel('Latitude', fontsize=12)
ax0.set_xlabel('Longitude', fontsize=12)
ax0.tick_params(axis='both', labelsize=12, rotation=30)
# ax0.set_title('(a) Basins', loc='left', fontsize=14, fontweight='bold')
pfun.dar(ax0)

# add wwtp locations
if WWTP_loc == True:
    edgecolor = 'crimson'
    facecolor = 'crimson'
    alpha = 0.5
    ax0.scatter(moh20_lon_wwtps,moh20_lat_wwtps,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                    linewidth=1, s=moh20_sizes_wwtps, label='WWTPs')
    ax0.scatter(was24_lon_wwtps,was24_lat_wwtps,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                    linewidth=1, s=was24_sizes_wwtps)
    leg_szs = [100, 1000, 10000]
    szs = [0.05*(leg_sz) for leg_sz in leg_szs]
    l0 = plt.scatter([],[], s=szs[0], color=facecolor, alpha=alpha, edgecolors=edgecolor, linewidth=1)
    l1 = plt.scatter([],[], s=szs[1], color=facecolor, alpha=alpha, edgecolors=edgecolor, linewidth=1)
    l2 = plt.scatter([],[], s=szs[2], color=facecolor, alpha=alpha, edgecolors=edgecolor, linewidth=1)
    labels = ['< 100', '1,000', '10,000']
    legend = ax0.legend([l0, l1, l2], labels, fontsize = 10, markerfirst=False,
        title='Loading \n'+r' (kg N d$^{-1}$)',loc='upper left', labelspacing=1, borderpad=0.8)
    plt.setp(legend.get_title(),fontsize=9)

    edgecolor = 'blue'
    facecolor = 'blue'
    alpha = 0.5
    ax0.scatter(triv_lon,triv_lat,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                    linewidth=1, s=triv_sizes, label='Rivers')
    ax0.scatter(LOriv_lon,LOriv_lat,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                    linewidth=1, s=LOriv_sizes)
    
    ax0.text(-123.2,48.05,'WWTPs',color='crimson',fontsize=12,fontweight='bold')
    ax0.text(-123.2,48.0,'Rivers',color='blue',fontsize=12,fontweight='bold')

# ##############################################################
# ##    Sub-basins and multiple years and percent change      ##
# ##############################################################

# # plot timeseries
#     for year in years:
#         # create time vector
#         startdate = year+'.01.01'
#         enddate = year+'.12.31'
#         dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
#         dates_local = [pfun.get_dt_local(x) for x in dates]
#         # loop through model runs
#         for k,region in enumerate(regions):
#             # get data for the basin and gtagex and year
#             gtagex = 'cas7_t1noDIN_x11ab'
#             avg_concentration = var_vol_norm[gtagex+region+year]
#             # pass through hanning window
#             avg_concentration_filtered = zfun.lowpass(avg_concentration,n=nwin)
#             # plot data
#             if region == 'All Puget Sound':
#                 ax1.plot(dates_local,avg_concentration_filtered,linewidth=1,
#                         color=colors[k], alpha=1,linestyle='--')
#             else:
#                 ax1.plot(dates_local,avg_concentration_filtered,linewidth=2,
#                     color=colors[k], alpha=0.8)
                
#     # format figure
#     ax1.grid(visible=True, axis='both', color='silver', linestyle='--')
#     ax1.set_xticklabels([])
#     ax1.tick_params(axis='y', labelsize=12, rotation=30)
#     ax1.set_ylabel(r'mmol/m3', fontsize=12)
#     # create time vector
#     startdate = years[0]+'.01.01'
#     enddate = years[-1]+'.12.31'
#     dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
#     dates_local = [pfun.get_dt_local(x) for x in dates]
#     ax1.set_xlim([dates_local[0],dates_local[-1]])
#     ax1.xaxis.set_major_locator(mdates.YearLocator())
#     ax1.set_title('(b) No-Loading Avg. Conc. ({}-day Hanning Window)'.format(nwin), loc='left', fontsize=14, fontweight='bold')

#     # set y-lims
#     if vars[j] == 'NO3':
#         ymax_conc = 35
#     elif vars[j] == 'NH4':
#         ymax_conc = 1.75
#     elif vars[j] == 'Phytoplankton':
#         ymax_conc = 1.75
#     elif vars[j] == 'Zooplankton':
#         ymax_conc = 0.175
#     elif vars[j] == 'Large Detritus':
#         ymax_conc = 0.04
#     elif vars[j] == 'Small Detritus':
#         ymax_conc = 1.2
#     elif vars[j] == 'Dissolved Oxygen':
#         ymax_conc = 300
#     ax1.set_ylim([0,ymax_conc])

#     # add difference plot -----------------------------------------------------------
#     divider = make_axes_locatable(ax1)
#     ax2 = divider.append_axes("bottom", size='50%', pad=0.4)
#     ax1.figure.add_axes(ax2)
#     ax2.set_xlim([dates_local[0],dates_local[-1]])

#     # plot difference
#     for year in years:
#         # create time vector
#         startdate = year+'.01.01'
#         enddate = year+'.12.31'
#         dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
#         dates_local = [pfun.get_dt_local(x) for x in dates]
#         # loop through model runs
#         for k,region in enumerate(regions):
#             # get data for the basin and gtagex and year
#             noloading = 'cas7_t1noDIN_x11ab'
#             avg_concentration_noloading = var_vol_norm[noloading+region+year]
#             loading = 'cas7_t1_x11ab'
#             avg_concentration_loading = var_vol_norm[loading+region+year]
#             # calculate difference
#             diff =  avg_concentration_loading - avg_concentration_noloading
#             # pass through hanning window
#             diff_filtered = zfun.lowpass(diff,n=nwin)
#             # plot data
#             if region == 'All Puget Sound':
#                 ax2.plot(dates_local,diff_filtered,linewidth=1,
#                         color=colors[k], alpha=1,linestyle='--')
#             else:
#                 ax2.plot(dates_local,diff_filtered,linewidth=2,
#                     color=colors[k], alpha=0.8)
                
#     # format figure
#     ax2.grid(visible=True, axis='x', color='silver', linestyle='--')
#     ax2.tick_params(axis='both', labelsize=12, rotation=30)
#     ax2.set_ylabel(r'mmol/m3', fontsize=12)
#     # create time vector
#     startdate = years[0]+'.01.01'
#     enddate = years[-1]+'.12.31'
#     dates = pd.date_range(start= startdate, end= enddate, freq= '1d')
#     dates_local = [pfun.get_dt_local(x) for x in dates]
#     ax2.set_xlim([dates_local[0],dates_local[-1]])
#     ax2.hlines(y=0, xmin=dates_local[0], xmax=dates_local[-1],
#                color='silver', linestyle='--', linewidth=0.75)
#     ax2.xaxis.set_major_locator(mdates.YearLocator())
#     ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
#     ax2.set_title('(c) Loading - No-Loading ({}-day Hanning Window)'.format(nwin), loc='left', fontsize=14, fontweight='bold')

#     # set y-lims
#     if vars[j] == 'NO3':
#         ymin_diff = -4
#     elif vars[j] == 'NH4':
#         ymin_diff = -0.4
#     elif vars[j] == 'Phytoplankton':
#         ymin_diff = -0.12
#     elif vars[j] == 'Zooplankton':
#         ymin_diff = -0.01
#     elif vars[j] == 'Large Detritus':
#         ymin_diff = -0.005
#     elif vars[j] == 'Small Detritus':
#         ymin_diff = -0.07
#     elif vars[j] == 'Dissolved Oxygen':
#         ymin_diff = -4
#     ymax_diff = ymin_diff*-1
#     ax2.set_ylim([ymin_diff,ymax_diff])

#     plt.tight_layout()
#     plt.show()