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
    nh4_dailyloaddf_triv = 86.4*nh4df_trivs*flowdf_trivs # kg/d = 86.4 * mg/L * m3/s
    no3_dailyloaddf_LOrivbio = 86.4*no3df_LOrivbio*flowdf_LOrivflo # kg/d = 86.4 * mg/L * m3/s
    nh4_dailyloaddf_LOrivbio = 86.4*nh4df_LOrivbio*flowdf_LOrivflo # kg/d = 86.4 * mg/L * m3/s

    # total_m3_s = (np.nansum(moh20_flowdf_wwtps.mean(axis=0)) +
    #                 np.nansum(was24_flowdf_wwtps.mean(axis=0)) +
    #                 np.nansum(flowdf_trivs.mean(axis=0)) +
    #                 np.nansum(flowdf_LOrivflo.mean(axis=0)))
    # print('Total river+WWTP flow: {} m3/s'.format(round(total_m3_s,2)))

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

    avgflow_moh20wwtp = moh20_flowdf_wwtps.mean(axis=0).to_frame(name='avg-daily-flow(m3/s)')
    moh20_avgload_wwtps = moh20_avgload_wwtps.join(avgflow_moh20wwtp['avg-daily-flow(m3/s)'])
    avgflow_was24wwtp = was24_flowdf_wwtps.mean(axis=0).to_frame(name='avg-daily-flow(m3/s)')
    was24_avgload_wwtps = was24_avgload_wwtps.join(avgflow_was24wwtp['avg-daily-flow(m3/s)'])

    avgload_trivs = no3_dailyloaddf_triv.mean(axis=0).to_frame(name='avg-daily-no3-load(kg/d)')
    nh4_avgload_triv = nh4_dailyloaddf_triv.mean(axis=0).to_frame(name='avg-daily-nh4-load(kg/d)')
    avgload_trivs = avgload_trivs.join(nh4_avgload_triv['avg-daily-nh4-load(kg/d)'])

    avgload_LOriv = no3_dailyloaddf_LOrivbio.mean(axis=0).to_frame(name='avg-daily-no3-load(kg/d)')
    nh4_avgload_LOrivbio = nh4_dailyloaddf_LOrivbio.mean(axis=0).to_frame(name='avg-daily-nh4-load(kg/d)')
    avgload_LOriv = avgload_LOriv.join(nh4_avgload_LOrivbio['avg-daily-nh4-load(kg/d)'])

    avgflow_triv = flowdf_trivs.mean(axis=0).to_frame(name='avg-daily-flow(m3/s)')
    avgload_trivs = avgload_trivs.join(avgflow_triv['avg-daily-flow(m3/s)'])
    avgflow_LOriv = flowdf_LOrivflo.mean(axis=0).to_frame(name='avg-daily-flow(m3/s)')
    avgload_LOriv = avgload_LOriv.join(avgflow_LOriv['avg-daily-flow(m3/s)'])

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
    # remove WWTPs outside of Puget Sound
    lat_min = 47
    moh20_avgload_wwtps = moh20_avgload_wwtps[~(moh20_avgload_wwtps['lat'] <= lat_min)]
    was24_avgload_wwtps = was24_avgload_wwtps[~(was24_avgload_wwtps['lat'] <= lat_min)]
    # remove rivers outside of Puget Sound
    avgload_trivs = avgload_trivs[~(avgload_trivs['lat'] <= lat_min)]
    avgload_LOriv = avgload_LOriv[~(avgload_LOriv['lat'] <= lat_min)]


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

    print('Land-based TN loads increased by: {} perc'.format( round( (loading_din-noloading_din)/noloading_din * 100 ,2)))
    print('Land-based NO3 loads increased by: {} perc'.format( round( (loading_no3-noloading_no3)/noloading_no3 * 100 ,2)))
    print('Land-based NH4 loads increased by: {} perc'.format( round( (loading_nh4-noloading_nh4)/noloading_nh4 * 100 ,2)))

    print('---------')
    total_river_flow = np.sum(avgload_trivs['avg-daily-flow(m3/s)']) + np.sum(avgload_LOriv['avg-daily-flow(m3/s)'])
    print('Total River Flow = {} m3/s'.format(total_river_flow))
    total_wwtp_flow = np.sum(moh20_avgload_wwtps['avg-daily-flow(m3/s)']) + np.sum(was24_avgload_wwtps['avg-daily-flow(m3/s)'])
    print('Total WWTP Flow = {} m3/s'.format(total_wwtp_flow))
    print('Total River+WWTP Flow = {} m3/s'.format(total_river_flow+total_wwtp_flow))

    print('---------')
    # concentrations
    C_loading = total_wwtp_DIN_load / total_wwtp_flow / 86.4 * 71.4 # kg/d / m3/s convert to mmol/m3
    print('Average WWTP DIN concentration = {} mmol/m3'.format(C_loading))

    # concentrations
    C_R_loading = loading_din / (total_river_flow+total_wwtp_flow) / 86.4 * 71.4 # kg/d / m3/s convert to mmol/m3
    C_R_noloading = noloading_din / (total_river_flow) / 86.4 * 71.4 # kg/d / m3/s convert to mmol/m3
    print('Average WWTP DIN concentration (loading) = {} mmol/m3'.format(C_R_loading))
    print('Average WWTP DIN concentration (noloading) = {} mmol/m3'.format(C_R_noloading))
    print('Difference = {} mmol/m3'.format(C_R_loading-C_R_noloading))




##############################################################
##                      PROCESS DATA                        ##
##############################################################

# read in masks
basin_mask_ds = grid_ds = xr.open_dataset('../../../LO_output/chapter_2/data/basin_masks_from_pugetsoundDObox.nc')
mask_rho = basin_mask_ds.mask_rho.values
mask_ps = basin_mask_ds.mask_pugetsound.values
h = basin_mask_ds['h'].values
lon = basin_mask_ds.lon_rho.values
lat = basin_mask_ds.lat_rho.values
plon, plat = pfun.get_plon_plat(lon,lat)

##############################################################
##                    Plot basin map                        ##
##############################################################

# Puget Sound bounds
xmin = -123.29
xmax = -122.1
ymin = 46.95
ymax = 48.50


# initialize figure
fig, axes = plt.subplots(1,3,figsize = (8,4), sharex = True, sharey=True)
ax = axes.ravel()
ax0 = ax[0]
ax1 = ax[1]
ax2 = ax[2]


# add wwtp locations
# Full region
ax0.pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho),
            vmin=0, vmax=45, cmap=cmocean.cm.ice_r)
# Puget Sound
ax0.pcolormesh(plon, plat, np.where(mask_ps == 0, np.nan, mask_ps),
            vmin=0, vmax=8, cmap=cmocean.cm.ice_r)
edgecolor = 'black'
facecolor = 'yellowgreen'
alpha = 0.5
ax0.scatter(moh20_lon_wwtps,moh20_lat_wwtps,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                linewidth=1, s=moh20_sizes_wwtps, label='WWTPs')
ax0.scatter(was24_lon_wwtps,was24_lat_wwtps,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                linewidth=1, s=was24_sizes_wwtps)
leg_szs = [100, 1000, 10000]
szs = [0.05*(leg_sz) for leg_sz in leg_szs]
l0 = plt.scatter([],[], s=szs[0], color='grey', alpha=alpha, edgecolors=edgecolor, linewidth=1)
l1 = plt.scatter([],[], s=szs[1], color='grey', alpha=alpha, edgecolors=edgecolor, linewidth=1)
l2 = plt.scatter([],[], s=szs[2], color='grey', alpha=alpha, edgecolors=edgecolor, linewidth=1)
labels = ['< 100', '1,000', '10,000']
legend = ax0.legend([l0, l1, l2], labels, fontsize = 10, markerfirst=False,
    title='Loading \n'+r' (kg N d$^{-1}$)',loc='upper left', labelspacing=1, borderpad=0.8)
plt.setp(legend.get_title(),fontsize=9)
# format figure
ax0.set_xlim([xmin,xmax])
ax0.set_ylim([ymin,ymax])
ax0.set_ylabel('Latitude', fontsize=12)
ax0.set_xlabel('Longitude', fontsize=12)
ax0.tick_params(axis='both', labelsize=12)
# ax0.set_title('(a) Basins', loc='left', fontsize=14, fontweight='bold')
pfun.dar(ax0)
# add load
ax0.text(-122.15,47.0, f'{int(total_wwtp_DIN_load):,d}' + r' kg d$^{-1}$',
         fontsize=12, fontweight='bold', ha='right')
ax1.set_title('(a) WWTPs',loc='left',fontsize=14,fontweight='bold')


# add river locations
# Full region
ax1.pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho),
            vmin=0, vmax=45, cmap=cmocean.cm.ice_r)
# Puget Sound
ax1.pcolormesh(plon, plat, np.where(mask_ps == 0, np.nan, mask_ps),
            vmin=0, vmax=8, cmap=cmocean.cm.ice_r)
edgecolor = 'black'
facecolor = 'mediumpurple'
alpha = 0.5
ax1.scatter(triv_lon,triv_lat,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                linewidth=1, s=triv_sizes, label='Rivers')
ax1.scatter(LOriv_lon,LOriv_lat,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                linewidth=1, s=LOriv_sizes)
# format figure
ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin,ymax])
ax1.set_xlabel('Longitude', fontsize=12)
ax1.tick_params(axis='both', labelsize=12)
# ax0.set_title('(a) Basins', loc='left', fontsize=14, fontweight='bold')
pfun.dar(ax1)
ax1.text(-122.15,47.0, f'{int(total_river_DIN_load):,d}' + r' kg d$^{-1}$',
         fontsize=12, fontweight='bold', ha='right')
ax1.set_title('(b) Rivers',loc='left',fontsize=14,fontweight='bold')


# add exchange flow at admiralty inlet
# Full region
ax2.pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho),
            vmin=0, vmax=45, cmap=cmocean.cm.ice_r)
# Puget Sound
ax2.pcolormesh(plon, plat, np.where(mask_ps == 0, np.nan, mask_ps),
            vmin=0, vmax=8, cmap=cmocean.cm.ice_r)
edgecolor = 'black'
facecolor = 'cornflowerblue'
alpha = 0.5
ocn_lon = -122.725
ocn_lat = 48.165
# Data from Mackas & Harrison (1997) through Strait of Juan de Fuca
ocn_load = 2600 * 1000 # (tonnes per day * 1000 kg/ton)
ocn_size = 0.05*ocn_load
ax2.scatter(ocn_lon,ocn_lat,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                linewidth=1, s=ocn_size)
# format figure
ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin,ymax])
# ax2.set_xlabel('Longitude', fontsize=12)
ax2.tick_params(axis='both', labelsize=12)
# ax0.set_title('(a) Basins', loc='left', fontsize=14, fontweight='bold')
pfun.dar(ax2)
ax2.text(-122.15,47.0, f'{int(ocn_load):,d}' + r' kg d$^{-1}$',
         fontsize=12, fontweight='bold', ha='right')
ax1.set_title('(c) Exchange Flow',loc='left',fontsize=14,fontweight='bold')

plt.tight_layout()
plt.show()
