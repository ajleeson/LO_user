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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

##############################################################
##                       USER INPUTS                        ##
##############################################################

vn = 'oxygen'

year =  '2017'

# which  model run to look at?
gtagex = 'cas7_t1_x11ab' # long hindcast

plt.close('all')

##############################################################
##                     INTIALIZE STUFF                      ##
##############################################################

# initialize empty dictionaries
avg_DO_bot_dict = {} # dictionary with average bottom DO in each inlet [mg/L]

# list 13 terminal inlets
inlets = ['dyes','sinclair','quartermaster','case','crescent','carr',
          'elliot','commencement','penn','portsusan','holmes','dabob','lynchcove']

# hypoxic season deep DO
hypminday = 242
hypmaxday = 302

# initialize dataframes
inlet_botDO_df = pd.DataFrame(columns=['Inlet', 'SepOctDeepDO[mg/L]', 'SepOctDeepDO_err[mg/L]'])
monthly_mean_df = pd.DataFrame(columns=['month','inlet','DOdeep(mg/L)'])

months = ['Jan','Feb','Mar','Apr','May','Jun',
        'Jul','Aug','Sep','Oct','Nov','Dec']

##############################################################
##          GET OFFSET BETWEEN FULL DOMAIN AND BOX          ##
##############################################################

ds_basin_mask = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / 'basin_masks_from_pugetsoundDObox.nc')
ds_his_file = xr.open_dataset(Ldir['roms_out'] / 'cas7_t1_x11ab' / 'f2017.02.12' / 'ocean_his_0001.nc')

# full grid lat/lon
lat_full = ds_his_file['lat_rho'].values
lon_full = ds_his_file['lon_rho'].values
# box extraction lat/lon
lat_box = ds_basin_mask['lat_rho'].values
lon_box = ds_basin_mask['lon_rho'].values

# get bot lat and lon at the bottom left
lat0 = lat_box[0, 0]
lon0 = lon_box[0, 0]

# get exact matching indices in full grid for reference point [0,0] in box extraction
zz_inds = np.where((lat_full == lat0) & (lon_full == lon0))
# check to make sure there is an exact match, or else throw error
if len(zz_inds[0]) == 0:
    raise ValueError('No exact matches for lat/lon in full domain')

# get value of i and j at from the box [0,0] in the full domain
j0_full = int(zz_inds[0][0])
i0_full = int(zz_inds[1][0])

# offsets so that: sub_index = full_index - offset
j_offset = j0_full
i_offset = i0_full

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# get bottom DO concentrations
ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / (gtagex + '_pugetsoundDO_' + year + '_DO_info.nc'))
DO_bot = ds['DO_bot'].values

# loop through inlets and get average bottom DO
for s,inlet in enumerate(inlets): # stations: 
    # get segment information
    seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
    seg_df = pd.read_pickle(seg_name)
    ji_list = seg_df[inlet+'_p']['ji_list']
    jj = [x[0] for x in ji_list]
    ii = [x[1] for x in ji_list]

    # convert ii and jj indices into the box sub-domain
    jj = [j - j_offset for j in jj]
    ii = [i - i_offset for i in ii]
    
    # get bottom DO at all points in the inlet
    inlet_bot_DO_allpoints = DO_bot[:,jj,ii]

    # TODO: get average bottom DO concentration in each inlet (volume-weighted average)

    # get average bottom DO within the inlet
    avg_DO_bot_dict[inlet] = np.nanmean(inlet_bot_DO_allpoints, axis = 1)

    # get inlet name
    if inlet == 'case':
            inlet_name = 'Case Inlet'
    elif inlet == 'portsusan':
            inlet_name = 'Port Susan'
    elif inlet == 'penn':
            inlet_name = 'Penn Cove'
    elif inlet == 'holmes':
            inlet_name = 'Holmes Harbor'
    elif inlet == 'dabob':
            inlet_name = 'Dabob Bay'
    elif inlet == 'lynchcove':
            inlet_name = 'Lynch Cove'
    elif inlet == 'crescent':
            inlet_name = 'Crescent Harbor'
    elif inlet == 'carr':
            inlet_name = 'Carr Inlet'
    elif inlet == 'dyes':
            inlet_name = 'Dyes Inlet'
    elif inlet == 'sinclair':
            inlet_name = 'Sinclair Inlet'
    elif inlet == 'elliot':
            inlet_name = 'Elliott Bay'
    elif inlet == 'commencement':
            inlet_name = 'Commencement Bay'
    elif inlet == 'quartermaster':
            inlet_name = 'Quartermaster Harbor'

    # average inlet bottom DO during hypoxic season
    DOinlet_avg = np.nanmean(avg_DO_bot_dict[inlet][hypminday:hypmaxday])
    DOinlet_err = np.nanstd(avg_DO_bot_dict[inlet][hypminday:hypmaxday])
    # add data to df
    new_data = {'Inlet': [inlet],
                'SepOctDeepDO[mg/L]': [DOinlet_avg],
                'SepOctDeepDO_err[mg/L]': [DOinlet_err]}
    df_new_rows = pd.DataFrame(new_data)
    inlet_botDO_df = pd.concat([inlet_botDO_df, df_new_rows],ignore_index=True)

    # get monthly means
    for m,month in enumerate(months):
        if m == 0:
            minday = 0 #1
            maxday = 30 #32
        elif m == 1:
            minday = 30# 32
            maxday = 58# 60
        elif m == 2:
            minday = 58# 60
            maxday = 89# 91
        elif m == 3:
            minday = 89# 91
            maxday = 119# 121
        elif m == 4:
            minday = 119# 121
            maxday = 150# 152
        elif m == 5:
            minday = 150# 152
            maxday = 180# 182
        elif m == 6:
            minday = 180# 182
            maxday = 211# 213
        elif m == 7:
            minday = 211# 213
            maxday = 242# 244
        elif m == 8:
            minday = 242# 244
            maxday = 272# 274
        elif m == 9:
            minday = 272# 274
            maxday = 303# 305
        elif m == 10:
            minday = 303# 305
            maxday = 333# 335
        elif m == 11:
            minday = 333# 335
            maxday = 363

        # get monthly  mean botto DO
        DOdeep_monthlymean = np.nanmean(avg_DO_bot_dict[inlet][minday:maxday])

        # add data to dataframe
        # save to dictionary
        if inlet == 'elliot':
            station_name = 'elliott' # correct typo
        else:
            station_name = inlet
        new_data = {'month': [months[m]], 
                    'inlet':   [station_name], 
                    'DOdeep(mg/L)':  [DOdeep_monthlymean]}
        df_new_rows = pd.DataFrame(new_data)
        monthly_mean_df = pd.concat([monthly_mean_df, df_new_rows],ignore_index=True)


# save to csv file
print(inlet_botDO_df)
inlet_botDO_df.to_csv('inlet_avg_bot_DO_mgL.csv', index=False)

print(monthly_mean_df)
monthly_mean_df.to_csv('monthly_mean_bot_DO_mgL.csv', index=False)
