"""
Calculate exchange flow terms for DO budget
Have shallow and deep layer
"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
import pandas as pd
import math
import csv
import cmocean
import matplotlib.pylab as plt
import gsw
import pickle

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
jobname = 'twentyoneinlets'
startdate = '2017.01.01'
# enddate = '2017.01.02'
enddate = '2017.12.31'
year = '2017' # for making a date label

dsf = Ldir['ds_fmt']
# dt0 = datetime.strptime(startdate,dsf)
# dt1 = datetime.strptime(enddate,dsf)
# # days = (dt0, dt1)
    
# # pandas Index objects
# dt_ind = pd.date_range(start=dt0, end=dt1)
# yd_ind = pd.Index(dt_ind.dayofyear)

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_'+startdate+'_'+enddate) / '2layer_traps'
Lfun.make_dir(out_dir)

# create dictionaries with interface depths
interface_dict = dict()
print('\n')

##########################################################
##get list of all traps in the segment and in LiveOcean ##
##########################################################

# prep info for rivers
seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
seg_dict = pd.read_pickle(seg_name)
# get list of all rivers and all wwtps in LiveOcean
LOrivs = []
LOwwtps = []
# get list of all WWTPs
wwtp_info = Ldir['data'] / 'grids' / 'cas7' / 'wwtp_info.csv'
with open(wwtp_info, 'r') as f:
    for line in f:
        rname,row_py,col_py = line.strip().split(',')
        LOwwtps = np.concatenate((LOwwtps,[str(rname)]))
# get list of all rivers
triv_info = Ldir['data'] / 'grids' / 'cas7' / 'triv_info.csv'
with open(triv_info, 'r') as f:
    for line in f:
        rname,row_py,col_py,idir,isign,uv = line.strip().split(',')
        LOrivs = np.concatenate((LOrivs,[str(rname)]))
river_info = Ldir['data'] / 'grids' / 'cas7' / 'river_info.csv'
with open(river_info, 'r') as f:
    for line in f:
        rname,row_py,col_py,idir,isign,uv = line.strip().split(',')
        LOrivs = np.concatenate((LOrivs,[str(rname)]))

# ##########################################################
# ##     calculate all river and wwtp DO transports       ##
# ##########################################################

frc = 'trapsF00'

# list of variables to extract
vn_list = ['transport']

dt0 = datetime.strptime(startdate, Lfun.ds_fmt)
dt1 = datetime.strptime(enddate, Lfun.ds_fmt)
ndays = (dt1-dt0).days + 1

# make mds_list: list of datestrings (e.g. 2021.01.01) to loop over
mds_list = []
mdt = dt0
while mdt <= dt1:
    mds_list.append(datetime.strftime(mdt, Lfun.ds_fmt))
    mdt = mdt + timedelta(days=1)

# get list of river names
# (this is a bit titchy because of NetCDF 3 limitations on strings, forcing them
# to be arrays of characters)
mds = mds_list[0]
# fn = Ldir['LOo'] / 'forcing' / gridname / ('f' + mds) / frc / 'rivers.nc'
fn = '/dat1/parker/LO_output/forcing/' + gridname + '/f' + mds + '/' +  frc + '/rivers.nc'
ds = xr.open_dataset(fn)
rn = ds['river_name'].values
NR = rn.shape[0]
riv_name_list = list(rn)

NT = len(mds_list)

# get state variable values
nanmat = np.nan * np.ones((NT, NR))
v_dict = dict()
for vn in vn_list:
    v_dict[vn] = nanmat.copy()
for tt,mds in enumerate(mds_list):
    this_dt = datetime.strptime(mds, Lfun.ds_fmt)
    if this_dt.day == 1 and this_dt.month == 1:
        print(' Year = %d' % (this_dt.year))
    # fn = Ldir['LOo'] / 'forcing' / gridname / ('f' + mds) / frc / 'rivers.nc'
    fn = '/dat1/parker/LO_output/forcing/' + gridname + '/f' + mds + '/' +  frc + '/rivers.nc'
    ds = xr.open_dataset(fn)
    # The river transport is given at noon of a number of days surrounding the forcing date.
    # Here we find the index of the time for the day "mds".
    RT = pd.to_datetime(ds['river_time'].values)
    mdt = this_dt + timedelta(hours=12)
    mask = RT == mdt
    for vn in vn_list:
        if vn == 'transport':
            v_dict[vn][tt,:] = ds['river_' + vn][mask,:]
        else:
            # the rest of the variables allow for depth variation, but we
            # don't use this, so, just use the bottom value
            v_dict[vn][tt,:] = ds['river_' + vn][mask,0,:]
    ds.close()

# make transport positive
v_dict['transport'] = np.abs(v_dict['transport'])

# store output in an xarray Dataset
mdt_list = [(datetime.strptime(item, Lfun.ds_fmt) + timedelta(hours=12)) for item in mds_list]
times = pd.Index(mdt_list)

# create dataset
x = xr.Dataset(coords={'time': times,'riv': riv_name_list})

# add state variables
for vn in vn_list:
    v = v_dict[vn]
    x[vn] = (('time','riv'), v)

x.to_netcdf('2017rivflow.nc')

