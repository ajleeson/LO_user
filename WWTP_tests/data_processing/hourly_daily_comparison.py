"""
2025.04.22

Artificially create daily averages from hourly data output
of cas7_t0_x4b.
Then save for use in downstream TEF code (i.e. process_sections and bulk_calc)
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
enddate = '2017.12.31'
enddate_hrly = '2018.01.01 00:00:00'
year = '2017' # for making a date label

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
# remove the second lynchcove extraction
del sta_dict['lynchcove2']

# get list of inlet names
inlets = list(sta_dict.keys())

# where to put output figures
out_dir = Ldir['LOo'] / 'loading_test' / 'artificial_daily_averages' / ('extractions_'+startdate+'_'+enddate)
Lfun.make_dir(out_dir)

# create time_vecotr
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]


##########################################################
##        Calculate shallow and deep exchange flow      ##
##########################################################

for i,station in enumerate(inlets):
    # print status
    print('({}/{}) Working on {}...'.format(i+1,len(sta_dict),station))

    # initialize empty dataframe for saving
    df = pd.DataFrame()

    # get section information
    ds = xr.open_dataset('../../../LO_output/extract/'+gtagex+
                         '/tef2/extractions_'+startdate+'_'+enddate+
                         '/'+station+'.nc')
    
    # calculate daily averages from hourly data
    # first, remove extra hour at the end of the year (first hour of 2018)
    ds_dropped_last_value = ds.isel(time=slice(0, 8760))
    # Then, resample the data to daily averages
    ds_daily_averages = ds_dropped_last_value.resample(time="1D").mean()

    # save new averaged data
    out_fn = out_dir / (station+'.nc')
    ds_daily_averages.to_netcdf(out_fn)
