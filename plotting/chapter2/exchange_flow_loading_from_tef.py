"""
Plot QinDINin at Admiralty Inlet
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
import matplotlib.patches as patches
import csv
import cmocean
import matplotlib.pylab as plt
import gsw
import pickle
import tef_fun as tef_fun

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

residual = True # recirculation
show_EU = False

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t1_x11b'
jobname = 'penn_hoodcanal'
startdate = '2017.01.01'
enddate = '2017.12.31'
enddate_hrly = '2018.01.01 00:00:00'
year = '2017' # for making a date label

# # parse gtagex
# gridname, tag, ex_name = gtagex.split('_')
# Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# # find job lists from the extract moor
# job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# # Get mooring stations:
# sta_dict = job_lists.get_sta_dict(jobname)

# # where to put output figures
# out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
# Lfun.make_dir(out_dir)

# create time vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]

print('\n')


# --------------------------- get TEF exchange flow terms ----------------------------------------
section = 'ai' # Admiralty Inlet

in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31')
tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir,section)
# get inflowing values
Q_in = tef_df['q_p'] # Qin [m3/s]
NO3_in = tef_df['NO3_p'] # NO3in [mmol/m3]
NH4_in = tef_df['NH4_p'] # NH4in [mmol/m3]
DIN_in = NO3_in + NH4_in # DINin [mmol/m3]
# determine Qin*DINin
QinDINin_mmol_s = Q_in * DIN_in # [mmol/m3]
# convert to kg/d
QinDINin_kg_d = QinDINin_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)

plt.close('all')
plt.plot(dates_local_daily[1:-1],QinDINin_kg_d)
print(np.nanmean(QinDINin_kg_d))
