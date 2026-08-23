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

section = 'ai' # Admiralty Inlet


# LOADING RUN
# --------------------------- get TEF exchange flow terms ---------------------------------------
print('Loading')
in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'cps' / ('bulk_'+year+'.01.01_'+year+'.12.31')
tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir,section)
# get inflowing values
Q_in = tef_df['q_p'] # Qin [m3/s]
NO3_in = tef_df['NO3_p'] # NO3in [mmol/m3]
NH4_in = tef_df['NH4_p'] # NH4in [mmol/m3]
phyt_in = tef_df['phytoplankton_p']  # mmol/m3]
zoop_in = tef_df['zooplankton_p']   # [mmol/m3]
Sdet_in = tef_df['SdetritusN_p']    # [mmol/m3]
Ldet_in = tef_df['LdetritusN_p']    # [mmol/m3]
# get outflowing values
Q_out = tef_df['q_m'] # Qin [m3/s]
NO3_out = tef_df['NO3_m'] # NO3in [mmol/m3]
NH4_out = tef_df['NH4_m'] # NH4in [mmol/m3]
phyt_out = tef_df['phytoplankton_m'] # [mmol/m3]
zoop_out = tef_df['zooplankton_m']   # [mmol/m3]
Sdet_out = tef_df['SdetritusN_m']    # [mmol/m3]
Ldet_out = tef_df['LdetritusN_m']    # [mmol/m3]
# combine terms
DIN_in = NO3_in + NH4_in # DINin [mmol/m3]
TN_in = DIN_in + phyt_in + zoop_in + Sdet_in + Ldet_in # TNin [mmol/m3]
DIN_out = NO3_out + NH4_out # DINin [mmol/m3]
TN_out = DIN_out + phyt_out + zoop_out + Sdet_out + Ldet_out # TNin [mmol/m3]
# determine Qin*DINin and Qin*TNin
QinDINin_mmol_s = Q_in * DIN_in # [mmol/m3]
QinTNin_mmol_s = Q_in * TN_in # [mmol/m3]
print('Loading')
print(np.nanmean(QinTNin_mmol_s))
print('------------------')
# determine Qout*DINout and Qout*TNout
QoutDINout_mmol_s = Q_out * DIN_out # [mmol/m3]
QoutTNout_mmol_s = Q_out * TN_out # [mmol/m3]
# convert to kg/d
QinDINin_kg_d = QinDINin_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QinTNin_kg_d = QinTNin_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QoutDINout_kg_d = QoutDINout_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QoutTNout_kg_d = QoutTNout_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)

plt.close('all')
plt.plot(dates_local_daily[1:-1],QoutDINout_kg_d+QinDINin_kg_d,linestyle='--',color='teal',label='Loading')
# plt.plot(dates_local_daily[1:-1],QinTNin_kg_d,linewidth=0.5,color='crimson',label='Loading')
# print(np.nanmean(QoutDINout_kg_d+QinDINin_kg_d))

# monthly_mean = QinDINin_kg_d.resample("M").mean()
# print(monthly_mean)

# NO-LOADING RUN
# --------------------------- get TEF exchange flow terms ----------------------------------------
print('No-Loading')
in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1noDIN_x11b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31')
tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir,section)
# get inflowing values
Q_in = tef_df['q_p'] # Qin [m3/s]
NO3_in = tef_df['NO3_p'] # NO3in [mmol/m3]
NH4_in = tef_df['NH4_p'] # NH4in [mmol/m3]
phyt_in = tef_df['phytoplankton_p']  # mmol/m3]
zoop_in = tef_df['zooplankton_p']   # [mmol/m3]
Sdet_in = tef_df['SdetritusN_p']    # [mmol/m3]
Ldet_in = tef_df['LdetritusN_p']    # [mmol/m3]
# get outflowing values
Q_out = tef_df['q_m'] # Qin [m3/s]
NO3_out = tef_df['NO3_m'] # NO3in [mmol/m3]
NH4_out = tef_df['NH4_m'] # NH4in [mmol/m3]
phyt_out = tef_df['phytoplankton_m'] # [mmol/m3]
zoop_out = tef_df['zooplankton_m']   # [mmol/m3]
Sdet_out = tef_df['SdetritusN_m']    # [mmol/m3]
Ldet_out = tef_df['LdetritusN_m']    # [mmol/m3]
# combine terms
DIN_in = NO3_in + NH4_in # DINin [mmol/m3]
TN_in = DIN_in + phyt_in + zoop_in + Sdet_in + Ldet_in # TNin [mmol/m3]
DIN_out = NO3_out + NH4_out # DINin [mmol/m3]
TN_out = DIN_out + phyt_out + zoop_out + Sdet_out + Ldet_out # TNin [mmol/m3]
# determine Qin*DINin and Qin*TNin
QinDINin_mmol_s = Q_in * DIN_in # [mmol/m3]
QinTNin_mmol_s = Q_in * TN_in # [mmol/m3]
print('No-loading')
print(np.nanmean(QinTNin_mmol_s))
print('------------------')
# determine Qout*DINout and Qout*TNout
QoutDINout_mmol_s = Q_out * DIN_out # [mmol/m3]
QoutTNout_mmol_s = Q_out * TN_out # [mmol/m3]
# convert to kg/d
QinDINin_kg_d = QinDINin_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QinTNin_kg_d = QinTNin_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QoutDINout_kg_d = QoutDINout_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)
QoutTNout_kg_d = QoutTNout_mmol_s / 71.4 * 86.4 # [kg/day] (71.4 gets from mmol/m3 to mg/L, and 86.4 gets to kg/d)

plt.plot(dates_local_daily[1:-1],QoutDINout_kg_d+QinDINin_kg_d,color='cadetblue',linewidth=2,alpha=0.5,label='No-Loading')
# plt.plot(dates_local_daily[1:-1],QinTNin_kg_d,color='red',linewidth=2,alpha=0.5,label='No-Loading')
# print(np.nanmean(QoutDINout_kg_d+QinDINin_kg_d))

# monthly_mean = QinDINin_kg_d.resample("M").mean()
# print(monthly_mean)
