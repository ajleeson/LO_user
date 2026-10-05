"""

Analyze change in DO for depths deeper than 10 m

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


# Hanning window length
nwin = 20

# years =  ['2014','2015']
years =  ['2015','2016','2017','2018','2019','2020']
# years =  ['2017']#['2015','2016','2017','2018','2019','2020']

plt.close('all')

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# initialize lists to store data for each basin
# loading
ps_avg_aug_DO_loading = []
ss_avg_aug_DO_loading = []
wb_avg_aug_DO_loading = []
mb_avg_aug_DO_loading = []
hc_avg_aug_DO_loading = []
# no-loading
ps_avg_aug_DO_noloading = []
ss_avg_aug_DO_noloading = []
wb_avg_aug_DO_noloading = []
mb_avg_aug_DO_noloading = []
hc_avg_aug_DO_noloading = []


# get aveage DO for all years (for the month of August)
for year in years:

    # read in data
    loading_avgDO   = xr.open_dataset('../../../LO_output/chapter_2/data/basins_shallow10m_deep_DO_cas7_t1_x11ab'+year+'.nc')
    noloading_avgDO = xr.open_dataset('../../../LO_output/chapter_2/data/basins_shallow10m_deep_DO_cas7_t1noDIN_x11ab'+year+'.nc')

    # crop to just august time period
    loading_avgDO_aug = loading_avgDO.sel(ocean_time=slice(str(year)+'-08-01', str(year)+'-09-01'))
    noloading_avgDO_aug = noloading_avgDO.sel(ocean_time=slice(str(year)+'-08-01', str(year)+'-09-01'))

    # loading_avgDO_aug = loading_avgDO
    # noloading_avgDO_aug = noloading_avgDO

    # concatenate to existing lists for each basin
    # loading run
    ps_avg_aug_DO_loading = np.concatenate((ps_avg_aug_DO_loading, loading_avgDO_aug['deep_DO_mgL_ps'].values))
    ss_avg_aug_DO_loading = np.concatenate((ss_avg_aug_DO_loading, loading_avgDO_aug['deep_DO_mgL_ss'].values))
    wb_avg_aug_DO_loading = np.concatenate((wb_avg_aug_DO_loading, loading_avgDO_aug['deep_DO_mgL_wb'].values))
    mb_avg_aug_DO_loading = np.concatenate((mb_avg_aug_DO_loading, loading_avgDO_aug['deep_DO_mgL_mb'].values))
    hc_avg_aug_DO_loading = np.concatenate((hc_avg_aug_DO_loading, loading_avgDO_aug['deep_DO_mgL_hc'].values))
    # no-loading run
    ps_avg_aug_DO_noloading = np.concatenate((ps_avg_aug_DO_noloading, noloading_avgDO_aug['deep_DO_mgL_ps'].values))
    ss_avg_aug_DO_noloading = np.concatenate((ss_avg_aug_DO_noloading, noloading_avgDO_aug['deep_DO_mgL_ss'].values))
    wb_avg_aug_DO_noloading = np.concatenate((wb_avg_aug_DO_noloading, noloading_avgDO_aug['deep_DO_mgL_wb'].values))
    mb_avg_aug_DO_noloading = np.concatenate((mb_avg_aug_DO_noloading, noloading_avgDO_aug['deep_DO_mgL_mb'].values))
    hc_avg_aug_DO_noloading = np.concatenate((hc_avg_aug_DO_noloading, noloading_avgDO_aug['deep_DO_mgL_hc'].values))


# get average concentration per basin
# loading
ps_mean_DO_sub10m_august_loading = np.nanmean(ps_avg_aug_DO_loading)
ss_mean_DO_sub10m_august_loading = np.nanmean(ss_avg_aug_DO_loading)
wb_mean_DO_sub10m_august_loading = np.nanmean(wb_avg_aug_DO_loading)
mb_mean_DO_sub10m_august_loading = np.nanmean(mb_avg_aug_DO_loading)
hc_mean_DO_sub10m_august_loading = np.nanmean(hc_avg_aug_DO_loading)
# no-loading
ps_mean_DO_sub10m_august_noloading = np.nanmean(ps_avg_aug_DO_noloading)
ss_mean_DO_sub10m_august_noloading = np.nanmean(ss_avg_aug_DO_noloading)
wb_mean_DO_sub10m_august_noloading = np.nanmean(wb_avg_aug_DO_noloading)
mb_mean_DO_sub10m_august_noloading = np.nanmean(mb_avg_aug_DO_noloading)
hc_mean_DO_sub10m_august_noloading = np.nanmean(hc_avg_aug_DO_noloading)

# get change in DO and percent change for each basin
# all puget sound
ps_diff = ps_mean_DO_sub10m_august_loading - ps_mean_DO_sub10m_august_noloading
ps_percent_change = (ps_diff / ps_mean_DO_sub10m_august_noloading) * 100
print('Puget Sound: change in DO = ' + str(round(ps_diff,3)) + ' mg/L, percent change = ' + str(round(ps_percent_change,2)) + '%')
# south sound
ss_diff = ss_mean_DO_sub10m_august_loading - ss_mean_DO_sub10m_august_noloading
ss_percent_change = (ss_diff / ss_mean_DO_sub10m_august_noloading) * 100
print('South Sound: change in DO = ' + str(round(ss_diff,3)) + ' mg/L, percent change = ' + str(round(ss_percent_change,2)) + '%')
# whidbey basin
wb_diff = wb_mean_DO_sub10m_august_loading - wb_mean_DO_sub10m_august_noloading
wb_percent_change = (wb_diff / wb_mean_DO_sub10m_august_noloading) * 100
print('Whidbey Basin: change in DO = ' + str(round(wb_diff,3)) + ' mg/L, percent change = ' + str(round(wb_percent_change,2)) + '%')
# main basin
mb_diff = mb_mean_DO_sub10m_august_loading - mb_mean_DO_sub10m_august_noloading
mb_percent_change = (mb_diff / mb_mean_DO_sub10m_august_noloading) * 100
print('Main Basin: change in DO = ' + str(round(mb_diff,3)) + ' mg/L, percent change = ' + str(round(mb_percent_change,2)) + '%')
# hood canal
hc_diff = hc_mean_DO_sub10m_august_loading - hc_mean_DO_sub10m_august_noloading
hc_percent_change = (hc_diff / hc_mean_DO_sub10m_august_noloading) * 100
print('Hood Canal: change in DO = ' + str(round(hc_diff,3)) + ' mg/L, percent change = ' + str(round(hc_percent_change,2)) + '%')
