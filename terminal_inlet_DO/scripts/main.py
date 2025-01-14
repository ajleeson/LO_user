"""
Main script to process data and generate figures for
studying drivers of low oxygen in Puget Sound terminal inlets.

Aurora Leeson
January 2025
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
import matplotlib.patches as mpatches
from scipy.linalg import lstsq
import math
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
import csv
import cmocean
from scipy.stats import pearsonr
from scipy.stats import ttest_ind
import pingouin as pg
import matplotlib.pylab as plt
import gsw
import pickle
# import get_two_layer

# import helper functions
import get_monthly_means
import budget_error
import budget_barchart
import QinDOin_correl_consumption
import plot_monthly_means
import dodeep_hypvol_timeseries
import net_decrease_boxplots

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

##########################################################
##                    Read in data                      ##
##########################################################

print('Reading data...')

# NOTE: data in these deeplay_dict and shallowlay_dict
# are tidally-averaged daily time series
# in units of kmol O2 per second
# Values have been passed through a 71-hour lowpass Godin filter
# (Thomson & Emery, 2014)

# terminal inlet deep layer values
with open('../data/deeplay_dict.pickle', 'rb') as handle:
    deeplay_dict = pickle.load(handle)

# terminal inlet shallow layer values
with open('../data/shallowlay_dict.pickle', 'rb') as handle:
    shallowlay_dict = pickle.load(handle)

# terminal inlet dimensions
with open('../data/dimensions_dict.pickle', 'rb') as handle:
    dimensions_dict = pickle.load(handle)

# terminal inlet DO concentrations [mg/L]
with open('../data/DOconcen_dict.pickle', 'rb') as handle:
    DOconcen_dict = pickle.load(handle)

# get inlet names
inlets = list(deeplay_dict.keys())

# list of hypoxic inlets
hyp_inlets = ['penn','case','holmes','portsusan','lynchcove','dabob']

##########################################################
##                 Key values                           ##
##########################################################

# convert from kmol O2 per m3 per second to mg/L per day
kmolm3sec_to_mgLday = 1000 * 32 * 60 * 60 * 24

# yearday of drawdown period (mid-July through mid-August)
minday = 194
maxday = 225

# minday = 164
# maxday = 225

##########################################################
##   Get dates for analysis (2017.01.02 to 2017.12.30)  ##
##########################################################

year = '2017'

# set up dates
startdate = year + '.01.01'
enddate = year + '.12.31'
enddate_hrly = str(int(year)+1)+'.01.01 00:00:00'

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local_hrly = [pfun.get_dt_local(x) for x in dates_hrly]
# crop time vector (because we only have jan 2 - dec 30)
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')[2::]
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]

##########################################################
##                 Get monthly means                    ## 
##########################################################

# MONTHLYmean_XXXX are arrays of monthly mean values
# for all inlets, compressed into a single array

# df_MONTHLY_mean_XXX are dataframes, where each column
# is an individual inlet. All columns contain monthly
# mean values corresponding to the inlet (ie., 12 rows)

[MONTHLYmean_DOdeep,
MONTHLYmean_DOin,
MONTHLYmean_Tflush,
MONTHLYmean_perchyp,
df_MONTHLYmean_DOdeep,
df_MONTHLYmean_DOin,
df_MONTHLYmean_Tflush,
df_MONTHLYmean_perchyp] = get_monthly_means.get_monthly_means(deeplay_dict,DOconcen_dict,
                                                                dimensions_dict,inlets)

##########################################################
##               Deep Budget Error Analysis             ##
##########################################################

# calculate and print error of budget
# expressed as a % of QinDOin and biological consumption
budget_error.budget_error(inlets,shallowlay_dict,deeplay_dict,
                          dimensions_dict,kmolm3sec_to_mgLday)

##########################################################
##                  Budget Bar Charts                   ##
##########################################################

budget_barchart.budget_barchart(inlets,shallowlay_dict,deeplay_dict,
                    dates_local_hrly,dates_local_daily,hyp_inlets,
                    minday,maxday,kmolm3sec_to_mgLday)

#########################################################
## Correlation of Cons & QinDOin (mid-jul to mid-aug) ##
#########################################################

QinDOin_correl_consumption.correl(inlets,deeplay_dict,minday,maxday,kmolm3sec_to_mgLday)

#########################################################
##Plot monthly mean DOdeep, DOin, Tflush, and % hyp vol##
#########################################################
        
plot_monthly_means.plot_monthly_means(MONTHLYmean_DOdeep,
                                        MONTHLYmean_DOin,
                                        MONTHLYmean_Tflush,
                                        MONTHLYmean_perchyp,
                                        df_MONTHLYmean_DOdeep,
                                        df_MONTHLYmean_DOin,
                                        df_MONTHLYmean_Tflush)

##########################################################
##   Mean DOdeep vs % hyp vol and  DOdeep time series   ## 
##########################################################

dodeep_hypvol_timeseries.dodeep_hypvol_timeseries(MONTHLYmean_DOdeep,
                                                MONTHLYmean_perchyp,
                                                DOconcen_dict,
                                                dates_local_daily,
                                                dates_local_hrly,
                                                inlets,minday,maxday)

##########################################################
##    Net decrease (mid-July to mid-August) boxplots    ## 
##########################################################

net_decrease_boxplots.net_decrease_boxplots(dimensions_dict,deeplay_dict,
                                            minday,maxday)

# ##########################################################
# ##               Multiple regression (1)                ## 
# ##########################################################

# # create array of predictors

# input_array = np.array([mean_DOin, mean_Tflush, [1]*len(mean_DOin)]).T


# B,a,b,c = lstsq(input_array,deep_lay_DO)
# slope_DOin = B[0]
# slope_Tflush = B[1]
# intercept = B[2]

# print('\nMean deep layer DO [mg/L] = {}*DOin + {}*Tflush + {}'.format(
#     round(slope_DOin,2),round(slope_Tflush,2),round(intercept,2)))


# # calculate r^2 and p value
# r,p = pearsonr(mean_DOin,deep_lay_DO)
# print('\n===============================================')
# print('DO_deep dependence on DO_in')
# print('   r = {}'.format(r))
# print('   R^2 = {}'.format(r**2))
# print('   p = {}'.format(p))
# print('===============================================\n')

# # calculate r^2 and p value
# predicted_DOdeep = slope_DOin * mean_DOin + slope_Tflush * mean_Tflush + intercept
# r,p = pearsonr(deep_lay_DO,predicted_DOdeep)
# print('\n===============================================')
# print('DO_deep dependence on DO_in and T_flush')
# print('   r = {}'.format(r))
# print('   R^2 = {}'.format(r**2))
# print('   p = {}'.format(p))
# print('===============================================\n')
