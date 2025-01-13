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
import budget_error
import budget_barchart
import QinDOin_correl_consumption

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

##########################################################
##                    Read in data                      ##
##########################################################

print('Reading data...')

# NOTE: data in these dicts are tidally-averaged daily time series
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

print(list(deeplay_dict['lynchcove'].columns))

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

# ##########################################################
# ##                   DO scatterplots                    ## 
# ##########################################################

# DO_analysis = True
# if DO_analysis == True:

#     # initialize arrays for plotting
#     intervals = 12
#     deep_lay_DO = np.zeros(len(sta_dict)*intervals)
#     bott_sig_DO = np.zeros(len(sta_dict)*intervals)
#     min_bott_sig_DO = np.zeros(len(sta_dict)*intervals)
#     annual_mean_DO = np.zeros(len(sta_dict)*intervals)
#     perc_hyp_vol = np.zeros(len(sta_dict)*intervals)
#     mean_DOin = np.zeros(len(sta_dict)*intervals)
#     mean_Tflush = np.zeros(len(sta_dict)*intervals)
#     mean_TEFin = np.zeros(len(sta_dict)*intervals)
#     mean_recirc = np.zeros(len(sta_dict)*intervals)
#     mean_cons = np.zeros(len(sta_dict)*intervals)
#     mean_depth = np.zeros(len(sta_dict)*intervals)
#     inlet_vol = np.zeros(len(sta_dict)*intervals)
#     aspect_ratio = np.zeros(len(sta_dict)*intervals)

#     annmean_DOin = np.zeros(len(sta_dict))
#     annmin_DOin = np.zeros(len(sta_dict))
#     colors_twentyone = []

#     colors = []

#     # get values for plotting and calculating r value
#     for i,inlet in enumerate(sta_dict):
#         # get interface depth from csv file
#         with open('interface_depths.csv', 'r') as f:
#             for line in f:
#                 inlet, interface_depth = line.strip().split(',')
#                 interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
#         z_interface = float(interface_dict[inlet])
#         fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + inlet + '_2014.01.01_2014.12.31.nc'
#         ds = xr.open_dataset(fn)
#         moor_depth = ds.h.values
        
#         for month in range(intervals):
#             if month == 0:
#                 minday = 0 #1
#                 maxday = 30 #32
#             elif month == 1:
#                 minday = 30# 32
#                 maxday = 58# 60
#             elif month == 2:
#                 minday = 58# 60
#                 maxday = 89# 91
#             elif month == 3:
#                 minday = 89# 91
#                 maxday = 119# 121
#             elif month == 4:
#                 minday = 119# 121
#                 maxday = 150# 152
#             elif month == 5:
#                 minday = 150# 152
#                 maxday = 180# 182
#             elif month == 6:
#                 minday = 180# 182
#                 maxday = 211# 213
#             elif month == 7:
#                 minday = 211# 213
#                 maxday = 242# 244
#             elif month == 8:
#                 minday = 242# 244
#                 maxday = 272# 274
#             elif month == 9:
#                 minday = 272# 274
#                 maxday = 303# 305
#             elif month == 10:
#                 minday = 303# 305
#                 maxday = 332# 335
#             elif month == 11:
#                 minday = 332# 335
#                 maxday = 363
#             # minday=month
#             # maxday=month+1
#             # save values
#             deep_lay_DO[i*intervals+month] =  np.nanmean(DOconcen_dict[inlet]['Deep Layer'][minday:maxday])
#             bott_sig_DO[i*intervals+month] =  np.nanmean(DOconcen_dict[inlet]['Bottom Sigma DO'][minday:maxday])
#             min_bott_sig_DO[i*intervals+month] =  np.nanmin(DOconcen_dict[inlet]['Minimum Bottom Layer DO'][minday:maxday])
#             annual_mean_DO[i*intervals+month] = np.nanmean(DOconcen_dict[inlet]['Deep Layer'])
#             perc_hyp_vol[i*intervals+month] = np.nanmean(DOconcen_dict[inlet]['percent hypoxic volume'][minday:maxday])
#             mean_DOin[i*intervals+month] = np.nanmean(DOconcen_dict[inlet]['Qin DO'][minday:maxday])
#             mean_Tflush[i*intervals+month] = np.nanmean(dimensions_dict[inlet]['Inlet volume'][0]/deeplay_dict[inlet]['Qin m3/s'][minday:maxday]) / (60*60*24)
#             mean_TEFin[i*intervals+month] = np.nanmean(deeplay_dict[inlet]['TEF Exchange Flow'][minday:maxday]/deeplay_dict[inlet]['Volume'][minday:maxday]) * (
#                                     32 * 1000) * (60*60*24)
#             mean_recirc[i*intervals+month] = np.nanmean(deeplay_dict[inlet]['TEF Recirculation'][minday:maxday]/deeplay_dict[inlet]['Volume'][minday:maxday]) * (
#                                     32 * 1000) * (60*60*24)
#             mean_cons[i*intervals+month] = np.nanmean(deeplay_dict[inlet]['Bio Consumption'][minday:maxday]/deeplay_dict[inlet]['Volume'][minday:maxday]) * (
#                                     32 * 1000) * (60*60*24)
#             mean_depth[i*intervals+month] = dimensions_dict[inlet]['Mean depth'][0]
#             inlet_vol[i*intervals+month] = dimensions_dict[inlet]['Inlet volume'][0]
#             aspect_ratio[i*intervals+month] = dimensions_dict[inlet]['L/W aspect ratio'][0]
#             colors.append(basin_color_dict[basin_dict[inlet]])

#         annmean_DOin[i] = np.nanmean(DOconcen_dict[inlet]['Qin DO'])
#         annmin_DOin[i] = np.nanmin(DOconcen_dict[inlet]['Qin DO'])
#         colors_twentyone.append(basin_color_dict[basin_dict[inlet]])
        


#     # DOin vs. Tflush colored by percent hypoxic volume ============== 4PART MONEY PLOT SCATTER ================================
#     percent_hypoxic = True
#     if percent_hypoxic == True:
#         # initialize figure
#         fig, ax = plt.subplots(2,2,figsize = (10,9))
#         ax = ax.ravel()

#         # format figure
#         # ax[0].set_title('(a) ' + year + r' monthly mean DO$_{deep}$ vs. DO$_{in}$' + '\n' + r'colored by mean T$_{flush}$',
#         #                 size=14, loc='left')
#         ax[0].set_title('(a) All Inlets', size=14, loc='left', fontweight='bold')
#         # ax[0].set_title(year + r' monthly mean DO$_{deep}$ vs. DO$_{in}$',
#         #                 size=16, loc='left')
#         # format grid
#         # ax[0].set_facecolor('#EEEEEE')
#         ax[0].tick_params(axis='x', labelrotation=30)
#         # ax[0].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#         # for border in ['top','right','bottom','left']:
#         #     ax[0].spines[border].set_visible(False)
#         ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#         ax[0].tick_params(axis='both', labelsize=12)
#         ax[0].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
#         ax[0].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
#         # plot
#         cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
#         cmap_oxy = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
#         ax[0].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, c='k', alpha=0.5)
#         ax[0].plot([0,12],[0,12],color='gray')
#         # ax[0].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4, fontsize=10)
#         ax[0].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
#         cs = ax[0].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
#         # create colorbarlegend
#         cbar = fig.colorbar(cs)
#         cbar.ax.tick_params(labelsize=12)
#         cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
#         cbar.outline.set_visible(False)
#         ax[0].set_xlim([0,10])
#         ax[0].set_ylim([0,10])

#         # format figure
#         # ax[1].set_title('(b) ' + year + ' monthly mean '+r'DO$_{in}$ vs. T$_{flush}$'+'\ncolored by % hypoxic volume',
#         #                 size=14, loc='left')
#         ax[1].set_title('(b) All Inlets', size=14, loc='left', fontweight='bold')
#         # format grid
#         # ax[1].set_facecolor('#EEEEEE')
#         ax[1].tick_params(axis='x', labelrotation=30)
#         # ax[1].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#         # for border in ['top','right','bottom','left']:
#         #     ax[1].spines[border].set_visible(False)
#         ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#         ax[1].tick_params(axis='both', labelsize=12)
#         ax[1].set_xlabel(r'Monthly mean T$_{flush}$ [days]', fontsize=12)
#         ax[1].set_ylabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
#         ax[1].set_ylim([0,10])
#         ax[1].set_xlim([0,85])
#         # plot
#         cmap_hyp = plt.cm.get_cmap('gist_heat_r')
#         cs_DO = ax[1].scatter(mean_Tflush,mean_DOin,s=60,zorder=5,edgecolor='gray',c=perc_hyp_vol,cmap=cmap_hyp)
#         # create colorbarlegend
#         cbar = fig.colorbar(cs_DO)
#         cbar.ax.tick_params(labelsize=12)
#         cbar.ax.set_ylabel('Monthly mean % hypoxic volume', rotation=90, fontsize=12)
#         cbar.outline.set_visible(False)

#         # crescent bay
#         # ax[2].set_title('(c) Crescent Bay 2017 monthly mean \n' + r'DO$_{deep}$ vs. DO$_{in}$ colored by T$_{flush}$', loc='left', size=14)
#         ax[2].set_title('(c) Crescent Bay', size=14, loc='left', fontweight='bold')
#         # format grid
#         # ax[2].set_facecolor('#EEEEEE')
#         ax[2].tick_params(axis='x', labelrotation=30)
#         # ax[2].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#         # for border in ['top','right','bottom','left']:
#         #     ax[2].spines[border].set_visible(False)
#         ax[2].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#         ax[2].tick_params(axis='both', labelsize=12)
#         ax[2].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
#         ax[2].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
#         # plot
#         cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
#         cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
#         ax[2].plot([0,11],[0,11],color='dimgray')
#         ax[2].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
#         # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
#         ax[2].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, color='gray',alpha=0.5, edgecolor='none')
#         for i,inlet in enumerate(sta_dict):
#             if inlet == 'crescent':
#                 cs = ax[2].scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
#                                 s=150, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_tflush,
#                             linewidth=2, vmin=0, vmax=40)
#             else:
#                 continue
#         # create colorbarlegend
#         cbar = fig.colorbar(cs)
#         cbar.ax.tick_params(labelsize=12)
#         cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
#         cbar.outline.set_visible(False)
#         ax[2].set_xlim([0,11])
#         ax[2].set_ylim([0,11])

#         # lynch cove
#         # ax[3].set_title('(d) Lynch Cove 2017 monthly mean \n' + r'DO$_{deep}$ vs. DO$_{in}$ colored by T$_{flush}$', loc='left', size=14)
#         ax[3].set_title('(d) Lynch Cove', size=14, loc='left', fontweight='bold')
#         # format grid
#         # ax[3].set_facecolor('#EEEEEE')
#         ax[3].tick_params(axis='x', labelrotation=30)
#         # ax[3].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
#         # for border in ['top','right','bottom','left']:
#         #     ax[3].spines[border].set_visible(False)
#         ax[3].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#         ax[3].tick_params(axis='both', labelsize=12)
#         ax[3].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
#         ax[3].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
#         # plot
#         cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
#         cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
#         ax[3].plot([0,11],[0,11],color='dimgray')
#         ax[3].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
#         # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
#         ax[3].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, color='gray',alpha=0.5, edgecolor='none')
#         for i,inlet in enumerate(sta_dict):
#             if inlet == 'lynchcove':
#                 cs = ax[3].scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
#                                 s=150, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_tflush,
#                             linewidth=2, vmin=0, vmax=40)
#             else:
#                 continue
#         # create colorbarlegend
#         cbar = fig.colorbar(cs)
#         cbar.ax.tick_params(labelsize=12)
#         cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
#         cbar.outline.set_visible(False)
#         ax[3].set_xlim([0,11])
#         ax[3].set_ylim([0,11])


#         plt.tight_layout()
#         # save figure
#         plt.show()

# ##########################################################
# ## MEAN DEEP DO vs % HYPOXIC VOLUME and  Deep DO time series ## 
# ##########################################################

# # mid Jul - mid Aug
# minday = 194
# maxday = 225

# # initialize figure
# fig, ax = plt.subplots(1,2,figsize = (12,5),gridspec_kw={'width_ratios': [1, 1.5]})


# # Deep DO vs. % hypoxic volume
# # ax[0].set_title('(a) ' + year + ' monthly mean \n     ' + r'% hypoxic volume vs. DO$_{deep}$',
# #                 size=14, loc='left')
# ax[0].set_title('(a) All inlets', size=14, loc='left', fontweight='bold')
# # format grid
# # ax[0].set_facecolor('#EEEEEE')
# ax[0].tick_params(axis='x', labelrotation=30)
# ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
# # ax[0].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
# # for border in ['top','right','bottom','left']:
# #     ax[0].spines[border].set_visible(False)
# ax[0].tick_params(axis='both', labelsize=12)
# ax[0].set_xlabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=14)
# ax[0].set_ylabel('Monthly mean % hypoxic volume', fontsize=14)
# # plot
# ax[0].scatter(deep_lay_DO,perc_hyp_vol,alpha=0.3,s=80,zorder=5,color='navy')
#         # color=colors)
# ax[0].set_xlim([0,10])
# ax[0].set_ylim([0,100])


# # Deep DO timeseries
# # ax[1].set_title('(b) ' + year + r' DO$_{deep}$ time series [mg/L]' +  '\n     (30-day Hanning Window)',
# #                 size=14, loc='left')
# ax[1].set_title('(b) All inlets', size=14, loc='left', fontweight='bold')
# # format grid
# # ax[1].set_facecolor('#EEEEEE')
# ax[1].tick_params(axis='x', labelrotation=30)
# loc = mdates.MonthLocator(interval=1)
# ax[1].xaxis.set_major_locator(loc)
# ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
# # ax[1].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
# # for border in ['top','right','bottom','left']:
# #     ax[1].spines[border].set_visible(False)
# ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
# ax[1].tick_params(axis='both', labelsize=12)
# # add drawdown period
# # ax[1].axvline(dates_local_daily[minday],0,12,color='pink')
# # ax[1].axvline(dates_local_daily[maxday],0,12,color='pink')
# ax[1].axvline(dates_local_daily[minday],0,12,color='grey')#,linestyle=':')
# ax[1].axvline(dates_local_daily[maxday],0,12,color='grey')#,linestyle=':')
# # loop through stations
# for i,inlet in enumerate(sta_dict):
#     # get average deep layer DO
#     deep_lay_DO_alltime = DOconcen_dict[inlet]['Deep Layer']
#     # 30-day hanning windwo
#     deep_lay_DO_alltime = zfun.lowpass(deep_lay_DO_alltime.values,n=30)
#     ax[1].plot(dates_local_daily,deep_lay_DO_alltime,linewidth=1,color='navy',alpha=0.5)

# # format labels
# ax[1].set_xlim([dates_local[0],dates_local[-2]])
# ax[1].set_ylim([0,10])
# ax[1].set_ylabel(r'DO$_{deep}$ [mg/L]',fontsize=14)
# plt.tight_layout()
# plt.show()

# ##########################################################
# ## net decrease (mid-July to mid-August) boxplots ## 
# ##########################################################

# # mid July to mid August
# minday = 194
# maxday = 225

# # initialize figure
# fig, ax = plt.subplots(1,1,figsize = (10,5))

# # format figure
# ax.set_title('d/dt(DO) (mid-Jul through mid-Aug)',size=14)
# ax.tick_params(axis='x', labelrotation=30)
# ax.grid(True,color='silver',linewidth=1,linestyle='--',axis='x')
# ax.tick_params(axis='both', labelsize=12)
# ax.set_ylabel('d/dt(DO) [mg/L per day]', fontsize=12)

# # # add line with slope -1
# # ax.plot([0,0.45], [0,-0.45], color='grey', linestyle='-')

# storage_all = []

# for i,inlet in enumerate(sta_dict):
    
#     # get daily net decrease rate
#     storage_daily =  deeplay_dict[inlet]['Storage'][minday:maxday]/(deeplay_dict[inlet]['Volume'][minday:maxday]) # kmol O2 /s /m3
    
#     # convert to mg/L per day
#     storage_daily = storage_daily.values * 1000 * 32 * 60 * 60 * 24

#     # add to array
#     storage_all.append(list(storage_daily))

# # create boxplot
# ax.axhline(y=0, xmin=-0.5, xmax=1.05,color='silver',linewidth=1,linestyle='--')
# bplot = plt.boxplot(storage_all, patch_artist=True, labels=sta_dict.keys(),
#             showmeans=True)

# for patch in bplot['boxes']:
#     patch.set_facecolor('honeydew')

# # condudct anova test
# # anova0 = f_oneway(storage_all[0],storage_all[1],storage_all[2],storage_all[3],storage_all[4],
# #          storage_all[5],storage_all[6],storage_all[7],storage_all[8],storage_all[9],
# #          storage_all[10],storage_all[11],storage_all[12])
# # anova1 = kruskal(storage_all[0],storage_all[1],storage_all[2],storage_all[3],storage_all[4],
# #          storage_all[5],storage_all[6],storage_all[7],storage_all[8],storage_all[9],
# #          storage_all[10],storage_all[11],storage_all[12])
# # anova2 = alexandergovern(storage_all[0],storage_all[1],storage_all[2],storage_all[3],storage_all[4],
# #          storage_all[5],storage_all[6],storage_all[7],storage_all[8],storage_all[9],
# #          storage_all[10],storage_all[11],storage_all[12])
# flat_storage = [x for xs in storage_all for x in xs]
# df = pd.DataFrame({'storage':flat_storage,
#                    'inlet': np.repeat(list(sta_dict.keys()), repeats=31)})
# # print(df)
# pingu = pg.welch_anova(data=df, dv='storage', between='inlet')


# print('\n===============================\n')
# # print(anova0)
# # print(anova1)
# # print(anova2)
# print(pingu)

# plt.tight_layout()
# plt.show()

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
