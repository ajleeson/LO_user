"""
Violin plots of average shallow (<= 10 m) and deep layer DO concentrations in Puget Sound
(the Puget Sound mask has already been applied in the extractions, so I do not need
to apply them in this script)

"""

# import things
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
from statsmodels.graphics.boxplots import violinplot

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
print(pth)
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

plt.close('all')

##############################################################
##                       USER INPUTS                        ##
##############################################################

years = ['2017']#['2015','2016','2017','2018','2019','2020']

# which  model run to look at?
gtagexes = ['cas7_t1noDIN_x11ab','cas7_t1_x11ab']

# # where to put output figures
# out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'figures'
# Lfun.make_dir(out_dir)

##############################################################
##                      PROCESS DATA                        ##
##############################################################

# open dataset for every year, and add to dictionary, with year as key

# open datasets
# initialize empty dictionary
ds_loading_dict = {}
ds_noloading_dict = {}
for year in years:
    # open datasets
    ds_loading = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / ('cas7_t1_x11ab_pugetsoundDO_' + year + '_shallow10m_deep_DO.nc'))
    ds_noloading = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / ('cas7_t1noDIN_x11ab_pugetsoundDO_' + year + '_shallow10m_deep_DO.nc'))
    # drop leap day if 2016 or 2020
    if year in ['2016','2020']:
        ds_loading = ds_loading.sel(ocean_time=~((ds_loading.ocean_time.dt.month == 2) & (ds_loading.ocean_time.dt.day == 29)))
        ds_noloading = ds_noloading.sel(ocean_time=~((ds_noloading.ocean_time.dt.month == 2) & (ds_noloading.ocean_time.dt.day == 29)))
    # add ds to dictionary
    ds_loading_dict[year] = ds_loading
    ds_noloading_dict[year] = ds_noloading

# get climatology of bottom DO for all seven years
shallowDO_loading_arrays = [ds['shallow_DO_mgL'].data for ds in ds_loading_dict.values()]  # [365]
deepDO_loading_arrays = [ds['deep_DO_mgL'].data for ds in ds_loading_dict.values()]  # [365]

shallowDO_noloading_arrays = [ds['shallow_DO_mgL'].data for ds in ds_noloading_dict.values()]  # [365]
deepDO_noloading_arrays = [ds['deep_DO_mgL'].data for ds in ds_noloading_dict.values()]  # [365]

# stack the arrays
shallowDO_loading_stacked = np.stack(shallowDO_loading_arrays, axis=0)  # [7, 365]
deepDO_loading_stacked = np.stack(deepDO_loading_arrays, axis=0)  # [7, 365]

shallowDO_noloading_stacked = np.stack(shallowDO_noloading_arrays, axis=0)  # [7, 365]
deepDO_noloading_stacked = np.stack(deepDO_noloading_arrays, axis=0)  # [7, 365]

# average over all seven years
shallowDO_clim_loading = shallowDO_loading_stacked.mean(axis=0)  # [365]
deepDO_clim_loading = deepDO_loading_stacked.mean(axis=0)  # [365]

shallowDO_clim_noloading = shallowDO_noloading_stacked.mean(axis=0)  # [365]
deepDO_clim_noloading = deepDO_noloading_stacked.mean(axis=0)  # [365]

# crop by season
jan = 0
feb = 31
mar = 59
apr = 90
may = 120
jun = 151
jul = 181
aug = 212
sep = 243
oct = 273
nov = 304
dec = 334

# print('January')
# print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[jan:feb])*1000/32,2)))
# print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[jan:feb])*1000/32,2)))
# print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[jan:feb])*1000/32,2)))
# print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[jan:feb])*1000/32,2)))

# print('February')
# print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[feb:mar])*1000/32,2)))
# print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[feb:mar])*1000/32,2)))
# print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[feb:mar])*1000/32,2)))
# print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[feb:mar])*1000/32,2)))

# print('March')
# print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[mar:apr])*1000/32,2)))
# print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[mar:apr])*1000/32,2)))
# print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[mar:apr])*1000/32,2)))
# print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[mar:apr])*1000/32,2)))

# print('April')
# print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[apr:may])*1000/32,2)))
# print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[apr:may])*1000/32,2)))
# print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[apr:may])*1000/32,2)))
# print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[apr:may])*1000/32,2)))

# print('May')
# print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[may:jun])*1000/32,2)))
# print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[may:jun])*1000/32,2)))
# print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[may:jun])*1000/32,2)))
# print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[may:jun])*1000/32,2)))

# print('June')
# print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[jun:jul])*1000/32,2)))
# print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[jun:jul])*1000/32,2)))
# print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[jun:jul])*1000/32,2)))
# print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[jun:jul])*1000/32,2)))

# print('July')
# print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[jul:aug])*1000/32,2)))
# print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[jul:aug])*1000/32,2)))
# print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[jul:aug])*1000/32,2)))
# print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[jul:aug])*1000/32,2)))


jan_loading_shallow = shallowDO_clim_loading[jan:feb]
jan_noloading_shallow = shallowDO_clim_noloading[jan:feb]
jan_loading_deep = deepDO_clim_loading[jan:feb]
jan_noloading_deep = deepDO_clim_noloading[jan:feb]

feb_loading_shallow = shallowDO_clim_loading[feb:mar]
feb_noloading_shallow = shallowDO_clim_noloading[feb:mar]
feb_loading_deep = deepDO_clim_loading[feb:mar]
feb_noloading_deep = deepDO_clim_noloading[feb:mar]

mar_loading_shallow = shallowDO_clim_loading[mar:apr]
mar_noloading_shallow = shallowDO_clim_noloading[mar:apr]
mar_loading_deep = deepDO_clim_loading[mar:apr]
mar_noloading_deep = deepDO_clim_noloading[mar:apr]

apr_loading_shallow = shallowDO_clim_loading[apr:may]
apr_noloading_shallow = shallowDO_clim_noloading[apr:may]
apr_loading_deep = deepDO_clim_loading[apr:may]
apr_noloading_deep = deepDO_clim_noloading[apr:may]

may_loading_shallow = shallowDO_clim_loading[may:jun]
may_noloading_shallow = shallowDO_clim_noloading[may:jun]
may_loading_deep = deepDO_clim_loading[may:jun]
may_noloading_deep = deepDO_clim_noloading[may:jun]

jun_loading_shallow = shallowDO_clim_loading[jun:jul]
jun_noloading_shallow = shallowDO_clim_noloading[jun:jul]
jun_loading_deep = deepDO_clim_loading[jun:jul]
jun_noloading_deep = deepDO_clim_noloading[jun:jul]

jul_loading_shallow = shallowDO_clim_loading[jul:aug]
jul_noloading_shallow = shallowDO_clim_noloading[jul:aug]
jul_loading_deep = deepDO_clim_loading[jul:aug]
jul_noloading_deep = deepDO_clim_noloading[jul:aug]

aug_loading_shallow = shallowDO_clim_loading[aug:sep]
aug_noloading_shallow = shallowDO_clim_noloading[aug:sep]
aug_loading_deep = deepDO_clim_loading[aug:sep]
aug_noloading_deep = deepDO_clim_noloading[aug:sep]

sep_loading_shallow = shallowDO_clim_loading[sep:oct]
sep_noloading_shallow = shallowDO_clim_noloading[sep:oct]
sep_loading_deep = deepDO_clim_loading[sep:oct]
sep_noloading_deep = deepDO_clim_noloading[sep:oct]

oct_loading_shallow = shallowDO_clim_loading[oct:nov]
oct_noloading_shallow = shallowDO_clim_noloading[oct:nov]
oct_loading_deep = deepDO_clim_loading[oct:nov]
oct_noloading_deep = deepDO_clim_noloading[oct:nov]

nov_loading_shallow = shallowDO_clim_loading[nov:dec]
nov_noloading_shallow = shallowDO_clim_noloading[nov:dec]
nov_loading_deep = deepDO_clim_loading[nov:dec]
nov_noloading_deep = deepDO_clim_noloading[nov:dec]

dec_loading_shallow = shallowDO_clim_loading[dec::]
dec_noloading_shallow = shallowDO_clim_noloading[dec::]
dec_loading_deep = deepDO_clim_loading[dec::]
dec_noloading_deep = deepDO_clim_noloading[dec::]

########################


winter_loading_shallow = shallowDO_clim_loading[jan:apr]
winter_noloading_shallow = shallowDO_clim_noloading[jan:apr]
winter_loading_deep = deepDO_clim_loading[jan:apr]
winter_noloading_deep = deepDO_clim_noloading[jan:apr]

spring_loading_shallow = shallowDO_clim_loading[apr:jul]
spring_noloading_shallow = shallowDO_clim_noloading[apr:jul]
spring_loading_deep = deepDO_clim_loading[apr:jul]
spring_noloading_deep = deepDO_clim_noloading[apr:jul]

summer_loading_shallow = shallowDO_clim_loading[jul:oct]
summer_noloading_shallow = shallowDO_clim_noloading[jul:oct]
summer_loading_deep = deepDO_clim_loading[jul:oct]
summer_noloading_deep = deepDO_clim_noloading[jul:oct]

fall_loading_shallow = shallowDO_clim_loading[oct::]
fall_noloading_shallow = shallowDO_clim_noloading[oct::]
fall_loading_deep = deepDO_clim_loading[oct::]
fall_noloading_deep = deepDO_clim_noloading[oct::]

print('Winter')
print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[jan:apr]),3)))
print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[jan:apr]),3)))
print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[jan:apr]),3)))
print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[jan:apr]),3)))

print('Spring')
print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[apr:jul]),3)))
print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[apr:jul]),3)))
print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[apr:jul]),3)))
print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[apr:jul]),3)))

print('Summer')
print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[jul:oct]),3)))
print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[jul:oct]),3)))
print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[jul:oct]),3)))
print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[jul:oct]),3)))

print('Fall')
print('   surf (loading): {}'.format(round(np.nanmean(shallowDO_clim_loading[oct::]),3)))
print('   surf (no-loading): {}'.format(round(np.nanmean(shallowDO_clim_noloading[oct::]),3)))
print('   deep (loading): {}'.format(round(np.nanmean(deepDO_clim_loading[oct::]),3)))
print('   deep (no-loading): {}'.format(round(np.nanmean(deepDO_clim_noloading[oct::]),3)))


fig,ax = plt.subplots(1,1, figsize=(9,4), sharey=True)
ax.plot(shallowDO_clim_loading)
ax.plot(shallowDO_clim_noloading)
ax.plot(deepDO_clim_loading)
ax.plot(deepDO_clim_noloading)



##############################################################
##                        PLOTTING                          ##
##############################################################

# Initialize figure
# fs = 10
# pfun.start_plot(fs=fs, figsize=(10,8))
fig,axes = plt.subplots(1,2, figsize=(8,2.4), sharey=True)
ax = axes.ravel()

# # plot shallow layer (depth <= 10 m)
# violinplot([jan_noloading_shallow,
#             feb_noloading_shallow,
#             mar_noloading_shallow,
#             apr_noloading_shallow,
#             may_noloading_shallow,
#             jun_noloading_shallow,
#             jul_noloading_shallow,
#             aug_noloading_shallow,
#             sep_noloading_shallow,
#             oct_noloading_shallow,
#             nov_noloading_shallow,
#             dec_noloading_shallow,
#             ],
#             positions=[1,2,3,4,5,6,7,8,9,10,11,12], show_boxplot=False,
#             side='left', ax=ax[0], plot_opts={'violin_fc':'black'})
# violinplot([jan_noloading_deep,
#             feb_noloading_deep,
#             mar_noloading_deep,
#             apr_noloading_deep,
#             may_noloading_deep,
#             jun_noloading_deep,
#             jul_noloading_deep,
#             aug_noloading_deep,
#             sep_noloading_deep,
#             oct_noloading_deep,
#             nov_noloading_deep,
#             dec_noloading_deep,
#             ],
#             positions=[1,2,3,4,5,6,7,8,9,10,11,12], show_boxplot=False,
#             side='left', ax=ax[1], plot_opts={'violin_fc':'black'})

# # plot deep layer (depth <= 10 m)
# violinplot([jan_loading_shallow,
#             feb_loading_shallow,
#             mar_loading_shallow,
#             apr_loading_shallow,
#             may_loading_shallow,
#             jun_loading_shallow,
#             jul_loading_shallow,
#             aug_loading_shallow,
#             sep_loading_shallow,
#             oct_loading_shallow,
#             nov_loading_shallow,
#             dec_loading_shallow,
#             ],
#             positions=[1,2,3,4,5,6,7,8,9,10,11,12], show_boxplot=False,
#             side='right', ax=ax[0], plot_opts={'violin_fc':'cornflowerblue'})
# violinplot([jan_loading_deep,
#             feb_loading_deep,
#             mar_loading_deep,
#             apr_loading_deep,
#             may_loading_deep,
#             jun_loading_deep,
#             jul_loading_deep,
#             aug_loading_deep,
#             sep_loading_deep,
#             oct_loading_deep,
#             nov_loading_deep,
#             dec_loading_deep,
#             ],
#             positions=[1,2,3,4,5,6,7,8,9,10,11,12], show_boxplot=False,
#             side='right', ax=ax[1], plot_opts={'violin_fc':'cornflowerblue'})

# # # format figure
# ax[0].set_ylabel('DO concentration [mg/L]\n(Surface 10 m)', fontsize=12)
# ax[1].set_ylabel('DO concentration [mg/L]\n(Deeper than 10 m)', fontsize=12)
# for axis in ax:
#     axis.grid(visible=True, axis='both', color='silver', linestyle='--')
#     axis.tick_params(axis='both', labelsize=12)
#     axis.tick_params(axis='x',rotation=30)
#     # axis.set_ylim([0,11])
#     xticks = np.arange(1, 13, 1)
#     axis.set_xticks(xticks, labels=['Jan','Feb','Mar','Apr','May','Jun',
#                                      'Jul','Aug','Sep','Oct','Nov','Dec'])
# ax[0].text(9.5,10,'No-Loading',fontweight='bold',color='black',alpha=0.6,fontsize=14)
# ax[0].text(9.5,9.3,'Loading',fontweight='bold',color='cornflowerblue',alpha=0.6,fontsize=14)

# plt.tight_layout()


# plot shallow layer (depth <= 10 m)
violinplot([winter_noloading_shallow,spring_noloading_shallow,summer_noloading_shallow,fall_noloading_shallow],
            positions=[1,2,3,4], show_boxplot=False,
            side='left', ax=ax[0], plot_opts={'violin_fc':'black'})
violinplot([winter_loading_shallow,spring_loading_shallow,summer_loading_shallow,fall_loading_shallow],
            positions=[1,2,3,4], show_boxplot=False,
            side='right', ax=ax[0], plot_opts={'violin_fc':'cornflowerblue'})

# plot deep layer (depth <= 10 m)
violinplot([winter_noloading_deep,spring_noloading_deep,summer_noloading_deep,fall_noloading_deep],
            positions=[1,2,3,4], show_boxplot=False,
            side='left', ax=ax[1], plot_opts={'violin_fc':'black'})
violinplot([winter_loading_deep,spring_loading_deep,summer_loading_deep,fall_loading_deep],
            positions=[1,2,3,4], show_boxplot=False,
            side='right', ax=ax[1], plot_opts={'violin_fc':'cornflowerblue'})

# format figure
ax[0].set_ylabel('DO [mg/L]', fontsize=12)
# ax[1].set_ylabel('DO [mg/L]', fontsize=12)
for axis in ax:
    axis.grid(visible=True, axis='both', color='silver', linestyle='--')
    axis.tick_params(axis='both', labelsize=12)
    # axis.tick_params(axis='x',rotation=30)
    # axis.set_ylim([0,11])
    xticks = np.arange(1, 5, 1)
    axis.set_xticks(xticks, labels=['Winter','Spring','Summer','Fall'])
ax[0].text(0.6,5.1,'No-Loading',fontweight='bold',color='black',alpha=0.6,fontsize=11)
ax[0].text(0.6,4.5,'Loading',fontweight='bold',color='cornflowerblue',alpha=0.6,fontsize=11)

ax[0].text(0.6,10.5,'(a) Surface 10 m',fontweight='bold',fontsize=12)
ax[1].text(0.6,10.5,'(b) Deep',fontweight='bold',fontsize=12)

plt.tight_layout()