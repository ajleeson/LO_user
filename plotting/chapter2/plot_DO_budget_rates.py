"""
Plot DO bgc rates for all sub-basins
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import get_two_layer
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t1_x11b'
year = '2017'

basins = ['hoodcanal','whidbey','mainbasin','southsound']
basin_seg_dict = {'hoodcanal':'hc_m',
                  'whidbey':'wb_p',
                  'mainbasin':'hc_p',
                  'southsound':'ss_m'}

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# set up dates
startdate = year + '.01.01'
# enddate = year + '.12.31'
# enddate_hrly = str(int(year)+1)+'.01.01 00:00:00'

enddate = year + '.01.02'
enddate_hrly = '2017.01.02 23:00:00'

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')[2::]
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]
# crop time vector (because we only have jan 2 - dec 30)
dates_no_crop = dates_local_daily
dates_local_daily = dates_local_daily

print('\n')

##########################################################
##            Get all variables for analysis            ##
##########################################################

print('Getting all data for analysis\n')

# # initialize dataframe for saving
# inlet_budget_df = pd.DataFrame(columns=['Inlet', 'QinDOin', 'QinDOin_err', 'QoutDOout',
#        'QoutDOout_err', 'Photosynthesis', 'Photosynthesis_err',
#        'Consumption', 'Consumption_err', 'd/dt(DO)', 'd/dt(DO)_err',
#        'PhysicalResupply', 'PhysicalResupply_err', 'NetEcosystemMetabolism',
#        'NetEcosystemMetabolism_err', 'SepOctInletDO[mg/L]', 'SepOctInletDO_err[mg/L]',
#        'MeanDepth[m]'])

# COLLAPSE
for i,basin in enumerate(basins):

    # initialize figure
    fig,ax = plt.subplots(1,1, figsize=(8,6))

    plt.suptitle(basin + '(10-day Hanning Window)',fontsize=14, fontweight='bold')

# ---------------------------------- get BGC rate terms --------------------------------------------
#     bgc_dir = Ldir['LOo'] / 'chapter_2' / 'data' / ('DO_budget_terms_' + startdate + '_' + enddate) / basin
#     # get months
#     months = [year+'.01.01_'+year+'.01.31',
#                 year+'.02.01_'+year+'.02.28',
#                 year+'.03.01_'+year+'.03.31',
#                 year+'.04.01_'+year+'.04.30',
#                 year+'.05.01_'+year+'.05.31',
#                 year+'.06.01_'+year+'.06.30',
#                 year+'.07.01_'+year+'.07.31',
#                 year+'.08.01_'+year+'.08.31',
#                 year+'.09.01_'+year+'.09.30',
#                 year+'.10.01_'+year+'.10.31',
#                 year+'.11.01_'+year+'.11.30',
#                 year+'.12.01_'+year+'.12.31',]
    bgc_dir = Ldir['LOo'] / 'chapter_2' / 'data' / ('DO_budget_terms_' + startdate + '_2017.12.31') / basin
    months = [year+'.01.01_'+year+'.01.02']

            
    # initialize arrays to save values
    # LOADING RUN
    photosynthesis_unfiltered_load = []
    nitrification_unfiltered_load = []
    respiration_unfiltered_load = [] # water column respiration
    sod_unfiltered_load = [] # sediment oxygen demand
    airsea_unfiltered_load = []
    o2vol_unfiltered_load = []
    # NO-LOADING RUN
    photosynthesis_unfiltered_noload = []
    nitrification_unfiltered_noload = []
    respiration_unfiltered_noload = [] # water column respiration
    sod_unfiltered_noload = [] # sediment oxygen demand
    airsea_unfiltered_noload = []
    o2vol_unfiltered_noload = []

    # combine all months
    for month in months:

        # LOADING RUN
        gtagex = 'cas7_t1_x11b'
        fn = gtagex + '_' + month + '.p'
        df_bgc = pd.read_pickle(bgc_dir/fn)
        # conversion factor to go from mmol O2/hr to kmol O2/s
        conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
        # get photosynthesis
        photosynthesis_unfiltered_load = np.concatenate((photosynthesis_unfiltered_load, df_bgc['photo [mmol/hr]'].values * conv)) # kmol/s
        # get nitrification
        nitrification_unfiltered_load = np.concatenate((nitrification_unfiltered_load, df_bgc['nitri [mmol/hr]'].values
                                                   * conv * -1)) # kmol/s; multiply by -1 b/c loss term
        # get water column respiration
        respiration_unfiltered_load = np.concatenate((respiration_unfiltered_load, df_bgc['respi [mmol/hr]'].values
                                                   * conv * -1)) # kmol/s; multiply by -1 b/c loss term
        # get sediment oxygen demand
        sod_unfiltered_load = np.concatenate((sod_unfiltered_load, df_bgc['SOD [mmol/hr]'].values
                                                   * conv * -1)) # kmol/s; multiply by -1 b/c loss term
        # get air-sea gas exchange
        airsea_unfiltered_load = np.concatenate((airsea_unfiltered_load, df_bgc['airsea [mmol/hr]'].values * conv)) # kmol/s
        # get (DO*V)
        o2vol_unfiltered_load = np.concatenate((o2vol_unfiltered_load, df_bgc['DO*V [mmol]'].values)) # mmol

        # NO-LOADING RUN
        gtagex = 'cas7_t1noDIN_x11b'
        fn = gtagex + '_' + month + '.p'
        df_bgc = pd.read_pickle(bgc_dir/fn)
        # conversion factor to go from mmol O2/hr to kmol O2/s
        conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
        # get photosynthesis
        photosynthesis_unfiltered_noload = np.concatenate((photosynthesis_unfiltered_noload, df_bgc['photo [mmol/hr]'].values * conv)) # kmol/s
        # get nitrification
        nitrification_unfiltered_noload = np.concatenate((nitrification_unfiltered_noload, df_bgc['nitri [mmol/hr]'].values
                                                   * conv * -1)) # kmol/s; multiply by -1 b/c loss term
        # get water column respiration
        respiration_unfiltered_noload = np.concatenate((respiration_unfiltered_noload, df_bgc['respi [mmol/hr]'].values
                                                   * conv * -1)) # kmol/s; multiply by -1 b/c loss term
        # get sediment oxygen demand
        sod_unfiltered_noload = np.concatenate((sod_unfiltered_noload, df_bgc['SOD [mmol/hr]'].values
                                                   * conv * -1)) # kmol/s; multiply by -1 b/c loss term
        # get air-sea gas exchange
        airsea_unfiltered_noload = np.concatenate((airsea_unfiltered_noload, df_bgc['airsea [mmol/hr]'].values * conv)) # kmol/s
        # get (DO*V)
        o2vol_unfiltered_noload = np.concatenate((o2vol_unfiltered_noload, df_bgc['DO*V [mmol]'].values)) # mmol

    # take time derivative of (DO*V) to get d/dt (DO*V)
    # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
    ddtDOV_unfiltered_load = np.diff((o2vol_unfiltered_load)) * conv
    ddtDOV_unfiltered_noload = np.diff((o2vol_unfiltered_noload)) * conv

#     # get DO concentration [mg/L]
#     DO_deep_unfiltered = o2vol_deep_unfiltered/vol_deep_unfiltered * 32/1000 # mg/L
#     DO_total_unfiltered = (o2vol_deep_unfiltered+o2vol_surf_unfiltered) / (vol_deep_unfiltered+vol_surf_unfiltered) * 32/1000 # mg/L

#     # apply Godin filter
#     photo  = zfun.lowpass(photosynthesis_unfiltered, f='godin')[36:-34:24]
#     nitri  = zfun.lowpass(nitrification_unfiltered, f='godin')[36:-34:24]
#     respi  = zfun.lowpass(respiration_unfiltered, f='godin')[36:-34:24]
#     sod    = zfun.lowpass(sod_unfiltered, f='godin')[36:-34:24]
#     airsea = zfun.lowpass(airsea_unfiltered, f='godin')[36:-34:24]
#     ddtDOV = zfun.lowpass(ddtDOV_unfiltered, f='godin')[36:-34:24]

    # plot budget time series
#     ax.plot(dates_local,zfun.lowpass(photosynthesis_unfiltered,n=10),color='#8F0445', label='Photosynthesis')
#     ax.plot(dates_local,zfun.lowpass(nitrification_unfiltered,n=10),color='#FCC2DD', label='Nitrification')
#     ax.plot(dates_local,zfun.lowpass(respiration_unfiltered,n=10),color='yellowgreen', label='Respiration') 
#     ax.plot(dates_local,zfun.lowpass(sod_unfiltered,n=10),color='#0D4B91', label='Sediment oxygen demand')
#     ax.plot(dates_local,zfun.lowpass(airsea_unfiltered,n=10),color='teal', label='Air-Sea') 

    # plot no-loading
    ax.plot(dates_local,photosynthesis_unfiltered_noload,color='#8F0445', label='Photosynthesis', linewidth=2, alpha=0.5)
    ax.plot(dates_local,nitrification_unfiltered_noload,color='#FCC2DD', label='Nitrification', linewidth=2, alpha=0.5)
    ax.plot(dates_local,respiration_unfiltered_noload,color='yellowgreen', label='Respiration', linewidth=2, alpha=0.5)
    ax.plot(dates_local,sod_unfiltered_noload,color='#0D4B91', label='Sediment oxygen demand', linewidth=2, alpha=0.5)
    ax.plot(dates_local,airsea_unfiltered_noload,color='teal', label='Air-Sea', linewidth=2, alpha=0.5)

    # plot loading
    ax.plot(dates_local,photosynthesis_unfiltered_load,color='#8F0445', linewidth=1, linestyle='--', alpha=1)
    ax.plot(dates_local,nitrification_unfiltered_load,color='#FCC2DD', linewidth=1, linestyle='--', alpha=1)
    ax.plot(dates_local,respiration_unfiltered_load,color='yellowgreen', linewidth=1, linestyle='--', alpha=1)
    ax.plot(dates_local,sod_unfiltered_load,color='#0D4B91', linewidth=1, linestyle='--', alpha=1)
    ax.plot(dates_local,airsea_unfiltered_load,color='teal', linewidth=1, linestyle='--', alpha=1)


    ax.legend(loc='upper right')



#     # get budget terms during decline period (volume normalize and convert to mg/L/day)
#     conversion = 1000 * 32 * 60 * 60 * 24 

#     deep_exchange_nonnormalized = np.nanmean(TEF_deep[minday:maxday])
#     photosynthesis_all_nonnormalized = np.nanmean(photo_total[minday:maxday])

#     deep_exchange_all = TEF_deep[minday:maxday]/vol_total[minday:maxday]
#     deep_exchange_avg = np.nanmean(deep_exchange_all) * conversion
#     deep_exchange_err = np.nanstd(deep_exchange_all) * conversion

#     surf_exchange_transport_all = TEF_surf[minday:maxday]/vol_total[minday:maxday]
#     surf_exchange_transport_avg = np.nanmean(surf_exchange_transport_all) * conversion
#     surf_exchange_transport_err = np.nanstd(surf_exchange_transport_all) * conversion

#     photosynthesis_all = photo_total[minday:maxday]/vol_total[minday:maxday]
#     photosynthesis_avg = np.nanmean(photosynthesis_all) * conversion
#     photosynthesis_err = np.nanstd(photosynthesis_all) * conversion

#     consumption_all = cons_total[minday:maxday]/vol_total[minday:maxday]
#     consumption_avg = np.nanmean(consumption_all) * conversion
#     consumption_err = np.nanstd(consumption_all) * conversion

#     # get percent of consumption occuring sub-oxycline
#     suboxycline_perc = np.nanmean(cons_deep[minday:maxday]/cons_total[minday:maxday]) * 100
#     print(f'{inlet_name} sub-oxycline consumption: {suboxycline_perc:.1f}%')
    

#     ddtDO_all = ddtDOV_total[minday:maxday]/vol_total[minday:maxday]
#     ddtDO_avg = np.nanmean(ddtDO_all) * conversion
#     ddtDO_err = np.nanstd(ddtDO_all) * conversion

#     airsea_all = airsea_surf[minday:maxday]/vol_total[minday:maxday]
#     airsea_avg = np.nanmean(airsea_all) * conversion
#     airsea_err = np.nanstd(airsea_all) * conversion

#     rivers_all = traps_total_DO[minday:maxday]/vol_total[minday:maxday]
#     rivers_avg = np.nanmean(rivers_all) * conversion
#     rivers_err = np.nanstd(rivers_all) * conversion

#     error_all = error_DO[minday:maxday]/vol_total[minday:maxday]
#     error_avg = np.nanmean(error_all) * conversion
#     error_err = np.nanstd(error_all) * conversion

#     physresup_avg = np.nanmean(deep_exchange_all+surf_exchange_transport_all) * conversion # Physical resupply = Exchange flow + Vertical transport
#     physresup_err = np.nanstd(deep_exchange_all+surf_exchange_transport_all) * conversion

#     NEM_avg = np.nanmean(photosynthesis_all+consumption_all) * conversion # NEM = Photosynthesis + Consumption
#     NEM_err = np.nanstd(photosynthesis_all+consumption_all) * conversion

#     # hypoxic season deep DO
#     hypminday = 242
#     hypmaxday = 302
#     DOinlet_avg = np.nanmean(DO_inlet[hypminday:hypmaxday])
#     DOinlet_err = np.nanstd(DO_inlet[hypminday:hypmaxday])

#     # decline period DO concentrations
#     DO_in = DO_p * 32/1000
#     DOin_DOinlet = np.nanmean(DO_in[minday:maxday] - DO_inlet[minday:maxday])
#     DOin_DOinlet_err = np.nanstd(DO_in[minday:maxday] - DO_inlet[minday:maxday])

#     # calculate mean depth
#     fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'vol_df_cas7_c21.p'
#     vol_df = pd.read_pickle(fn)
#     inlet_vol = vol_df['volume m3'].loc[station+'_p']
#     inlet_area = vol_df['area m2'].loc[station+'_p']
#     mean_depth = inlet_vol / inlet_area

#     # print fraction of inlet mean depth that is sub-oxycline
#     # open interface depths
#     # get interface depth from csv file
#     with open('interface_depths_oxycline.csv', 'r') as f:
#         for line in f:
#             inlet, interface_depth = line.strip().split(',')
#             interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
#     z_interface = float(interface_dict[station])
#     suboxy_depth = mean_depth + z_interface # interface is measured from the surface and is a negative number
#     print('fraction of inlet mean depth that is sub-oxycline: {:.1f}%'.format((suboxy_depth/mean_depth)*100))


#     # add data to df
#     new_data = {'Inlet': [inlet_name],
#                 'QinDOin': [deep_exchange_avg],
#                 'QinDOin_err': [deep_exchange_err],
#                 'QoutDOout': [surf_exchange_transport_avg],
#                 'QoutDOout_err': [surf_exchange_transport_err],
#                 'Photosynthesis': [photosynthesis_avg],
#                 'Photosynthesis_err': [photosynthesis_err],
#                 'Consumption': [consumption_avg],
#                 'Consumption_err': [consumption_err],
#                 'd/dt(DO)': [ddtDO_avg],
#                 'd/dt(DO)_err': [ddtDO_err],
#                 'PhysicalResupply': [physresup_avg],
#                 'PhysicalResupply_err': [physresup_err],
#                 'NetEcosystemMetabolism': [NEM_avg],
#                 'NetEcosystemMetabolism_err': [NEM_err],
#                 'SepOctInletDO[mg/L]': [DOinlet_avg],
#                 'SepOctInletDO_err[mg/L]': [DOinlet_err],
#                 'MeanDepth[m]': [mean_depth],
#                 'Error': [error_avg],
#                 'Error_err': [error_err],
#                 'AirSea': [airsea_avg],
#                 'AirSea_err': [airsea_err],
#                 'Rivers': [rivers_avg],
#                 'Rivers_err': [rivers_err]}
#                 # 'QinDOin_nonorm': [deep_exchange_nonnormalized],
#                 # 'Photo_nonorm': [photosynthesis_all_nonnormalized]}
#     df_new_rows = pd.DataFrame(new_data)
#     inlet_budget_df = pd.concat([inlet_budget_df, df_new_rows],ignore_index=True)

#     # save values to dictionary (with 30-day Hanning Window filter applied)
#     # Note that these values have already been Godin-filters, and we are applying
#     # a Hanning window filter on top of that.
#     if station == 'elliot':
#         station = 'elliott' # correct typo
#     DOTI_timeseries[station] = zfun.lowpass(DO_inlet,n=30) # mg/L
#     DOin_DOout_timeseries[station] = zfun.lowpass((DO_p.values-DO_m.values)* 32/1000,n=30) # mg/L
#     Qin_Qout_timeseries[station] = zfun.lowpass(Q_p.values+Q_m.values,n=30) # m3/s

# ######################################################################
# # BUDGET ERROR

# # calculate budget error (mg/L per day) ------------------------------
#     conversion = (1000 * 32 * 60 * 60 * 24)
#     error_budget = (error_DO/inlet_vol) * conversion # [mg/L/day]
#     inlet_error_ann_avg = np.nanmean(error_budget)
#     # calculate QinDOin (mg/L per day) 
#     QinDOin = (TEF_deep/inlet_vol) * conversion # [mg/L/day]
#     inlet_QinDOin_ann_avg = np.nanmean(QinDOin)
#     # calculate biological consumption in deep layer (mg/L per day)
#     consumption = (cons_deep/inlet_vol) * conversion # [mg/L/day]
#     consumption_1lay = ((cons_deep+cons_surf)/inlet_vol) * conversion # [mg/L/day]
#     inlet_consumption_ann_avg = np.nanmean(consumption)
#     # calculate d/dt(DO) (mg/L per day)
#     ddtDO = (ddtDOV_deep/inlet_vol) * conversion # [mg/L/day]
#     inlet_ddtDO_ann_avg = np.nanmean(ddtDO)

#     # calculating division before annual averaging
#     # # full year
#     # err_minday = 0
#     # err_maxday = 363
#     # winter
#     err_minday = 0
#     err_maxday = 90
#     # # spring
#     # err_minday = 90
#     # err_maxday = 181
#     # # summer
#     # err_minday = 181
#     # err_maxday = 272
#     # # fall
#     # err_minday = 272
#     # err_maxday = 363

#     error_mgL_ann_avg.append(np.nanmean(error_budget[err_minday:err_maxday])) # [mg/L/day]

#     # print(station)
#     # print(error_mgL_ann_avg[i])

#     error_QinDOin_ann_avg.append(np.abs(np.nanmean(error_DO[err_minday:err_maxday])/np.nanmean(TEF_deep[err_minday:err_maxday])))
#     error_consumption_ann_avg.append(np.abs(np.nanmean(error_DO[err_minday:err_maxday])/
#                                             np.nanmean(cons_deep[err_minday:err_maxday]+cons_surf[err_minday:err_maxday])))

#     decline_per_vol_norm_error = np.nanmean((error_DO[err_minday:err_maxday]/vol_deep[err_minday:err_maxday]))* conversion

#     # if station == 'quartermaster':
#     #     print(station + ' annual mean error [mg/L/day]')
#     #     print(error_mgL_ann_avg[i])
#     #     print(np.nanstd(error_budget))
#     #     print(station + '(annual mean error)/(annual mean QinDOin) [expressed as percentage]')
#     #     print('    {}%'.format(round(np.abs(error_QinDOin_ann_avg[i]) * 100,2)))
#     #     print(station + '(annual mean error)/(annual mean consumption) [expressed as percentage]')
#     #     print('    {}%'.format(round(np.abs(error_consumption_ann_avg[i]) * 100,2)))

# # calculate bulk statistics
# error_QinDOin = np.abs(np.nanmean(error_QinDOin_ann_avg)) * 100
# error_consumption = np.abs(np.nanmean(error_consumption_ann_avg)) * 100
# error_ddtDO = np.abs(np.nanmean(error_ddtDO_ann_avg)) * 100
# error_ddtDO_onelayer = np.abs(np.nanmean(error_ddtDO_onelayer_ann_avg)) * 100

# print('-----------------------------')
# print('max percent of consumption error error')
# print(np.nanmax(np.abs(error_consumption_ann_avg)))
# print('-----------------------------')


# # print bulk statistics
# print('(annual mean error)/(annual mean QinDOin) [expressed as percentage]')
# print('    {}%'.format(round(error_QinDOin,2)))
# print('\n')
# print('(annual mean error)/(annual mean consumption) [expressed as percentage]')
# print('    {}%'.format(round(error_consumption,4)))


# print('\n')
# print('annual mean error [mg/L/day]')
# print(np.nanmean(error_mgL_ann_avg))
# print(np.nanstd(error_mgL_ann_avg))


# # save dictsto csv file
# # dates
# dates = pd.date_range(start='2017-01-02', end='2017-12-30', freq='D')
# date_list = dates.strftime('%Y-%m-%d').tolist()
# # data
# DOTI_df = pd.DataFrame.from_dict(DOTI_timeseries)
# DOTI_df.insert(0, 'date', date_list)
# DOin_DOout_df = pd.DataFrame.from_dict(DOin_DOout_timeseries)
# DOin_DOout_df.insert(0, 'date', date_list)
# Qin_Qout_df = pd.DataFrame.from_dict(Qin_Qout_timeseries)
# Qin_Qout_df.insert(0, 'date', date_list)
# # save
# DOTI_df.to_csv('../../../../terminal_inlet_DO_rev3/terminletDO_mgL_30dayHanning.csv', index=False)
# DOin_DOout_df.to_csv('../../../../terminal_inlet_DO_rev3/DOin_DOout_mgL_30dayHanning.csv', index=False)


# # exchange flow vs. inlet volume
# fig, ax = plt.subplots(1,1,figsize=(6,6))
# # ax.scatter(inlet_budget_df['QinDOin_nonorm'],inlet_budget_df['Photo_nonorm'],
# #                     s=50, edgecolor='white',linewidth=0.5,zorder=5)
# # for inlet in inlet_budget_df['Inlet']:
# #       ax.text(inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'QinDOin_nonorm'],
# #               inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'Photo_nonorm']+0.005,
# #                 inlet,color='silver')
# ax.scatter(inlet_budget_df['QinDOin'],inlet_budget_df['Photosynthesis'],
#                     s=50, edgecolor='white',linewidth=0.5,zorder=5)
# for inlet in inlet_budget_df['Inlet']:
#       ax.text(inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'QinDOin'],
#               inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'Photosynthesis']+0.005,
#                 inlet,color='silver')
# plt.show()

# # get average of all inlets and add data to df
# # error of averages is error of all values added in quadrature, divided by number of items
# new_data = {'Inlet': ['All inlets'],
#             'QinDOin': [np.nanmean(inlet_budget_df['QinDOin'])],
#             'QinDOin_err': [np.nanstd(inlet_budget_df['QinDOin_err'])],
#             'QoutDOout': [np.nanmean(inlet_budget_df['QoutDOout'])],
#             'QoutDOout_err': [np.nanstd(inlet_budget_df['QoutDOout_err'])],
#             'Photosynthesis': [np.nanmean(inlet_budget_df['Photosynthesis'])],
#             'Photosynthesis_err': [np.nanstd(inlet_budget_df['Photosynthesis'])],
#             'Consumption': [np.nanmean(inlet_budget_df['Consumption'])],
#             'Consumption_err': [np.nanstd(inlet_budget_df['Consumption'])],
#             'd/dt(DO)': [np.nanmean(inlet_budget_df['d/dt(DO)'])],
#             'd/dt(DO)_err': [np.nanstd(inlet_budget_df['d/dt(DO)'])],
#             'PhysicalResupply': [np.nanmean(inlet_budget_df['PhysicalResupply'])],
#             'PhysicalResupply_err': [np.nanstd(inlet_budget_df['PhysicalResupply'])],
#             'NetEcosystemMetabolism': [np.nanmean(inlet_budget_df['NetEcosystemMetabolism'])],
#             'NetEcosystemMetabolism_err': [np.nanstd(inlet_budget_df['NetEcosystemMetabolism'])],
#             'SepOctInletDO[mg/L]': [np.nanmean(inlet_budget_df['SepOctInletDO[mg/L]'])],
#             'SepOctInletDO_err[mg/L]': [np.nanstd(inlet_budget_df['SepOctInletDO[mg/L]'])],
#             'MeanDepth[m]': [np.nanmean(inlet_budget_df['MeanDepth[m]'])],
#             'Error': [np.nanmean(inlet_budget_df['Error'])],
#             'Error_err': [np.nanstd(inlet_budget_df['Error_err'])],
#             'AirSea': [np.nanmean(inlet_budget_df['AirSea'])],
#             'AirSea_err': [np.nanstd(inlet_budget_df['AirSea'])],
#             'Rivers': [np.nanmean(inlet_budget_df['Rivers'])],
#             'Rivers_err': [np.nanstd(inlet_budget_df['Rivers'])]
#             ,}
# df_new_rows = pd.DataFrame(new_data)
# inlet_budget_df = pd.concat([inlet_budget_df, df_new_rows],ignore_index=True)

# # # save to csv file
# # print(inlet_budget_df)
# # inlet_budget_df.to_csv('inlet_budgets_decline_period_mgL_day_onelayer.csv', index=False)


# # # save lynch cove budget to csv file
# # lynchcove_budget_df = pd.DataFrame.from_dict(lynchcove_dict_10dayhanning)
# # dates = pd.date_range(start='2017-01-02', end='2017-12-30', freq='D')
# # date_list = dates.strftime('%Y-%m-%d').tolist()
# # lynchcove_budget_df.insert(0, 'date', date_list)
# # lynchcove_budget_df.to_csv('../../../../terminal_inlet_DO_rev3/lynchcove_2017_budget_kmolO2_s_10dayHanning_onelayer.csv', index=False)
# # print(lynchcove_budget_df)
