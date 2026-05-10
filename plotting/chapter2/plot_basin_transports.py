"""
Plot transport of NPZD+O terms all basins
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import get_two_layer
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from matplotlib import colormaps
from matplotlib.colors import ListedColormap
from scipy.stats import pearsonr
from scipy.linalg import lstsq
from matplotlib.ticker import AutoMinorLocator

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t1_x11b'
year = '2017'

stations = 'all'

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# set up dates
startdate = year + '.01.01'
enddate = year + '.12.31'
enddate_hrly = str(int(year)+1)+'.01.01 00:00:00'

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# Get basin boundaries:
# admiralty inlet, deception pass, hood canal, south sound, whidbey basin
boundaries = ['ai','dp','hc','ss','wb']
boundary_dict = {'ai':'Admiralty Inlet', 'dp':'Deception Pass', 'hc':'Hood Canal',
                  'ss':'South Sound', 'wb':'Whidbey Basin'}

# state variables
vars = ['salt','oxygen', 'NO3', 'NH4', 'phytoplankton', 'zooplankton', 'SdetritusN', 'LdetritusN']


# # where to put output figures
# out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
# Lfun.make_dir(out_dir)

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

# initialize dataframe for saving
inlet_budget_df = pd.DataFrame(columns=['Inlet', 'ExchangeFlow', 'ExchangeFlow_err', 'VerticalTransport',
       'VerticalTransport_err', 'Photosynthesis', 'Photosynthesis_err',
       'Consumption', 'Consumption_err', 'd/dt(DO)', 'd/dt(DO)_err',
       'PhysicalResupply', 'PhysicalResupply_err', 'NetEcosystemMetabolism',
       'NetEcosystemMetabolism_err', 'SepOctDeepDO[mg/L]', 'SepOctDeepDO_err[mg/L]',
       'MeanDepth[m]'])

# COLLAPSE
for i,boundary in enumerate(boundaries):
        
    # initialize figure
    fig,axes = plt.subplots(4,2, figsize=(13,8), sharex=True)
    ax = axes.ravel()

    plt.suptitle(boundary_dict[boundary] + ' (10-day Hanning Window)',fontsize=14, fontweight='bold')

# --------------------------- get TEF terms ----------------------------------------
    in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'cps' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (boundary + '.nc')
    bulk_loading = xr.open_dataset(in_dir)
    tef_df_loading, vn_list, vec_list = get_two_layer.get_two_layer(bulk_loading)

    # in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1noDIN_x11b' / 'tef2' / 'cps' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (boundary + '.nc')
    # bulk_NOloading = xr.open_dataset(in_dir)
    # tef_df_NOloading, vn_list, vec_list = get_two_layer.get_two_layer(bulk_NOloading)

    # get transports
    Q_p_loading   = tef_df_loading['q_p'] # Qin [m3/s]
    Q_m_loading   = tef_df_loading['q_m'] # Qout [m3/s]
    # Q_p_NOloading = tef_df_NOloading['q_p'] # Qin [m3/s]
    # Q_m_NOloading = tef_df_NOloading['q_m'] # Qout [m3/s]

    # loop through variables and plot the transports
    for v, var in enumerate(vars):
        ax[v].set_title(var, fontsize=12, fontweight='bold', loc='left')

        transport_p_loading   = tef_df_loading[var+'_p'].values * Q_p_loading.values
        transport_m_loading   = tef_df_loading[var+'_m'].values * Q_m_loading.values
        # transport_p_NOloading = tef_df_NOloading[var+'_p'].values * Q_p_NOloading.values
        # transport_m_NOloading = tef_df_NOloading[var+'_m'].values * Q_m_NOloading.values

        ax[v].plot(dates_local_daily,zfun.lowpass(transport_p_loading + transport_m_loading,n=10),
                   color='royalblue',label='Loading')
        # ax[v].plot(dates_local_daily,zfun.lowpass(transport_p_NOloading + transport_m_NOloading,n=10),
        #            color='crimson',label='No-Loading')


    # # plot budget time series
    # ax[0].plot(dates_local_daily,zfun.lowpass(Q_m.values,n=10),color='#0D4B91',label='TEF')
    # ax[2].plot(dates_local_daily,zfun.lowpass(Q_p.values,n=10),color='#0D4B91',label='TEF')
    # ax[1].plot(dates_local_daily,zfun.lowpass(TEF_surf,n=10),color='#0D4B91',label='TEF')
    # ax[3].plot(dates_local_daily,zfun.lowpass(TEF_deep,n=10),color='#0D4B91',label='TEF')
    # # ax[3].plot(dates_local_daily,TEF_deep,color='#0D4B91',label='TEF')

    # # format budget time series figures
    # for axnum in [ax[0],ax[1],ax[2],ax[3]]:
    #     axnum.set_xlim([dates_hrly[0],dates_hrly[-2]])
    #     axnum.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')
    #     axnum.tick_params(axis='x', labelrotation=30, labelsize=12)
    #     axnum.tick_params(axis='y', labelsize=12)
    #     loc = mdates.MonthLocator(interval=1)
    #     axnum.xaxis.set_major_locator(loc)
    # ax[3].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    # ax[0].set_ylabel(r'Volume transport [m$^3$ s$^{-1}$]', fontsize=12)
    # ax[2].set_ylabel(r'Volume transport [m$^3$ s$^{-1}$]', fontsize=12)
    # ax[1].set_ylabel(r'DO transport [kmol O$_2$ s$^{-1}$]', fontsize=12)
    # ax[3].set_ylabel(r'DO transport [kmol O$_2$ s$^{-1}$]', fontsize=12)

    # # add data to df
    # new_data = {'Inlet': [inlet_name],
    #             'ExchangeFlow': [deep_exchange_avg],
    #             'ExchangeFlow_err': [deep_exchange_err],
    #             'VerticalTransport': [deep_vert_transport_avg],
    #             'VerticalTransport_err': [deep_vert_transport_err],
    #             'Photosynthesis': [deep_photosynthesis_avg],
    #             'Photosynthesis_err': [deep_photosynthesis_err],
    #             'Consumption': [deep_consumption_avg],
    #             'Consumption_err': [deep_consumption_err],
    #             'd/dt(DO)': [deep_ddtDO_avg],
    #             'd/dt(DO)_err': [deep_ddtDO_err],
    #             'PhysicalResupply': [deep_physresup_avg],
    #             'PhysicalResupply_err': [deep_physresup_err],
    #             'NetEcosystemMetabolism': [deep_NEM_avg],
    #             'NetEcosystemMetabolism_err': [deep_NEM_err],
    #             'SepOctDeepDO[mg/L]': [deep_DO_avg],
    #             'SepOctDeepDO_err[mg/L]': [deep_DO_err],
    #             'MeanDepth[m]': [mean_depth],
    #             'DOin-DOdeep[mg/L]': [DOin_DOdeep],
    #             'DOin-DOdeep_err[mg/L]': [DOin_DOdeep_err],
    #             'Inflow no normalization [kmol/s]': [deep_exchange_nonnormalized],
    #             'Inlet volume [m3]': [total_vol]}
    # df_new_rows = pd.DataFrame(new_data)
    # inlet_budget_df = pd.concat([inlet_budget_df, df_new_rows],ignore_index=True)

    plt.show()
