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

# enddate = year + '.01.02'
# enddate_hrly = '2017.01.02 23:00:00'
enddate = year + '.01.31'
enddate_hrly = '2017.01.31 23:00:00'

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
    months = [year+'.01.01_'+year+'.01.31']

            
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

    # apply Godin filter
    # loading
    photo_load  = zfun.lowpass(photosynthesis_unfiltered_load, f='godin')[36:-34:24]
    nitri_load  = zfun.lowpass(nitrification_unfiltered_load, f='godin')[36:-34:24]
    respi_load  = zfun.lowpass(respiration_unfiltered_load, f='godin')[36:-34:24]
    sod_load    = zfun.lowpass(sod_unfiltered_load, f='godin')[36:-34:24]
    airsea_load = zfun.lowpass(airsea_unfiltered_load, f='godin')[36:-34:24]
    ddtDOV_load = zfun.lowpass(ddtDOV_unfiltered_load, f='godin')[36:-34:24]
    # no loading
    photo_noload  = zfun.lowpass(photosynthesis_unfiltered_noload, f='godin')[36:-34:24]
    nitri_noload  = zfun.lowpass(nitrification_unfiltered_noload, f='godin')[36:-34:24]
    respi_noload  = zfun.lowpass(respiration_unfiltered_noload, f='godin')[36:-34:24]
    sod_noload    = zfun.lowpass(sod_unfiltered_noload, f='godin')[36:-34:24]
    airsea_noload = zfun.lowpass(airsea_unfiltered_noload, f='godin')[36:-34:24]
    ddtDOV_noload = zfun.lowpass(ddtDOV_unfiltered_noload, f='godin')[36:-34:24]

    # # UNFILETERED!!!!!!!!!!!!!!!!!
    # # plot no-loading
    # ax.plot(dates_local,photosynthesis_unfiltered_noload,color='#8F0445', label='Photosynthesis', linewidth=2, alpha=0.5)
    # ax.plot(dates_local,nitrification_unfiltered_noload,color='#FCC2DD', label='Nitrification', linewidth=2, alpha=0.5)
    # ax.plot(dates_local,respiration_unfiltered_noload,color='yellowgreen', label='Respiration', linewidth=2, alpha=0.5)
    # ax.plot(dates_local,sod_unfiltered_noload,color='#0D4B91', label='Sediment oxygen demand', linewidth=2, alpha=0.5)
    # ax.plot(dates_local,airsea_unfiltered_noload,color='teal', label='Air-Sea', linewidth=2, alpha=0.5)
    # # plot loading
    # ax.plot(dates_local,photosynthesis_unfiltered_load,color='#8F0445', linewidth=1, linestyle='--', alpha=1)
    # ax.plot(dates_local,nitrification_unfiltered_load,color='#FCC2DD', linewidth=1, linestyle='--', alpha=1)
    # ax.plot(dates_local,respiration_unfiltered_load,color='yellowgreen', linewidth=1, linestyle='--', alpha=1)
    # ax.plot(dates_local,sod_unfiltered_load,color='#0D4B91', linewidth=1, linestyle='--', alpha=1)
    # ax.plot(dates_local,airsea_unfiltered_load,color='teal', linewidth=1, linestyle='--', alpha=1)

    # GODIN-FILTERED!!!!!!!!!!!!!!!!!!
    # plot no-loading
    ax.plot(dates_local_daily,photo_noload,color='#8F0445', label='Photosynthesis', linewidth=2, alpha=0.5)
    ax.plot(dates_local_daily,nitri_noload,color='#FCC2DD', label='Nitrification', linewidth=2, alpha=0.5)
    ax.plot(dates_local_daily,respi_noload,color='yellowgreen', label='Respiration', linewidth=2, alpha=0.5)
    ax.plot(dates_local_daily,sod_noload,color='#0D4B91', label='Sediment oxygen demand', linewidth=2, alpha=0.5)
    ax.plot(dates_local_daily,airsea_noload,color='teal', label='Air-Sea', linewidth=2, alpha=0.5)
    # plot loading
    ax.plot(dates_local_daily,photo_load,color='#8F0445', linewidth=1, linestyle='--', alpha=1)
    ax.plot(dates_local_daily,nitri_load,color='#FCC2DD', linewidth=1, linestyle='--', alpha=1)
    ax.plot(dates_local_daily,respi_load,color='yellowgreen', linewidth=1, linestyle='--', alpha=1)
    ax.plot(dates_local_daily,sod_load,color='#0D4B91', linewidth=1, linestyle='--', alpha=1)
    ax.plot(dates_local_daily,airsea_load,color='teal', linewidth=1, linestyle='--', alpha=1)


    ax.legend(loc='upper right')

