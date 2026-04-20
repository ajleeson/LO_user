"""
Plot surface and deep budget for all 21 inlets
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
jobname = 'twentyoneinlets'
year = '2017'

stations = 'all'

interface_types = ['tef']#,'og','drdz','tef','halocline','oxycline'] # how to define dividing depth

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

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get stations:
if stations == 'all':
    sta_dict = job_lists.get_sta_dict(jobname)
    # remove lynchcove2
    del sta_dict['lynchcove2']
    # remove shallow inlets (< 10 m deep)
    del sta_dict['hammersley']
    del sta_dict['henderson']
    del sta_dict['oak']
    del sta_dict['totten']
    del sta_dict['similk']
    del sta_dict['budd']
    del sta_dict['eld']
    del sta_dict['killsut']
    # del sta_dict['dabob']
else:
    sta_dict = stations

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


######################
# wind data

# Add daily average wind speed
col_names = [ "YY","MM","DD","hh","mm",
    "WDIR","WSPD","GST","WVHT","DPD","APD","MWD",
    "PRES","ATMP","WTMP","DEWP","VIS","TIDE"]
df = pd.read_csv("46121h2017.txt",
    sep=r"\s+",
    comment="#",
    names=col_names,
    na_values=[99, 99.0, 999, 999.0, 9999, 9999.0])
# rename for datetime parsing
df = df.rename(columns={"YY": "year",
    "MM": "month",
    "DD": "day",
    "hh": "hour",
    "mm": "minute"})
# create datetime index
df["datetime"] = pd.to_datetime(df[["year","month","day","hour","minute"]])
df = df.set_index("datetime").drop(columns=["year","month","day","hour","minute"])
df_daily = df.resample("D").mean()
full_index = pd.date_range(
start="2017-01-02",
end="2017-12-30",
freq="D")
df_daily_aligned = df_daily.reindex(full_index)
df_daily_aligned.index.name = "date"

##########################################################
##            Get all variables for analysis            ##
##########################################################

print('Getting all data for analysis\n')

# get lat and lon of grid
Ldir['ds0'] = startdate
in_dir = Ldir['roms_out'] / Ldir['gtagex']
# G, S, T = zrfun.get_basic_info(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
fn0 = xr.open_dataset(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
lonr = fn0.lon_rho.values
latr = fn0.lat_rho.values
# open box extraction
box_fn = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
ds_box = xr.open_dataset(box_fn)
DX = (ds_box.pm.values)**-1
DY = (ds_box.pn.values)**-1
DA = DX*DY # get area of each grid cell in m^2

# initialize dataframe for saving
inlet_budget_df = pd.DataFrame(columns=['Inlet', 'ExchangeFlow', 'ExchangeFlow_err', 'VerticalTransport',
       'VerticalTransport_err', 'Photosynthesis', 'Photosynthesis_err',
       'Consumption', 'Consumption_err', 'd/dt(DO)', 'd/dt(DO)_err',
       'PhysicalResupply', 'PhysicalResupply_err', 'NetEcosystemMetabolism',
       'NetEcosystemMetabolism_err', 'SepOctDeepDO[mg/L]', 'SepOctDeepDO_err[mg/L]',
       'MeanDepth[m]'])
# initialize lynchcove budget dict
lynchcove_dict_10dayhanning = {}

# initialize error statistics
error_QinDOin_ann_avg = []
error_consumption_ann_avg = []
error_ddtDO_ann_avg = []
error_ddtDO_onelayer_ann_avg = []
error_over_ddtDO = []
error_kmolO2s_ann_avg = []
error_Qinvol_ann_avg = []

# COLLAPSE
for i,station in enumerate(sta_dict):
        
    # initialize figure
    fig,axes = plt.subplots(2,2, figsize=(13,8), sharex=True)
    ax = axes.ravel()

    plt.suptitle(station + '(10-day Hanning Window)',fontsize=14, fontweight='bold')

# --------------------------- get TEF exchange flow terms ----------------------------------------
    in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'c21' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
    bulk = xr.open_dataset(in_dir)
    tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
    Q_p = tef_df['q_p'] # Qin [m3/s]
    Q_m = tef_df['q_m'] # Qout [m3/s]
    DO_p = tef_df['oxygen_p'] # DOin [mmol/m3]
    DO_m = tef_df['oxygen_m'] # DOout [mmol/m3]
    # convert from mmol/s to kmol/s
    TEF_surf = (Q_m.values * DO_m.values) * (1/1000) * (1/1000) # Qout * DOout
    TEF_deep = (Q_p.values * DO_p.values) * (1/1000) * (1/1000) # Qout * DOout

    # plot budget time series
    ax[0].plot(dates_local_daily,zfun.lowpass(Q_m.values,n=10),color='#0D4B91',label='TEF')
    ax[2].plot(dates_local_daily,zfun.lowpass(Q_p.values,n=10),color='#0D4B91',label='TEF')
    ax[1].plot(dates_local_daily,zfun.lowpass(TEF_surf,n=10),color='#0D4B91',label='TEF')
    ax[3].plot(dates_local_daily,zfun.lowpass(TEF_deep,n=10),color='#0D4B91',label='TEF')
    # ax[3].plot(dates_local_daily,TEF_deep,color='#0D4B91',label='TEF')

    # format budget time series figures
    for axnum in [ax[0],ax[1],ax[2],ax[3]]:
        axnum.set_xlim([dates_hrly[0],dates_hrly[-2]])
        axnum.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')
        axnum.tick_params(axis='x', labelrotation=30, labelsize=12)
        axnum.tick_params(axis='y', labelsize=12)
        loc = mdates.MonthLocator(interval=1)
        axnum.xaxis.set_major_locator(loc)
    ax[3].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax[0].set_ylabel(r'Volume transport [m$^3$ s$^{-1}$]', fontsize=12)
    ax[2].set_ylabel(r'Volume transport [m$^3$ s$^{-1}$]', fontsize=12)
    ax[1].set_ylabel(r'DO transport [kmol O$_2$ s$^{-1}$]', fontsize=12)
    ax[3].set_ylabel(r'DO transport [kmol O$_2$ s$^{-1}$]', fontsize=12)

# ---------------------------------- get BGC terms --------------------------------------------
    bgc_dir = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_bgc' / station
    # get months
    months = [year+'.01.01_'+year+'.01.31',
                year+'.02.01_'+year+'.02.28',
                year+'.03.01_'+year+'.03.31',
                year+'.04.01_'+year+'.04.30',
                year+'.05.01_'+year+'.05.31',
                year+'.06.01_'+year+'.06.30',
                year+'.07.01_'+year+'.07.31',
                year+'.08.01_'+year+'.08.31',
                year+'.09.01_'+year+'.09.30',
                year+'.10.01_'+year+'.10.31',
                year+'.11.01_'+year+'.11.30',
                year+'.12.01_'+year+'.12.31',]
    
    # loop through different interface types
    for t,type in enumerate(interface_types):
            
        # initialize arrays to save values
        photo_surf_unfiltered = []
        photo_deep_unfiltered = []
        cons_surf_unfiltered = [] # nitrification, respiration
        cons_deep_unfiltered = [] # nitrification, respiration, sediment oxygen demand
        airsea_surf_unfiltered = []
        o2vol_surf_unfiltered = []
        o2vol_deep_unfiltered = []
        vol_surf_unfiltered = []
        vol_deep_unfiltered = []
        vol_total_unfiltered = []

        # combine all months
        for month in months:
            fn = type + '_' + station + '_' + month + '.p'
            df_bgc = pd.read_pickle(bgc_dir/fn)
            # conversion factor to go from mmol O2/hr to kmol O2/s
            conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
            # get photosynthesis
            photo_surf_unfiltered = np.concatenate((photo_surf_unfiltered, df_bgc['surf photo [mmol/hr]'].values * conv)) # kmol/s
            photo_deep_unfiltered = np.concatenate((photo_deep_unfiltered, df_bgc['deep photo [mmol/hr]'].values * conv)) # kmol/s
            # get consumption
            surf_cons_terms = df_bgc['surf nitri [mmol/hr]'].values + df_bgc['surf respi [mmol/hr]'].values
            deep_cons_terms = df_bgc['deep nitri [mmol/hr]'].values + df_bgc['deep respi [mmol/hr]'].values + df_bgc['SOD [mmol/hr]'].values
            cons_surf_unfiltered = np.concatenate((cons_surf_unfiltered, surf_cons_terms * conv * -1)) # kmol/s; multiply by -1 b/c loss term
            cons_deep_unfiltered = np.concatenate((cons_deep_unfiltered, deep_cons_terms * conv * -1)) # kmol/s; multiply by -1 b/c loss term
            # get air-sea gas exchange
            airsea_surf_unfiltered = np.concatenate((airsea_surf_unfiltered, df_bgc['airsea [mmol/hr]'].values * conv)) # kmol/s
            # get (DO*V)
            o2vol_surf_unfiltered = np.concatenate((o2vol_surf_unfiltered, df_bgc['surf DO*V [mmol]'].values)) # mmol
            o2vol_deep_unfiltered = np.concatenate((o2vol_deep_unfiltered, df_bgc['deep DO*V [mmol]'].values)) # mmol
            # get volume
            vol_surf_unfiltered = np.concatenate((vol_surf_unfiltered, df_bgc['surf vol [m3]'].values)) # m3, and replace zeros with nan
            vol_deep_unfiltered = np.concatenate((vol_deep_unfiltered, df_bgc['deep vol [m3]'].values)) # m3
            vol_total_unfiltered = np.concatenate((vol_total_unfiltered, df_bgc['deep vol [m3]'].values + df_bgc['surf vol [m3]'].values)) # m3

        # take time derivative of (DO*V) to get d/dt (DO*V)
        ddtDOV_surf_unfiltered = np.diff(o2vol_surf_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
        ddtDOV_deep_unfiltered = np.diff(o2vol_deep_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
        ddtDOV_total_unfiltered = np.diff((o2vol_deep_unfiltered+o2vol_surf_unfiltered)) * conv
        # take time derivative of V to get d/dt (V)
        ddtV_surf_unfiltered = np.diff(vol_surf_unfiltered) /60 /60 # diff gets us d(V) dt, where t=1 hr (m3/hr). Then /3600 to get m3/s
        ddtV_deep_unfiltered = np.diff(vol_deep_unfiltered) /60 /60 # diff gets us d(V) dt, where t=1 hr (m3/hr). Then /3600 to get m3/s
        # get DO concentration [mg/L]
        DO_deep_unfiltered = o2vol_deep_unfiltered/vol_deep_unfiltered * 32/1000 # mg/L

        # apply Godin filter
        photo_surf  = zfun.lowpass(photo_surf_unfiltered, f='godin')[36:-34:24]
        photo_deep  = zfun.lowpass(photo_deep_unfiltered, f='godin')[36:-34:24]
        cons_surf   = zfun.lowpass(cons_surf_unfiltered, f='godin')[36:-34:24]
        cons_deep   = zfun.lowpass(cons_deep_unfiltered, f='godin')[36:-34:24]
        airsea_surf = zfun.lowpass(airsea_surf_unfiltered, f='godin')[36:-34:24]
        ddtDOV_surf = zfun.lowpass(ddtDOV_surf_unfiltered, f='godin')[36:-34:24]
        ddtDOV_deep = zfun.lowpass(ddtDOV_deep_unfiltered, f='godin')[36:-34:24]
        ddtDOV_onelayer = zfun.lowpass(ddtDOV_total_unfiltered, f='godin')[36:-34:24]
        ddtvol_surf = zfun.lowpass(ddtV_surf_unfiltered, f='godin')[36:-34:24]
        ddtvol_deep = zfun.lowpass(ddtV_deep_unfiltered, f='godin')[36:-34:24]
        vol_deep = zfun.lowpass(vol_deep_unfiltered, f='godin')[36:-34:24]
        vol_total = zfun.lowpass(vol_total_unfiltered, f='godin')[36:-34:24]
        DO_deep = zfun.lowpass(DO_deep_unfiltered, f='godin')[36:-34:24]

        # plot budget time series
        if t == 0:
            ax[0].plot(dates_local_daily,zfun.lowpass(ddtvol_surf,n=10),color='black',label='d/dt(Volume)')
            ax[2].plot(dates_local_daily,zfun.lowpass(ddtvol_deep,n=10),color='black',label='d/dt(Volume)')

            ax[1].plot(dates_local_daily,zfun.lowpass(photo_surf,n=10),color='#8F0445',label='Photosynthesis')
            ax[3].plot(dates_local_daily,zfun.lowpass(photo_deep,n=10),color='#8F0445',label='Photosynthesis')
            # ax[3].plot(dates_local_daily,photo_deep,color='#8F0445',label='Photosynthesis')

            ax[1].plot(dates_local_daily,zfun.lowpass(cons_surf,n=10),color='#FCC2DD',label='Consumption')
            ax[3].plot(dates_local_daily,zfun.lowpass(cons_deep,n=10),color='#FCC2DD',label='Consumption')
            # ax[3].plot(dates_local_daily,cons_deep,color='#FCC2DD',label='Consumption')

            ax[1].plot(dates_local_daily,zfun.lowpass(ddtDOV_surf,n=10),color='black',label='d/dt(DO)')
            ax[3].plot(dates_local_daily,zfun.lowpass(ddtDOV_deep,n=10),color='black',label='d/dt(DO)')
            # ax[3].plot(dates_local_daily,ddtDOV_deep,color='black',label='d/dt(DO)')

            ax[1].plot(dates_local_daily,zfun.lowpass(airsea_surf,n=10),color='yellowgreen',label='Air-Sea')
        else:
            ax[0].plot(dates_local_daily,zfun.lowpass(ddtvol_surf,n=10),color='black',label='d/dt(Volume)')
            ax[2].plot(dates_local_daily,zfun.lowpass(ddtvol_deep,n=10),color='black',label='d/dt(Volume)')

            ax[1].plot(dates_local_daily,zfun.lowpass(photo_surf,n=10),color='#8F0445')
            ax[3].plot(dates_local_daily,zfun.lowpass(photo_deep,n=10),color='#8F0445')

            ax[1].plot(dates_local_daily,zfun.lowpass(cons_surf,n=10),color='#FCC2DD')
            ax[3].plot(dates_local_daily,zfun.lowpass(cons_deep,n=10),color='#FCC2DD')

            ax[1].plot(dates_local_daily,zfun.lowpass(ddtDOV_surf,n=10),color='black')
            ax[3].plot(dates_local_daily,zfun.lowpass(ddtDOV_deep,n=10),color='black')

            ax[1].plot(dates_local_daily,zfun.lowpass(airsea_surf,n=10),color='yellowgreen') 

# ------------------------------- get rivers and WWTPs ----------------------------------------
    fn = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_traps' / (station + '.p')
    df_traps = pd.read_pickle(fn)
    rivers_surf_DO_unfiltered = df_traps['surface trapsDO [kmol/s]'].values
    wwtps_deep_DO_unfiltered = df_traps['deep trapsDO [kmol/s]'].values
    rivers_surf_flow_unfiltered = df_traps['surface flowtraps [m3/s]'].values
    wwtps_deep_flow_unfiltered = df_traps['deep flowtraps [m3/s]'].values

    # Pad with zeros if no rivers or WWTPs
    if rivers_surf_DO_unfiltered.size == 0 or np.isnan(rivers_surf_DO_unfiltered[0]):
            rivers_surf_DO_unfiltered = np.zeros(8761)
    if wwtps_deep_DO_unfiltered.size == 0:
            wwtps_deep_DO_unfiltered = np.zeros(8761)
    if rivers_surf_flow_unfiltered.size == 0 or np.isnan(rivers_surf_flow_unfiltered[0]):
            rivers_surf_flow_unfiltered = np.zeros(8761)
    if wwtps_deep_flow_unfiltered.size == 0:
            wwtps_deep_flow_unfiltered = np.zeros(8761)
    # Godin filter
    traps_surf_DO = zfun.lowpass(rivers_surf_DO_unfiltered, f='godin')[36:-34:24]
    traps_deep_DO =  zfun.lowpass(wwtps_deep_DO_unfiltered, f='godin')[36:-34:24]
    traps_surf_flow = zfun.lowpass(rivers_surf_flow_unfiltered, f='godin')[36:-34:24]
    traps_deep_flow =  zfun.lowpass(wwtps_deep_flow_unfiltered, f='godin')[36:-34:24]
    traps_color = 'teal'

    # plot budget time series
    ax[0].plot(dates_local_daily,zfun.lowpass(traps_surf_flow,n=10),color=traps_color,label='TRAPS')
    ax[2].plot(dates_local_daily,zfun.lowpass(traps_deep_flow,n=10),color=traps_color,label='TRAPS')
    ax[1].plot(dates_local_daily,zfun.lowpass(traps_surf_DO,n=10),color=traps_color,label='TRAPS')
    ax[3].plot(dates_local_daily,zfun.lowpass(traps_deep_DO,n=10),color=traps_color,label='TRAPS')
    # ax[3].plot(dates_local_daily,traps_deep_DO,color=traps_color,label='TRAPS')

# ------------------------------- get vertical exchange ----------------------------------------

    # pick color
    vertX_color = '#99C5F7'

    # calculate residual term, and ascribe to the vertical exchange
    vertX_surf_DO_TEF = ddtDOV_surf - (TEF_surf 
                                + photo_surf
                                + cons_surf
                                + airsea_surf
                                + traps_surf_DO)
    
    vertX_deep_DO_TEF = ddtDOV_deep - (TEF_deep 
                                + photo_deep
                                + cons_deep
                                + traps_deep_DO)
    
    # calculate residual term, and ascribe to the vertical exchange
    vertX_surf_flow_TEF = ddtvol_surf - (Q_m.values + traps_surf_flow)
    
    vertX_deep_flow_TEF = ddtvol_deep - (Q_p.values + traps_deep_flow)
    
    # plot budget time series
    ax[0].plot(dates_local_daily,zfun.lowpass(vertX_surf_flow_TEF,n=10),color=vertX_color,label='Vertical')
    ax[2].plot(dates_local_daily,zfun.lowpass(vertX_deep_flow_TEF,n=10),color=vertX_color,label='Vertical')
    ax[1].plot(dates_local_daily,zfun.lowpass(vertX_surf_DO_TEF,n=10),color=vertX_color,label='Vertical')
    ax[3].plot(dates_local_daily,zfun.lowpass(vertX_deep_DO_TEF,n=10),color=vertX_color,label='Vertical')
    # ax[3].plot(dates_local_daily,vertX_deep_DO_TEF,color=vertX_color,label='Vertical')

# ------------------------------- get budget error ----------------------------------------

    # pick color
    error_color = 'darkorange'

    # calculate error
    error_DO = vertX_surf_DO_TEF + vertX_deep_DO_TEF
    error_flow = vertX_surf_flow_TEF + vertX_deep_flow_TEF
    
    # plot budget time series
    ax[0].plot(dates_local_daily,zfun.lowpass(error_flow,n=10),color=error_color,label='Error')
    ax[2].plot(dates_local_daily,zfun.lowpass(error_flow,n=10),color=error_color,label='Error')
    ax[1].plot(dates_local_daily,zfun.lowpass(error_DO,n=10),color=error_color,label='Error')
    ax[3].plot(dates_local_daily,zfun.lowpass(error_DO,n=10),color=error_color,label='Error')
    
    ax[0].legend(loc='lower right',ncol=5, fontsize=10, handletextpad=0.15)
    ax[1].legend(loc='lower right',ncol=4, fontsize=10, handletextpad=0.15)

    # # add wind
    # ax2 = ax[3].twinx()
    # ax2.plot(dates_local_daily,df_daily_aligned['WSPD'],color='violet',
    #         linewidth=2,label='Daily average wind')
    # ax2.set_ylim([0,10])


    # format figure
    ax[0].text(0.02,0.92,'Surface Volume',fontsize=12,fontweight='bold',ha='left',transform=ax[0].transAxes)
    ax[1].text(0.02,0.92,'Surface DO',fontsize=12,fontweight='bold',ha='left',transform=ax[1].transAxes)
    ax[2].text(0.02,0.92,'Deep Volume',fontsize=12,fontweight='bold',ha='left',transform=ax[2].transAxes)
    ax[3].text(0.02,0.92,'Deep DO',fontsize=12,fontweight='bold',ha='left',transform=ax[3].transAxes)
    plt.tight_layout()
    plt.show()


# ------------------------------- save budget terms to df ----------------------------------------

    # Save Lynch Cove budget terms
    if station == 'lynchcove':
          lynchcove_dict_10dayhanning['d/dt(DO)'] = zfun.lowpass(ddtDOV_deep,n=10)
          lynchcove_dict_10dayhanning['Error'] = zfun.lowpass(error_DO,n=10)
          lynchcove_dict_10dayhanning['Exchange Flow'] = zfun.lowpass(TEF_deep,n=10)
          lynchcove_dict_10dayhanning['Vertical Transport'] = zfun.lowpass(vertX_deep_DO_TEF,n=10)
          lynchcove_dict_10dayhanning['Photosynthesis'] = zfun.lowpass(photo_deep,n=10)
          lynchcove_dict_10dayhanning['Consumption'] = zfun.lowpass(cons_deep,n=10)
        #   print(lynchcove_dict_10dayhanning['Error'])

     # get inlet name
    if station == 'case':
            inlet_name = 'Case Inlet'
    elif station == 'portsusan':
            inlet_name = 'Port Susan'
    elif station == 'penn':
            inlet_name = 'Penn Cove'
    elif station == 'holmes':
            inlet_name = 'Holmes Harbor'
    elif station == 'dabob':
            inlet_name = 'Dabob Bay'
    elif station == 'lynchcove':
            inlet_name = 'Lynch Cove'
    elif station == 'crescent':
            inlet_name = 'Crescent Harbor'
    elif station == 'carr':
            inlet_name = 'Carr Inlet'
    elif station == 'dyes':
            inlet_name = 'Dyes Inlet'
    elif station == 'sinclair':
            inlet_name = 'Sinclair Inlet'
    elif station == 'elliot':
            inlet_name = 'Elliott Bay'
    elif station == 'commencement':
            inlet_name = 'Commencement Bay'
    elif station == 'quartermaster':
            inlet_name = 'Quartermaster Harbor'

    # get budget terms during decline period (volume normalize and convert to mg/L/day)
    conversion = 1000 * 32 * 60 * 60 * 24 
    minday = 194
    maxday = 256

    deep_exchange_nonnormalized = np.nanmean(TEF_deep[minday:maxday])
    total_vol = np.nanmean(vol_total[minday:maxday])

    deep_exchange_all = TEF_deep[minday:maxday]/vol_deep[minday:maxday]
    deep_exchange_avg = np.nanmean(deep_exchange_all) * conversion
    deep_exchange_err = np.nanstd(deep_exchange_all) * conversion

    deep_vert_transport_all = vertX_deep_DO_TEF[minday:maxday]/vol_deep[minday:maxday]
    deep_vert_transport_avg = np.nanmean(deep_vert_transport_all) * conversion
    deep_vert_transport_err = np.nanstd(deep_vert_transport_all) * conversion

    deep_photosynthesis_all = photo_deep[minday:maxday]/vol_deep[minday:maxday]
    deep_photosynthesis_avg = np.nanmean(deep_photosynthesis_all) * conversion
    deep_photosynthesis_err = np.nanstd(deep_photosynthesis_all) * conversion

    deep_consumption_all = cons_deep[minday:maxday]/vol_deep[minday:maxday]
    deep_consumption_avg = np.nanmean(deep_consumption_all) * conversion
    deep_consumption_err = np.nanstd(deep_consumption_all) * conversion

    deep_ddtDO_all = ddtDOV_deep[minday:maxday]/vol_deep[minday:maxday]
    deep_ddtDO_avg = np.nanmean(deep_ddtDO_all) * conversion
    deep_ddtDO_err = np.nanstd(deep_ddtDO_all) * conversion

    deep_physresup_avg = np.nanmean(deep_exchange_all+deep_vert_transport_all) * conversion # Physical resupply = Exchange flow + Vertical transport
    deep_physresup_err = np.nanstd(deep_exchange_all+deep_vert_transport_all) * conversion

    deep_NEM_avg = np.nanmean(deep_photosynthesis_all+deep_consumption_all) * conversion # NEM = Photosynthesis + Consumption
    deep_NEM_err = np.nanstd(deep_photosynthesis_all+deep_consumption_all) * conversion

    # hypoxic season deep DO
    hypminday = 242
    hypmaxday = 302
    deep_DO_avg = np.nanmean(DO_deep[hypminday:hypmaxday])
    deep_DO_err = np.nanstd(DO_deep[hypminday:hypmaxday])

    # decline period DO concentrations
    DO_in = DO_p * 32/1000
    DOin_DOdeep = np.nanmean(DO_in[minday:maxday] - DO_deep[minday:maxday])
    DOin_DOdeep_err = np.nanstd(DO_in[minday:maxday] - DO_deep[minday:maxday])

    # calculate mean depth
    fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'vol_df_cas7_c21.p'
    vol_df = pd.read_pickle(fn)
    inlet_vol = vol_df['volume m3'].loc[station+'_p']
    inlet_area = vol_df['area m2'].loc[station+'_p']
    mean_depth = inlet_vol / inlet_area

    # add data to df
    new_data = {'Inlet': [inlet_name],
                'ExchangeFlow': [deep_exchange_avg],
                'ExchangeFlow_err': [deep_exchange_err],
                'VerticalTransport': [deep_vert_transport_avg],
                'VerticalTransport_err': [deep_vert_transport_err],
                'Photosynthesis': [deep_photosynthesis_avg],
                'Photosynthesis_err': [deep_photosynthesis_err],
                'Consumption': [deep_consumption_avg],
                'Consumption_err': [deep_consumption_err],
                'd/dt(DO)': [deep_ddtDO_avg],
                'd/dt(DO)_err': [deep_ddtDO_err],
                'PhysicalResupply': [deep_physresup_avg],
                'PhysicalResupply_err': [deep_physresup_err],
                'NetEcosystemMetabolism': [deep_NEM_avg],
                'NetEcosystemMetabolism_err': [deep_NEM_err],
                'SepOctDeepDO[mg/L]': [deep_DO_avg],
                'SepOctDeepDO_err[mg/L]': [deep_DO_err],
                'MeanDepth[m]': [mean_depth],
                'DOin-DOdeep[mg/L]': [DOin_DOdeep],
                'DOin-DOdeep_err[mg/L]': [DOin_DOdeep_err],
                'Inflow no normalization [kmol/s]': [deep_exchange_nonnormalized],
                'Inlet volume [m3]': [total_vol]}
    df_new_rows = pd.DataFrame(new_data)
    inlet_budget_df = pd.concat([inlet_budget_df, df_new_rows],ignore_index=True)


    # # Error analysis
    # print('===========================\n'+station)
    # QinDOin = (TEF_deep/inlet_vol) * (1000 * 32 * 60 * 60 * 24)
    # error = (error_DO/inlet_vol) * (1000 * 32 * 60 * 60 * 24)
    # print('     Error (ann avg.) = ' + str(round(np.nanmean(np.abs(error)),3)))
    # print('     QinDOin (ann avg.) = ' + str(round(np.nanmean(np.abs(QinDOin)),3)))
    # print('     perc of QinDOin (ann avg.) = ' + str(round((np.nanmean(np.abs(error))/np.nanmean(np.abs(QinDOin)))*100,2)) + '%')
    # initialize lists

    # calculate budget error (mg/L per day) ------------------------------
    conversion = (1000 * 32 * 60 * 60 * 24)
    error_TEF = (error_DO/inlet_vol) * conversion # [mg/L/day]
    inlet_error_ann_avg = np.nanmean(error_TEF)
    # calculate QinDOin (mg/L per day) 
    QinDOin = (TEF_deep/inlet_vol) * conversion # [mg/L/day]
    inlet_QinDOin_ann_avg = np.nanmean(QinDOin)
    # calculate biological consumption in deep layer (mg/L per day)
    consumption = (cons_deep/inlet_vol) * conversion # [mg/L/day]
    consumption_1lay = ((cons_deep+cons_surf)/inlet_vol) * conversion # [mg/L/day]
    inlet_consumption_ann_avg = np.nanmean(consumption)
    # calculate d/dt(DO) (mg/L per day)
    ddtDO = (ddtDOV_deep/inlet_vol) * conversion # [mg/L/day]
    inlet_ddtDO_ann_avg = np.nanmean(ddtDO)
    # print(inlet_ddtDO_ann_avg)
    ddtDO_onelayer = (ddtDOV_onelayer/inlet_vol) * conversion # [mg/L/day]
    inlet_ddtDO_onelayer_ann_avg = np.nanmean(ddtDO_onelayer)

    # square_error = np.sum(x**2 for x in error_TEF[:-1])
    # square_QinDOin = np.sum(x**2 for x in QinDOin[:-1])
    # square_consumption = np.sum(x**2 for x in consumption[:-1])
    # rms_error = np.sqrt(square_error/(len(error_TEF)-1))
    # rms_QinDOin = np.sqrt(square_QinDOin/(len(QinDOin)-1))
    # rms_consumption = np.sqrt(square_consumption/(len(consumption)-1))

    # error_QinDOin_ann_avg.append(rms_error/rms_QinDOin)
    # error_consumption_ann_avg.append(rms_error/rms_consumption)

    error_kmolO2s_ann_avg.append(np.nanmean(error_DO))

    # # calculating annual average before diving
    # error_QinDOin_ann_avg.append(inlet_error_ann_avg/inlet_QinDOin_ann_avg)
    # error_consumption_ann_avg.append(inlet_error_ann_avg/inlet_consumption_ann_avg)
    # error_ddtDO_ann_avg.append(inlet_error_ann_avg/inlet_ddtDO_ann_avg)
    # error_ddtDO_onelayer_ann_avg.append(inlet_error_ann_avg/inlet_ddtDO_onelayer_ann_avg)

    # calculating division before annual averaging
    # 0 90 181 272 363 # decline period: 194/256
    err_minday = 0
    err_maxday = 363

    error_QinDOin_ann_avg.append(np.abs(np.nanmean(error_DO[err_minday:err_maxday])/np.nanmean(TEF_deep[err_minday:err_maxday])))
    error_consumption_ann_avg.append(np.abs(np.nanmean(error_DO[err_minday:err_maxday])/np.nanmean(cons_deep[err_minday:err_maxday])))
    error_Qinvol_ann_avg.append(np.abs(np.nanmean(error_flow[err_minday:err_maxday])/np.nanmean(Q_p[err_minday:err_maxday])))

    # error_ddtDO_ann_avg.append(np.nanmean(error_TEF/ddtDO))
    # error_ddtDO_onelayer_ann_avg.append(np.nanmean(error_TEF/ddtDO_onelayer))

    # print('----------------')
    # print(station)
    # print(np.abs(np.nanmean(error_TEF[err_minday:err_maxday])/np.nanmean(QinDOin[err_minday:err_maxday])))
    # print(rms_error)
    # print(np.abs(np.nanmean(error_TEF[err_minday:err_maxday])/np.nanmean(consumption[err_minday:err_maxday]))*100)
    decline_per_vol_norm_error = np.nanmean((error_DO[err_minday:err_maxday]/vol_deep[err_minday:err_maxday]))* conversion
    error_over_ddtDO.append(np.abs(decline_per_vol_norm_error/deep_ddtDO_avg))
    # print(np.round(decline_per_vol_norm_error/deep_ddtDO_avg*100,2))
    # print(inlet_error_ann_avg/inlet_consumption_ann_avg)
    # print(np.abs(np.nanmean(error_TEF[err_minday:err_maxday])/np.nanmean(consumption[err_minday:err_maxday])))
    # print(np.nanmean(np.abs(error_TEF[err_minday:err_maxday]/consumption_1lay[err_minday:err_maxday]))*100)
    # print(round(np.abs(np.nanmean(error_TEF/consumption))*100,1))
    # print(round((inlet_error_ann_avg/inlet_consumption_ann_avg)*100,1))

# calculate bulk statistics
error_QinDOin = np.abs(np.nanmean(error_QinDOin_ann_avg)) * 100
error_consumption = np.abs(np.nanmean(error_consumption_ann_avg)) * 100
error_ddtDO = np.abs(np.nanmean(error_ddtDO_ann_avg)) * 100
error_ddtDO_onelayer = np.abs(np.nanmean(error_ddtDO_onelayer_ann_avg)) * 100
error_Qinvol = np.abs(np.nanmean(error_Qinvol_ann_avg)) * 100
#     # add values to list
#     error_QinDOin_ann_avg.append(np.nanmean(np.abs(error_TEF)/np.abs(QinDOin)))
#     error_consumption_ann_avg.append(np.nanmean(error_TEF/consumption))
#     error_ddtDO_ann_avg.append(np.nanmean(error_TEF/ddtDO))
#     error_ddtDO_onelayer_ann_avg.append(np.nanmean(np.abs(error_TEF)/np.abs(ddtDO_onelayer)))
#     print(ddtDO_onelayer[100:105])
#     print(np.abs(ddtDO_onelayer)[100:105])
# # calculate bulk statistics
# error_QinDOin = np.abs(np.nanmean(error_QinDOin_ann_avg)) * 100
# error_consumption = np.abs(np.nanmean(error_consumption_ann_avg)) * 100
# error_ddtDO = np.abs(np.nanmean(error_ddtDO_ann_avg)) * 100
# error_ddtDO_onelayer = np.abs(np.nanmean(error_ddtDO_onelayer_ann_avg)) * 100

# print bulk statistics
print('(annual mean error)/(annual mean QinDOin) [expressed as percentage]')
print('    {}%'.format(round(error_QinDOin,2)))
print('\n')
print('(annual mean error)/(annual mean deep consumption) [expressed as percentage]')
print('    {}%'.format(round(error_consumption,4)))
print('\n')
print('decline period err/ddt(DO) % [which had been volume normalized]')
print('    {}%'.format(round(np.nanmean(error_over_ddtDO)*100,2)))
# print('\n')
# print('(annual mean error)/(annual mean deep d/dt DO) [expressed as percentage]')
# print('    {}%'.format(round(error_ddtDO,2)))

print('VOLUME BUDGET (annual mean error)/(annual mean Qin) [expressed as percentage]')
print('    {}%'.format(round(error_Qinvol,2)))
print('\n')

print('\n')
print(np.nanstd(error_kmolO2s_ann_avg))



# # exchange flow vs. inlet volume
# fig, ax = plt.subplots(1,1,figsize=(6,6))
# ax.scatter(inlet_budget_df['MeanDepth[m]'],inlet_budget_df['Inflow no normalization [kmol/s]'],
#                     s=50, edgecolor='white',linewidth=0.5,zorder=5)
# for inlet in inlet_budget_df['Inlet']:
#       ax.text(inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'MeanDepth[m]'],
#               inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'Inflow no normalization [kmol/s]']+0.005,
#                 inlet)
# plt.show()

# get average of all inlets and add data to df
# error of averages is error of all values added in quadrature, divided by number of items
new_data = {'Inlet': ['All inlets'],
            'ExchangeFlow': [np.nanmean(inlet_budget_df['ExchangeFlow'])],
            'ExchangeFlow_err': [np.nanstd(inlet_budget_df['ExchangeFlow'])],
            'VerticalTransport': [np.nanmean(inlet_budget_df['VerticalTransport'])],
            'VerticalTransport_err': [np.nanstd(inlet_budget_df['VerticalTransport'])],
            'Photosynthesis': [np.nanmean(inlet_budget_df['Photosynthesis'])],
            'Photosynthesis_err': [np.nanstd(inlet_budget_df['Photosynthesis'])],
            'Consumption': [np.nanmean(inlet_budget_df['Consumption'])],
            'Consumption_err': [np.nanstd(inlet_budget_df['Consumption'])],
            'd/dt(DO)': [np.nanmean(inlet_budget_df['d/dt(DO)'])],
            'd/dt(DO)_err': [np.nanstd(inlet_budget_df['d/dt(DO)'])],
            'PhysicalResupply': [np.nanmean(inlet_budget_df['PhysicalResupply'])],
            'PhysicalResupply_err': [np.nanstd(inlet_budget_df['PhysicalResupply'])],
            'NetEcosystemMetabolism': [np.nanmean(inlet_budget_df['NetEcosystemMetabolism'])],
            'NetEcosystemMetabolism_err': [np.nanstd(inlet_budget_df['NetEcosystemMetabolism'])],
            'SepOctDeepDO[mg/L]': [np.nanmean(inlet_budget_df['SepOctDeepDO[mg/L]'])],
            'SepOctDeepDO_err[mg/L]': [np.nanstd(inlet_budget_df['SepOctDeepDO[mg/L]'])],
            'MeanDepth[m]': [np.nanmean(inlet_budget_df['MeanDepth[m]'])],
            'DOin-DOdeep[mg/L]': [np.nanmean(inlet_budget_df['DOin-DOdeep[mg/L]'])],
            'DOin-DOdeep_err[mg/L]': [np.nanstd(inlet_budget_df['DOin-DOdeep[mg/L]'])],}
df_new_rows = pd.DataFrame(new_data)
inlet_budget_df = pd.concat([inlet_budget_df, df_new_rows],ignore_index=True)

# # save to csv file
# print(inlet_budget_df)
# inlet_budget_df.to_csv('../../../../terminal_inlet_DO_rev2/inlet_budgets_decline_period_mgL_day_' + interface_types[0] + '.csv', index=False)
# inlet_budget_df.to_csv('inlet_budgets_decline_period_mgL_day_' + interface_types[0] + '.csv', index=False)

# print(inlet_budget_df[['Inlet', 'ExchangeFlow', 'ExchangeFlow_err']])
# print(inlet_budget_df[['Inlet', 'VerticalTransport', 'VerticalTransport_err']])
# print(inlet_budget_df[['Inlet', 'Photosynthesis', 'Photosynthesis_err']])
# print(inlet_budget_df[['Inlet', 'Consumption', 'Consumption_err']])
# print(inlet_budget_df[['Inlet', 'd/dt(DO)', 'd/dt(DO)_err']])
# print(inlet_budget_df[['Inlet', 'PhysicalResupply', 'PhysicalResupply_err']])
# print(inlet_budget_df[['Inlet', 'NetEcosystemMetabolism', 'NetEcosystemMetabolism_err']])
# print(inlet_budget_df[['Inlet', 'SepOctDeepDO[mg/L]', 'SepOctDeepDO_err[mg/L]']])
# print(inlet_budget_df[['Inlet', 'MeanDepth[m]']])
# print(inlet_budget_df[['Inlet', 'DOin-DOdeep[mg/L]', 'DOin-DOdeep_err[mg/L]']])


# # save lynch cove budget to csv file
# lynchcove_budget_df = pd.DataFrame.from_dict(lynchcove_dict_10dayhanning)
# dates = pd.date_range(start='2017-01-02', end='2017-12-30', freq='D')
# date_list = dates.strftime('%Y-%m-%d').tolist()
# lynchcove_budget_df.insert(0, 'date', date_list)
# lynchcove_budget_df.to_csv('../../../../terminal_inlet_DO_rev2/lynchcove_2017_budget_kmolO2_s_10dayHanning.csv', index=False)
# print(lynchcove_budget_df)

# # ----------------------- test new scatter plot (budget terms vs DOdeep) -------------------------------

# # initialize figure
# fig = plt.figure(figsize=(9.5, 7))
# gs = GridSpec(2, 12, figure=fig, height_ratios=[3, 4], wspace=2.8)
# fig.subplots_adjust(top=0.80, right=0.99, left=0.12)  # leave room for colorbar
# axes = []
# for i in range(4):
#     if i == 0:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3])
#     else:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3], sharex=axes[0])
#     axes.append(ax)
# for i in range(3):
#     ax = fig.add_subplot(gs[1, i*4:(i+1)*4], sharex=axes[0])
#     axes.append(ax)
    
# # plot scatter points
# vars = ['ExchangeFlow','VerticalTransport',
#        'Photosynthesis','Consumption', 'd/dt(DO)',
#        'PhysicalResupply', 'NetEcosystemMetabolism']
# colors = ['#0D4B91','#99C5F7','#8F0445','#FCC2DD','black','#488DDB','#F069A8']
# letters = ['(a) Exchange Flow','(b) Vertical','(c) Photosynthesis','(d) Consumption',
#            '(e) d/dt(DO)','(f) Physical Resupply', '(g) Net Ecosystem\nMetabolism']
# ylims = [ [-1,4], [-4.5,1], [-0.1,0.4], [-0.15,0.05],
#          [-0.25,0.15], [-0.4,0.1], [-0.06,0.25] ]

# for i,var in enumerate(vars):
#     # add error bars
#     axes[i].errorbar(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df[var][:-1],
#                     xerr=inlet_budget_df['SepOctDeepDO_err[mg/L]'][:-1],
#                     yerr=inlet_budget_df[var+'_err'][:-1],
#                     fmt='o',color='black')
#     # plot points
#     inlet_crosssection_df = pd.read_csv('inlet_crosssections.csv')
#     cmap_temp = colormaps['gnuplot2_r'].resampled(256)
#     cmap_depth = ListedColormap(cmap_temp(np.linspace(0.1, 0.8, 256)))# get range of colormap
#     cs = axes[i].scatter(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df[var][:-1],
#                     s=50, zorder=5,c=inlet_budget_df['MeanDepth[m]'][:-1], cmap=cmap_depth, vmin=0, vmax=110)
    
#     # create colorbarlegend
#     if i == 0 and t == 0:
#         cbar_ax = fig.add_axes([0.1, 0.92, 0.8, 0.03])
#         cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal')
#         cbar.ax.tick_params(labelsize=12)
#         cbar.set_label(r'Mean depth [m]', fontsize=12)
#         cbar.outline.set_visible(False)

#     # add zero line
#     axes[i].axhline(0,0,8, color='gray',linestyle=':')
#     # format panel
#     axes[i].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
#                  transform=axes[i].transAxes, zorder=6, va='top')
#     axes[i].set_xlim([0,8])
#     axes[i].set_ylim(ylims[i])

#     # # add linear fit
#     # input_array = np.array([inlet_budget_df['SepOctDeepDO[mg/L]'][:-1], [1]*len(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1])]).T
#     # B,a,b,c = lstsq(input_array,inlet_budget_df[var][:-1])
#     # slope = B[0]
#     # intercept = B[1]
#     # x=[0,8]
#     # y=[intercept+slope*xval for xval in x]
#     # axes[i].plot(x,y,color='grey')
#     # calculate r^2 and p value
#     r,p = pearsonr(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df[var][:-1])
#     axes[i].text(0.03,0.16,'R = {}\np = {}'.format(round(r,2),round(p,2)),
#                  transform=axes[i].transAxes, zorder=6, va='top')
#     # print('\n===============================================')
#     # print('{}: d/dt(DO) dependence on DO_deep'.format(type))
#     # print('   r = {}'.format(r))
#     # print('   R^2 = {}'.format(r**2))
#     # print('   p = {}'.format(p))
#     # print('===============================================\n')

#     axes[i].tick_params(axis='both', labelsize=12)

#     if i in [0,4]:
#           axes[i].set_ylabel('Decline period rates\n' + r'[mg L$^{-1}$ d$^{-1}$]',fontsize=12)
#     if i >=4 :
#           axes[i].set_xlabel(r'Sep-Oct DO$_{deep}$ [mg L$^{-1}$]',fontsize=12)

# # plt.subplots_adjust(top=0.88)
# # plt.tight_layout()
# plt.show()


# ----------------------- test new scatter plot (budget terms vs DOdeep with inlet labels) -------------------------------

# initialize figure
fig, axes = plt.subplots(2,4,figsize=(11, 7), sharex=True)
fig.subplots_adjust(top=0.80, right=0.99, left=0.08, wspace=0.3)  # leave room for colorbar
axes = axes.ravel()
    
# plot scatter points
vars = ['ExchangeFlow','VerticalTransport',
       'Photosynthesis','Consumption',
        'd/dt(DO)','PhysicalResupply', 'NetEcosystemMetabolism', 'Inlet']
# colors = ['#0D4B91','#99C5F7','#8F0445','#FCC2DD','black','#488DDB','#F069A8']
letters = ['(a) Exchange Flow','(b) Vertical','(c) Photosynthesis','(d) Consumption',
           '(e) d/dt(DO)','(f) Physical Resupply', '(g) Net Ecosystem\nMetabolism','(h)']
ylims = [[-1,4], [-4.5,1], [-0.1,0.4], [-0.15,0.05],
         [-0.25,0.15], [-0.4,0.12], [-0.06,0.25],  [-0.5,12.5]]

for i,var in enumerate(vars):

    # define colormap based on inlet mean depth
    cmap_temp = colormaps['gnuplot2_r'].resampled(256)
    cmap_depth = ListedColormap(cmap_temp(np.linspace(0.1, 0.8, 256)))# get range of colormap

    # first panel lists inlets from lowest to highest DO
    if i == 7:
        inlets_only_df = inlet_budget_df.iloc[:-1]
        df_sorted = inlets_only_df.sort_values(by=['SepOctDeepDO[mg/L]'])
        axes[i].errorbar(df_sorted['SepOctDeepDO[mg/L]'],df_sorted['Inlet'],
                    xerr=df_sorted['SepOctDeepDO_err[mg/L]'],
                    fmt='o',color='black')
        axes[i].scatter(df_sorted['SepOctDeepDO[mg/L]'], df_sorted['Inlet'],
                        s=50, zorder=5,c=df_sorted['MeanDepth[m]'], cmap=cmap_depth, vmin=0, vmax=110) 
        # add inlet labels
        inlets = ['Lynch Cove','Dabob','Holmes','Port Susan',
                  'Commencement','Elliott Bay', 'Carr Inlet',
                  'Penn Cove','Crescent','Case Inlet',
                  'Quartermster','Sinclair','Dyes Inlet']
        for idx in range(13):
            inlet = inlets[idx] #df_sorted['Inlet'].values[idx]
            DOdeep = df_sorted['SepOctDeepDO[mg/L]'].values[idx]
            if idx < 4:
                axes[i].text(DOdeep+0.4, idx, inlet, fontsize=10, va='center', ha = 'left', color='gray') 
            else:
                axes[i].text(DOdeep-0.4, idx, inlet, fontsize=10, va='center', ha = 'right', color='gray') 
        minor_locator = AutoMinorLocator(2)
        axes[i].yaxis.set_minor_locator(minor_locator)
        axes[i].grid(which='minor', color='gray', linestyle=':')
        axes[i].set_yticklabels([])
        # continue
    else:
        # add error bars
        axes[i].errorbar(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df[var][:-1],
                        xerr=inlet_budget_df['SepOctDeepDO_err[mg/L]'][:-1],
                        yerr=inlet_budget_df[var+'_err'][:-1],
                        fmt='o',color='black')
        # plot points
        cs = axes[i].scatter(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df[var][:-1],
                        s=50, zorder=5,c=inlet_budget_df['MeanDepth[m]'][:-1], cmap=cmap_depth, vmin=0, vmax=110)
        
        # calculate r^2 and p value
        r,p = pearsonr(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df[var][:-1])
        axes[i].text(0.03,0.16,'R = {}\np = {}'.format(round(r,2),round(p,2)),
                    transform=axes[i].transAxes, zorder=6, va='top', fontsize=12, color='gray')
        
        # add zero line
        axes[i].axhline(0,0,8, color='gray',linestyle=':')
    
    # create colorbarlegend
    if i == 0:
        cbar_ax = fig.add_axes([0.1, 0.92, 0.87, 0.03])
        cbar = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=14)
        cbar.set_label(r'Mean depth [m]', fontsize=14)
        cbar.outline.set_visible(False)

    
    # format panel
    axes[i].text(0.03,0.96,letters[i],fontsize=12,fontweight='bold',
                 transform=axes[i].transAxes, zorder=6, va='top')
    axes[i].set_xlim([0,8])
    axes[i].set_ylim(ylims[i])


    axes[i].tick_params(axis='both', labelsize=12)
    axes[i].tick_params(axis='y', labelsize=12, rotation=45)

    if i in [0,4]:
          axes[i].set_ylabel('Rates ' + r'[mg L$^{-1}$ d$^{-1}$]',fontsize=12)
    if i >=4 :
          axes[i].set_xlabel(r'DO$_{deep}$ [mg L$^{-1}$]',fontsize=12)

# plt.subplots_adjust(top=0.88)
# plt.tight_layout()
plt.show()






















































# # ----------------------- test new scatter plot (budget terms vs mean depth) -------------------------------

# # initialize figure
# fig = plt.figure(figsize=(10, 5.5))
# gs = GridSpec(2, 12, figure=fig, height_ratios=[3, 4])#, wspace=2.4)
# axes = []  # list to store axes
# for i in range(4):
#     if i == 0:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3])
#     else:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3], sharex=axes[0])
#     axes.append(ax)
# for i in range(3):
#     ax = fig.add_subplot(gs[1, i*4:(i+1)*4], sharex=axes[0])
#     axes.append(ax)
# axes = [ax for ax in axes]

# # plot scatter points
# vars = ['ExchangeFlow','VerticalTransport',
#        'Photosynthesis','Consumption', 'd/dt(DO)',
#        'PhysicalResupply', 'NetEcosystemMetabolism']
# colors = ['#0D4B91','#99C5F7','#8F0445','#FCC2DD','black','#488DDB','#F069A8']
# letters = ['(a) Exchange Flow','(b) Vertical','(c) Photosynthesis','(d) Consumption',
#            '(e) d/dt(DO)','(f) Physical Resupply', '(g) Net Ecosystem\nMetabolism']
# ylims = [ [-0.5,4], [-4.5,1], [-0.05,0.4], [-0.15,0.05],
#          [-0.2,0.1], [-0.5,0.2], [-0.06,0.25] ]

# for i,var in enumerate(vars):
#     # add error bars
#     axes[i].errorbar(inlet_budget_df['MeanDepth[m]'][:-1],inlet_budget_df[var][:-1],
#                     yerr=inlet_budget_df[var+'_err'][:-1],
#                     fmt='o',color='black')
#     # plot points
#     axes[i].scatter(inlet_budget_df['MeanDepth[m]'][:-1],inlet_budget_df[var][:-1], color=colors[i],
#                     s=50, edgecolor='white',linewidth=0.5,zorder=5)
#     # add mean line
#     axes[i].axhline(inlet_budget_df[var].values[-1],0,110, color=colors[i])
#     minval = inlet_budget_df[var].values[-1] - inlet_budget_df[var+'_err'].values[-1]
#     maxval = inlet_budget_df[var].values[-1] + inlet_budget_df[var+'_err'].values[-1]
#     axes[i].fill_between([0,110], [minval,minval], [maxval,maxval],
#                          color=colors[i],alpha=0.3)
#     # add zero line
#     axes[i].axhline(0,0,110, color='gray',linestyle=':')
#     # format panel
#     axes[i].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
#                  transform=axes[i].transAxes, zorder=6, va='top')
#     axes[i].set_xlim([0,110])
#     axes[i].set_ylim(ylims[i])

#     axes[i].tick_params(axis='both', labelsize=12)

#     if i in [0,4]:
#           axes[i].set_ylabel('Decline period rates\n' + r'[mg L$^{-1}$ d$^{-1}$]',fontsize=12)
#     if i >=4 :
#           axes[i].set_xlabel('Inlet mean depth [m]',fontsize=12)

# plt.tight_layout()
# plt.show()

# # ----------------------- just d/dt(DO) -------------------------------

# # initialize figure
# fig, ax = plt.subplots(1,1,figsize=(5,4))

# # plot scatter points
# var =  'd/dt(DO)'

# # add error bars
# ax.errorbar(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df[var][:-1],
#                 xerr=inlet_budget_df['SepOctDeepDO_err[mg/L]'][:-1],
#                 yerr=inlet_budget_df[var+'_err'][:-1],
#                 fmt='o',color='black',elinewidth=0.5)
# # plot points
# ax.scatter(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df[var][:-1], color='black',
#                 s=50, edgecolor='black',linewidth=0.5,zorder=5)
# # add mean line
# ax.axhline(inlet_budget_df[var].values[-1],0,8, color='teal')
# minval = inlet_budget_df[var].values[-1] - inlet_budget_df[var+'_err'].values[-1]
# maxval = inlet_budget_df[var].values[-1] + inlet_budget_df[var+'_err'].values[-1]
# ax.fill_between([0,9], [minval,minval], [maxval,maxval],
#                         color='lightseagreen',alpha=0.3)
# # add zero line
# ax.axhline(0,0,8, color='gray',linestyle=':')

# ax.set_xlim([0,8])

# ax.tick_params(axis='both', labelsize=12)
# ax.set_ylabel(r'Decline period d/dt(DO) [mg L$^{-1}$ d$^{-1}$]',fontsize=12)
# ax.set_xlabel(r'Hypoxic season DO$_{deep}$ [mg L$^{-1}$]',fontsize=12)

# plt.tight_layout()
# plt.show()

# # ----------------------- test new scatter plot (budget terms vs DOin-DOdeep) -------------------------------

# # initialize figure
# fig = plt.figure(figsize=(10, 5.5))
# gs = GridSpec(2, 12, figure=fig, height_ratios=[3, 4])#, wspace=2.4)
# axes = []  # list to store axes
# for i in range(4):
#     if i == 0:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3])
#     else:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3], sharex=axes[0])
#     axes.append(ax)
# for i in range(3):
#     ax = fig.add_subplot(gs[1, i*4:(i+1)*4], sharex=axes[0])
#     axes.append(ax)
# axes = [ax for ax in axes]

# # plot scatter points
# vars = ['ExchangeFlow','VerticalTransport',
#        'Photosynthesis','Consumption', 'd/dt(DO)',
#        'PhysicalResupply', 'NetEcosystemMetabolism']
# colors = ['#0D4B91','#99C5F7','#8F0445','#FCC2DD','black','#488DDB','#F069A8']
# letters = ['(a) Exchange Flow','(b) Vertical','(c) Photosynthesis','(d) Consumption',
#            '(e) d/dt(DO)','(f) Physical Resupply', '(g) Net Ecosystem\nMetabolism']
# ylims = [ [-0.5,4], [-4.5,1], [-0.05,0.4], [-0.15,0.05],
#          [-0.2,0.1], [-0.5,0.2], [-0.06,0.25] ]

# for i,var in enumerate(vars):
#     # add error bars
#     axes[i].errorbar(inlet_budget_df['DOin-DOdeep[mg/L]'][:-1],inlet_budget_df[var][:-1],
#                     xerr=inlet_budget_df['DOin-DOdeep_err[mg/L]'][:-1],
#                     yerr=inlet_budget_df[var+'_err'][:-1],
#                     fmt='o',color='black')
#     # plot points
#     axes[i].scatter(inlet_budget_df['DOin-DOdeep[mg/L]'][:-1],inlet_budget_df[var][:-1], color=colors[i],
#                     s=50, edgecolor='white',linewidth=0.5,zorder=5)
#     # add mean line
#     axes[i].axhline(inlet_budget_df[var].values[-1],-1,3, color=colors[i])
#     minval = inlet_budget_df[var].values[-1] - inlet_budget_df[var+'_err'].values[-1]
#     maxval = inlet_budget_df[var].values[-1] + inlet_budget_df[var+'_err'].values[-1]
#     axes[i].fill_between([-1,3], [minval,minval], [maxval,maxval],
#                          color=colors[i],alpha=0.3)
#     # add zero line
#     axes[i].axhline(0,0,8, color='gray',linestyle=':')
#     axes[i].axvline(0,0,8, color='gray',linestyle=':')
#     # format panel
#     axes[i].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
#                  transform=axes[i].transAxes, zorder=6, va='top')
#     axes[i].set_xlim([-1,3])
#     axes[i].set_ylim(ylims[i])

#     axes[i].tick_params(axis='both', labelsize=12)

#     if i in [0,4]:
#           axes[i].set_ylabel('Decline period rates\n' + r'[mg L$^{-1}$ d$^{-1}$]',fontsize=12)
#     if i >=4 :
#           axes[i].set_xlabel(r'Decline period DO$_{in}$ - DO$_{deep}$ [mg L$^{-1}$]',fontsize=12)

# plt.tight_layout()
# plt.show()