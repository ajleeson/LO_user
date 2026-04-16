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
inlet_budget_df = pd.DataFrame(columns=['Inlet', 'QinDOin', 'QinDOin_err', 'QoutDOout',
       'QoutDOout_err', 'Photosynthesis', 'Photosynthesis_err',
       'Consumption', 'Consumption_err', 'd/dt(DO)', 'd/dt(DO)_err',
       'PhysicalResupply', 'PhysicalResupply_err', 'NetEcosystemMetabolism',
       'NetEcosystemMetabolism_err', 'SepOctInletDO[mg/L]', 'SepOctInletDO_err[mg/L]',
       'MeanDepth[m]'])
# initialize lynchcove budget dict
lynchcove_dict_10dayhanning = {}

# COLLAPSE
for i,station in enumerate(sta_dict):
        
    # initialize figure
    fig,ax = plt.subplots(1,1, figsize=(8,6))

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
    # ax.plot(dates_local_daily,zfun.lowpass(TEF_surf,n=10),color='#99C5F7',label='QoutDOout')
    # ax.plot(dates_local_daily,zfun.lowpass(TEF_deep,n=10),color='#0D4B91',label='QinDOin')
    ax.plot(dates_local_daily,zfun.lowpass(TEF_deep+TEF_surf,n=10),color='#0D4B91',label='Exchange Flow')

    # format budget time series figures
    ax.set_xlim([dates_hrly[0],dates_hrly[-2]])
    ax.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')
    ax.tick_params(axis='x', labelrotation=30, labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    loc = mdates.MonthLocator(interval=1)
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.set_ylabel(r'DO transport [kmol O$_2$ s$^{-1}$]', fontsize=12)

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

    # combine all months
    for month in months:
        fn = 'tef_' + station + '_' + month + '.p'
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

    # take time derivative of (DO*V) to get d/dt (DO*V)
    ddtDOV_surf_unfiltered = np.diff(o2vol_surf_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
    ddtDOV_deep_unfiltered = np.diff(o2vol_deep_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
    ddtDOV_total_unfiltered = np.diff((o2vol_deep_unfiltered+o2vol_surf_unfiltered)) * conv

    # get DO concentration [mg/L]
    DO_deep_unfiltered = o2vol_deep_unfiltered/vol_deep_unfiltered * 32/1000 # mg/L
    DO_total_unfiltered = (o2vol_deep_unfiltered+o2vol_surf_unfiltered) / (vol_deep_unfiltered+vol_surf_unfiltered) * 32/1000 # mg/L

    # apply Godin filter
    photo_surf  = zfun.lowpass(photo_surf_unfiltered, f='godin')[36:-34:24]
    photo_deep  = zfun.lowpass(photo_deep_unfiltered, f='godin')[36:-34:24]
    photo_total = zfun.lowpass((photo_surf_unfiltered+photo_deep_unfiltered), f='godin')[36:-34:24]
    cons_surf   = zfun.lowpass(cons_surf_unfiltered, f='godin')[36:-34:24]
    cons_deep   = zfun.lowpass(cons_deep_unfiltered, f='godin')[36:-34:24]
    cons_total   = zfun.lowpass((cons_surf_unfiltered+cons_deep_unfiltered), f='godin')[36:-34:24]
    airsea_surf = zfun.lowpass(airsea_surf_unfiltered, f='godin')[36:-34:24]
    ddtDOV_surf = zfun.lowpass(ddtDOV_surf_unfiltered, f='godin')[36:-34:24]
    ddtDOV_deep = zfun.lowpass(ddtDOV_deep_unfiltered, f='godin')[36:-34:24]
    ddtDOV_total = zfun.lowpass(ddtDOV_total_unfiltered, f='godin')[36:-34:24]
    vol_surf = zfun.lowpass(vol_surf_unfiltered, f='godin')[36:-34:24]
    vol_deep = zfun.lowpass(vol_deep_unfiltered, f='godin')[36:-34:24]
    vol_total = zfun.lowpass((vol_surf_unfiltered+vol_deep_unfiltered), f='godin')[36:-34:24]
    DO_deep = zfun.lowpass(DO_deep_unfiltered, f='godin')[36:-34:24]
    DO_inlet = zfun.lowpass(DO_total_unfiltered, f='godin')[36:-34:24]

    # plot budget time series
    ax.plot(dates_local_daily,zfun.lowpass(ddtDOV_total,n=10),color='black',label='d/dt(DO)')

    ax.plot(dates_local_daily,zfun.lowpass(photo_total,n=10),color='#8F0445')

    ax.plot(dates_local_daily,zfun.lowpass(cons_total,n=10),color='#FCC2DD')

    ax.plot(dates_local_daily,zfun.lowpass(airsea_surf,n=10),color='yellowgreen') 

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
    traps_total_DO = zfun.lowpass((rivers_surf_DO_unfiltered+wwtps_deep_DO_unfiltered), f='godin')[36:-34:24]
    traps_surf_flow = zfun.lowpass(rivers_surf_flow_unfiltered, f='godin')[36:-34:24]
    traps_deep_flow =  zfun.lowpass(wwtps_deep_flow_unfiltered, f='godin')[36:-34:24]
    traps_color = 'teal'

    # plot budget time series
    ax.plot(dates_local_daily,zfun.lowpass(traps_total_DO,n=10),color=traps_color,label='TRAPS')

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

# ------------------------------- get budget error ----------------------------------------

    # pick color
    error_color = 'darkorange'

    # calculate error
    error_DO = vertX_surf_DO_TEF + vertX_deep_DO_TEF
    
    # plot budget time series
    ax.plot(dates_local_daily,zfun.lowpass(error_DO,n=10),color=error_color,label='Error')
    
    ax.legend(loc='lower right',ncol=5, fontsize=10, handletextpad=0.15)


    # format figure
    ax.text(0.02,0.92,'One-layer DO budget',fontsize=12,fontweight='bold',ha='left',transform=ax.transAxes)
    plt.tight_layout()
    plt.show()


# # ------------------------------- save budget terms to df ----------------------------------------

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
    photosynthesis_all_nonnormalized = np.nanmean(photo_total[minday:maxday])

    deep_exchange_all = TEF_deep[minday:maxday]/vol_total[minday:maxday]
    deep_exchange_avg = np.nanmean(deep_exchange_all) * conversion
    deep_exchange_err = np.nanstd(deep_exchange_all) * conversion

    surf_exchange_transport_all = TEF_surf[minday:maxday]/vol_total[minday:maxday]
    surf_exchange_transport_avg = np.nanmean(surf_exchange_transport_all) * conversion
    surf_exchange_transport_err = np.nanstd(surf_exchange_transport_all) * conversion

    photosynthesis_all = photo_total[minday:maxday]/vol_total[minday:maxday]
    photosynthesis_avg = np.nanmean(photosynthesis_all) * conversion
    photosynthesis_err = np.nanstd(photosynthesis_all) * conversion

    consumption_all = cons_total[minday:maxday]/vol_total[minday:maxday]
    consumption_avg = np.nanmean(consumption_all) * conversion
    consumption_err = np.nanstd(consumption_all) * conversion

    ddtDO_all = ddtDOV_total[minday:maxday]/vol_total[minday:maxday]
    ddtDO_avg = np.nanmean(ddtDO_all) * conversion
    ddtDO_err = np.nanstd(ddtDO_all) * conversion

    physresup_avg = np.nanmean(deep_exchange_all+surf_exchange_transport_all) * conversion # Physical resupply = Exchange flow + Vertical transport
    physresup_err = np.nanstd(deep_exchange_all+surf_exchange_transport_all) * conversion

    NEM_avg = np.nanmean(photosynthesis_all+consumption_all) * conversion # NEM = Photosynthesis + Consumption
    NEM_err = np.nanstd(photosynthesis_all+consumption_all) * conversion

    # hypoxic season deep DO
    hypminday = 242
    hypmaxday = 302
    DOinlet_avg = np.nanmean(DO_inlet[hypminday:hypmaxday])
    DOinlet_err = np.nanstd(DO_inlet[hypminday:hypmaxday])

    # decline period DO concentrations
    DO_in = DO_p * 32/1000
    DOin_DOinlet = np.nanmean(DO_in[minday:maxday] - DO_inlet[minday:maxday])
    DOin_DOinlet_err = np.nanstd(DO_in[minday:maxday] - DO_inlet[minday:maxday])

    # calculate mean depth
    fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'vol_df_cas7_c21.p'
    vol_df = pd.read_pickle(fn)
    inlet_vol = vol_df['volume m3'].loc[station+'_p']
    inlet_area = vol_df['area m2'].loc[station+'_p']
    mean_depth = inlet_vol / inlet_area

    # add data to df
    new_data = {'Inlet': [inlet_name],
                'QinDOin': [deep_exchange_avg],
                'QinDOin_err': [deep_exchange_err],
                'QoutDOout': [surf_exchange_transport_avg],
                'QoutDOout_err': [surf_exchange_transport_err],
                'Photosynthesis': [photosynthesis_avg],
                'Photosynthesis_err': [photosynthesis_err],
                'Consumption': [consumption_avg],
                'Consumption_err': [consumption_err],
                'd/dt(DO)': [ddtDO_avg],
                'd/dt(DO)_err': [ddtDO_err],
                'PhysicalResupply': [physresup_avg],
                'PhysicalResupply_err': [physresup_err],
                'NetEcosystemMetabolism': [NEM_avg],
                'NetEcosystemMetabolism_err': [NEM_err],
                'SepOctInletDO[mg/L]': [DOinlet_avg],
                'SepOctInletDO_err[mg/L]': [DOinlet_err],
                'MeanDepth[m]': [mean_depth],
                'DOin-DOinlet[mg/L]': [DOin_DOinlet],
                'DOin-DOinlet_err[mg/L]': [DOin_DOinlet_err],
                'QinDOin_nonorm': [deep_exchange_nonnormalized],
                'Photo_nonorm': [photosynthesis_all_nonnormalized]}
    df_new_rows = pd.DataFrame(new_data)
    inlet_budget_df = pd.concat([inlet_budget_df, df_new_rows],ignore_index=True)

# exchange flow vs. inlet volume
fig, ax = plt.subplots(1,1,figsize=(6,6))
# ax.scatter(inlet_budget_df['QinDOin_nonorm'],inlet_budget_df['Photo_nonorm'],
#                     s=50, edgecolor='white',linewidth=0.5,zorder=5)
# for inlet in inlet_budget_df['Inlet']:
#       ax.text(inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'QinDOin_nonorm'],
#               inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'Photo_nonorm']+0.005,
#                 inlet,color='silver')
ax.scatter(inlet_budget_df['QinDOin'],inlet_budget_df['Photosynthesis'],
                    s=50, edgecolor='white',linewidth=0.5,zorder=5)
for inlet in inlet_budget_df['Inlet']:
      ax.text(inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'QinDOin'],
              inlet_budget_df.loc[inlet_budget_df['Inlet'] == inlet, 'Photosynthesis']+0.005,
                inlet,color='silver')
plt.show()

# get average of all inlets and add data to df
# error of averages is error of all values added in quadrature, divided by number of items
new_data = {'Inlet': ['All inlets'],
            'QinDOin': [np.nanmean(inlet_budget_df['QinDOin'])],
            'QinDOin_err': [np.nanstd(inlet_budget_df['QinDOin_err'])],
            'QoutDOout': [np.nanmean(inlet_budget_df['QoutDOout'])],
            'QoutDOout_err': [np.nanstd(inlet_budget_df['QoutDOout_err'])],
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
            'SepOctInletDO[mg/L]': [np.nanmean(inlet_budget_df['SepOctInletDO[mg/L]'])],
            'SepOctInletDO_err[mg/L]': [np.nanstd(inlet_budget_df['SepOctInletDO[mg/L]'])],
            'MeanDepth[m]': [np.nanmean(inlet_budget_df['MeanDepth[m]'])],
            'DOin-DOinlet[mg/L]': [np.nanmean(inlet_budget_df['DOin-DOinlet[mg/L]'])],
            'DOin-DOinlet_err[mg/L]': [np.nanstd(inlet_budget_df['DOin-DOinlet[mg/L]'])],}
df_new_rows = pd.DataFrame(new_data)
inlet_budget_df = pd.concat([inlet_budget_df, df_new_rows],ignore_index=True)

# save to csv file
print(inlet_budget_df)
inlet_budget_df.to_csv('inlet_budgets_decline_period_mgL_day_onelaye.csv', index=False)

# # print(inlet_budget_df[['Inlet', 'ExchangeFlow', 'ExchangeFlow_err']])
# # print(inlet_budget_df[['Inlet', 'VerticalTransport', 'VerticalTransport_err']])
# # print(inlet_budget_df[['Inlet', 'Photosynthesis', 'Photosynthesis_err']])
# # print(inlet_budget_df[['Inlet', 'Consumption', 'Consumption_err']])
# # print(inlet_budget_df[['Inlet', 'd/dt(DO)', 'd/dt(DO)_err']])
# # print(inlet_budget_df[['Inlet', 'PhysicalResupply', 'PhysicalResupply_err']])
# # print(inlet_budget_df[['Inlet', 'NetEcosystemMetabolism', 'NetEcosystemMetabolism_err']])
# # print(inlet_budget_df[['Inlet', 'SepOctDeepDO[mg/L]', 'SepOctDeepDO_err[mg/L]']])
# # print(inlet_budget_df[['Inlet', 'MeanDepth[m]']])
# print(inlet_budget_df[['Inlet', 'DOin-DOdeep[mg/L]', 'DOin-DOdeep_err[mg/L]']])


# ----------------------- test new scatter plot (budget terms vs DOinlet) -------------------------------

# initialize figure
fig = plt.figure(figsize=(10, 5.5))
gs = GridSpec(2, 12, figure=fig, height_ratios=[3, 4])#, wspace=2.4)
axes = []  # list to store axes
for i in range(4):
    if i == 0:
        ax = fig.add_subplot(gs[0, i*3:(i+1)*3])
    else:
        ax = fig.add_subplot(gs[0, i*3:(i+1)*3], sharex=axes[0])
    axes.append(ax)
for i in range(3):
    ax = fig.add_subplot(gs[1, i*4:(i+1)*4], sharex=axes[0])
    axes.append(ax)
axes = [ax for ax in axes]

# plot scatter points
vars = ['QinDOin','QoutDOout',
       'Photosynthesis','Consumption', 'd/dt(DO)',
       'PhysicalResupply', 'NetEcosystemMetabolism']
colors = ['#0D4B91','#99C5F7','#8F0445','#FCC2DD','black','#488DDB','#F069A8']
letters = ['(a) Inflow','(b) Outflow','(c) Photosynthesis','(d) Consumption',
           '(e) d/dt(DO)','(f) Physical Resupply', '(g) Net Ecosystem\nMetabolism']
ylims = [ [-0.5,4], [-4.5,1], [-0.05,0.4], [-0.15,0.05],
         [-0.2,0.1], [-0.5,0.2], [-0.06,0.25] ]

for i,var in enumerate(vars):
    # add error bars
    axes[i].errorbar(inlet_budget_df['SepOctInletDO[mg/L]'][:-1],inlet_budget_df[var][:-1],
                    xerr=inlet_budget_df['SepOctInletDO_err[mg/L]'][:-1],
                    yerr=inlet_budget_df[var+'_err'][:-1],
                    fmt='o',color='black')
    # plot points
    axes[i].scatter(inlet_budget_df['SepOctInletDO[mg/L]'][:-1],inlet_budget_df[var][:-1], color=colors[i],
                    s=50, edgecolor='white',linewidth=0.5,zorder=5)
    # add mean line
    axes[i].axhline(inlet_budget_df[var].values[-1],0,8, color=colors[i])
    minval = inlet_budget_df[var].values[-1] - inlet_budget_df[var+'_err'].values[-1]
    maxval = inlet_budget_df[var].values[-1] + inlet_budget_df[var+'_err'].values[-1]
    axes[i].fill_between([0,9], [minval,minval], [maxval,maxval],
                         color=colors[i],alpha=0.3)
    # add zero line
    axes[i].axhline(0,0,8, color='gray',linestyle=':')
    # format panel
    axes[i].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
                 transform=axes[i].transAxes, zorder=6, va='top')
    axes[i].set_xlim([0,8])
    axes[i].set_ylim(ylims[i])

    axes[i].tick_params(axis='both', labelsize=12)

    if i in [0,4]:
          axes[i].set_ylabel('Decline period rates\n' + r'[mg L$^{-1}$ d$^{-1}$]',fontsize=12)
    if i >=4 :
          axes[i].set_xlabel(r'Sep-Oct DO$_{inlet}$ [mg L$^{-1}$]',fontsize=12)

plt.tight_layout()
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

# ----------------------- test new scatter plot (budget terms vs DOin-DOinlet) -------------------------------

# initialize figure
fig = plt.figure(figsize=(10, 5.5))
gs = GridSpec(2, 12, figure=fig, height_ratios=[3, 4])#, wspace=2.4)
axes = []  # list to store axes
for i in range(4):
    if i == 0:
        ax = fig.add_subplot(gs[0, i*3:(i+1)*3])
    else:
        ax = fig.add_subplot(gs[0, i*3:(i+1)*3], sharex=axes[0])
    axes.append(ax)
for i in range(3):
    ax = fig.add_subplot(gs[1, i*4:(i+1)*4], sharex=axes[0])
    axes.append(ax)
axes = [ax for ax in axes]

# plot scatter points
vars = ['QinDOin','QoutDOout',
       'Photosynthesis','Consumption', 'd/dt(DO)',
       'PhysicalResupply', 'NetEcosystemMetabolism']
colors = ['#0D4B91','#99C5F7','#8F0445','#FCC2DD','black','#488DDB','#F069A8']
letters = ['(a) Exchange Flow','(b) Vertical','(c) Photosynthesis','(d) Consumption',
           '(e) d/dt(DO)','(f) Physical Resupply', '(g) Net Ecosystem\nMetabolism']
ylims = [ [-0.5,4], [-4.5,1], [-0.05,0.4], [-0.15,0.05],
         [-0.2,0.1], [-0.5,0.2], [-0.06,0.25] ]

for i,var in enumerate(vars):
    # add error bars
    # axes[i].errorbar(inlet_budget_df['DOin-DOinlet[mg/L]'][:-1],inlet_budget_df[var][:-1],
    #                 xerr=inlet_budget_df['DOin-DOinlet_err[mg/L]'][:-1],
    #                 yerr=inlet_budget_df[var+'_err'][:-1],
    #                 fmt='o',color='black')
    # plot points
    axes[i].scatter(inlet_budget_df['DOin-DOinlet[mg/L]'][:-1],inlet_budget_df[var][:-1], color=colors[i],
                    s=50, edgecolor='white',linewidth=0.5,zorder=5)
    # add mean line
    axes[i].axhline(inlet_budget_df[var].values[-1],-1,3, color=colors[i])
    minval = inlet_budget_df[var].values[-1] - inlet_budget_df[var+'_err'].values[-1]
    maxval = inlet_budget_df[var].values[-1] + inlet_budget_df[var+'_err'].values[-1]
    axes[i].fill_between([-1,3], [minval,minval], [maxval,maxval],
                         color=colors[i],alpha=0.3)
    # add zero line
    axes[i].axhline(0,0,8, color='gray',linestyle=':')
    axes[i].axvline(0,0,8, color='gray',linestyle=':')
    # format panel
    axes[i].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
                 transform=axes[i].transAxes, zorder=6, va='top')
    axes[i].set_xlim([-1,3])
    axes[i].set_ylim(ylims[i])

    axes[i].tick_params(axis='both', labelsize=12)

    if i in [0,4]:
          axes[i].set_ylabel('Decline period rates\n' + r'[mg L$^{-1}$ d$^{-1}$]',fontsize=12)
    if i >=4 :
          axes[i].set_xlabel(r'Decline period DO$_{in}$ - DO$_{inlet}$ [mg L$^{-1}$]',fontsize=12)

plt.tight_layout()
plt.show()