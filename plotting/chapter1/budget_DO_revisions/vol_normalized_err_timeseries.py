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
import matplotlib.colors as mcolors
import matplotlib.cm as cm

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

# initialize figure
fig,axes = plt.subplots(2,1, figsize=(10,6), sharex=True)
ax = axes.ravel()
# format budget time series figures
for axis in ax:
    axis.set_xlim([dates_hrly[0],dates_hrly[-2]])
    axis.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')
    axis.tick_params(axis='x', labelrotation=30, labelsize=12)
    axis.tick_params(axis='y', labelsize=12)
    loc = mdates.MonthLocator(interval=1)
    axis.xaxis.set_major_locator(loc)
ax[0].set_title('2017 error terms (10-day Hanning Window)',fontsize=14, fontweight='bold', loc='left')
ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax[0].set_ylabel(r'DO transport [kmol O$_2$ s$^{-1}$]', fontsize=12)
ax[1].set_ylabel(r'DO transport [mg L$^{-1}$ d$^{-1}$]', fontsize=12)
ax[0].text(0.05,0.9,'(a) Volume-integrated', transform=ax[0].transAxes, fontsize=12, fontweight='bold')
ax[1].text(0.05,0.9,'(b) Volume-normalized', transform=ax[1].transAxes, fontsize=12, fontweight='bold')


# COLLAPSE
for i,station in enumerate(sta_dict):
        

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

    # # calculate correlation between surface and deep vertical transports (this should be -1)
    # # calculate r^2 and p value
    # r,p = pearsonr(vertX_surf_DO_TEF[:-1],vertX_deep_DO_TEF[:-1])
    # R2 = round(r**2,2)
    # print('=================')
    # print(station)
    # print('r = {}'.format(round(r,4)))
    # print('p = {}'.format(p))

# ------------------------------- get budget error ----------------------------------------

    # calculate error
    error_DO = vertX_surf_DO_TEF + vertX_deep_DO_TEF

    # calculate mean depth
    fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'vol_df_cas7_c21.p'
    vol_df = pd.read_pickle(fn)
    inlet_vol = vol_df['volume m3'].loc[station+'_p']
    inlet_area = vol_df['area m2'].loc[station+'_p']
    mean_depth = inlet_vol / inlet_area

    norm = mcolors.Normalize(vmin=0, vmax=110)
    cmap_temp = plt.get_cmap('gnuplot2_r')
    cmap = ListedColormap(cmap_temp(np.linspace(0.1, 0.8, 256)))
    # cmap = plt.get_cmap('rainbow')
    error_color = cmap(norm(mean_depth))
    
    # plot volume-integrated error
    conversion = (1000 * 32 * 60 * 60 * 24)
    # integrated terms
    # ax[0].plot(dates_local_daily,zfun.lowpass(error_DO,n=10),color='white',
    #         alpha=0.8, linewidth=2)
    ax[0].plot(dates_local_daily,zfun.lowpass(error_DO,n=10),color=error_color,
            alpha=0.8, linewidth=1.5)
    print('-------------------')
    print(station)
    # print('Volume-integrated: {}'.format(np.nanmean(error_DO)))
    print('Volume-integrated: {}'.format(np.nanmean(np.abs(error_DO))))
    # volume-normalized
    # ax[1].plot(dates_local_daily,zfun.lowpass(error_DO/vol_deep*conversion,n=10),color='white',
    #         alpha=0.8, linewidth=2)
    ax[1].plot(dates_local_daily,zfun.lowpass(error_DO/vol_deep*conversion,n=10),color=error_color,
            alpha=0.8, linewidth=1.5)
    # print('Volume-normalized: {}'.format(np.nanmean(error_DO/vol_deep*conversion)))
    
    
# create colorbar
# This object "connects" the values to the colors
sm = cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([]) # Required for some versions of Matplotlib
cbar_ax = fig.add_axes([0.87, 0.1, 0.03, 0.8])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'Mean depth [m]', fontsize=12)
cbar.outline.set_visible(False)

plt.subplots_adjust(left=0.17,right=0.85)
# plt.tight_layout()
plt.show()

