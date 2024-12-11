"""
Plot surface and deep budget for all 21 inlets
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
import matplotlib.pylab as plt
import gsw
import pickle
import get_two_layer

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
jobname = 'twentyoneinlets'
year = '2017'

# Pick ONE cluster below to analyze

group = 'All inlets'
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
else:
    sta_dict = stations

    
# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
Lfun.make_dir(out_dir)

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')[2::]
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]
# crop time vector (because we only have jan 2 - dec 30)
dates_no_crop = dates_local_daily
dates_local_daily = dates_local_daily

# create color dictionary for stations
basin_dict = {'similk': 'Whidbey',
              'oak': 'Whidbey',
              'crescent': 'Whidbey',
              'penn': 'Whidbey',
              'portsusan': 'Whidbey',
              'holmes': 'Whidbey',
              'dabob': 'Hood Canal',
              'lynchcove': 'Hood Canal',
              'dyes': 'Main Basin',
              'sinclair': 'Main Basin',
              'elliot': 'Main Basin',
              'quartermaster': 'Main Basin',
              'commencement': 'Main Basin',
              'case': 'South Sound',
              'carr': 'South Sound',
              'hammersley': 'South Sound', 
              'totten': 'South Sound', 
              'eld': 'South Sound', 
              'budd': 'South Sound', 
              'henderson': 'South Sound',
              'killsut': 'Admiralty Inlet'}
basin_color_dict = {'Whidbey': 'limegreen',
                    'Hood Canal': 'hotpink',
                    'Main Basin': 'deepskyblue',
                    'South Sound': 'blueviolet',
                    'Admiralty Inlet': 'black'}

print('\n')

##########################################################
##            Get all variables for analysis            ##
##########################################################

print('Getting all data for analysis\n')

# create dictionaries with interface depths
interface_dict = dict()

# create dictionaries to save dataframe of budget terms for each station
bottomlay_dict = {}
surfacelay_dict = {}
DOconcen_dict = {}
dimensions_dict = {}
for i,station in enumerate(sta_dict):
    bottomlay_dict[station] = pd.DataFrame()
    surfacelay_dict[station] = pd.DataFrame()
    DOconcen_dict[station] = pd.DataFrame()
    dimensions_dict[station] = pd.DataFrame()

# get lat and lon of grid
Ldir['ds0'] = startdate
in_dir = Ldir['roms_out'] / Ldir['gtagex']
# G, S, T = zrfun.get_basic_info(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
fn0 = xr.open_dataset(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
lonr = fn0.lon_rho.values
latr = fn0.lat_rho.values
# open box extraction
box_fn = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
ds_box = xr.open_dataset(box_fn)
DX = (ds_box.pm.values)**-1
DY = (ds_box.pn.values)**-1
DA = DX*DY # get area of each grid cell in m^2

# COLLAPSE
for i,station in enumerate(sta_dict):

        # get interface depth from csv file
        with open('interface_depths.csv', 'r') as f:
            for line in f:
                inlet, interface_depth = line.strip().split(',')
                interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
        z_interface = float(interface_dict[station])

# --------------------------- get EU exchange flow terms ----------------------------------------
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_EU_exchange' / (station + '.p')
        df_exchange = pd.read_pickle(fn)
        exchange_surf_unfiltered = df_exchange['surface [kmol/s]'].values 
        exchange_deep_unfiltered = df_exchange['deep [kmol/s]'].values
        # Godin filter
        EU_surf = zfun.lowpass(exchange_surf_unfiltered, f='godin')[36:-34:24] 
        EU_deep = zfun.lowpass(exchange_deep_unfiltered, f='godin')[36:-34:24]
        exchange_color = 'cornflowerblue'

# --------------------------- get TEF exchange flow terms ----------------------------------------
        in_dir = Ldir['LOo'] / 'extract' / 'cas7_t0_x4b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
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
        bgc_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_bgc' / station
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
        photo_color = 'turquoise'
        cons_surf_unfiltered = [] # nitrification, respiration
        cons_deep_unfiltered = [] # nitrification, respiration, sediment oxygen demand
        cons_color = 'deeppink'
        airsea_surf_unfiltered = []
        airsea_color = 'plum'
        o2vol_surf_unfiltered = []
        o2vol_deep_unfiltered = []
        ddtDOV_color = 'gray'

        # combine all months
        for month in months:
            fn = station + '_' + month + '.p'
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

        # take time derivative of (DO*V) to get d/dt (DO*V)
        ddtDOV_surf_unfiltered = np.diff(o2vol_surf_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
        ddtDOV_deep_unfiltered = np.diff(o2vol_deep_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s

        # apply Godin filter
        photo_surf = zfun.lowpass(photo_surf_unfiltered, f='godin')[36:-34:24]
        photo_deep = zfun.lowpass(photo_deep_unfiltered, f='godin')[36:-34:24]
        cons_surf   = zfun.lowpass(cons_surf_unfiltered, f='godin')[36:-34:24]
        cons_deep   = zfun.lowpass(cons_deep_unfiltered, f='godin')[36:-34:24]
        airsea_surf = zfun.lowpass(airsea_surf_unfiltered, f='godin')[36:-34:24]
        ddtDOV_surf = zfun.lowpass(ddtDOV_surf_unfiltered, f='godin')[36:-34:24]
        ddtDOV_deep = zfun.lowpass(ddtDOV_deep_unfiltered, f='godin')[36:-34:24]

# ------------------------------- get rivers and WWTPs ----------------------------------------
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_traps' / (station + '.p')
        df_traps = pd.read_pickle(fn)
        rivers_surf_unfiltered = df_traps['surface [kmol/s]'].values
        wwtps_deep_unfiltered = df_traps['deep [kmol/s]'].values

        # Pad with zeros if no rivers or WWTPs
        if rivers_surf_unfiltered.size == 0 or np.isnan(rivers_surf_unfiltered[0]):
                rivers_surf_unfiltered = np.zeros(8761)
        if wwtps_deep_unfiltered.size == 0:
                wwtps_deep_unfiltered = np.zeros(8761)
        # Godin filter
        traps_surf = zfun.lowpass(rivers_surf_unfiltered, f='godin')[36:-34:24]
        traps_deep =  zfun.lowpass(wwtps_deep_unfiltered, f='godin')[36:-34:24]
        traps_color = 'black'

# ------------------------------- get vertical exchange ----------------------------------------

        # pick color
        vertX_color = 'mediumorchid'

        # calculate error term, and ascribe to the vertical exchange
        vertX_surf_TEF = ddtDOV_surf - (TEF_surf 
                                    + photo_surf
                                    + cons_surf
                                    + airsea_surf
                                    + traps_surf)
        
        vertX_deep_TEF = ddtDOV_deep - (TEF_deep 
                                    + photo_deep
                                    + cons_deep
                                    + traps_deep)
        
        # calculate error term, and ascribe to the vertical exchange
        vertX_surf_EU = ddtDOV_surf - (EU_surf 
                                    + photo_surf
                                    + cons_surf
                                    + airsea_surf
                                    + traps_surf)
        
        vertX_deep_EU = ddtDOV_deep - (EU_deep 
                                    + photo_deep
                                    + cons_deep
                                    + traps_deep)
        
# ---------------------------------- DO concentrations --------------------------------------------

        # get average V of each layer
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
        df_V = pd.read_pickle(fn)
        # Godin filter already applied earlier in workflow
        surf_V = df_V['surface [m3]'].values[1:-1]
        deep_V = df_V['deep [m3]'].values[1:-1]

        # apply Godin filter to DO*V (mmol)
        o2vol_surf = zfun.lowpass(o2vol_surf_unfiltered, f='godin')[36:-34:24]
        o2vol_deep = zfun.lowpass(o2vol_deep_unfiltered, f='godin')[36:-34:24]

        # get average DO of each layer by calculating (DO*V)/V
        # and convert from mmol/m3 to mg/L
        o2_surf = o2vol_surf / surf_V * 32/1000
        o2_deep = o2vol_deep / deep_V * 32/1000

        # get average inflowing DO concentration (Qin concentration)
        DO_in = DO_p.values * 32/1000

        # get bottom DO
        # get section information
        # ds = xr.open_dataset('../../../../LO_output/extract/'+gtagex+
        #                     '/tef2/extractions_'+startdate+'_'+enddate+
        #                     '/'+station+'.nc')
        
        #% Get indices of inlet
        seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
        seg_df = pd.read_pickle(seg_name)
        ji_list = seg_df[station+'_p']['ji_list']
        jj_LO = [x[0] for x in ji_list] # y; lat; jj
        ii_LO = [x[1] for x in ji_list] # x; lon; ii
        # get lat and lon corresponding to ii and jj indices
        lat_LO = latr[jj_LO,0]
        lon_LO = lonr[0,ii_LO]
        # get corresponding ii and jj indices in box extraction
        lat_box_all = ds_box['lat_rho'].values[:,0]
        lon_box_all = ds_box['lon_rho'].values[0,:]
        jj = np.zeros(len(jj_LO))
        ii = np.zeros(len(ii_LO))
        for j_ind,lat in enumerate(lat_LO):
            jj[j_ind] = np.where(lat_box_all==lat)[0][0]
        for i_ind,lon in enumerate(lon_LO):
            ii[i_ind] = np.where(lon_box_all==lon)[0][0]
        # convert to array of ints
        jj = jj.astype(int)
        ii = ii.astype(int)

        # get bottom DO from my bottom DO get_DO script extraction file
        ds_oxy = xr.open_dataset(Ldir['LOo'] / 'pugetsound_DO' / 'data' / (year + '_DO_info_' + 'withStraits' + '.nc'))

        # calculate average bottom oxygen throughout section
        bottom_oxygen_alldays = np.nanmean(ds_oxy['DO_bot'].values[:,jj,ii],axis=1) # mg/L
        bottom_oxygen = bottom_oxygen_alldays[1:-1]
        # calculate minimum bottom DO throughout section
        oxygen_min_alldays = np.nanmin(ds_oxy['DO_bot'].values[:,jj,ii],axis=1) # mg/L
        oxygen_min = oxygen_min_alldays[1:-1]

        # calculate percent hypoxic volume
        hyp_thick = ds_oxy['hyp_thick'].values # m

        # get hypoxic thickness and cell area corresponding to the inlet
        hyp_thick = hyp_thick[:,jj,ii] # m
        DA_inlet = DA[jj,ii]
        # calculate hypoxic volume
        hyp_vol = np.sum(hyp_thick * DA_inlet, axis=(1)) # m^3
        # crop to 363 days
        hyp_vol = hyp_vol[1:-1]


# -------------------------- get volume and mean depth of layers -------------------------

        # get volume of layers
        fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
        df_V = pd.read_pickle(fn)
        # Godin filter already applied earlier in workflow
        surf_V = df_V['surface [m3]'].values[1:-1]
        deep_V = df_V['deep [m3]'].values[1:-1]

# --------------------------- get other physical dimensions ------------------------------

        # get full volume of inlet, and mean depth
        fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'vol_df_cas7_c21.p'
        vol_df = pd.read_pickle(fn)
        inlet_vol = vol_df['volume m3'].loc[station+'_p']
        inlet_area = vol_df['area m2'].loc[station+'_p']
        mean_depth = inlet_vol / inlet_area
        # estimate aspect ratio from inlet area and width at mouth
        # aspect ratio = length/width (so longer inlet has larger ratio)
        # thus, aspect ratio = area / width^2
        # since this is an approximation, we use the nominal cell width of 500 m
        fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'sect_df_cas7_c21.p'
        mouth_df = pd.read_pickle(fn)
        mouth_df = mouth_df.loc[mouth_df['sn'] == station]
        # get number of horizontal and vertical entrance points
        uv_list = mouth_df['uv'].values.tolist()
        u_count = uv_list.count('u')
        v_count = uv_list.count('v')
        mouth_width = 500 * np.sqrt(u_count**2 + v_count**2)
        est_aspect_ratio = inlet_area / mouth_width**2
        # get mouth area, for calculating fluxes 

        
# ------------------------- save data in dataframe dict -----------------------------------
        # Note: everything is in units of kmol O2 /s

        # surface layer
        surfacelay_dict[station]['EU Exchange Flow'] = EU_surf
        surfacelay_dict[station]['TEF Exchange Flow'] = TEF_surf
        surfacelay_dict[station]['EU Recirculation'] = EU_surf + vertX_surf_EU
        surfacelay_dict[station]['TEF Recirculation'] = TEF_surf + vertX_surf_TEF
        surfacelay_dict[station]['TRAPS'] = traps_surf
        surfacelay_dict[station]['Photosynthesis'] = photo_surf
        surfacelay_dict[station]['Bio Consumption'] = cons_surf
        surfacelay_dict[station]['Air-Sea Transfer'] = airsea_surf
        surfacelay_dict[station]['EU Vertical'] = vertX_surf_EU
        surfacelay_dict[station]['TEF Vertical'] = vertX_surf_TEF
        surfacelay_dict[station]['Storage'] = ddtDOV_surf
        surfacelay_dict[station]['Volume'] = surf_V

        # bottom layer
        bottomlay_dict[station]['Storage'] = ddtDOV_deep
        bottomlay_dict[station]['EU Exchange Flow'] = EU_deep
        bottomlay_dict[station]['TEF Exchange Flow'] = TEF_deep
        bottomlay_dict[station]['EU Recirculation'] = EU_deep + vertX_deep_EU
        bottomlay_dict[station]['TEF Recirculation'] = TEF_deep + vertX_deep_TEF
        bottomlay_dict[station]['EU Vertical'] = vertX_deep_EU
        bottomlay_dict[station]['TEF Vertical'] = vertX_deep_TEF
        bottomlay_dict[station]['TRAPS'] = traps_deep
        bottomlay_dict[station]['Photosynthesis'] = photo_deep
        bottomlay_dict[station]['Bio Consumption'] = cons_deep
        bottomlay_dict[station]['Photosynthesis & Consumption'] = photo_deep + cons_deep
        bottomlay_dict[station]['Volume'] = deep_V
        bottomlay_dict[station]['Qin m3/s'] = Q_p.values 

# ------------------------- save DO concentrations in dataframe dict -----------------------------------

        # mg/L units

        # DO concentrations
        DOconcen_dict[station]['Surface Layer'] = o2_surf
        DOconcen_dict[station]['Deep Layer'] = o2_deep
        DOconcen_dict[station]['Bottom Sigma DO'] = bottom_oxygen
        DOconcen_dict[station]['Minimum Bottom Layer DO'] = oxygen_min
        DOconcen_dict[station]['Qin DO'] = DO_in
        DOconcen_dict[station]['percent hypoxic volume'] = hyp_vol/[inlet_vol] * 100

# ------------------------ save inlet dimensions in dataframe dict ----------------------------------------

        dimensions_dict[station]['Inlet volume'] = [inlet_vol] # m3
        dimensions_dict[station]['Mean depth'] = [mean_depth] # m
        dimensions_dict[station]['Mouth width'] = mouth_width
        dimensions_dict[station]['L/W aspect ratio'] = [est_aspect_ratio] # dimensionless

        # print keys
        if i == 0:
            print(list(surfacelay_dict[station].keys()))
            print(list(DOconcen_dict[station].keys()))
            print(list(dimensions_dict[station].keys()))


##########################################################
##                   DO scatterplots                    ## 
##########################################################

# print(dates_local_daily)

DO_analysis = True
if DO_analysis == True:

    # initialize arrays for plotting
    intervals = 12
    deep_lay_DO = np.zeros(len(sta_dict)*intervals)
    bott_sig_DO = np.zeros(len(sta_dict)*intervals)
    min_bott_sig_DO = np.zeros(len(sta_dict)*intervals)
    annual_mean_DO = np.zeros(len(sta_dict)*intervals)
    perc_hyp_vol = np.zeros(len(sta_dict)*intervals)
    mean_DOin = np.zeros(len(sta_dict)*intervals)
    mean_Tflush = np.zeros(len(sta_dict)*intervals)
    mean_TEFin = np.zeros(len(sta_dict)*intervals)
    mean_recirc = np.zeros(len(sta_dict)*intervals)
    mean_cons = np.zeros(len(sta_dict)*intervals)
    mean_depth = np.zeros(len(sta_dict)*intervals)
    inlet_vol = np.zeros(len(sta_dict)*intervals)
    aspect_ratio = np.zeros(len(sta_dict)*intervals)

    annmean_DOin = np.zeros(len(sta_dict))
    annmin_DOin = np.zeros(len(sta_dict))
    colors_twentyone = []

    colors = []

    # get values for plotting and calculating r value
    for i,station in enumerate(sta_dict):
        # get interface depth from csv file
        with open('interface_depths.csv', 'r') as f:
            for line in f:
                inlet, interface_depth = line.strip().split(',')
                interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
        z_interface = float(interface_dict[station])
        fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_2014.01.01_2014.12.31.nc'
        ds = xr.open_dataset(fn)
        moor_depth = ds.h.values
        
        for month in range(intervals):
            if month == 0:
                minday = 0 #1
                maxday = 30 #32
            elif month == 1:
                minday = 30# 32
                maxday = 58# 60
            elif month == 2:
                minday = 58# 60
                maxday = 89# 91
            elif month == 3:
                minday = 89# 91
                maxday = 119# 121
            elif month == 4:
                minday = 119# 121
                maxday = 150# 152
            elif month == 5:
                minday = 150# 152
                maxday = 180# 182
            elif month == 6:
                minday = 180# 182
                maxday = 211# 213
            elif month == 7:
                minday = 211# 213
                maxday = 242# 244
            elif month == 8:
                minday = 242# 244
                maxday = 272# 274
            elif month == 9:
                minday = 272# 274
                maxday = 303# 305
            elif month == 10:
                minday = 303# 305
                maxday = 332# 335
            elif month == 11:
                minday = 332# 335
                maxday = 363
            # minday=month
            # maxday=month+1
            # save values
            deep_lay_DO[i*intervals+month] =  np.nanmean(DOconcen_dict[station]['Deep Layer'][minday:maxday])
            bott_sig_DO[i*intervals+month] =  np.nanmean(DOconcen_dict[station]['Bottom Sigma DO'][minday:maxday])
            min_bott_sig_DO[i*intervals+month] =  np.nanmin(DOconcen_dict[station]['Minimum Bottom Layer DO'][minday:maxday])
            annual_mean_DO[i*intervals+month] = np.nanmean(DOconcen_dict[station]['Deep Layer'])
            perc_hyp_vol[i*intervals+month] = np.nanmean(DOconcen_dict[station]['percent hypoxic volume'][minday:maxday])
            mean_DOin[i*intervals+month] = np.nanmean(DOconcen_dict[station]['Qin DO'][minday:maxday])
            mean_Tflush[i*intervals+month] = np.nanmean(dimensions_dict[station]['Inlet volume'][0]/bottomlay_dict[station]['Qin m3/s'][minday:maxday]) / (60*60*24)
            mean_TEFin[i*intervals+month] = np.nanmean(bottomlay_dict[station]['TEF Exchange Flow'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]) * (
                                    32 * 1000) * (60*60*24)
            mean_recirc[i*intervals+month] = np.nanmean(bottomlay_dict[station]['TEF Recirculation'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]) * (
                                    32 * 1000) * (60*60*24)
            mean_cons[i*intervals+month] = np.nanmean(bottomlay_dict[station]['Bio Consumption'][minday:maxday]/bottomlay_dict[station]['Volume'][minday:maxday]) * (
                                    32 * 1000) * (60*60*24)
            mean_depth[i*intervals+month] = dimensions_dict[station]['Mean depth'][0]
            inlet_vol[i*intervals+month] = dimensions_dict[station]['Inlet volume'][0]
            aspect_ratio[i*intervals+month] = dimensions_dict[station]['L/W aspect ratio'][0]
            colors.append(basin_color_dict[basin_dict[station]])

        annmean_DOin[i] = np.nanmean(DOconcen_dict[station]['Qin DO'])
        annmin_DOin[i] = np.nanmin(DOconcen_dict[station]['Qin DO'])
        colors_twentyone.append(basin_color_dict[basin_dict[station]])
        
# create dataframe
df = pd.DataFrame()

df['DOdeep [mg/L]'] = deep_lay_DO
df['DOin [mg/L]'] = mean_DOin
df['WRT [days]'] = mean_Tflush

# Save the DataFrame to an Excel file
df.to_excel('monthly_mean_DO_13inlets.xlsx', index=False)


    # # DOin vs. Tflush colored by percent hypoxic volume ============== 4PART MONEY PLOT SCATTER ================================
    # percent_hypoxic = True
    # if percent_hypoxic == True:
    #     # initialize figure
    #     fig, ax = plt.subplots(2,2,figsize = (10,9))
    #     ax = ax.ravel()

    #     # format figure
    #     ax[0].set_title('(a) ' + year + r' monthly mean DO$_{deep}$ vs. DO$_{in}$' + '\n' + r'colored by mean T$_{flush}$',
    #                     size=14, loc='left')
    #     # ax[0].set_title(year + r' monthly mean DO$_{deep}$ vs. DO$_{in}$',
    #     #                 size=16, loc='left')
    #     # format grid
    #     # ax[0].set_facecolor('#EEEEEE')
    #     ax[0].tick_params(axis='x', labelrotation=30)
    #     # ax[0].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    #     # for border in ['top','right','bottom','left']:
    #     #     ax[0].spines[border].set_visible(False)
    #     ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    #     ax[0].tick_params(axis='both', labelsize=12)
    #     ax[0].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    #     ax[0].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
    #     # plot
    #     cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
    #     cmap_oxy = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
    #     ax[0].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, c='k', alpha=0.5)
    #     ax[0].plot([0,12],[0,12],color='gray')
    #     # ax[0].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='#EEEEEE',zorder=4, fontsize=10)
    #     ax[0].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
    #     cs = ax[0].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
    #     # create colorbarlegend
    #     cbar = fig.colorbar(cs)
    #     cbar.ax.tick_params(labelsize=12)
    #     cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
    #     cbar.outline.set_visible(False)
    #     ax[0].set_xlim([0,10])
    #     ax[0].set_ylim([0,10])


    #     # format figure
    #     ax[1].set_title('(b) ' + year + ' monthly mean '+r'DO$_{in}$ vs. T$_{flush}$'+'\ncolored by % hypoxic volume',
    #                     size=14, loc='left')
    #     # format grid
    #     # ax[1].set_facecolor('#EEEEEE')
    #     ax[1].tick_params(axis='x', labelrotation=30)
    #     # ax[1].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    #     # for border in ['top','right','bottom','left']:
    #     #     ax[1].spines[border].set_visible(False)
    #     ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    #     ax[1].tick_params(axis='both', labelsize=12)
    #     ax[1].set_xlabel(r'Monthly mean T$_{flush}$ [days]', fontsize=12)
    #     ax[1].set_ylabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    #     ax[1].set_ylim([0,10])
    #     ax[1].set_xlim([0,85])
    #     # plot
    #     cmap_hyp = plt.cm.get_cmap('gist_heat_r')
    #     cs_DO = ax[1].scatter(mean_Tflush,mean_DOin,s=60,zorder=5,edgecolor='gray',c=perc_hyp_vol,cmap=cmap_hyp)
    #     # create colorbarlegend
    #     cbar = fig.colorbar(cs_DO)
    #     cbar.ax.tick_params(labelsize=12)
    #     cbar.ax.set_ylabel('Monthly mean % hypoxic volume', rotation=90, fontsize=12)
    #     cbar.outline.set_visible(False)

    #     # crescent bay
    #     ax[2].set_title('(c) Crescent Bay 2017 monthly mean \n' + r'DO$_{deep}$ vs. DO$_{in}$ colored by T$_{flush}$', loc='left', size=14)
    #     # format grid
    #     # ax[2].set_facecolor('#EEEEEE')
    #     ax[2].tick_params(axis='x', labelrotation=30)
    #     # ax[2].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    #     # for border in ['top','right','bottom','left']:
    #     #     ax[2].spines[border].set_visible(False)
    #     ax[2].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    #     ax[2].tick_params(axis='both', labelsize=12)
    #     ax[2].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    #     ax[2].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
    #     # plot
    #     cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
    #     cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
    #     ax[2].plot([0,11],[0,11],color='dimgray')
    #     ax[2].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
    #     # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
    #     ax[2].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, color='gray',alpha=0.5, edgecolor='none')
    #     for i,station in enumerate(sta_dict):
    #         if station == 'crescent':
    #             cs = ax[2].scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
    #                             s=150, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_tflush,
    #                         linewidth=2, vmin=0, vmax=40)
    #         else:
    #             continue
    #     # create colorbarlegend
    #     cbar = fig.colorbar(cs)
    #     cbar.ax.tick_params(labelsize=12)
    #     cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
    #     cbar.outline.set_visible(False)
    #     ax[2].set_xlim([0,11])
    #     ax[2].set_ylim([0,11])

    #     # lynch cove
    #     ax[3].set_title('(d) Lynch Cove 2017 monthly mean \n' + r'DO$_{deep}$ vs. DO$_{in}$ colored by T$_{flush}$', loc='left', size=14)
    #     # format grid
    #     # ax[3].set_facecolor('#EEEEEE')
    #     ax[3].tick_params(axis='x', labelrotation=30)
    #     # ax[3].grid(True,color='w',linewidth=1,linestyle='-',axis='both')
    #     # for border in ['top','right','bottom','left']:
    #     #     ax[3].spines[border].set_visible(False)
    #     ax[3].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    #     ax[3].tick_params(axis='both', labelsize=12)
    #     ax[3].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    #     ax[3].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
    #     # plot
    #     cmap_temp = plt.cm.get_cmap('cubehelix_r', 256)
    #     cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
    #     ax[3].plot([0,11],[0,11],color='dimgray')
    #     ax[3].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
    #     # cs = ax.scatter(mean_DOin,deep_lay_DO,s=80, zorder=5, c=mean_Tflush, cmap=cmap_oxy)
    #     ax[3].scatter(mean_DOin,deep_lay_DO,s=60, zorder=5, color='gray',alpha=0.5, edgecolor='none')
    #     for i,station in enumerate(sta_dict):
    #         if station == 'lynchcove':
    #             cs = ax[3].scatter(mean_DOin[i*intervals:(i+1)*intervals],deep_lay_DO[i*intervals:(i+1)*intervals],marker='s',
    #                             s=150, zorder=6, c=mean_Tflush[i*intervals:(i+1)*intervals], edgecolor='black',cmap=cmap_tflush,
    #                         linewidth=2, vmin=0, vmax=40)
    #         else:
    #             continue
    #     # create colorbarlegend
    #     cbar = fig.colorbar(cs)
    #     cbar.ax.tick_params(labelsize=12)
    #     cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
    #     cbar.outline.set_visible(False)
    #     ax[3].set_xlim([0,11])
    #     ax[3].set_ylim([0,11])


    #     plt.tight_layout()
    #     # save figure
    #     plt.show()