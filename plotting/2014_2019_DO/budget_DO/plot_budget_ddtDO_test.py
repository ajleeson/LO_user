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
from scipy import stats
import math
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap
import csv
import cmocean
from scipy.stats import pearsonr
from scipy.stats import ttest_ind
from scipy.stats import f_oneway
from scipy.stats import kruskal
from scipy.stats import alexandergovern
import pingouin as pg
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

DIN_in_direct = [None] * len(sta_dict)

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
        TEF_deep = (Q_p.values * DO_p.values) * (1/1000) * (1/1000) # Qin * DOin

        # nutrients (NO3 + NH4)
        DIN_p = tef_df['NO3_p']+tef_df['NH4_p'] # NH4in [mmol/m3] # NO3in [mmol/m3]
        DIN_m = tef_df['NO3_m']+tef_df['NH4_m'] # NH4out [mmol/m3] # NO3out [mmol/m3]
        # print(np.nanmean(DIN_p))
        # convert from mmol/s to kmol/s
        QoutDINout = (Q_m.values * DIN_m.values) * (1/1000) * (1/1000) # Qout * DINout
        QinDINin = (Q_p.values * DIN_p.values) * (1/1000) * (1/1000) # Qin * DINin

        DIN_in_direct[i] = np.nanmean(DIN_p[164:225])

        # phytoplankton
        P_p = tef_df['phytoplankton_p'] # Pin
        P_m = tef_df['phytoplankton_m'] # Pout
        # convert from mmol/s to kmol/s
        QoutPout = (Q_m.values * P_m.values) * (1/1000) * (1/1000) # Qout * Pout
        QinPin = (Q_p.values * P_p.values) * (1/1000) * (1/1000) # Qin * Pin

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
            cons_surf_unfiltered = np.concatenate((cons_surf_unfiltered, surf_cons_terms * conv * -1)) # kmol/s; multiply by -1 QtoP/c loss term
            cons_deep_unfiltered = np.concatenate((cons_deep_unfiltered, deep_cons_terms * conv * -1)) # kmol/s; multiply by -1 QtoP/c loss term
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
        surfacelay_dict[station]['TEF Exchange Flow'] = TEF_surf

        surfacelay_dict[station]['QouDINout'] = QoutDINout

        surfacelay_dict[station]['TEF Recirculation'] = TEF_surf + vertX_surf_TEF
        surfacelay_dict[station]['TRAPS'] = traps_surf
        surfacelay_dict[station]['Photosynthesis'] = photo_surf
        surfacelay_dict[station]['Bio Consumption'] = cons_surf
        surfacelay_dict[station]['Air-Sea Transfer'] = airsea_surf
        surfacelay_dict[station]['TEF Vertical'] = vertX_surf_TEF
        surfacelay_dict[station]['Storage'] = ddtDOV_surf
        surfacelay_dict[station]['Volume'] = surf_V

        # bottom layer
        bottomlay_dict[station]['Storage'] = ddtDOV_deep

        bottomlay_dict[station]['QinDINin'] = QinDINin
        bottomlay_dict[station]['QinPin'] = QinPin

        bottomlay_dict[station]['TEF Exchange Flow'] = TEF_deep
        bottomlay_dict[station]['TEF Recirculation'] = TEF_deep + vertX_deep_TEF
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

print(np.nanmean(DIN_in_direct)/1000)

# initialize arrays for plotting
deep_lay_DO = np.zeros(len(sta_dict))
mean_DOin = np.zeros(len(sta_dict))
mean_Tflush = np.zeros(len(sta_dict))
mean_TEFin = np.zeros(len(sta_dict))
mean_recirc = np.zeros(len(sta_dict))
mean_cons = np.zeros(len(sta_dict))
mean_depth = np.zeros(len(sta_dict))

mean_Qin = np.zeros(len(sta_dict))
mean_QinDINin = np.zeros(len(sta_dict))
mean_QinDOin = np.zeros(len(sta_dict))
mean_QoutDOout = np.zeros(len(sta_dict))
mean_QinPin = np.zeros(len(sta_dict))
mean_photo = np.zeros(len(sta_dict))
mean_airsea = np.zeros(len(sta_dict))
mean_rivers = np.zeros(len(sta_dict))

inlet_vol = np.zeros(len(sta_dict))

colors = []

# get values for plotting and calculating r value
for i,station in enumerate(sta_dict):

    # get individual inlet volume
    ind_inlet_vol = dimensions_dict[station]['Inlet volume'][0]

    # drawdown period
    minday = 164
    maxday = 225 

    deep_lay_DO[i] =  np.nanmean(DOconcen_dict[station]['Deep Layer'][minday:maxday])
    mean_DOin[i] = np.nanmean(DOconcen_dict[station]['Qin DO'][minday:maxday])

    mean_Qin[i] = np.nanmean(bottomlay_dict[station]['Qin m3/s'][minday:maxday]/ind_inlet_vol) * (60*60*24)  #1/Tflush (days-1)
    mean_QinDINin[i] = np.nanmean(bottomlay_dict[station]['QinDINin'][minday:maxday]/ind_inlet_vol) * (1000) * (60*60*24) # mmol/L per day
    mean_QinPin[i] = np.nanmean(bottomlay_dict[station]['QinPin'][minday:maxday]/ind_inlet_vol) * (1000) * (60*60*24) # mmol/L per day
    mean_QinDOin[i] = np.nanmean(bottomlay_dict[station]['TEF Exchange Flow'][minday:maxday]/ind_inlet_vol) * (32 * 1000) * (60*60*24) # QinDOin [mg/L per day]
    mean_QoutDOout[i] = np.nanmean(surfacelay_dict[station]['TEF Exchange Flow'][minday:maxday]/ind_inlet_vol) * (32 * 1000) * (60*60*24) # QinDOin [mg/L per day]
                            
    
    mean_photo[i] = np.nanmean((bottomlay_dict[station]['Photosynthesis'][minday:maxday]+
                                surfacelay_dict[station]['Photosynthesis'][minday:maxday])
                                /ind_inlet_vol) * (32 * 1000) * (60*60*24) # mg/L per day
    
    mean_airsea[i] = np.nanmean(surfacelay_dict[station]['Air-Sea Transfer'][minday:maxday]
                                /ind_inlet_vol) * (32 * 1000) * (60*60*24) # mg/L per day
    mean_rivers[i] = np.nanmean((surfacelay_dict[station]['TRAPS'][minday:maxday]+
                                bottomlay_dict[station]['TRAPS'][minday:maxday])
                                /ind_inlet_vol) * (32 * 1000) * (60*60*24) # mg/L per day
    
    mean_Tflush[i] = np.nanmean(ind_inlet_vol/bottomlay_dict[station]['Qin m3/s'][minday:maxday]) / (60*60*24)
    
    mean_cons[i] = np.nanmean((bottomlay_dict[station]['Bio Consumption'][minday:maxday]+
                                surfacelay_dict[station]['Bio Consumption'][minday:maxday])
                                /ind_inlet_vol) * (32 * 1000) * (60*60*24) # mg/L per day

    mean_depth[i] = dimensions_dict[station]['Mean depth'][0]
    inlet_vol[i] = ind_inlet_vol
    

# Plot scatter plots
# initialize figure
fig, ax = plt.subplots(2,3,figsize = (14,9))
ax = ax.ravel()

# DIN
ax[0].set_title('(a) DIN fluxes\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
ax[0].tick_params(axis='x', labelrotation=30)
ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[0].tick_params(axis='both', labelsize=12)
ax[0].set_xlabel(r'Mean $\frac{Q_{in}}{V_{inlet}}=\frac{1}{T_{flush}}$ [day-1]', fontsize=12)
ax[0].set_ylabel(r'Mean $\frac{Q_{in}DIN_{in}}{V_{inlet}}$ [mmol/L per day]', fontsize=12)
# plot
ax[0].scatter(mean_Tflush**-1,mean_QinDINin,s=60, zorder=5, c='navy', alpha=0.5)
#     ax[0].set_ylim([0,0.018])
ax[0].set_xlim([0,0.35])
#     print(mean_QinDINin)
#     print(mean_QinDINin/mean_Tflush**-1)
# conduct linear regession
x = mean_Tflush**-1
res = stats.linregress(x,mean_QinDINin)
print(res.slope)
ax[0].plot(x, res.intercept + res.slope*x, 'orchid')
ax[0].text(0.1,0.9,f'R-squared: {res.rvalue**2:.3f}',transform=ax[0].transAxes,
            fontsize=12,fontweight='bold')
ax[0].text(0.1,0.82,f'Slope: {res.slope:.4f}',transform=ax[0].transAxes,
            fontsize=12,fontweight='bold')
QtoN = res.slope * 1000 

# Consumption vs. Photosynthesis
ax[1].set_title('(b) consumpton vs. photosynthesis\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
ax[1].tick_params(axis='x', labelrotation=30)
ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[1].tick_params(axis='both', labelsize=12)
ax[1].set_ylabel(r'Mean $\frac{Consumption}{V_{inlet}}$ [mg/L per day]', fontsize=12)
ax[1].set_xlabel(r'Mean $\frac{Photosynthesis}{V_{inlet}}$ [mg/L per day]', fontsize=12)
# plot
ax[1].scatter(mean_photo,mean_cons,s=60, zorder=5, c='navy', alpha=0.5)
ax[1].set_ylim([-0.35,0])
ax[1].set_xlim([0,0.6])
# conduct linear regession
x = mean_photo
res = stats.linregress(x,mean_cons)
ax[1].plot(x, res.intercept + res.slope*x, 'orchid')
ax[1].text(0.5,0.9,f'R-squared: {res.rvalue**2:.3f}',transform=ax[1].transAxes,
            fontsize=12,fontweight='bold')
ax[1].text(0.5,0.82,f'Slope: {res.slope:.3f}',transform=ax[1].transAxes,
            fontsize=12,fontweight='bold')
ax[1].text(0.5,0.74,f'intercept: {res.intercept:.3f}',transform=ax[1].transAxes,
            fontsize=12,fontweight='bold')

# # Phytoplankton 
# ax[2].set_title('(c) Phytoplankton fluxes\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
# ax[2].tick_params(axis='x', labelrotation=30)
# ax[2].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
# ax[2].tick_params(axis='both', labelsize=12)
# ax[2].set_xlabel(r'Mean $\frac{Q_{in}}{V_{inlet}}=\frac{1}{T_{flush}}$ [day-1]', fontsize=12)
# ax[2].set_ylabel(r'Mean $\frac{Q_{in}P_{in}}{V_{inlet}}$ [mmol/L per day]', fontsize=12)
# # plot
# ax[2].scatter(mean_Tflush**-1,mean_QinPin,s=60, zorder=5, c='navy', alpha=0.5)
# ax[2].set_ylim([0,0.0025])
# ax[2].set_xlim([0,0.35])
# # conduct linear regession
# x = mean_Tflush**-1
# res = stats.linregress(mean_Tflush**-1,mean_QinPin)
# ax[2].plot(x, res.intercept + res.slope*x, 'orchid')
# ax[2].text(0.1,0.9,f'R-squared: {res.rvalue**2:.3f}',transform=ax[2].transAxes,
#            fontsize=12,fontweight='bold')
# ax[2].text(0.1,0.82,f'Slope: {res.slope:.6f}',transform=ax[2].transAxes,
#            fontsize=12,fontweight='bold')
# QtoP = res.slope * 1000

# Photosynthesis growth rate
ax[2].set_title('(c) volume-normalized photosynthesis\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
ax[2].tick_params(axis='x', labelrotation=30)
ax[2].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[2].tick_params(axis='both', labelsize=12)
ax[2].set_xlabel(r'Mean $\frac{Q_{in}}{V_{inlet}}=\frac{1}{T_{flush}}$ [day-1]', fontsize=12)
ax[2].set_ylabel(r'Mean $\frac{Photosynthesis}{V_{inlet}}$ [mg/L per day]', fontsize=12)
# plot
ax[2].scatter(mean_Tflush**-1,mean_photo,s=60, zorder=5, c='navy', alpha=0.5)
#     # create growth curve
#     Tflush_inv = np.linspace(0,0.35,100)
#     ks = 0.1
#     Cox = 138/16 # mmol O2 : mmol N
#     mu0 = 1.7
#     alpha = 0.07
#     E = 5 # what is a typical value for this????
#     sunlight = mu0 * alpha * E /(np.sqrt(mu0**2 + alpha**2*E**2))
# #     print(sunlight)
#     N_Tflush = np.nanmean(mean_QinDINin/mean_Qin) # QinDINin/Qin = DINin. Need to multiply by Tflush (or divide by 1/Tflush)
#     P_Tflush = np.nanmean(mean_QinPin/mean_Qin)
#     growth_prediction = (N_Tflush/Tflush_inv / (ks + 2*np.sqrt(ks*N_Tflush/Tflush_inv) + N_Tflush/Tflush_inv)) * (P_Tflush/Tflush_inv) * 32 * Cox * sunlight
#     ax[3].plot(Tflush_inv, growth_prediction, 'orchid')
#     # create growth curve
#     Q = np.linspace(0,0.35,100)
#     ks = 0.1
#     Cox = 138/16 # mmol O2 : mmol N
#     mu0 = 1.7
#     alpha = 0.07
#     E = 5 # what is a typical value for this????
#     sunlight = mu0 * alpha * E /(np.sqrt(mu0**2 + alpha**2*E**2))
# #     print(sunlight)
#     growth_prediction = (QtoN*Q / (ks + 2*np.sqrt(ks*QtoN*Q) + QtoN*Q)) * (QtoP*Q) * 32 * Cox * sunlight
#     ax[3].plot(Q, growth_prediction, 'orchid')
ax[2].set_ylim([0,0.6])
ax[2].set_xlim([0,0.35])

# QoutDOout vs. QinDOin
ax[3].set_title('(d) QoutDOout vs. QinDOin\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
ax[3].tick_params(axis='x', labelrotation=30)
ax[3].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[3].tick_params(axis='both', labelsize=12)
ax[3].set_xlabel(r'Mean $\frac{Q_{in}DO_{in}}{V_{inlet}}$[mg/L per day]', fontsize=12)
ax[3].set_ylabel(r'Mean $\frac{Q_{out}DO_{out}}{V_{inlet}}$ [mg/L per day]', fontsize=12)
# plot
ax[3].scatter(mean_QinDOin,mean_QoutDOout,s=60, zorder=5, c='navy', alpha=0.5)
#     ax[1].set_ylim([-0.35,0])
# ax[4].set_xlim([0,0.35])
# conduct linear regession
x = mean_QinDOin
res = stats.linregress(x,mean_QoutDOout)
ax[3].plot(x, res.intercept + res.slope*x, 'orchid')
ax[3].text(0.5,0.9,f'R-squared: {res.rvalue**2:.5f}',transform=ax[3].transAxes,
            fontsize=12,fontweight='bold')
ax[3].text(0.5,0.82,f'Slope: {res.slope:.5f}',transform=ax[3].transAxes,
            fontsize=12,fontweight='bold')
ax[3].text(0.5,0.74,f'Intercept: {res.intercept:.5f}',transform=ax[3].transAxes,
            fontsize=12,fontweight='bold')

# air-sea+rivers vs. Photosynthesis
ax[4].set_title('(e) Air-sea+Rivers vs. photosynthesis\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
ax[4].tick_params(axis='x', labelrotation=30)
ax[4].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[4].tick_params(axis='both', labelsize=12)
ax[4].set_ylabel(r'Mean $\frac{Air-sea +Rivers}{V_{inlet}}$ [mg/L per day]', fontsize=12)
ax[4].set_xlabel(r'Mean $\frac{Photosynthesis}{V_{inlet}}$ [mg/L per day]', fontsize=12)
# plot
ax[4].scatter(mean_photo,mean_airsea+mean_rivers,s=60, zorder=5, c='navy', alpha=0.5)
#     ax[1].set_ylim([-0.35,0])
ax[4].set_xlim([0,0.6])
# conduct linear regession
x = mean_photo
res = stats.linregress(x,mean_airsea+mean_rivers)
ax[4].plot(x, res.intercept + res.slope*x, 'orchid')
ax[4].text(0.5,0.9,f'R-squared: {res.rvalue**2:.3f}',transform=ax[4].transAxes,
            fontsize=12,fontweight='bold')
ax[4].text(0.5,0.82,f'Slope: {res.slope:.3f}',transform=ax[4].transAxes,
            fontsize=12,fontweight='bold')
ax[4].text(0.5,0.74,f'Intercept: {res.intercept:.3f}',transform=ax[4].transAxes,
            fontsize=12,fontweight='bold')
# plot expected relationship
photo_values = np.linspace(0,0.6,5)
predicted_values = photo_values*-0.441 + 0.029
ax[4].plot(photo_values,predicted_values,color='darkorange',linestyle=':')

# checking d/dt(DO)-- intercept should be -0.03!!
ax[5].set_title('(f) photo+cons+airsea vs. exchange+TRAPS\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
ax[5].tick_params(axis='x', labelrotation=30)
ax[5].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
ax[5].tick_params(axis='both', labelsize=12)
ax[5].set_ylabel(r'Mean $\frac{photo+cons+airsea}{V_{inlet}}$ [mg/L per day]', fontsize=12)
ax[5].set_xlabel(r'Mean $\frac{exchange+traps}{V_{inlet}}$ [mg/L per day]', fontsize=12)
# plot
# ax[5].scatter(mean_QinDOin+mean_QoutDOout+mean_rivers,mean_photo+mean_cons+mean_airsea,s=60, zorder=5, c='navy', alpha=0.5)
mean_exchange = mean_QinDOin+mean_QoutDOout
ax[5].scatter(mean_QinDOin+mean_QoutDOout+mean_rivers,mean_photo+mean_cons+mean_airsea,s=60, zorder=5, c='navy', alpha=0.5)
# conduct linear regession
x = mean_QinDOin+mean_QoutDOout+mean_rivers
res = stats.linregress(x,mean_photo+mean_cons+mean_airsea)
ax[5].plot(x, res.intercept + res.slope*x, 'orchid')
ax[5].text(0.5,0.9,f'R-squared: {res.rvalue**2:.3f}',transform=ax[5].transAxes,
            fontsize=12,fontweight='bold')
ax[5].text(0.5,0.82,f'Slope: {res.slope:.3f}',transform=ax[5].transAxes,
            fontsize=12,fontweight='bold')
ax[5].text(0.5,0.74,f'Intercept: {res.intercept:.3f}',transform=ax[5].transAxes,
            fontsize=12,fontweight='bold')

# # air-sea vs. exchange flow strength
#     ax[4].set_title('(e) Air-sea vs. exchange flow strength\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
#     ax[4].tick_params(axis='x', labelrotation=30)
#     ax[4].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#     ax[4].tick_params(axis='both', labelsize=12)
#     ax[4].set_ylabel(r'Mean $\frac{Air-sea}{V_{inlet}}$ [mg/L per day]', fontsize=12)
#     ax[4].set_xlabel(r'Mean $\frac{Q_{in}DO_{in}}{V_{inlet}}=\frac{1}{T_{flush}}$ [mg/L per day]', fontsize=12)
#     # plot
#     ax[4].scatter(mean_QinDOin,mean_airsea,s=60, zorder=5, c='navy', alpha=0.5)
# #     ax[1].set_ylim([-0.35,0])
#     # ax[4].set_xlim([0,0.35])
#     # conduct linear regession
#     x = mean_QinDOin
#     res = stats.linregress(mean_QinDOin,mean_airsea)
#     ax[4].plot(x, res.intercept + res.slope*x, 'orchid')
#     ax[4].text(0.5,0.9,f'R-squared: {res.rvalue**2:.3f}',transform=ax[4].transAxes,
#                fontsize=12,fontweight='bold')
#     ax[4].text(0.5,0.82,f'Slope: {res.slope:.3f}',transform=ax[4].transAxes,
#                fontsize=12,fontweight='bold')

#     # Photosynthesis growth rate
#     ax[2].set_title('(c) volume-normalized photosynthesis\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
#     ax[2].tick_params(axis='x', labelrotation=30)
#     ax[2].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#     ax[2].tick_params(axis='both', labelsize=12)
#     ax[2].set_xlabel(r'Mean 1/T$_{flush}$ [day-1]', fontsize=12)
#     ax[2].set_ylabel(r'Mean Photosynthesis [mg/L per day]', fontsize=12)
#     # plot
#     ax[2].scatter(mean_Tflush**-1,mean_photo,s=60, zorder=5, c='navy', alpha=0.5)
#     ax[2].set_ylim([0,0.6])
#     ax[2].set_xlim([0,0.35])

#     # Consumption
#     ax[3].set_title('(d) volume-normalized consumption\nmid-Jun to mid-Aug', size=12, loc='left', fontweight='bold')
#     ax[3].tick_params(axis='x', labelrotation=30)
#     ax[3].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
#     ax[3].tick_params(axis='both', labelsize=12)
#     ax[3].set_xlabel(r'Mean 1/T$_{flush}$ [day-1]', fontsize=12)
#     ax[3].set_ylabel(r'Mean Consumption [mg/L per day]', fontsize=12)
#     # plot
#     ax[3].scatter(mean_Tflush**-1,mean_cons,s=60, zorder=5, c='navy', alpha=0.5)
#     ax[3].set_ylim([-0.35,0])
#     ax[3].set_xlim([0,0.35])

plt.tight_layout()

