"""
Plot surface and deep budget for all 21 inlets
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import get_two_layer
import matplotlib.dates as mdates

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
box_fn = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
ds_box = xr.open_dataset(box_fn)
DX = (ds_box.pm.values)**-1
DY = (ds_box.pn.values)**-1
DA = DX*DY # get area of each grid cell in m^2

# COLLAPSE
for i,station in enumerate(sta_dict):
        
        # initialize figure
        fig,axes = plt.subplots(2,1, figsize=(10,7), sharex=True)
        ax = axes.ravel()

        plt.suptitle(station + '(10-day Hanning Window)',fontsize=14)

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
        ax[0].plot(dates_local_daily,zfun.lowpass(TEF_surf,n=10),color='#0D4B91',label='TEF')
        ax[1].plot(dates_local_daily,zfun.lowpass(TEF_deep,n=10),color='#0D4B91',label='TEF')
        # format budget time series figures
        for axnum in [ax[0],ax[1]]:
            axnum.set_xlim([dates_hrly[0],dates_hrly[-2]])
            axnum.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')
            axnum.tick_params(axis='x', labelrotation=30, labelsize=12)
            axnum.tick_params(axis='y', labelsize=12)
            axnum.set_ylabel(r'DO transport [kmol O$_2$ s$^{-1}$]', fontsize=12)
            loc = mdates.MonthLocator(interval=1)
            axnum.xaxis.set_major_locator(loc)
        ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

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
        

        interface_types = ['og','drdz','tef','halocline','oxycline'] # how to define dividing depth
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

            # take time derivative of (DO*V) to get d/dt (DO*V)
            ddtDOV_surf_unfiltered = np.diff(o2vol_surf_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
            ddtDOV_deep_unfiltered = np.diff(o2vol_deep_unfiltered) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s

            # apply Godin filter
            photo_surf  = zfun.lowpass(photo_surf_unfiltered, f='godin')[36:-34:24]
            photo_deep  = zfun.lowpass(photo_deep_unfiltered, f='godin')[36:-34:24]
            cons_surf   = zfun.lowpass(cons_surf_unfiltered, f='godin')[36:-34:24]
            cons_deep   = zfun.lowpass(cons_deep_unfiltered, f='godin')[36:-34:24]
            airsea_surf = zfun.lowpass(airsea_surf_unfiltered, f='godin')[36:-34:24]
            ddtDOV_surf = zfun.lowpass(ddtDOV_surf_unfiltered, f='godin')[36:-34:24]
            ddtDOV_deep = zfun.lowpass(ddtDOV_deep_unfiltered, f='godin')[36:-34:24]

            # plot budget time series
            if t == 0:
                ax[0].plot(dates_local_daily,zfun.lowpass(photo_surf,n=10),color='#8F0445',alpha=0.7,label='Photosynthesis')
                ax[1].plot(dates_local_daily,zfun.lowpass(photo_deep,n=10),color='#8F0445',alpha=0.7,label='Photosynthesis')

                ax[0].plot(dates_local_daily,zfun.lowpass(cons_surf,n=10),color='#FCC2DD',alpha=0.7,label='Consumption')
                ax[1].plot(dates_local_daily,zfun.lowpass(cons_deep,n=10),color='#FCC2DD',alpha=0.7,label='Consumption')

                ax[0].plot(dates_local_daily,zfun.lowpass(ddtDOV_surf,n=10),color='black',alpha=0.7,label='d/dt(DO)')
                ax[1].plot(dates_local_daily,zfun.lowpass(ddtDOV_deep,n=10),color='black',alpha=0.7,label='d/dt(DO)')

                ax[0].plot(dates_local_daily,zfun.lowpass(airsea_surf,n=10),color='yellowgreen',label='Air-Sea')
            else:
                ax[0].plot(dates_local_daily,zfun.lowpass(photo_surf,n=10),color='#8F0445',alpha=0.7)
                ax[1].plot(dates_local_daily,zfun.lowpass(photo_deep,n=10),color='#8F0445',alpha=0.7)

                ax[0].plot(dates_local_daily,zfun.lowpass(cons_surf,n=10),color='#FCC2DD',alpha=0.7)
                ax[1].plot(dates_local_daily,zfun.lowpass(cons_deep,n=10),color='#FCC2DD',alpha=0.7)

                ax[0].plot(dates_local_daily,zfun.lowpass(ddtDOV_surf,n=10),color='black',alpha=0.7)
                ax[1].plot(dates_local_daily,zfun.lowpass(ddtDOV_deep,n=10),color='black',alpha=0.7)

                ax[0].plot(dates_local_daily,zfun.lowpass(airsea_surf,n=10),color='yellowgreen') 
            ax[0].legend(loc='lower right',ncol=6, fontsize=12, handletextpad=0.15)

            # # plot box plots
            # boxprops = dict(facecolor='black', alpha=0.5)
            # medianprops = dict(color='silver', linewidth=1)
            # ddt_bplot = ax[2].boxplot(ddtDOV_deep[~np.isnan(ddtDOV_deep)], patch_artist=True, widths=0.5,
            #                           positions=[t], boxprops=boxprops, medianprops=medianprops)
            # ax[2].text(t,np.nanmax(ddtDOV_deep),str(round(np.nanmean(ddtDOV_deep),3)),ha='center',va='bottom',
            #            fontweight='bold')

            # boxprops = dict(facecolor='#FCC2DD', alpha=0.5)
            # medianprops = dict(color='hotpink', linewidth=1)
            # ddt_bplot = ax[3].boxplot(cons_deep[~np.isnan(cons_deep)], patch_artist=True, widths=0.5,
            #                           positions=[t], boxprops=boxprops, medianprops=medianprops)
            # ax[3].text(t,np.nanmax(cons_deep),str(round(np.nanmean(cons_deep),3)),ha='center',va='bottom',
            #            fontweight='bold')

            # boxprops = dict(facecolor='#8F0445', alpha=0.5)
            # medianprops = dict(color='pink', linewidth=1)
            # ddt_bplot = ax[4].boxplot(photo_deep[~np.isnan(photo_deep)], patch_artist=True, widths=0.5,
            #                           positions=[t], boxprops=boxprops, medianprops=medianprops)
            # ax[4].text(t,np.nanmax(photo_deep),str(round(np.nanmean(photo_deep),3)),ha='center',va='bottom',
            #            fontweight='bold')
            # ax[4].set_xticks([0,1,2,3,4], interface_types, rotation=45, fontsize=12)

# # ------------------------------- get rivers and WWTPs ----------------------------------------
#         fn = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('DO_budget_' + startdate + '_' + enddate) / '2layer_traps' / (station + '.p')
#         df_traps = pd.read_pickle(fn)
#         rivers_surf_unfiltered = df_traps['surface [kmol/s]'].values
#         wwtps_deep_unfiltered = df_traps['deep [kmol/s]'].values

#         # Pad with zeros if no rivers or WWTPs
#         if rivers_surf_unfiltered.size == 0 or np.isnan(rivers_surf_unfiltered[0]):
#                 rivers_surf_unfiltered = np.zeros(8761)
#         if wwtps_deep_unfiltered.size == 0:
#                 wwtps_deep_unfiltered = np.zeros(8761)
#         # Godin filter
#         traps_surf = zfun.lowpass(rivers_surf_unfiltered, f='godin')[36:-34:24]
#         traps_deep =  zfun.lowpass(wwtps_deep_unfiltered, f='godin')[36:-34:24]
#         traps_color = 'black'

# # ------------------------------- get vertical exchange ----------------------------------------

#         # pick color
#         vertX_color = 'mediumorchid'

#         # calculate error term, and ascribe to the vertical exchange
#         vertX_surf_TEF = ddtDOV_surf - (TEF_surf 
#                                     + photo_surf
#                                     + cons_surf
#                                     + airsea_surf
#                                     + traps_surf)
        
#         vertX_deep_TEF = ddtDOV_deep - (TEF_deep 
#                                     + photo_deep
#                                     + cons_deep
#                                     + traps_deep)
        
# # ---------------------------------- DO concentrations --------------------------------------------

#         # get average V of each layer
#         fn = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
#         df_V = pd.read_pickle(fn)
#         # Godin filter already applied earlier in workflow
#         surf_V = df_V['surface [m3]'].values[1:-1]
#         deep_V = df_V['deep [m3]'].values[1:-1]

#         # apply Godin filter to DO*V (mmol)
#         o2vol_surf = zfun.lowpass(o2vol_surf_unfiltered, f='godin')[36:-34:24]
#         o2vol_deep = zfun.lowpass(o2vol_deep_unfiltered, f='godin')[36:-34:24]

#         # get average DO of each layer by calculating (DO*V)/V
#         # and convert from mmol/m3 to mg/L
#         o2_surf = o2vol_surf / surf_V * 32/1000
#         o2_deep = o2vol_deep / deep_V * 32/1000

#         # get average inflowing DO concentration (Qin concentration)
#         DO_in = DO_p.values * 32/1000
        
#         #% Get indices of inlet
#         seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
#         seg_df = pd.read_pickle(seg_name)
#         ji_list = seg_df[station+'_p']['ji_list']
#         jj_LO = [x[0] for x in ji_list] # y; lat; jj
#         ii_LO = [x[1] for x in ji_list] # x; lon; ii
#         # get lat and lon corresponding to ii and jj indices
#         lat_LO = latr[jj_LO,0]
#         lon_LO = lonr[0,ii_LO]
#         # get corresponding ii and jj indices in box extraction
#         lat_box_all = ds_box['lat_rho'].values[:,0]
#         lon_box_all = ds_box['lon_rho'].values[0,:]
#         jj = np.zeros(len(jj_LO))
#         ii = np.zeros(len(ii_LO))
#         for j_ind,lat in enumerate(lat_LO):
#             jj[j_ind] = np.where(lat_box_all==lat)[0][0]
#         for i_ind,lon in enumerate(lon_LO):
#             ii[i_ind] = np.where(lon_box_all==lon)[0][0]
#         # convert to array of ints
#         jj = jj.astype(int)
#         ii = ii.astype(int)


# # -------------------------- get volume and mean depth of layers -------------------------

#         # get volume of layers
#         fn = Ldir['LOo'] / 'pugetsound_DO' / ('VOLUME_budget_' + startdate + '_' + enddate) / '2layer_volume_storage' / (station + '.p')
#         df_V = pd.read_pickle(fn)
#         # Godin filter already applied earlier in workflow
#         surf_V = df_V['surface [m3]'].values[1:-1]
#         deep_V = df_V['deep [m3]'].values[1:-1]

# # --------------------------- get other physical dimensions ------------------------------

#         # get full volume of inlet, and mean depth
#         fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'vol_df_cas7_c21.p'
#         vol_df = pd.read_pickle(fn)
#         inlet_vol = vol_df['volume m3'].loc[station+'_p']
#         inlet_area = vol_df['area m2'].loc[station+'_p']
#         mean_depth = inlet_vol / inlet_area
#         # estimate aspect ratio from inlet area and width at mouth
#         # aspect ratio = length/width (so longer inlet has larger ratio)
#         # thus, aspect ratio = area / width^2
#         # since this is an approximation, we use the nominal cell width of 500 m
#         fn =  Ldir['LOo'] / 'extract' / 'tef2' / 'sect_df_cas7_c21.p'
#         mouth_df = pd.read_pickle(fn)
#         mouth_df = mouth_df.loc[mouth_df['sn'] == station]
#         # get number of horizontal and vertical entrance points
#         uv_list = mouth_df['uv'].values.tolist()
#         u_count = uv_list.count('u')
#         v_count = uv_list.count('v')
#         mouth_width = 500 * np.sqrt(u_count**2 + v_count**2)
#         est_aspect_ratio = inlet_area / mouth_width**2
#         # get mouth area, for calculating fluxes 

        
# # ------------------------- save data in dataframe dict -----------------------------------
#         # Note: everything is in units of kmol O2 /s

#         # surface layer
#         surfacelay_dict[station]['TEF'] = TEF_surf
#         surfacelay_dict[station]['Physical Resupply'] = TEF_surf + vertX_surf_TEF
#         surfacelay_dict[station]['TRAPS'] = traps_surf
#         surfacelay_dict[station]['Photosynthesis'] = photo_surf
#         surfacelay_dict[station]['Consumption'] = cons_surf
#         surfacelay_dict[station]['Air-Sea Transfer'] = airsea_surf
#         surfacelay_dict[station]['Vertical'] = vertX_surf_TEF
#         surfacelay_dict[station]['Storage'] = ddtDOV_surf
#         surfacelay_dict[station]['Volume'] = surf_V

#         # bottom layer
#         bottomlay_dict[station]['Storage'] = ddtDOV_deep
#         bottomlay_dict[station]['TEF'] = TEF_deep
#         bottomlay_dict[station]['Physical Resupply'] = TEF_deep + vertX_deep_TEF
#         bottomlay_dict[station]['Vertical'] = vertX_deep_TEF
#         bottomlay_dict[station]['TRAPS'] = traps_deep
#         bottomlay_dict[station]['Photosynthesis'] = photo_deep
#         bottomlay_dict[station]['Consumption'] = cons_deep
#         bottomlay_dict[station]['NEM'] = photo_deep + cons_deep
#         bottomlay_dict[station]['Volume'] = deep_V
#         bottomlay_dict[station]['Qin m3/s'] = Q_p.values 

# # ------------------------- save DO concentrations in dataframe dict -----------------------------------

#         # mg/L units

#         # DO concentrations
#         DOconcen_dict[station]['Surface Layer'] = o2_surf
#         DOconcen_dict[station]['Deep Layer'] = o2_deep
#         # DOconcen_dict[station]['Bottom Sigma DO'] = bottom_oxygen
#         # DOconcen_dict[station]['Minimum Bottom Layer DO'] = oxygen_min
#         DOconcen_dict[station]['Qin DO'] = DO_in
#         # DOconcen_dict[station]['percent hypoxic volume'] = hyp_vol/[inlet_vol] * 100

# # ------------------------ save inlet dimensions in dataframe dict ----------------------------------------

#         dimensions_dict[station]['Inlet volume'] = [inlet_vol] # m3
#         dimensions_dict[station]['Mean depth'] = [mean_depth] # m
#         dimensions_dict[station]['Mouth width'] = mouth_width
#         dimensions_dict[station]['L/W aspect ratio'] = [est_aspect_ratio] # dimensionless

#         # print keys
#         if i == 0:
#             print(list(surfacelay_dict[station].keys()))
#             print(list(DOconcen_dict[station].keys()))
#             print(list(dimensions_dict[station].keys()))


        plt.tight_layout()
        plt.show()

