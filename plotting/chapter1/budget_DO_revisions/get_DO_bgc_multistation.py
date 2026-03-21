# O2 production from new (NO3) and regenrated production(NH4)
# O2 consumption in water column (nitrification + remi) + sediment (denitrification)

# Original file from Jilian: https://github.com/Jilian0717/LO_user/blob/main/tracer_budget/two_layer/get_DO_bgc_air_sea_shallow_deep.py
# modified by me

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from lo_tools import zfun, zrfun, Lfun
from datetime import datetime, timedelta
import scipy
import pickle, sys
import pandas as pd
from time import time
tt0 = time()

#-----------------------------------------------

# stations = ['lynchcove','penn','budd','case','carr']
# create dictionaries with interface depths
interface_dict = dict()

jobname = 'twentyoneinlets'

burial=50 # 50% burial of sinking detritus in the Salish Sea

interface_types = ['og','drdz','tef','halocline','oxycline'] # how to define dividing depth

#%------------------------------------------------
Ldir = Lfun.Lstart()
# Ldir['roms_out'] = Ldir['roms_out2']
# Ldir['roms_out'] = Ldir['roms_out1']
# Ldir['roms_out'] = Ldir['roms_out5'] # for apogee
Ldir['roms_out'] = Ldir['roms_out'] # testing on local pc
Ldir['gtagex'] = 'cas7_t1_x11b'

year = '2017'

# testing
ds0 = '2017.01.01'
ds1 = '2017.01.02'

# jan
# ds0 = '2017.01.01'
# ds1 = '2017.01.31'
# feb
# ds0 = '2017.02.01'
# ds1 = '2017.02.28'
# mar
# ds0 = '2017.03.01'
# ds1 = '2017.03.31'
# apr
# ds0 = '2017.04.01'
# ds1 = '2017.04.30'
# may
# ds0 = '2017.05.01'
# ds1 = '2017.05.31'
# jun
# ds0 = '2017.06.01'
# ds1 = '2017.06.30'

# jul (running)
# ds0 = '2017.07.01'
# ds1 = '2017.07.31'
# aug (running)
# ds0 = '2017.08.01'
# ds1 = '2017.08.31'
# sep (running)
# ds0 = '2017.09.01'
# ds1 = '2017.09.30'
# oct (running)
# ds0 = '2017.10.01'
# ds1 = '2017.10.31'
# nov (running)
# ds0 = '2017.11.01'
# ds1 = '2017.11.30'
# dec
ds0 = '2017.12.01'
ds1 = '2017.12.31'

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / 'budget_revisons' / ('DO_budget_'+year+'.01.01_'+year+'.12.31') / '2layer_bgc'
Lfun.make_dir(out_dir)


Ldir['ds0'] = ds0
in_dir = Ldir['roms_out'] / Ldir['gtagex']
G, S, T = zrfun.get_basic_info(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')

fn0 = xr.open_dataset(in_dir / ('f' + Ldir['ds0']) / 'ocean_his_0002.nc')
dx = 1/fn0.pm.values
dy = 1/fn0.pn.values
lonr = fn0.lon_rho.values
latr = fn0.lat_rho.values
area = dx * dy
NX, NY = dx.shape

AttSW = np.zeros((NX, NY)) + 0.05
AttSW[(lonr > -123.89) & (latr < 50.29) & (latr > 47.02)] = 0.15
AttSW[(lonr > -125.31) & (lonr < -123.89) & (latr < 51.02) & (latr > 49.13)] = 0.15

dt0 = datetime.strptime(ds0, Lfun.ds_fmt)
dt1 = datetime.strptime(ds1, Lfun.ds_fmt)
dt00 = dt0

#% air-sea flux related parameters
dtdays = 3600 * (1.0/86400) / 1
cff2_air = dtdays * 0.31 * 24 / 100
# formulation for Schmidt number coefficient
A_O2 = 1953.4; B_O2 = 128.0; C_O2 = 3.9918
D_O2 = 0.050091; E_O2 = 0.0
# Calculate O2 saturation concentration using Garcia and Gordon
#  L and O (1992) formula, (EXP(AA) is in ml/l).
OA0 = 2.00907       # Oxygen saturation coefficients
OA1 = 3.22014;      OA2 = 4.05010;       OA3 = 4.94457
OA4 = -0.256847;    OA5 = 3.88767;       OB0 = -0.00624523
OB1 = -0.00737614;  OB2 = -0.0103410;    OB3 = -0.00817083
OC0 = -0.000000488682    

#%%

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
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

# create dictionary of empty dataframes
# df_dict = {'lynchcove': pd.DataFrame(),
#            'penn': pd.DataFrame(),
#         #    'budd': pd.DataFrame(),
#            'carr': pd.DataFrame(),
#            'case': pd.DataFrame(),
#            'commencement': pd.DataFrame(),
#            'crescent': pd.DataFrame(),
#            'dabob': pd.DataFrame(),
#            'dyes': pd.DataFrame(),
#         #    'eld': pd.DataFrame(),
#            'elliot': pd.DataFrame(),
#         #    'hammersley': pd.DataFrame(),
#         #    'henderson': pd.DataFrame(),
#            'holmes': pd.DataFrame(),
#         #    'killsut': pd.DataFrame(),
#         #    'oak': pd.DataFrame(),
#            'portsusan': pd.DataFrame(),
#            'quartermaster': pd.DataFrame(),
#         #    'similk': pd.DataFrame(),
#            'sinclair': pd.DataFrame()}
#         #    'totten': pd.DataFrame()}

df_dict = {'lynchcove_og': pd.DataFrame(),
           'lynchcove_drdz': pd.DataFrame(),
           'lynchcove_tef': pd.DataFrame(),
           'lynchcove_halocline': pd.DataFrame(),
           'lynchcove_oxycline': pd.DataFrame(),

           'penn_og': pd.DataFrame(),
           'penn_drdz': pd.DataFrame(),
           'penn_tef': pd.DataFrame(),
           'penn_halocline': pd.DataFrame(),
           'penn_oxycline': pd.DataFrame(),

           'carr_og': pd.DataFrame(),
           'carr_drdz': pd.DataFrame(),
           'carr_tef': pd.DataFrame(),
           'carr_halocline': pd.DataFrame(),
           'carr_oxycline': pd.DataFrame(),

           'case_og': pd.DataFrame(),
           'case_drdz': pd.DataFrame(),
           'case_tef': pd.DataFrame(),
           'case_halocline': pd.DataFrame(),
           'case_oxycline': pd.DataFrame(),

           'commencement_og': pd.DataFrame(),
           'commencement_drdz': pd.DataFrame(),
           'commencement_tef': pd.DataFrame(),
           'commencement_halocline': pd.DataFrame(),
           'commencement_oxycline': pd.DataFrame(),

           'crescent_og': pd.DataFrame(),
           'crescent_drdz': pd.DataFrame(),
           'crescent_tef': pd.DataFrame(),
           'crescent_halocline': pd.DataFrame(),
           'crescent_oxycline': pd.DataFrame(),

           'dabob_og': pd.DataFrame(),
           'dabob_drdz': pd.DataFrame(),
           'dabob_tef': pd.DataFrame(),
           'dabob_halocline': pd.DataFrame(),
           'dabob_oxycline': pd.DataFrame(),

           'dyes_og': pd.DataFrame(),
           'dyes_drdz': pd.DataFrame(),
           'dyes_tef': pd.DataFrame(),
           'dyes_halocline': pd.DataFrame(),
           'dyes_oxycline': pd.DataFrame(),

           'elliot_og': pd.DataFrame(),
           'elliot_drdz': pd.DataFrame(),
           'elliot_tef': pd.DataFrame(),
           'elliot_halocline': pd.DataFrame(),
           'elliot_oxycline': pd.DataFrame(),

           'holmes_og': pd.DataFrame(),
           'holmes_drdz': pd.DataFrame(),
           'holmes_tef': pd.DataFrame(),
           'holmes_halocline': pd.DataFrame(),
           'holmes_oxycline': pd.DataFrame(),

           'portsusan_og': pd.DataFrame(),
           'portsusan_drdz': pd.DataFrame(),
           'portsusan_tef': pd.DataFrame(),
           'portsusan_halocline': pd.DataFrame(),
           'portsusan_oxycline': pd.DataFrame(),

           'quartermaster_og': pd.DataFrame(),
           'quartermaster_drdz': pd.DataFrame(),
           'quartermaster_tef': pd.DataFrame(),
           'quartermaster_halocline': pd.DataFrame(),
           'quartermaster_oxycline': pd.DataFrame(),

           'sinclair_og': pd.DataFrame(),
           'sinclair_drdz': pd.DataFrame(),
           'sinclair_tef': pd.DataFrame(),
           'sinclair_halocline': pd.DataFrame(),
           'sinclair_oxycline': pd.DataFrame()}

cnt = 0
#%%
while dt00 <= dt1:  # loop each day and every history file
    print(dt00)
    sys.stdout.flush()
    ds00 = dt00.strftime(Lfun.ds_fmt)
    fn_list = Lfun.get_fn_list('hourly', Ldir, ds00, ds00)
    #%%
    for hr,fn in enumerate(fn_list[0:-1]): 
        print('    hour {}'.format(hr+1))
        sys.stdout.flush()
        # print(fn)

        ds = xr.open_dataset(fn)
        swrad = ds.swrad.values.squeeze()
        chl = ds.chlorophyll.values.squeeze()
        zeta = ds.zeta.values.squeeze()
        h = ds.h.values
        zw = zrfun.get_z(h, zeta, S, only_w=True)
        zrho = zrfun.get_z(h, zeta, S, only_rho=True)
        dz = np.diff(zw, axis=0)
        salt = ds.salt.values.squeeze()
        NH4 = ds.NH4.values.squeeze()
        NO3 = ds.NO3.values.squeeze()
        phy = ds.phytoplankton.values.squeeze()
        SDeN = ds.SdetritusN.values.squeeze()
        LDeN = ds.LdetritusN.values.squeeze()
        Oxy = ds.oxygen.values.squeeze()

       # PARsur = PARfrac * swrad #* rho0 * Cp # surface PAR, watts/m2
        Att = np.zeros(salt.shape)
        #ExpAtt = np.zeros(Att.shape)
        nk,ni,nj = Att.shape
        PAR = np.zeros([nk+1, ni, nj])
        #PAR[-1,:,:] = PARsur  
        PAR[-1,:,:] = 0.43 * swrad  
#
        Oxy_pro = np.zeros(Att.shape) # O2 production
        Oxy_nitri = np.zeros(Att.shape) # O2 consumption by nitrification in water column
        Oxy_remi = np.zeros(Att.shape) # O2 consumption by remineralization in water column
        Oxy_sed = np.zeros([ni,nj]) # O2 consumption by SOD
        
        z_w = zrfun.get_z(h, zeta, S, only_rho=False, only_w=True)
        vol = np.diff(z_w,axis=0) * area # grid cell volume
        
        #if np.nanmin(PARsur) > 0: # doing photosysnthesis
        if np.nanmin(0.43*swrad) > 0: # doing photosysnthesis
            for k in np.arange(29,-1,-1): # surface to bottom, i.e. water column                
                # O2 production from photosynthesis, mmol O2/m3/hr ??????????????
                Att[k,:,:] = (AttSW + 0.012*chl[k,:,:] - 0.0065*(salt[k,:,:]-32))* (z_w[k+1,:,:] - z_w[k,:,:])                       
                PAR[k,:,:] = PAR[k+1,:,:] * (1-np.exp(-Att[k,:,:]))/Att[k,:,:] # average at cell center
                fac1 = PAR[k,:,:] * 0.07 
                Epp = 1.7/np.sqrt(1.7*1.7 + fac1*fac1)
                t_PPmax = Epp*fac1
                cff1 = NH4[k,:,:] * 10; cff1[cff1<0] = 0
                cff2 = NO3[k,:,:] * 10; cff2[cff2<0] = 0

                inhNH4 = 1.0/(1.0+cff1)


                cff4 = dtdays * t_PPmax * 10 * inhNH4/(1.0 + cff2 + 2.0*np.sqrt(cff2)) * phy[k,:,:]
                cff5 = dtdays * t_PPmax * 10 / (1.0 + cff1 + 2.0 * np.sqrt(cff1)) * phy[k,:,:]
                       
                #Oxy_pro[k,:,:] = N_Flux_NewProd*rOxNO3 + N_Flux_RegProd*rOxNH4
                Oxy_pro[k,:,:] = NO3[k,:,:] * cff4 * 138/16 + NH4[k,:,:] * cff5 * 106/16 

                # O2 comsumption by nitrification in water column, mmol O2/m3/hr ??????????????            
                fac2 = Oxy[k,:,:]
                fac2[fac2<0] = 0
                fac3 = fac2/(3+fac2)
                fac3[fac3<0] = 0
                Oxy_nitri[k,:,:] = 2.0 * NH4[k,:,:] * dtdays * 0.05 * fac3
              
                PAR[k,:,:] = PAR[k+1,:,:] * np.exp(-Att[k,:,:])  # !!! light attenuation at the bottom of grid cell !!!!
            
        else: # PARsur = 0, nitrification occurs at the maximum rate (NitriR)
            for k in np.arange(29,-1,-1):
                Oxy_nitri[k,:,:] = 2.0 * NH4[k,:,:] * dtdays * 0.05
            
        
        #O2 consumption by remineralization in water column, mmol O2/m3/hr ??????????????
        for k in np.arange(0,30):
            Oxy_remi[k,:,:] = (SDeN[k,:,:]*dtdays*0.1*1 + LDeN[k,:,:]*dtdays*0.1*1) * 106/16

        # mmol O2/m3/hr to mmol O2/hr
        Oxy_pro = Oxy_pro * vol 
        Oxy_nitri = Oxy_nitri * vol # 
        Oxy_remi = Oxy_remi * vol #   

        #---------- sediment SOD: a more strict way ----------
        # LdetritusN decomposition in sediment
        cff1_L = (LDeN[0,:,:] *(1-burial/100) * 80)/24 / dz[0,:,:] * 1  # mmol/m2/hr --> mmol/m3     
        NH4_gain_L = np.zeros([NX,NY])        
        for i in range(NX): # loop the whole domain
            for j in range(NY):
                if cff1_L[i,j]*106/16 <= Oxy[0,:,:][i,j]:  # use O2: PN + O2 --> NH4; Oxy[0,:,:] is bottom oxygen                                        
                    NH4_gain_L[i,j] = cff1_L[i,j] #mmol/m3, NH4 gain from decomposing LdetN                                                                             
        NH4_gain_flux_L = NH4_gain_L * area * dz[0,:,:] # mmol but the actual unit is mmol/hr
        # SdetritusN decomposition in sediment
        cff1_S = (SDeN[0,:,:] *(1-burial/100) * 8)/24 / dz[0,:,:] * 1  # mmol/m2/hr --> mmol/m3     
        NH4_gain_S = np.zeros([NX,NY])        
        for i in range(NX): # loop the whole domain
            for j in range(NY):
                if cff1_S[i,j]*106/16 <= Oxy[0,:,:][i,j]:  # use O2: PN + O2 --> NH4; Oxy[0,:,:] is bottom oxygen                                        
                    NH4_gain_S[i,j] = cff1_S[i,j] #mmol/m3, NH4 gain from decomposing LdetN                                                                             
        NH4_gain_flux_S = NH4_gain_S * area * dz[0,:,:] # mmol but the actual unit is mmol/hr
        
        NH4_gain_flux = NH4_gain_flux_L + NH4_gain_flux_S
                
        #---------- air-sea flux ----------
        Uwind = ds.Uwind.values.squeeze()
        Vwind = ds.Vwind.values.squeeze()
        temp_surf = ds.temp.values[0,-1,:,:] # surface temp
        salt_surf = ds.salt.values[0,-1,:,:] # surface salt
        Oxy_surf = ds.oxygen.values[0,-1,:,:] # surface O2
        # Compute O2 transfer velocity: u10squared (u10 in m/s)
        u10squ = Uwind * Uwind + Vwind * Vwind  # ifdef BULK_FLUXES
        
        SchmidtN_Ox = A_O2 - temp_surf*(B_O2 - temp_surf*(C_O2 - temp_surf*(D_O2 - temp_surf*E_O2)))

        cff3 = cff2_air * u10squ * np.sqrt(660.0/SchmidtN_Ox)  # m
        TS = np.log((298.15-temp_surf)/(273.15+temp_surf))
        AA = OA0 + TS*(OA1+TS*(OA2+TS*(OA3+TS*(OA4+TS*OA5)))) + salt_surf*(OB0+TS*(OB1+TS*(OB2+TS*OB3))) + OC0*salt_surf*salt_surf
         # Convert from ml/l to mmol/m3
        O2satu = 1000./22.3916 * np.exp(AA) # mmol/m3

        ###################################################################
        ## GET VALUES IN EACH TERMINAL INLET AND SAVE IN INDIVIDUAL FILE ##
        ###################################################################

        for type in interface_types:
            for station in sta_dict: # stations: 

                # get interface depth from csv file
                with open('interface_depths_' + type + '.csv', 'r') as f:
                    for line in f:
                        inlet, interface_depth = line.strip().split(',')
                        interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
                z_interface = float(interface_dict[station])

                # get segment information
                seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
                seg_df = pd.read_pickle(seg_name)
                ji_list = seg_df[station+'_p']['ji_list']
                jj = [x[0] for x in ji_list]
                ii = [x[1] for x in ji_list]

                # get storage term
                tmp_zrho = zrho[:,jj,ii] # in domain
                ix_shallow = tmp_zrho>=z_interface # shallower than interface
                ix_deep = tmp_zrho<z_interface  # deeper than interface
                tmp_DOV = Oxy[:,jj,ii] * vol[:,jj,ii]
                ret_DOvol_surf = np.nansum(tmp_DOV[ix_shallow])
                ret_DOvol_deep = np.nansum(tmp_DOV[ix_deep])

                # get volume term
                tmp_V = vol[:,jj,ii]
                ret_vol_surf = np.nansum(tmp_V[ix_shallow])
                ret_vol_deep = np.nansum(tmp_V[ix_deep])

                # get photosynthesis, nitrification, and respiration terms
                ret_photo_surf = np.nansum(Oxy_pro[:,jj,ii][ix_shallow])
                ret_photo_deep = np.nansum(Oxy_pro[:,jj,ii][ix_deep])
                ret_nitri_surf = np.nansum(Oxy_nitri[:,jj,ii][ix_shallow])
                ret_nitri_deep = np.nansum(Oxy_nitri[:,jj,ii][ix_deep])
                ret_respi_surf = np.nansum(Oxy_remi[:,jj,ii][ix_shallow])
                ret_respi_deep = np.nansum(Oxy_remi[:,jj,ii][ix_deep])

                # get sediment oxygen demand 
                ret_sod = np.nansum(NH4_gain_flux[jj,ii] * 106/16)    # mmol O2/hr  

                # get O2 gas exchange
                ret_airsea = np.nansum(cff3[jj,ii] * (O2satu[jj,ii]-Oxy_surf[jj,ii]) * area[jj,ii])

                # get dataframe for saving
                if cnt == 0:
                    # start data
                    df_dict[station+'_'+type]['surf DO*V [mmol]'] = [ret_DOvol_surf]
                    df_dict[station+'_'+type]['deep DO*V [mmol]'] = [ret_DOvol_deep]
                    df_dict[station+'_'+type]['surf vol [m3]'] = [ret_vol_surf]
                    df_dict[station+'_'+type]['deep vol [m3]'] = [ret_vol_deep]
                    df_dict[station+'_'+type]['surf photo [mmol/hr]'] = [ret_photo_surf]
                    df_dict[station+'_'+type]['deep photo [mmol/hr]'] = [ret_photo_deep]
                    df_dict[station+'_'+type]['surf nitri [mmol/hr]'] = [ret_nitri_surf]
                    df_dict[station+'_'+type]['deep nitri [mmol/hr]'] = [ret_nitri_deep]
                    df_dict[station+'_'+type]['surf respi [mmol/hr]'] = [ret_respi_surf]
                    df_dict[station+'_'+type]['deep respi [mmol/hr]'] = [ret_respi_deep]
                    df_dict[station+'_'+type]['SOD [mmol/hr]'] = [ret_sod]
                    df_dict[station+'_'+type]['airsea [mmol/hr]'] = [ret_airsea]
                else:
                    # get temp dataframe
                    df_tmp = pd.DataFrame()
                    df_tmp['surf DO*V [mmol]'] = [ret_DOvol_surf]
                    df_tmp['deep DO*V [mmol]'] = [ret_DOvol_deep]
                    df_tmp['surf vol [m3]'] = [ret_vol_surf]
                    df_tmp['deep vol [m3]'] = [ret_vol_deep]
                    df_tmp['surf photo [mmol/hr]'] = [ret_photo_surf]
                    df_tmp['deep photo [mmol/hr]'] = [ret_photo_deep]
                    df_tmp['surf nitri [mmol/hr]'] = [ret_nitri_surf]
                    df_tmp['deep nitri [mmol/hr]'] = [ret_nitri_deep]
                    df_tmp['surf respi [mmol/hr]'] = [ret_respi_surf]
                    df_tmp['deep respi [mmol/hr]'] = [ret_respi_deep]
                    df_tmp['SOD [mmol/hr]'] = [ret_sod]
                    df_tmp['airsea [mmol/hr]'] = [ret_airsea]
                    # append data
                    df_dict[station+'_'+type] = pd.concat([df_dict[station+'_'+type],df_tmp])
                    # reset index
                    df_dict[station+'_'+type].reset_index(drop=True, inplace=True)

        cnt += 1
        ds.close()       
    dt00 = dt00 + timedelta(days=1)

for type in interface_types:
    for station in sta_dict: #for station in stations:
        # get dataframe for saving
        df = df_dict[station+'_'+type]
        # save to pickle file
        Lfun.make_dir(out_dir/station)
        df.to_pickle(out_dir / station / (type + '_' + station + '_' + ds0 + '_' + ds1 + '.p'))
print('Done')