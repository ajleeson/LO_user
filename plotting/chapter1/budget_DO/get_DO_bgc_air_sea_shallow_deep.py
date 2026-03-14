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
station = 'budd'
z_interface = -6

#%------------------------------------------------
Ldir = Lfun.Lstart()
# Ldir['roms_out'] = Ldir['roms_out2']
# Ldir['roms_out'] = Ldir['roms_out1']
Ldir['roms_out'] = Ldir['roms_out5'] # for apogee
Ldir['gtagex'] = 'cas7_t0_x4b'

# jan
# ds0 = '2014.01.01'
# ds1 = '2014.01.31'
# feb
# ds0 = '2014.02.01'
# ds1 = '2014.02.28'
# mar
# ds0 = '2014.03.01'
# ds1 = '2014.03.31'
# apr
# ds0 = '2014.04.01'
# ds1 = '2014.04.30'
# may
# ds0 = '2014.05.01'
# ds1 = '2014.05.31'
# jun
ds0 = '2014.06.01'
ds1 = '2014.06.30'


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

#% load salish sea j,i
seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
seg_df = pd.read_pickle(seg_name)
ji_list = seg_df[station+'_p']['ji_list']
jj = [x[0] for x in ji_list]
ii = [x[1] for x in ji_list]

inDomain = np.zeros([NX, NY])
inDomain[jj,ii] = 1 # inside domain, the index=1, outside domian, the index =0, this step is to reduce the size of spatial variable (hopefully)

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

Oxy_sed_sum = []
Oxy_sed_sum2 = []
Oxy_pro_sum_shallow = []
Oxy_pro_sum_deep = []
Oxy_nitri_sum_shallow = []
Oxy_nitri_sum_deep = []
Oxy_remi_sum_shallow = []
Oxy_remi_sum_deep = []
Oxy_vol_sum_shallow = []
Oxy_vol_sum_deep = []
Oxy_air_flux_sum = []
t = []

cnt = 0
#%%
while dt00 <= dt1:  # loop each day and every history file
    print(dt00)
    sys.stdout.flush()
    ds00 = dt00.strftime(Lfun.ds_fmt)
    fn_list = Lfun.get_fn_list('hourly', Ldir, ds00, ds00)
    #%%
    for fn in fn_list[0:-1]: 
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
        
        tmp_zrho = zrho[:,jj,ii] # in domain
        ix_shallow = tmp_zrho>=z_interface # shallower than interface
        ix_deep = tmp_zrho<z_interface  # deeper than interface
        
        tmp_DOV = Oxy[:,jj,ii] * vol[:,jj,ii]
        Oxy_vol_sum_shallow.append(np.nansum(tmp_DOV[ix_shallow]))
        Oxy_vol_sum_deep.append(   np.nansum(tmp_DOV[ix_deep]))
        
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
        
        Oxy_pro_sum_shallow.append(np.nansum(Oxy_pro[:,jj,ii][ix_shallow]))
        Oxy_pro_sum_deep.append(np.nansum(Oxy_pro[:,jj,ii][ix_deep]))
        Oxy_nitri_sum_shallow.append(np.nansum(Oxy_nitri[:,jj,ii][ix_shallow]))
        Oxy_nitri_sum_deep.append(np.nansum(Oxy_nitri[:,jj,ii][ix_deep]))
        Oxy_remi_sum_shallow.append(np.nansum(Oxy_remi[:,jj,ii][ix_shallow]))
        Oxy_remi_sum_deep.append(np.nansum(Oxy_remi[:,jj,ii][ix_deep]))
        
        
        #---------- sediment SOD ----------
        # refer to Parker's code https://github.com/parkermac/LPM/blob/main/extract/moor/benthic_flux.py
        F_Det = LDeN[0,:,:] * 80 + SDeN[0,:,:] * 8 # mmol N/m2/d
        F_NO3 = F_Det.copy()
        F_NO3[F_Det > 1.2] = 1.2 # 1.2 mmol N/m2/d
        F_NH4 = F_Det - F_NO3
        F_NH4[F_NH4<0] = 0
        
        Oxy_sed_sum.append(np.nansum(F_NH4[jj,ii] * 106/16 * area[jj,ii]))       

        #---------- sediment SOD: a more strict way ----------
        # LdetritusN decomposition in sediment
        cff1_L = (LDeN[0,:,:] * 80)/24 / dz[0,:,:] * 1  # mmol/m2/hr --> mmol/m3     
        NH4_gain_L = np.zeros([NX,NY])        
        for i in range(NX): # loop the whole domain
            for j in range(NY):
                if cff1_L[i,j]*106/16 <= Oxy[0,:,:][i,j]:  # use O2: PN + O2 --> NH4; Oxy[0,:,:] is bottom oxygen                                        
                    NH4_gain_L[i,j] = cff1_L[i,j] #mmol/m3, NH4 gain from decomposing LdetN                                                                             
        NH4_gain_flux_L = NH4_gain_L * area * dz[0,:,:] # mmol but the actual unit is mmol/hr
        # SdetritusN decomposition in sediment
        cff1_S = (SDeN[0,:,:] * 8)/24 / dz[0,:,:] * 1  # mmol/m2/hr --> mmol/m3     
        NH4_gain_S = np.zeros([NX,NY])        
        for i in range(NX): # loop the whole domain
            for j in range(NY):
                if cff1_S[i,j]*106/16 <= Oxy[0,:,:][i,j]:  # use O2: PN + O2 --> NH4; Oxy[0,:,:] is bottom oxygen                                        
                    NH4_gain_S[i,j] = cff1_S[i,j] #mmol/m3, NH4 gain from decomposing LdetN                                                                             
        NH4_gain_flux_S = NH4_gain_S * area * dz[0,:,:] # mmol but the actual unit is mmol/hr
        
        NH4_gain_flux = NH4_gain_flux_L + NH4_gain_flux_S
        Oxy_sed_sum2.append(np.nansum(NH4_gain_flux[jj,ii] * 106/16))    # mmol O2/hr  
                
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
        
        #  O2 gas exchange
        Oxy_air_flux_sum.append(np.nansum(cff3[jj,ii] * (O2satu[jj,ii]-Oxy_surf[jj,ii]) * area[jj,ii])) #
        
        cnt += 1
        t.append(ds.ocean_time.values)
        ds.close()       
    dt00 = dt00 + timedelta(days=1)
        
# save netcdf
from netCDF4 import Dataset
nc = Dataset('O2_bgc_shallow_deep_'+ds0+'_'+ds1+'.nc','w')
time = nc.createDimension('time', len(t))
eta_rho = nc.createDimension('eta_rho', NX)
xi_rho = nc.createDimension('xi_rho', NY)
s_rho = nc.createDimension('s_rho', 30)

times = nc.createVariable('time','f8',('time',))
times.units = 'seconds*1e9 since 1970-01-01 00:00:00'
Oxy_pro_sum_shallow_tmp = nc.createVariable('Oxy_pro_sum_shallow','f4',('time',),compression='zlib',complevel=9)
Oxy_pro_sum_shallow_tmp.units = 'mmol O2/hr'
Oxy_pro_sum_deep_tmp = nc.createVariable('Oxy_pro_sum_deep','f4',('time',),compression='zlib',complevel=9)
Oxy_pro_sum_deep_tmp.units = 'mmol O2/hr'
Oxy_nitri_sum_shallow_tmp = nc.createVariable('Oxy_nitri_sum_shallow','f4',('time',),compression='zlib',complevel=9)
Oxy_nitri_sum_shallow_tmp.units = 'mmol O2/hr'
Oxy_nitri_sum_deep_tmp = nc.createVariable('Oxy_nitri_sum_deep','f4',('time',),compression='zlib',complevel=9)
Oxy_nitri_sum_deep_tmp.units = 'mmol O2/hr'
Oxy_remi_sum_shallow_tmp = nc.createVariable('Oxy_remi_sum_shallow','f4',('time',),compression='zlib',complevel=9)
Oxy_remi_sum_shallow_tmp.units = 'mmol O2/hr'
Oxy_remi_sum_deep_tmp = nc.createVariable('Oxy_remi_sum_deep','f4',('time',),compression='zlib',complevel=9)
Oxy_remi_sum_deep_tmp.units = 'mmol O2/hr'
Oxy_sed_sum_tmp = nc.createVariable('Oxy_sed_sum','f4',('time',),compression='zlib',complevel=9)
Oxy_sed_sum_tmp.units = 'mmol O2/day'
Oxy_sed_sum_tmp2 = nc.createVariable('Oxy_sed_sum2','f4',('time',),compression='zlib',complevel=9)
Oxy_sed_sum_tmp2.units = 'mmol O2/hr'
Oxy_vol_sum_shallow_tmp = nc.createVariable('Oxy_vol_sum_shallow','f4', ('time',),compression='zlib',complevel=9)
Oxy_vol_sum_shallow_tmp.units = 'mmol O2'
Oxy_vol_sum_deep_tmp = nc.createVariable('Oxy_vol_sum_deep','f4', ('time',),compression='zlib',complevel=9)
Oxy_vol_sum_deep_tmp.units = 'mmol O2'
Oxy_air_flux_sum_tmp = nc.createVariable('Oxy_air_flux_sum','f4', ('time',),compression='zlib',complevel=9)
Oxy_air_flux_sum_tmp.units = 'mmol O2/hr'


times[:] = t
Oxy_pro_sum_shallow_tmp[:] = Oxy_pro_sum_shallow
Oxy_pro_sum_deep_tmp[:] = Oxy_pro_sum_deep
Oxy_nitri_sum_shallow_tmp[:] = Oxy_nitri_sum_shallow
Oxy_nitri_sum_deep_tmp[:] = Oxy_nitri_sum_deep
Oxy_remi_sum_shallow_tmp[:] = Oxy_remi_sum_shallow
Oxy_remi_sum_deep_tmp[:] = Oxy_remi_sum_deep
Oxy_sed_sum_tmp[:] = Oxy_sed_sum
Oxy_sed_sum_tmp2[:] = Oxy_sed_sum2
Oxy_vol_sum_shallow_tmp[:] = Oxy_vol_sum_shallow
Oxy_vol_sum_deep_tmp[:] = Oxy_vol_sum_deep
Oxy_air_flux_sum_tmp[:] = Oxy_air_flux_sum

nc.close()
