# calculate SOD
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

burial=50 # 50% burial of sinking detritus in the Salish Sea

#%------------------------------------------------
Ldir = Lfun.Lstart()
Ldir['roms_out'] = Ldir['roms_out']
Ldir['gtagex'] = 'cas7_t0_x9b'

ds0 = '2018.05.10'
ds1 = '2018.05.20'
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
seg_name = 'seg_info_dict_cas7_c2_noriv.p'
seg_df = pd.read_pickle(seg_name)
ji_list = seg_df['sog6_m']['ji_list']
jj = [x[0] for x in ji_list]
ii = [x[1] for x in ji_list]

inDomain = np.zeros([NX, NY])
inDomain[jj,ii] = 1 # inside domain, the index=1, outside domian, the index =0, this step is to reduce the size of spatial variable (hopefully)

#% air-sea flux related parameters
#rho0 = 1025.0 # ROMS/Modules/mod_scalars.F
#sec2day = 1.0/86400
#dt = 3600 #40 # sec  ???????????????????????
#BioIter = 1
#dtdays = dt * sec2day / BioIter
dtdays = 3600 * (1.0/86400) / 1
#cff1 = rho0 * 550.0
#cff1 = 1025.0 * 550.0
#RW14_OXYGEN_SC = False
#if RW14_OXYGEN_SC:   # cff2: s2/m
#    cff2 = dtdays * 0.251 * 24 / 100  # 0.251: (cm/h)(m/s)^(-2), 
#else:
#    cff2 = dtdays * 0.31 * 24 / 100
cff2_air = dtdays * 0.31 * 24 / 100
# formulation for Schmidt number coefficient
#OCMIP_OXYGEN_SC = False
#RW14_OXYGEN_SC = False
#if OCMIP_OXYGEN_SC: # Keeling et al., 1998
#    A_O2 = 1638.0; B_O2 = 81.83; C_O2 = 1.483
#    D_O2 = 0.008004; E_O2 = 0.0
#elif RW14_OXYGEN_SC: # Wanninkhof, 2014
#    A_O2 = 1920.4; B_O2 = 135.6; C_O2 = 5.2122
#    D_O2 = 0.10939; E_O2 = 0.00093777
#else: # Wanninkhof, 1992
#    A_O2 = 1953.4; B_O2 = 128.0; C_O2 = 3.9918
#    D_O2 = 0.050091; E_O2 = 0.0
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
# light attenuation
#PARfrac = 0.43
#Cp = 3985 # Joules/kg/degC
#rho0 = 1025 # kg/m3
#AttSW = 0.05
#AttChl = 0.012
#Vp = 1.7;
#PhyIS = 0.07 # initial slope of P-I curve [1/(Watts m-2 day)]
#rOxNO3 = 138/16
#rOxNH4 = 106/16
#K_NO3 = 10 # [1/(millimole_N m-3)]
#K_NH4 = 10
#SDeRRN = 0.1 #Small detritus remineralization rate N-fraction [1/day]
#LDeRRN = 0.1 
#NitriR = 0.05 # Nitrification rate: oxidation of NH4 to NO3 [1/day]
 
#sec2day = 1.0/86400 # ROMS/Modules/mod_scalars.F
#dt = 3600 #40 # sec ???????????????????????
#BioIter = 1
#dtdays = dt * sec2day / BioIter

#Ws_L = 80 # m/d
#Ws_S = 8 # m/d

t = []
Oxy_sed_spatial2 = np.zeros((24*32, NX, NY))

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
        zeta = ds.zeta.values.squeeze()
        h = ds.h.values
        zw = zrfun.get_z(h, zeta, S, only_w=True)
        dz = np.diff(zw, axis=0)
        SDeN = ds.SdetritusN.values.squeeze()
        LDeN = ds.LdetritusN.values.squeeze()
        Oxy = ds.oxygen.values.squeeze()
   
        #---------- sediment SOD: a more strict way ----------
        # LdetritusN decomposition in sediment
        # FC_L = (LdetN_bot * Ws_L)/24 # mmol/m2/hr
        # cff1_L = FC_L / dz[0,:,:] * 1  # mmol/m2/hr --> mmol/m3 
        cff1_L = (LDeN[0,:,:] *(1-burial/100) * 80)/24 / dz[0,:,:] * 1  # mmol/m2/hr --> mmol/m3
        NH4_gain_L = np.zeros([NX,NY])
        for i in range(NX): # loop the whole domain
            for j in range(NY):
                if cff1_L[i,j]*106/16 <= Oxy[0,:,:][i,j]:  # use O2: PN + O2 --> NH4; Oxy[0,:,:] is bottom oxygen                                        
                    NH4_gain_L[i,j] = cff1_L[i,j] #mmol/m3, NH4 gain from decomposing LdetN                                                                             
        NH4_gain_flux_L = NH4_gain_L * area * dz[0,:,:] # mmol but the actual unit is mmol/hr
        # SdetritusN decomposition in sediment
        # FC_S = (SdetN_bot * Ws_S)/24 # mmol/m2/hr
        # cff1_S = FC_S / dz[0,:,:] * 1  # mmol/m2/hr --> mmol/m3 
        cff1_S = (SDeN[0,:,:] *(1-burial/100) * 8)/24 / dz[0,:,:] * 1  # mmol/m2/hr --> mmol/m3
        NH4_gain_S = np.zeros([NX,NY])
        for i in range(NX): # loop the whole domain
            for j in range(NY):
                if cff1_S[i,j]*106/16 <= Oxy[0,:,:][i,j]:  # use O2: PN + O2 --> NH4; Oxy[0,:,:] is bottom oxygen                                        
                    NH4_gain_S[i,j] = cff1_S[i,j] #mmol/m3, NH4 gain from decomposing LdetN                                                                             
        NH4_gain_flux_S = NH4_gain_S * area * dz[0,:,:] # mmol but the actual unit is mmol/hr
        
        NH4_gain_flux = NH4_gain_flux_L + NH4_gain_flux_S
        Oxy_sed_spatial2[cnt,:,:] = NH4_gain_flux * inDomain * 106/16
                
        cnt += 1
        t.append(ds.ocean_time.values)
        ds.close()       
    dt00 = dt00 + timedelta(days=1)
        
# save netcdf
from netCDF4 import Dataset
nc = Dataset('O2_SOD_burial_'+str(burial)+'_'+ds0+'.nc','w')
time = nc.createDimension('time', len(t))
eta_rho = nc.createDimension('eta_rho', NX)
xi_rho = nc.createDimension('xi_rho', NY)
s_rho = nc.createDimension('s_rho', 30)

times = nc.createVariable('time','f8',('time',))
times.units = 'seconds*1e9 since 1970-01-01 00:00:00'
Oxy_sed_spatial_tmp2 = nc.createVariable('Oxy_sed_spatial2','f4',('time','eta_rho','xi_rho'),compression='zlib',complevel=9)
Oxy_sed_spatial_tmp2.units = 'mmol O2/hr'

times[:] = t

Oxy_sed_spatial_tmp2[:] = Oxy_sed_spatial2[0:len(t),:,:]

nc.close()

#print('total time = %0.1f sec' % (time()-tt0))

