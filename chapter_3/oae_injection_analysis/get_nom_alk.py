'''
mooring plots to help w/ 
oae calcs

'''

import argparse
import sys 
import matplotlib.pyplot as plt
import xarray as xr 
import numpy as np 
from time import time
import PyCO2SYS as pyco2

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
from matplotlib.patches import Polygon
from pathlib import PosixPath, PurePosixPath

import gsw 
from scipy.optimize import brentq

# command line arugments
parser = argparse.ArgumentParser()
# which run was used:
# select years 
parser.add_argument('-gtx', '--gtagex', type=str, default = 'cas7_t1_x11ab')   
parser.add_argument('-d0', '--ds0', type=str, default = '2015.01.01') 
parser.add_argument('-d1', '--ds1', type=str, default = '2024.12.31') 
parser.add_argument('-job', '--jobname', type=str, default = 'whidbey_basin') # random job name default to whidbey basin
parser.add_argument('-moor', '--moorname', type=str, default = 'wb')

Ldir = Lfun.Lstart()

# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided
argsd = args.__dict__
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]

data_path_in = Ldir['LOo'] / 'extract' / args.gtagex / 'moor' / args.jobname
fn_in = data_path_in / (args.moorname + '_' + args.ds0 + '_' + args.ds1 + '.nc')

ds = xr.open_dataset(fn_in, decode_times = True)

'''if np.any(ds.wetdry_mask_rho.values==0):
    print('dry conditions present, exiting')
    sys.exit()'''

#dye = ds.dye_01.values
z_rho = ds.z_rho.values
z_w = ds.z_w.values
dz = np.diff(ds.z_w.values)
zeta = ds.zeta.values
ocean_time = ds.ocean_time.values

tempC = ds.temp.values 
SP = ds.salt.values
ALK = ds.alkalinity.values 
TIC = ds.TIC.values 

lat = ds.lat_rho.values
lon = ds.lon_rho.values

# will be an issue in co2sys and is just strange
if np.nanmin(SP)<0:
    print('Weirdness in salinity, negative salt')
    SP[SP<0] = 0


##########################################################
# STEP 1. Calc SA and SIG0
##########################################################
tt0 = time()

P = gsw.p_from_z(z_rho,lat)            # pressure [dbar]
SA = gsw.SA_from_SP(SP, P, lon, lat)   # absolute salinity [g kg-1]
CT = gsw.CT_from_pt(SA, tempC)         # conservative temperature [degC]
SIG0 = gsw.sigma0(SA,CT)   

print('Time to apply GSW = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

##########################################################
# STEP 3. cchem calculations
##########################################################
tt0 = time()

ti = gsw.t_from_CT(SA, CT, P) # in situ temperature [degC] 
rho = gsw.rho(SA, CT, P) # in situ density [kg m-3]

# Convert from micromol/L to micromol/kg using in situ dentity because these are the
# units expected by pyco2.
ALK1 = 1000 * ALK / rho
TIC1 = 1000 * TIC / rho

CO2dict = pyco2.sys(par1=ALK1, par1_type=1, par2=TIC1, par2_type=2,
        salinity=SP, temperature=ti, pressure=P,
        total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)

#ARAG = CO2dict['saturation_aragonite']
pH_total = CO2dict['pH_total']

print('Time to calculate pH for all layers = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

##########################################################
# STEP 4. test dalk calculations
##########################################################
tt0 = time()
pH_max = 8.5 
max_delta_pH = 0.5# 0.2     # I think we can do 0.5
NZ = np.shape(ALK1)[1] #30

# take the minimum of the two WQ thresholds for alk calculations
pH_plusdelta = pH_total + max_delta_pH
# set values above 8.5 to 8.5
pH_plusdelta[pH_plusdelta>pH_max] = pH_max
pH_target = pH_plusdelta.copy()

# This is assuming a DIC change of 0; so would be for something where you're adding
# sodium hydoxide, if add NAOH can assume delta dic = 0 
# but delta dic = 0.5*TA if adding something like soda ash na2co3
target_CO2dict = pyco2.sys(par1=TIC1, par1_type=2, par2=pH_target, par2_type=3,
        salinity=SP, temperature=ti, pressure=P,
        total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)

ALK_target = target_CO2dict['alkalinity']
dALK = ALK_target - ALK1

print('Time to calculate alk adjustment for all layers one time = %0.2f sec' % (time()-tt0))
sys.stdout.flush()


##########################################################
# STEP 5. plotting vertical profiles 
##########################################################
plt.close('all')
fs=14
plt.rc('font', size=fs)
fig = plt.figure(figsize=(12,8))
fig.set_size_inches(12,8, forward=False)

ax1 = plt.subplot2grid((4,4), (0,0), colspan=1,rowspan=4)
plt.plot(dALK,z_rho,'.',color='pink',alpha=0.05)
ax1.set_title('dALk umol/kg')

ax2 = plt.subplot2grid((4,4), (0,1), colspan=1,rowspan=4)
plt.plot(ALK,z_rho,'.',color='LightGrey')
ax2.set_title('ALK umol/kg')


ax34 = plt.subplot2grid((4,4), (0,2), colspan=2,rowspan=4)
plt.plot(pH_total,z_rho,'.',color='LightGrey')
plt.axvline(x=9, color='red', linestyle='-')
plt.axvline(x=8.5, color='Violet', linestyle=':')
plt.axvline(x=8.2, color='Black', linestyle=':')
plt.axvline(x=8, color='DimGrey', linestyle=':')
plt.axvline(x=7.6, color='Green', alpha = 0.1, linestyle='-')
ax34.set_title('pH total')

plt.plot(pH_target,z_rho,'D',color='Navy', alpha=0.05, markeredgecolor='none')

# p25 = np.percentile(pH_total, 25, axis=0)
# p50 = np.percentile(pH_total, 50, axis=0)
# p75 = np.percentile(pH_total, 75, axis=0)
# zm = np.nanmean(z_rho,axis=0)

# plt.plot(p25,zm,'|',color = 'Navy', markersize = 5)
# plt.plot(p75,zm,'|',color = 'Navy', markersize = 5)
# plt.plot(p50,zm,'|',color = 'DodgerBlue', markersize = 5)

##########################################################
# STEP 6. Get weighted averaged of dAlk in the top 8 layers 
##########################################################

# crop dALK to be the top 8 sigma layers
dALK_top8 = dALK[:,-8:]
# and density too
rho_top8 = rho[:,-8:]

# get dz of top 8 sigma layers (last 9 values)
z_w = ds.z_w.values
dz = np.diff(z_w, axis=1) # [m]
dz_top8 = dz[:,-8:]

# get thickness of top 8 sigma layers
thickness_top8 = np.sum(dz_top8, axis=1)

# get weights for averaging
weights_masked = dz_top8 / thickness_top8[:, np.newaxis]

# inside takes a weighted average of each var, outside extracts 
nom_dALK_perday = np.ma.getdata(np.ma.average(dALK_top8, weights = weights_masked, axis=1, keepdims=True))
nom_rho_perday = np.ma.getdata(np.ma.average(rho_top8, weights = weights_masked, axis=1, keepdims=True))

# get mean nominal alkalinity concentration
nom_dALK_mean = np.nanmean(nom_dALK_perday)
ax1.axvline(x=nom_dALK_mean, color='hotpink', linestyle='-')

# get mean rho per day (mean rho in top 8 sigma layers)
nom_rho_mean = np.nanmean(nom_rho_perday)

# convert back to mmol/m3 from umol/kg
nom_dALK_mmol_m3 = nom_dALK_mean * nom_rho_mean / 1000

print('Nominal alkalinity addition per day: {} mmol/m3'.format(nom_dALK_mmol_m3))
# 190.41386851288794 mmol/m3
