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
max_delta_pH = 0.2     # I think we can do 0.5
NZ = np.shape(ALK1)[1] #30

TIC0 = TIC1[0,:]
ALK0 = ALK1[0,:]
pH0 = pH_total[0,:]

# max across 11 years of average output files and add 0.2 (max_delta_pH)
pH_upper = np.nanmax(pH_total,axis=0)+max_delta_pH 

# take the minimum of the two WQ thresholds for alk calculations - no need for solver 
pH_a = np.ones(NZ)*pH_max
pH_target = np.minimum(pH_a,pH_upper)

# This is assuming a DIC change of 0; so would be for something where you're adding
# sodium hydoxide, if add NAOH can assume delta dic = 0 
# but delta dic = 0.5*TA if adding something like soda ash na2co3
target_CO2dict = pyco2.sys(par1=TIC0, par1_type=2, par2=pH_target, par2_type=3,
        salinity=SP[0,:], temperature=ti[0,:], pressure=P[0,:],
        total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)

ALK_target = target_CO2dict['alkalinity']
dALK = ALK_target - ALK0
print('Alkalinity addition: {}'.format(dALK))

print('Time to calculate alk adjustment for all layers one time = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

'''AA = ALK0+dALK

CO2dict2 = pyco2.sys(par1=AA, par1_type=1, par2=TIC0, par2_type=2,
        salinity=SP[0,:], temperature=ti[0,:], pressure=P[0,:],
        total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
PHC = CO2dict2['pH_total']'''

##########################################################
# STEP 5. use ARPA-E calcs 
##########################################################
tt0 = time()
arpae_dalk = (2000*rho[0,:])/1000  #2000mmol m-3, same as ROMS units
ae_dalk = np.ones(NZ)*arpae_dalk
ae_alk = ALK0+ae_dalk 
ae_CO2dict = pyco2.sys(par1=ae_alk, par1_type=1, par2=TIC0, par2_type=2,
        salinity=SP[0,:], temperature=ti[0,:], pressure=P[0,:],
        total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
ae_pH = ae_CO2dict['pH_total']

##########################################################
# STEP X. plotting vertical profiles 
##########################################################
plt.close('all')
fs=14
plt.rc('font', size=fs)
fig = plt.figure(figsize=(12,8))
fig.set_size_inches(12,8, forward=False)

'''ax1 = plt.subplot2grid((4,4), (0,0), colspan=1,rowspan=4)
plt.plot(SIG0,z_rho,'.',color='LightGrey')
ax1.set_title('SIG0 kg/m3')
ax1.set_ylabel('z m')'''

ax1 = plt.subplot2grid((4,4), (0,0), colspan=1,rowspan=4)
plt.plot(TIC,z_rho,'.',color='LightGrey')
ax1.set_title('TIC umol/kg')

ax2 = plt.subplot2grid((4,4), (0,1), colspan=1,rowspan=4)
plt.plot(ALK,z_rho,'.',color='LightGrey')
ax2.set_title('ALK umol/kg')

'''ax3 = plt.subplot2grid((4,4), (0,2), colspan=1,rowspan=4)
plt.plot(TIC,z_rho,'.',color='LightGrey')
ax3.set_title('TIC umol/kg')'''

ax34 = plt.subplot2grid((4,4), (0,2), colspan=2,rowspan=4)
plt.plot(pH_total,z_rho,'.',color='LightGrey')
plt.axvline(x=9, color='red', linestyle='-')
plt.axvline(x=8.5, color='Violet', linestyle=':')
plt.axvline(x=8.2, color='Black', linestyle=':')
plt.axvline(x=8, color='DimGrey', linestyle=':')
plt.axvline(x=7.6, color='Green', alpha = 0.1, linestyle='-')
ax34.set_title('pH total')

plt.plot(pH_total[0,:],z_rho[0,:],'*',color='Violet')
plt.plot(pH_target,z_rho[0,:],'*',color='Navy')
plt.plot(ae_pH,z_rho[0,:],'x',color='Black')
#plt.plot(PHC,z_rho[0,:],'o',color='CYAN')

p25 = np.percentile(pH_total, 25, axis=0)
p50 = np.percentile(pH_total, 50, axis=0)
p75 = np.percentile(pH_total, 75, axis=0)
zm = np.nanmean(z_rho,axis=0)

'''IQR = p75-p25
lower = p25 - 1.5*IQR 
upper = p75 + 1.5*IQR '''

plt.plot(p25,zm,'|',color = 'Navy', markersize = 5)
plt.plot(p75,zm,'|',color = 'Navy', markersize = 5)
plt.plot(p50,zm,'|',color = 'DodgerBlue', markersize = 5)

'''plt.plot(p50+0.2,zm,'o',color = 'pink', markersize = 5)
plt.plot(p50+0.5,zm,'o',color = 'red', markersize = 5)
plt.plot(p50+1.5,zm,'o',color = 'black', markersize = 5)'''

fig2 = plt.figure(figsize=(12,8))
fig2.set_size_inches(12,8, forward=False)

drho = np.diff(SIG0)
dz = np.diff(z_rho)
dzrho = -drho / dz 
zd = z_w[:,1:-1]

ax21 = plt.subplot2grid((4,4), (0,0), colspan=1,rowspan=4)
plt.plot(SIG0,z_rho,'.',color='LightGrey')
ax21.set_title('SIG0 kg/m3')
ax21.set_ylabel('z m')

Sp25 = np.percentile(dzrho, 25, axis=0)
Sp50 = np.percentile(dzrho, 50, axis=0)
Sp75 = np.percentile(dzrho, 75, axis=0)
zdm = np.nanmean(zd,axis=0)

ax22 = plt.subplot2grid((4,4), (0,1), colspan=3,rowspan=4)
plt.plot(dzrho,zd,'.',color='LightGrey')
#plt.plot(Sp25,zdm,'|',color = 'Navy', markersize = 5)
#plt.plot(Sp75,zdm,'|',color = 'Navy', markersize = 5)
plt.axvline(x=0.015, color='Blue', alpha = 0.5, linestyle='-')
plt.plot(Sp50,zdm,'|',color = 'Red', markersize = 5)

ax22.set_title('drho dz kg/m3')
ax22.set_ylabel('z m')

##########################################################
# STEP X. plotting dic ta space
##########################################################
fs=14
plt.rc('font', size=fs)
fig3 = plt.figure(figsize=(12,8))
fig3.set_size_inches(12,8, forward=False)

# make a range of values and take avg of salt and temp 
# will calc at P = 0 here 
dic_range = np.linspace(1200,2800,50)
alk_range = np.linspace(1200,2800,50)
gDIC, gALK = np.meshgrid(dic_range, alk_range)
gSALT = np.nanmean(SP)
gTEMP = np.nanmean(ti)

# pass flattened array to pyco2sys 
results = pyco2.sys(
    par1=gALK.flatten(), par1_type=1, # Total Alkalinity
    par2=gDIC.flatten(), par2_type=2, # Total Inorganic Carbon
    salinity=gSALT, temperature=gTEMP, pressure=0
)

# reshape back to 2D grid shape 
gPH = results["pH"].reshape(gDIC.shape)
#gOMEGA = results["saturation_aragonite"].reshape(gDIC.shape)

#contours = plt.contour(gDIC, gALK, gPH, levels=15) #, cmap='viridis')
contours = plt.contour(gDIC, gALK, gPH, colors = 'dimGrey', levels=15, linewidths=0.8, linestyles = '-') 
plt.clabel(contours, inline=True, fontsize=8, fmt='%.2f')

'''contours2 = plt.contour(gDIC, gALK, gOMEGA, colors = 'dimGrey', linewidths=0.8, linestyles = ':') 
plt.clabel(contours2, inline=True, fontsize=8, fmt='%.2f')'''

plt.plot(TIC1,ALK1,'.',color='LightGrey',alpha=0.6)
plt.plot(TIC1[0,:],ALK1[0,:],'.',color='Violet')
plt.plot(TIC1[0,:],ALK_target,'.',color='Navy')

plt.xlabel('DIC ($\mu mol/kg$)')
plt.ylabel('Alkalinity ($\mu mol/kg$)')
plt.title('test at '+ args.moorname)
plt.show()