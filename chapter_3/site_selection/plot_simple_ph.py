'''
Plot simple plot for Aurora, 
calculations of pH using pyco2sys

This is goign to make a plot for 1 April 2020
but I'm also showing you how to use indexing to speed up pyco2sys for the whole year 

NOTE: PYCO2 SYS wants practical salinity, but I only saved absolute salinity in the SML extraction. 
I can rerun the extractions and save SP - but shouldn't make a huge difference here. For this 
example I used gsw to go back to SP, but we should just be using the value SP from LO (or obs values of SP)

Showing you and example of how to run a loop thru the slices because pyco2sys can get bogged down 
by big arrays. The instance, here, where we have the whole LO domain (1 layer and 365(6)days) 
I step thru each day 'slice', 
but a similar approach would be used for one history file... where we have the LO domain and one time but 
30 layers. 
When you have a big file, it's faster to step thru the slices and use our land mask. 
When have smaller files you can pass the whole thing to pyco2sys (like year of average files for one mooring). 

TIME!!
It takes ~4-5seconds to calc pH per day for the LO domain. 
So one year would be ~30 mins. 
Similarly it'd take ~2-3 minutes to do all 30 layers of the LO domain from a history file

My testing is set to default True, which will put nlay to 1 and only calc 1 Jan 2020 , from 
my defaults. This is because I didn't want to wait 30 mins for the calculation ...

If this were me, and I were making plots, I'd probably calc the pH once and save a smaller file
and then do my calcs and plotting from it so that I wouldn't wait 30 mins for each plot.
'''

# imports
import sys
import argparse
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import calendar
import xarray as xr
from time import time
import numpy as np
import pandas as pd
from dateutil import parser
import gsw
import sys
import PyCO2SYS as pyco2

import matplotlib.pyplot as plt
#import matplotlib.cm as cm
import cmcrameri.cm as cm
import pickle 

# command line arugments
parser = argparse.ArgumentParser()
# which run to use
parser.add_argument('-gtx', '--gtagex', type=str, default = 'cas7_t1_x11ab')   # e.g. cas6_v3_lo8b
# select time period and frequency
parser.add_argument('-0', '--ds0', type=str, default = '2020.05.22') # e.g. 2015.01.01
parser.add_argument('-thresh', '--threshold', type=str, default = '0.125') # e.g. 2025.12.31
parser.add_argument('-lt', '--list_type', type=str, default = 'average') # list type: hourly, daily, weekly, lowpass
# select job name
parser.add_argument('-job', '--jobname', type=str, default = 'LO_domain') # job name
# Optional: set max number of subprocesses to run at any time
parser.add_argument('-Nproc', type=int, default=10)
# Optional: for testing
parser.add_argument('-test', '--testing', default=True, type=Lfun.boolean_string)
# get the args and put into Ldir
args = parser.parse_args()
# test that main required arguments were provided
argsd = args.__dict__
for a in ['gtagex']:
    if argsd[a] == None:
        print('*** Missing required argument: ' + a)
        sys.exit()
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]

# get the args and put into Ldir
args = parser.parse_args()

Ldir = Lfun.Lstart()

# this is kate code specific -- i calc'd the sml using two diff thresholds
if args.threshold == '0.125':
    folder_in = 'threshold_p125'
elif args.threshold == '0.03':
    folder_in = 'threshold_p03'
else:
    print('exiting')
    sys.exit()

# set your paths and open the dataset  
fin = Ldir['LOo'] / 'extract' / args.gtagex / 'sml_plus' / folder_in 
fout = Ldir['parent'] / 'crest26_output' / 'sml_evaluation' / args.gtagex / args.jobname / 'dye_test_dates' / 'pH' 
Lfun.make_dir(fout, clean=False)

thisyr = args.ds0.split('.')[0]
fn = args.jobname+'_sml_plus_'+str(thisyr)+'.01.01_'+str(thisyr) + '.12.31'
fn_in = fin / (args.jobname+'_'+str(thisyr)+'.01.01_'+str(thisyr) + '.12.31') / (fn + '.nc')

##############################################################
# STEP 1. open data set; use gsw and prep for pyco2sys calcs
##############################################################
# this takes about a minute for the whole year on my laptop
# could speedup by just doing one day 
tt0 = time()

# we are just going to be plotting one example in this script
ds = xr.open_dataset(fn_in, decode_times=True)  
ot = pd.to_datetime(ds.ocean_time.values)
SA = ds.SA.values
CT = ds.CT.values
Z = ds.zSML.values
ALK = ds.ALK.values
TIC = ds.TIC.values

mask_rho = ds.mask_rho.values       # 0 land / 1 water
lon = np.nanmean(ds.lon_rho)
lat = np.nanmean(ds.lat_rho)

# enter data as requested by pyco2sys
P = gsw.p_from_z(Z, lon, lat)     # using the base of the SML as z, usually you'd use z_rho here
ti = gsw.t_from_CT(SA, CT, P)                  # in situ temperature [degC] 
rho = gsw.rho(SA, CT, P)                       # in situ density [kg m-3]
# this was my bad not saving density I only saved CT and SA, and pyco2sys needs SP 
# This next step is a little sloppy to recalc here, but is okay for this one plot 
SP = gsw.SP_from_SA(SA, P, lon, lat) 

# Convert from micromol/L to micromol/kg using in situ dentity because these are the
# units expected by pyco2.
ALK1 = 1000 * ALK / rho
TIC1 = 1000 * TIC / rho

print('Time to set up variables = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

##############################################################
# STEP 2. calc pH 
##############################################################
# We do this one day at a time because pyco2 can get bogged down if we fed it
# the whole 3-D array for the whole LO domain at once. It can handle 3D arrays, 
# but large ones slow down the calculations
tt0 = time()
PH = np.nan * np.ones(SP.shape)       # Initialize array to hold results
nt, nr, nc = SP.shape                 # handy dimension sizes
pmat = np.nan * np.ones((nr,nc))      # Initialize array for single time
if args.testing:
    nlay = 1     # number of NT to process
else:
    nlay = nt
for ii in range(nlay):
    # ii = nz-1 # hand override to test surface layer
    tt00 = time()
    # Note that by using the [mask_rho==1] operator on each day (slice)
    # we go from a 2-D array with nan's to a 1-D vector with no nan's. We
    # do this to speed up the calculation. Note that the LO grid is about half
    # landmask.
    aALK = ALK1[ii,:,:].squeeze()[mask_rho==1]
    aTIC = TIC1[ii,:,:].squeeze()[mask_rho==1]
    aTemp = ti[ii,:,:].squeeze()[mask_rho==1]
    aPres = P[ii,:,:].squeeze()[mask_rho==1]
    aSalt = SP[ii,:,:].squeeze()[mask_rho==1]
    # Note: here is where to get info on the inputs, outputs, and units:
    # https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/
    CO2dict = pyco2.sys(par1=aALK, par1_type=1, par2=aTIC, par2_type=2,
        salinity=aSalt, temperature=aTemp, pressure=aPres,
        total_silicate=50, total_phosphate=2, opt_k_carbonic=10, opt_buffers_mode=0)
    aPH = CO2dict['pH_total']
    # Then write the pH total field into the PH array, indexing with [mask_rho==1].
    ppmat = pmat.copy()
    ppmat[mask_rho==1] = aPH
    PH[ii,:,:] = ppmat 
    print('  ii = %d' % (ii))
    print('  Time to get one slice = %0.2f sec' % (time()-tt00))
    sys.stdout.flush()
print('Time to calculate PH for all slices = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

# Save pH to a new netCDF file
# where to put output figures
out_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'sml_plus' / 'LO_domain_2020.01.01_2020.12.31' / 'pH'
Lfun.make_dir(out_dir)
ds_ph = xr.DataArray(PH, dims=('time','y','x'), name='pH')
# Save to .nc file
ds_ph.to_netcdf(out_dir / 'pH_2020_sml_p125.nc')

##############################################################
# STEP 2. plot pH 
##############################################################

if args.testing:
    jj = 0
else:
    jj = 142     # May 22, 2020

# PLOTTING
plt.close('all')
fs=12
plt.rc('font', size=fs)
fig = plt.figure(figsize=(8,10))
fig.set_size_inches(8,10, forward=False)

xmin = -126 
xmax = -122
xxticks = [-122, -123, -124, -125, -126]
xxticklabels = ['-122', '-123', '-124', '-125', '-126']

ymin = 45.5
ymax = 50.5
yyticks = [45.5,46,47,48,49,50,50.5]

vmax = 8.5
vmin = 7

axw = plt.subplot2grid((2,2), (0,0), colspan=2,rowspan=2)
pfun.add_coast(axw)
pfun.dar(axw)
axw.set_xlabel('Longitude')
axw.set_ylabel('Latitude')
axw.set_xlim([xmin, xmax])
axw.set_xticks(xxticks)
axw.set_ylim([ymin,ymax])
axw.set_yticks(yyticks)

cmap = plt.get_cmap(cm.roma)
cmap.set_extremes(under = 'pink', over = 'darkviolet')

cpw = axw.pcolormesh(ds['lon_rho'].values, ds['lat_rho'].values, PH[jj,:,:].squeeze(),vmin=vmin,vmax=vmax,cmap=cmap)
cbaxes = inset_axes(axw, width="5%", height="30%", loc='upper right')
fig.colorbar(cpw, cax=cbaxes, location='left', orientation='vertical',label = 'pH total', extend = 'both')

axw.contour(ds['lon_rho'], ds['lat_rho'], ds['h'], [400],
colors=['white'], linewidths=0.5, linestyles='solid',alpha=0.7)
axw.text(0.75, 0.05,'400m',color='DimGrey',weight='light',fontstyle='italic',transform=axw.transAxes,ha='right')

axw.set_title(f"pH:  {ot[jj]:%d-%m-%Y}")

#figname = 'pH_'+bb+'.png'
#fig.savefig(fout / figname)