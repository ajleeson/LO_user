'''
Modify ocean history file to add 
dTIC and dTA, with values of 0 everywhere.

run add_deltas_to_hisfile -gtx cas7_t1_x11ab -0 2020.05.30 -ro 5

Note that this script can only be run for one day at a time

'''

import sys
import shutil
import argparse
from datetime import datetime, timedelta
from time import time
import xarray as xr
import numpy as np
import matplotlib.pylab as plt
import csv
import sys
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

####################################################
# argument parsing

parser = argparse.ArgumentParser()
# arguments without defaults are required
parser.add_argument('-gtx', '--gtagex', default='cas7_t1_x11ab', type=str) # e.g. cas7_t1_x11b
parser.add_argument('-0', '--ds0', type=str)        # e.g. 2019.07.04
parser.add_argument('-ro', '--roms_out_num', type=int) # 2 = Ldir['roms_out2'], etc.

args = parser.parse_args()

gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

###################################################

gtagex_new = 'cas7_t1twin_x11aboaetwin'

##################################################

# get all information from arguments
argsd = args.__dict__
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]


ds0 = args.ds0
dt0 = datetime.strptime(ds0, Lfun.ds_fmt)

date = 'f' + dt0.strftime(Lfun.ds_fmt)

# set output location
out_dir = ('../../../LO_roms/' + gtagex_new + '/' + date)
Lfun.make_dir(out_dir)

# get original history file (use the prior day's history file)

roms_out_dir = Ldir['roms_out'] / Ldir['gtagex'] / date
ds_og_his = xr.open_dataset(roms_out_dir / 'ocean_his_0002.nc')

# make a copy of the original dataset to modify
ds_new = ds_og_his.copy()

# add delta variable to dataset, using a copy of TIC as 
# a reference of a 3D variable mapped to rho grid.
# and set all values to be zero
# dTIC
dTIC_conc = xr.zeros_like(ds_new['TIC'])
ds_new['dTIC'] = dTIC_conc
ds_new['dTIC'].attrs['long_name'] = 'Difference in TIC due to OAE'
ds_new['dTIC'].attrs['standard_name'] = 'molar_concentration_of_dTIC_'
ds_new['dTIC'].attrs['units'] = 'mmol m-3' 
ds_new['dTIC'].attrs['field'] = 'dTIC_' 
# dTA
dTA_conc = xr.zeros_like(ds_new['alkalinity'])
ds_new['dTA'] = dTA_conc
ds_new['dTA'].attrs['long_name'] = 'Difference in TA due to OAE'
ds_new['dTA'].attrs['standard_name'] = 'molar_concentration_of_dTA_'
ds_new['dTA'].attrs['units'] = 'meq m-3' 
ds_new['dTA'].attrs['field'] = 'dTA_' 

################################################################
# Save .nc files
print('Saving {}'.format(date))
ds_new.to_netcdf(str(out_dir) + '/ocean_his_0002.nc')

print('Done')
