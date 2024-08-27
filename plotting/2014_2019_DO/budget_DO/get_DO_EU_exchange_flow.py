"""
Calculate exchange flow terms for DO budget
Have shallow and deep layer
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
import math
import csv
import cmocean
import matplotlib.pylab as plt
import gsw
import pickle

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t0_x4b'
jobname = 'twentyoneinlets'
startdate = '2014.01.01'
enddate = '2014.12.31'
enddate_hrly = '2015.01.01 00:00:00'
year = '2014' # for making a date label

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)

# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_EU_exchange'
Lfun.make_dir(out_dir)

# create time_vecotr
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]

# create dictionaries with interface depths
interface_dict = dict()
print('\n')

##########################################################
##        Calculate shallow and deep exchange flow      ##
##########################################################

stations = ['lynchcove','penn','budd','case','carr']

for i,station in enumerate(stations): # enumerate(sta_dict):
    # print status
    print('({}/{}) Working on {}...'.format(i+1,len(sta_dict),station))

    # initialize empty dataframe for saving
    df = pd.DataFrame()

    # get interface depth from csv file
    with open('interface_depths.csv', 'r') as f:
        for line in f:
            inlet, interface_depth = line.strip().split(',')
            interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
    z_interface = float(interface_dict[station])

    # get section information
    ds = xr.open_dataset('../../../../LO_output/extract/'+gtagex+
                         '/tef2/extractions_'+startdate+'_'+enddate+
                         '/'+station+'.nc')
    
    # calculate exchange flow DO budget at cross section
    oxygen = ds['oxygen'] # convert to mmol/m3 (time,z,p)
    velocity = ds['vel'] # m/s (time,z,p)
    area = ds['DZ'] * ds['dd'] # m^2 (time,z,p)
    exchange = oxygen * velocity * area * (1/1000) * (1/1000) # mmol/s * 1/1000 * 1/1000 = kmol/s (time,z,p)


    # manage one or two layers
    if math.isnan(z_interface):
        print('    One-layer...add code!!!')

    else:
        print('    Two-layers...Interface depth [m]: {}'.format(z_interface))

        # Get z_rho (time,z,p), the depth of all z_rho layers
        DZ = ds['DZ'].values
        # append a set of zeros at the very bottom
        DZ_w0 = np.insert(DZ, 0, np.zeros([8761,len(ds['p'].values)]), axis=1)
        # calculate cumulative zum
        DZ_cumsum = np.cumsum(DZ_w0, axis=1) # this is effectively z_w recreated, but measured from bottom
        DZ_cumsum = np.delete(DZ_cumsum, -1, axis=1) # drop the surface level
        # divide DZ by two to get the distance between z_w and z_rho
        zw_2_zrho = DZ/2
        # get z_rho (and move datum to surface)
        z_rho = DZ_cumsum + zw_2_zrho - ds['h'].values

        # mask surface and deep layers
        surf_layer = np.where(z_rho >= z_interface, 1, np.nan)
        deep_layer = np.where(z_rho < z_interface, 1, np.nan)

        # separate exchange into deep and surface layer
        exchange_surf = exchange * surf_layer
        exchange_deep = exchange * deep_layer

        # print('\nLayer thickness, with zero appended to beginning')
        # print(DZ_w0[0,:,0])
        # print('\nCumulative sum of layer thicknesses (starting from bottom layer)')
        # print(DZ_cumsum[0,:,0])
        # print('\n1/2 of DZ-- so the distance between a z_w value and a z-rho value')
        # print(zw_2_zrho[0,:,0])
        # print('\nz_rho, after shifting the datum from bottom back to surface, using h')
        # print(z_rho[0,:,0])

        # sum up exchange flow over surface and deep areas
        exchange_surf_total = exchange_surf.sum(axis=1).sum(axis=1)
        exchange_deep_total = exchange_deep.sum(axis=1).sum(axis=1)
        exchange_total = exchange.sum(axis=1).sum(axis=1)

        # create dataframe of values
        df['surface [kmol/s]'] = exchange_surf_total
        df['deep [kmol/s]'] = exchange_deep_total
        df['total [kmol/s]'] = exchange_total
        # save to pickle file
        df.to_pickle(out_dir / (station + '.p'))
