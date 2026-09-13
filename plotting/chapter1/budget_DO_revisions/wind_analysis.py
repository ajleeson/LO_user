"""
Plot surface and deep budget for all 21 inlets
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import get_two_layer
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from matplotlib import colormaps
from matplotlib.colors import ListedColormap
from scipy.stats import pearsonr
from scipy.linalg import lstsq
from matplotlib.ticker import AutoMinorLocator

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

interface_type = 'tef' # doesn't matter for one-layer

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

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')[2::]
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]
# crop time vector (because we only have jan 2 - dec 30)
dates_no_crop = dates_local_daily
dates_local_daily = dates_local_daily

print('\n')


######################
# wind data

# Add daily average wind speed
col_names = [ "YY","MM","DD","hh","mm",
    "WDIR","WSPD","GST","WVHT","DPD","APD","MWD",
    "PRES","ATMP","WTMP","DEWP","VIS","TIDE"]
df = pd.read_csv("46121h2017.txt",
    sep=r"\s+",
    comment="#",
    names=col_names,
    na_values=[99, 99.0, 999, 999.0, 9999, 9999.0])
# rename for datetime parsing
df = df.rename(columns={"YY": "year",
    "MM": "month",
    "DD": "day",
    "hh": "hour",
    "mm": "minute"})
# create datetime index
df["datetime"] = pd.to_datetime(df[["year","month","day","hour","minute"]])
print(df)
df = df.set_index("datetime").drop(columns=["year","month","day","hour","minute"])
df_daily = df.resample("D").mean()
full_index = pd.date_range(
start="2017-01-02",
end="2017-12-30",
freq="D")
df_daily_aligned = df_daily.reindex(full_index)
df_daily_aligned.index.name = "date"

##########################################################
##            Get all variables for analysis            ##
##########################################################

print('Getting all data for analysis\n')

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
    fig,ax = plt.subplots(1,1, figsize=(10,5))

    plt.suptitle(station,fontsize=14, fontweight='bold')

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
    
            
    # initialize arrays to save values
    o2vol_surf_unfiltered = []
    o2vol_deep_unfiltered = []
    vol_surf_unfiltered = []
    vol_deep_unfiltered = []

    # combine all months
    for month in months:
        fn = interface_type + '_' + station + '_' + month + '.p'
        df_bgc = pd.read_pickle(bgc_dir/fn)
        # conversion factor to go from mmol O2/hr to kmol O2/s
        conv = (1/1000) * (1/1000) * (1/60) * (1/60) # 1 mol/1000 mmol and 1 kmol/1000 mol and 1 hr/3600 sec
        # get (DO*V)
        o2vol_surf_unfiltered = np.concatenate((o2vol_surf_unfiltered, df_bgc['surf DO*V [mmol]'].values)) # mmol
        o2vol_deep_unfiltered = np.concatenate((o2vol_deep_unfiltered, df_bgc['deep DO*V [mmol]'].values)) # mmol
        # get volume
        vol_surf_unfiltered = np.concatenate((vol_surf_unfiltered, df_bgc['surf vol [m3]'].values)) # m3, and replace zeros with nan
        vol_deep_unfiltered = np.concatenate((vol_deep_unfiltered, df_bgc['deep vol [m3]'].values)) # m3


    # take time derivative of (DO*V) to get d/dt (DO*V)
    ddtDOV_total_unfiltered = np.diff((o2vol_deep_unfiltered+o2vol_surf_unfiltered)) * conv # diff gets us d(DO*V) dt, where t=1 hr (mmol/hr). Then * conv to get kmol/s
    # Get volume
    vol_total_unfiltered = vol_surf_unfiltered+vol_deep_unfiltered
    # Godin filter
    ddtDOV_onelayer_godin = zfun.lowpass(ddtDOV_total_unfiltered, f='godin')[36:-34:24]
    vol_total_godin = zfun.lowpass((vol_total_unfiltered), f='godin')[36:-34:24]

    # get volume-normalized d/dt(DO)
    conversion = (1000 * 32 * 60 * 60 * 24)
    ddtDO_volnorm_unfiltered = ddtDOV_total_unfiltered / vol_total_unfiltered[1::] * conversion # mg/L/day
    ddtDO_volnorm_godin = ddtDOV_onelayer_godin/vol_total_godin * conversion # mg/L/day

    # plot budget time series
    ax.plot(dates_local[2::],ddtDO_volnorm_unfiltered,color='pink',
            linewidth=0.5,label='d/dt(DO) unfiltered')
    ax.plot(dates_local_daily,ddtDO_volnorm_godin,color='deeppink',
            linewidth=2,label='d/dt(DO) Godin-filtered', zorder=5)
        
    # format budget time series figures
    ax.set_xlim([dates_hrly[0],dates_hrly[-2]])
    # zero line
    ax.plot([dates_hrly[0],dates_hrly[-2]],[0,0],color='black',linewidth=1)
    ax.grid(True,color='gainsboro',linewidth=1,linestyle='--',axis='both')
    ax.tick_params(axis='x', labelrotation=30, labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    loc = mdates.MonthLocator(interval=1)
    ax.xaxis.set_major_locator(loc)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.set_ylabel(r'DO change [mg L$^{-1}$ d$^{-1}$]', fontsize=14)
    ax.legend(loc='upper right', fontsize=12)