"""
Compare average bottom DO between multiple years
(Set up to run for 6 years)

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator
import matplotlib.patches as patches
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import pandas as pd
import cmocean
import matplotlib.pylab as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

##############################################################
##                       USER INPUTS                        ##
##############################################################

remove_straits = True

WWTP_loc = False

# Hanning window length
nwin = 20

# years =  ['2015']
years =  ['2015','2016','2017','2018','2019','2020']

# which  model run to look at?
gtagexes = ['cas7_t1_x11ab','cas7_t1noDIN_x11ab'] 

# where to put output figures
out_dir = Ldir['LOo'] / 'chapter_2' / 'figures'
Lfun.make_dir(out_dir)

regions = ['Hood Canal', 'South Sound', 'Whidbey Basin', 'Main Basin', 'All Puget Sound']
plt.close('all')

##############################################################
##                    HELPER FUNCTIONS                      ##
##############################################################

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = '../../../LO_data/trapsD00/LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_LO = repeatrivs_df.loc[repeatrivs_df['SSM_rname'] == rname, 'LO_rname'].values[0]
    return rname_LO

def LO2SSM_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = Ldir['data'] / 'trapsD00' / 'LiveOcean_SSM_rivers.xlsx'
    repeatrivs_df = pd.read_excel(repeatrivs_fn)
    rname_SSM = repeatrivs_df.loc[repeatrivs_df['LO_rname'] == rname, 'SSM_rname'].values[0]
    return rname_SSM


if WWTP_loc == True:
    # set up the time index for the record
    Ldir = Lfun.Lstart()
    dsf = Ldir['ds_fmt']
    dt0 = datetime.strptime('2020.01.01',dsf)
    dt1 = datetime.strptime('2020.12.31',dsf)
    days = (dt0, dt1)
        
    # pandas Index objects
    dt_ind = pd.date_range(start=dt0, end=dt1)
    yd_ind = pd.Index(dt_ind.dayofyear)

    # Get LiveOcean grid info --------------------------------------------------

    # get the grid data
    ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
    z = -ds.h.values
    mask_rho = np.transpose(ds.mask_rho.values)
    lon = ds.lon_rho.values
    lat = ds.lat_rho.values
    X = lon[0,:] # grid cell X values
    Y = lat[:,0] # grid cell Y values
    plon, plat = pfun.get_plon_plat(lon,lat)
    # make a version of z with nans where masked
    zm = z.copy()
    zm[np.transpose(mask_rho) == 0] = np.nan
    zm[np.transpose(mask_rho) != 0] = -1

    # get flow, nitrate, and ammonium values
    fp_wwtps = '../../../LO_output/pre/trapsP01/moh20_wwtps/lo_base/Data_historical/'
    moh20_flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow.p')    # m3/s
    moh20_no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3.p')      # mmol/m3
    moh20_nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4.p')      # mmol/m3

    fp_wwtps = '../../../LO_output/pre/trapsP01/was24_wwtps/lo_base/Data_historical/'
    was24_flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow.p')    # m3/s
    was24_no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3.p')      # mmol/m3
    was24_nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4.p')      # mmol/m3

    # calculate total DIN concentration in mg/L
    moh20_dindf_wwtps = (moh20_no3df_wwtps + moh20_nh4df_wwtps)/71.4    # mg/L
    was24_dindf_wwtps = (was24_no3df_wwtps + was24_nh4df_wwtps)/71.4    # mg/L

    # calculate daily loading timeseries in kg/d
    moh20_dailyloaddf_wwtps = 86.4*moh20_dindf_wwtps*moh20_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s
    was24_dailyloaddf_wwtps = 86.4*was24_dindf_wwtps*was24_flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s

    # calculate average daily load over the year (kg/d)
    moh20_avgload_wwtps = moh20_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')
    was24_avgload_wwtps = was24_dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

    # add row and col index for plotting on LiveOcean grid
    griddf0_wwtps = pd.read_csv('../../../LO_data/grids/cas7/moh20_wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    moh20_avgload_wwtps = moh20_avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    moh20_avgload_wwtps = moh20_avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    griddf0_wwtps = pd.read_csv('../../../LO_data/grids/cas7/was24_wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    was24_avgload_wwtps = was24_avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    was24_avgload_wwtps = was24_avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    # get point source lat and lon
    moh20_lon_wwtps = [X[int(col)] for col in moh20_avgload_wwtps['col_py']]
    moh20_lat_wwtps = [Y[int(row)] for row in moh20_avgload_wwtps['row_py']]
    was24_lon_wwtps = [X[int(col)] for col in was24_avgload_wwtps['col_py']]
    was24_lat_wwtps = [Y[int(row)] for row in was24_avgload_wwtps['row_py']]
    
    # define marker sizes (minimum size is 10 so dots don't get too small)
    moh20_sizes_wwtps = [max(0.05*load,5) for load in moh20_avgload_wwtps['avg-daily-load(kg/d)']]
    was24_sizes_wwtps = [max(0.05*load,5) for load in was24_avgload_wwtps['avg-daily-load(kg/d)']]


##############################################################
##                      PROCESS DATA                        ##
##############################################################

# read in masks
basin_mask_ds = xr.open_dataset('../../../LO_output/chapter_2/data/basin_masks_from_pugetsoundDObox.nc')
mask_rho = basin_mask_ds.mask_rho.values
mask_hc = basin_mask_ds.mask_hoodcanal.values
mask_ss = basin_mask_ds.mask_southsound.values
mask_wb = basin_mask_ds.mask_whidbeybasin.values
mask_mb = basin_mask_ds.mask_mainbasin.values
mask_ps = basin_mask_ds.mask_pugetsound.values
lon = basin_mask_ds['lon_rho'].values
lat = basin_mask_ds['lat_rho'].values
h = basin_mask_ds['h'].values
plon, plat = pfun.get_plon_plat(lon,lat)

##############################################################
# get average concentration per basin

# initialize empty dictionaries and fill with vertical integrals
NO3_vert_dict = {}
phyto_vert_dict = {}
zoop_vert_dict = {}
NH4_vert_dict = {}
Ldet_vert_dict = {}
Sdet_vert_dict = {}

for year in years:
    for gtagex in gtagexes:
        ds = xr.open_dataset(Ldir['LOo'] / 'chapter_2' / 'data' / (gtagex + '_pugetsoundDO_' + year + '_NPZD_vert_ints.nc'))
        NO3_vert_int = ds['NO3_vert_int'].values
        phyto_vert_int = ds['phyto_vert_int'].values
        zoop_vert_int = ds['zoop_vert_int'].values
        NH4_vert_int = ds['NH4_vert_int'].values
        LdetritusN_vert_int = ds['LdetritusN_vert_int'].values
        SdetritusN_vert_int = ds['SdetritusN_vert_int'].values
        # add data to dictionaries
        NO3_vert_dict[gtagex+year] = NO3_vert_int
        phyto_vert_dict[gtagex+year] = phyto_vert_int
        zoop_vert_dict[gtagex+year] = zoop_vert_int
        NH4_vert_dict[gtagex+year] = NH4_vert_int
        Ldet_vert_dict[gtagex+year] = LdetritusN_vert_int
        Sdet_vert_dict[gtagex+year] = SdetritusN_vert_int

# grid cell areas
fp = Ldir['LOo'] / 'extract' / 'cas7_t1_x11ab' / 'box' / ('pugetsoundDO_2014.01.01_2014.12.31.nc')
box_ds = xr.open_dataset(fp)
DX = (box_ds.pm.values)**-1
DY = (box_ds.pn.values)**-1
DA = DX*DY # get area in m2

# initialize dictionary for average concentration (volume integrals [mol], normalized by volume)
NO3_vol_int = {}
phyto_vol_int = {}
zoop_vol_int = {}
NH4_vol_int = {}
Ldet_vol_int = {}
Sdet_vol_int = {}
TN_vol_int = {}

for year in years:
    for region in regions:

        # get mask for the region
        if region == 'Hood Canal':
            mask = mask_hc
        elif region == 'South Sound':
            mask = mask_ss
        elif region == 'Whidbey Basin':
            mask = mask_wb
        elif region == 'Main Basin':
            mask = mask_mb
        elif region == 'All Puget Sound':
            mask = mask_ps

        for gtagex in gtagexes:

            NO3_vert_int = NO3_vert_dict[gtagex+year]
            NO3_vert_int_masked = NO3_vert_int * mask
            NO3_vol_timeseries = np.sum(NO3_vert_int_masked * DA, axis=(1, 2)) # [mol]
            NO3_vol_int[gtagex+region+year] = NO3_vol_timeseries # [mol]

            NH4_vert_int = NH4_vert_dict[gtagex+year]
            NH4_vert_int_masked = NH4_vert_int * mask
            NH4_vol_timeseries = np.sum(NH4_vert_int_masked * DA, axis=(1, 2)) # [mol]
            NH4_vol_int[gtagex+region+year] = NH4_vol_timeseries # [mol]

            phyto_vert_int = phyto_vert_dict[gtagex+year]
            phyto_vert_int_masked = phyto_vert_int * mask
            phyto_vol_timeseries = np.sum(phyto_vert_int_masked * DA, axis=(1, 2)) # [mol]
            phyto_vol_int[gtagex+region+year] = phyto_vol_timeseries # [mol]

            zoop_vert_int = zoop_vert_dict[gtagex+year]
            zoop_vert_int_masked = zoop_vert_int * mask
            zoop_vol_timeseries = np.sum(zoop_vert_int_masked * DA, axis=(1, 2)) # [mol]
            zoop_vol_int[gtagex+region+year] = zoop_vol_timeseries # [mol]


            Ldet_vert_int = Ldet_vert_dict[gtagex+year]
            Ldet_vert_int_masked = Ldet_vert_int * mask
            Ldet_vol_timeseries = np.sum(Ldet_vert_int_masked * DA, axis=(1, 2)) # [mol]
            Ldet_vol_int[gtagex+region+year] = Ldet_vol_timeseries # [mol]

            Sdet_vert_int = Sdet_vert_dict[gtagex+year]
            Sdet_vert_int_masked = Sdet_vert_int * mask
            Sdet_vol_timeseries = np.sum(Sdet_vert_int_masked * DA, axis=(1, 2)) # [mol]
            Sdet_vol_int[gtagex+region+year] = Sdet_vol_timeseries # [mol]

            # get total nitrogen [mol]
            TN_vol_int_timeseries = NO3_vol_timeseries + NH4_vol_timeseries + phyto_vol_timeseries + zoop_vol_timeseries + Ldet_vol_timeseries + Sdet_vol_timeseries
            TN_vol_int[gtagex+region+year] = TN_vol_int_timeseries

# initialize dictionary for ratio of TN
NO3_TN_ratio = {}
phyto_TN_ratio = {}
zoop_TN_ratio = {}
NH4_TN_ratio = {}
Ldet_TN_ratio = {}
Sdet_TN_ratio = {}

for year in years:
    for region in regions:
        for gtagex in gtagexes:
                NO3_TN_ratio[gtagex+region+year] = NO3_vol_int[gtagex+region+year]     / TN_vol_int[gtagex+region+year]
                NH4_TN_ratio[gtagex+region+year] = NH4_vol_int[gtagex+region+year]     / TN_vol_int[gtagex+region+year]
                phyto_TN_ratio[gtagex+region+year] = phyto_vol_int[gtagex+region+year] / TN_vol_int[gtagex+region+year]
                zoop_TN_ratio[gtagex+region+year] = zoop_vol_int[gtagex+region+year]   / TN_vol_int[gtagex+region+year]
                Ldet_TN_ratio[gtagex+region+year] = Ldet_vol_int[gtagex+region+year]   / TN_vol_int[gtagex+region+year]
                Sdet_TN_ratio[gtagex+region+year] = Sdet_vol_int[gtagex+region+year]   / TN_vol_int[gtagex+region+year]


##############################################################
##          Plot fraction of TN of each species             ##
##############################################################

# variables
vars = ['NO3','NH4','Phytoplankton','Zooplankton','Large Detritus','Small Detritus']
# ratios
ratios = [NO3_TN_ratio, NH4_TN_ratio, phyto_TN_ratio, zoop_TN_ratio, Ldet_TN_ratio, Sdet_TN_ratio]

# define function for generating subplots
def generate_axes(fig):
    gridspec = fig.add_gridspec(nrows=6, ncols=12, wspace=4, hspace=0)
    ax = {}
    ax['map']           = fig.add_subplot(gridspec[0:6, 0:4])
    ax['NO3']           = fig.add_subplot(gridspec[0:2, 4:8])
    ax['NH4']           = fig.add_subplot(gridspec[0:2, 8:12], sharex=ax['NO3'])
    ax['Phytoplankton'] = fig.add_subplot(gridspec[2:4, 4:8],  sharex=ax['NO3'])
    ax['Zooplankton']   = fig.add_subplot(gridspec[2:4, 8:12], sharex=ax['NO3'])
    ax['Small Detritus']    = fig.add_subplot(gridspec[4:6, 4:8],  sharex=ax['NO3'])
    ax['Large Detritus']    = fig.add_subplot(gridspec[4:6, 8:12], sharex=ax['NO3'])
    return ax

# read in masks
basin_mask_ds = xr.open_dataset('../../../LO_output/chapter_2/data/basin_masks_from_pugetsoundDObox.nc')
mask_rho = basin_mask_ds.mask_rho.values
mask_hc = basin_mask_ds.mask_hoodcanal.values
mask_ss = basin_mask_ds.mask_southsound.values
mask_wb = basin_mask_ds.mask_whidbeybasin.values
mask_mb = basin_mask_ds.mask_mainbasin.values
mask_ps = basin_mask_ds.mask_pugetsound.values
lon = basin_mask_ds['lon_rho'].values
lat = basin_mask_ds['lat_rho'].values
h = basin_mask_ds['h'].values
plon, plat = pfun.get_plon_plat(lon,lat)

plt.close('all')

# LOADING AND NO-LOADING TRANSPORTS ON SAME PLOT
for b,basin in enumerate(regions):
        
    # Initialize figure
    fig = plt.figure(figsize=(12,8))
    ax = generate_axes(fig)

    # hanning window size
    nwin = 10

    plt.suptitle('2017 '+basin+' NPZD:TN ratio ('+str(nwin)+'-day Hanning Window)',fontsize=14)

    # add map
    # All Puget Sound
    ax['map'].pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho), vmin=0, vmax=1.1, cmap='bone' )
    # basins
    if basin == 'Hood Canal':
        mask = mask_hc
    elif basin == 'South Sound':    
        mask = mask_ss
    elif basin == 'Whidbey Basin':
        mask = mask_wb
    elif basin == 'Main Basin':
        mask = mask_mb
    elif basin == 'All Puget Sound':
        mask = mask_ps
    ax['map'].pcolormesh(plon, plat, np.where(mask == 0, np.nan, mask), vmin=0, vmax=1.5, cmap='bone' )
    # fix aspect ratio
    ax['map'].set_ylabel('Latitude', fontsize=12)
    ax['map'].set_xlabel('Longitude', fontsize=12)
    ax['map'].tick_params(axis='both', labelsize=12)
    ax['map'].set_title('Basin', loc='left', fontsize=12, fontweight='bold')
    pfun.dar(ax['map'])

    # add wwtp locations
    if WWTP_loc == True:
        edgecolor = 'black'
        facecolor = 'none'
        alpha = 1
        ax['map'].scatter(moh20_lon_wwtps,moh20_lat_wwtps,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                     linewidth=1, s=moh20_sizes_wwtps, label='WWTPs')
        ax['map'].scatter(was24_lon_wwtps,was24_lat_wwtps,color=facecolor, edgecolors=edgecolor, alpha=alpha,
                     linewidth=1, s=was24_sizes_wwtps)
        leg_szs = [100, 1000, 10000]
        szs = [0.05*(leg_sz) for leg_sz in leg_szs]
        l0 = plt.scatter([],[], s=szs[0], color=facecolor, alpha=alpha, edgecolors=edgecolor, linewidth=1)
        l1 = plt.scatter([],[], s=szs[1], color=facecolor, alpha=alpha, edgecolors=edgecolor, linewidth=1)
        l2 = plt.scatter([],[], s=szs[2], color=facecolor, alpha=alpha, edgecolors=edgecolor, linewidth=1)
        labels = ['< 100', '1,000', '10,000']
        legend = ax['map'].legend([l0, l1, l2], labels, fontsize = 10, markerfirst=False,
            title='WWTP loading \n'+r' (kg N d$^{-1}$)',loc='upper left', labelspacing=1, borderpad=0.8)
        plt.setp(legend.get_title(),fontsize=9)

    # loop through variables and plot the transports
    for v, var in enumerate(vars):
        ax[var].text(0.05,0.85,var, fontsize=12, fontweight='bold',transform=ax[var].transAxes)

        for y,year in enumerate(years):

            # create time vector
            startdate = year+'.01.01'
            enddate = year+'.12.31'
            enddate_hrly = str(int(year)+1)+'.01.01 00:00:00'
            # create time_vector
            dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
            dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
            dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')
            dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]

            # ratio of NPZD variable to TN
            ax[var].plot(dates_local_daily,zfun.lowpass(ratios[v]['cas7_t1noDIN_x11ab'+basin+year],n=nwin),
                    color='cadetblue',alpha=0.4,linewidth=3,label='No-Loading')
            ax[var].plot(dates_local_daily,zfun.lowpass(ratios[v]['cas7_t1_x11ab'+basin+year],n=nwin),
                    color='teal',alpha=1,linewidth=1,linestyle='--',label='Loading')
            # ax[var].axhline(y=0, color='black', alpha=0.5, linewidth=1, linestyle=':')
        
        if v < len(vars)-2:
            ax[var].tick_params(axis='x', which='both', labelbottom=False)
        else:
            startdate = years[0]+'.01.01'
            enddate = years[-1]+'.12.31'
            enddate_hrly = str(int(year)+1)+'.01.01 00:00:00'
            # create time_vector
            dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
            dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
            ax[var].set_xlim([dates_local[0], dates_local[-1]])
            ax[var].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

        # add a legend
        if y==0 and v == 1:
            ax[var].legend(loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=2)

    plt.show()
