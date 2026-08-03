"""
Plot transport of NPZD+O terms all basins
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
year = '2017'

stations = 'all'

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

# Get basin boundaries:
# admiralty inlet, deception pass, hood canal, south sound, whidbey basin
boundaries = ['ss']#['ai','dp','hc','ss','wb']
boundary_dict = {'ai':'Admiralty Inlet', 'dp':'Deception Pass', 'hc':'Hood Canal',
                  'ss':'South Sound', 'wb':'Whidbey Basin'}

# state variables
vars = ['salt','oxygen', 'NO3', 'NH4', 'phytoplankton', 'zooplankton', 'SdetritusN', 'LdetritusN']


# # where to put output figures
# out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
# Lfun.make_dir(out_dir)

# create time_vector
dates_hrly = pd.date_range(start= startdate, end=enddate_hrly, freq= 'h')
dates_local = [pfun.get_dt_local(x) for x in dates_hrly]
dates_daily = pd.date_range(start= startdate, end=enddate, freq= 'd')[2::]
dates_local_daily = [pfun.get_dt_local(x) for x in dates_daily]
# crop time vector (because we only have jan 2 - dec 30)
dates_no_crop = dates_local_daily
dates_local_daily = dates_local_daily

print('\n')


##########################################################
##            Get all variables for analysis            ##
##########################################################

# basins
basins = ['South Sound','Whidbey Basin','Hood Canal','Main Basin']
basin_boundaries = {'South Sound':['ss'],
                    'Whidbey Basin':['wb','dp'],
                    'Hood Canal':['hc'],
                    'Main Basin':['ai','hc','ss','wb']}

# define function for generating subplots
def generate_axes(fig):
    gridspec = fig.add_gridspec(nrows=8, ncols=12, wspace=4, hspace=0)
    ax = {}
    ax['map']           = fig.add_subplot(gridspec[0:8, 0:4]) 
    ax['salt']          = fig.add_subplot(gridspec[0:2, 4:8])
    ax['NO3']           = fig.add_subplot(gridspec[2:4, 4:8],  sharex=ax['salt'])
    ax['NH4']           = fig.add_subplot(gridspec[2:4, 8:12], sharex=ax['salt'])
    ax['phytoplankton'] = fig.add_subplot(gridspec[4:6, 4:8],  sharex=ax['salt'])
    ax['oxygen']        = fig.add_subplot(gridspec[0:2, 8:12], sharex=ax['salt'])
    ax['zooplankton']   = fig.add_subplot(gridspec[4:6, 8:12], sharex=ax['salt'])
    ax['SdetritusN']    = fig.add_subplot(gridspec[6:8, 4:8],  sharex=ax['salt'])
    ax['LdetritusN']    = fig.add_subplot(gridspec[6:8, 8:12], sharex=ax['salt'])
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

# colors for basin boundaries
boundary_colors = ['black','red','orange','dodgerblue']
noloading_color = ['gray','red','orange','dodgerblue']
loading_colors = ['black','crimson','darkorange','blue']

plt.close('all')

# LOADING AND NO-LOADING TRANSPORTS ON SAME PLOT
for b,basin in enumerate(basins):
        
    # Initialize figure
    fig = plt.figure(figsize=(12,8))
    ax = generate_axes(fig)

    # hanning window size
    nwin = 10

    plt.suptitle('2017 '+basin+' fluxes ('+str(nwin)+'-day Hanning Window)\n+ is into basin, - is out of basin',fontsize=14)

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
    ax['map'].pcolormesh(plon, plat, np.where(mask == 0, np.nan, mask), vmin=0, vmax=1.5, cmap='bone' )
    # fix aspect ratio
    ax['map'].set_ylabel('Latitude', fontsize=12)
    ax['map'].set_xlabel('Longitude', fontsize=12)
    ax['map'].tick_params(axis='both', labelsize=12)
    ax['map'].set_title('Basin', loc='left', fontsize=12, fontweight='bold')
    pfun.dar(ax['map'])

    # get boundaries
    basin_bounds = basin_boundaries[basin]

    # loop through boundaries
    for i, boundary in enumerate(basin_bounds):

    # --------------------------- get TEF terms ----------------------------------------
        in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'cps' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (boundary + '.nc')
        bulk_loading = xr.open_dataset(in_dir)
        tef_df_loading, vn_list, vec_list = get_two_layer.get_two_layer(bulk_loading)

        in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1noDIN_x11b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (boundary + '.nc')
        bulk_NOloading = xr.open_dataset(in_dir) 
        tef_df_NOloading, vn_list, vec_list = get_two_layer.get_two_layer(bulk_NOloading)

        # determine which direction goes into the basin (positive) and which goes out (negative)
        flip = 1 # default is positive goes into the basin
        if basin == 'South Sound':
            if boundary == 'ss':
                flip = -1 # positive is going from South Sound to Main Basin, so we need to flip
        if basin == 'Hood Canal':
            if boundary == 'hc':
                flip = -1 # positive is going from Hood Canal to Main Basin, so we need to flip
        if basin == 'Whidbey Basin':
            if boundary == 'wb':
                flip = 1 # positive is going from Main Basin to Whidbey Basin, so NO flip
            elif boundary == 'dp':
                flip = -1 # positive is going from Whidbey to Straits, so we need to flip
        if basin == 'Main Basin':
            if boundary == 'ai':
                flip = 1 # positive is going from Straits into Main Basin, so NO flip
            elif boundary == 'hc':
                flip = 1 # positive is going from Hood Canal into Main Basin, so NO flip
            elif boundary == 'ss':
                flip = 1 # positive is going from South Sound into Main Basin, so NO flip
            elif boundary == 'wb':
                flip = -1 # positive is going from Whidbey into Main Basin, so we need to flip

        # get transports (positive goes into the basin)
        Qin_loading    = flip * tef_df_loading['q_p'] # Qin [m3/s]
        Qout_loading   = flip * tef_df_loading['q_m'] # Qout [m3/s]
        Qin_NOloading  = flip * tef_df_NOloading['q_p'] # Qin [m3/s]
        Qout_NOloading = flip * tef_df_NOloading['q_m'] # Qout [m3/s]

        # plot basin boundaries
        df = pd.read_pickle( Ldir['LOo'] / 'extract' / 'tef2' / 'sections_cas7_cps' / (boundary + '.p'))
        ax['map'].plot(df['x'], df['y'], color=boundary_colors[i], linewidth=2, alpha=1)

        # loop through variables and plot the transports
        for v, var in enumerate(vars):
            if var == 'salt':
                units = r'$[10^{3}\ m^{3}s^{-1}\ g\ kg^{-1}]$'
            else:
                units = '[mol/s]'
            ax[var].text(0.05,0.85,var + ' ' + units, fontsize=12, fontweight='bold',transform=ax[var].transAxes)

            # TEF flux
            # [mol/s] = [mmol/m3] * [m3/s] * 1/1000 mol/mmmol
            QinCin_loading     = tef_df_loading[var+'_p'].values * Qin_loading.values / 1000
            QoutCout_loading   = tef_df_loading[var+'_m'].values * Qout_loading.values / 1000
            QinCin_NOloading   = tef_df_NOloading[var+'_p'].values * Qin_NOloading.values / 1000
            QoutCout_NOloading = tef_df_NOloading[var+'_m'].values * Qout_NOloading.values / 1000

            # total TEF (inflow + outflow)
            ax[var].plot(dates_local_daily,zfun.lowpass(QinCin_NOloading + QoutCout_NOloading,n=nwin),
                    color=noloading_color[i],alpha=0.4,linewidth=3,label='No-Loading')
            ax[var].plot(dates_local_daily,zfun.lowpass(QinCin_loading + QoutCout_loading,n=nwin),
                    color=loading_colors[i],alpha=1,linewidth=1,linestyle='--',label='Loading')
            ax[var].axhline(y=0, color='black', alpha=0.5, linewidth=1, linestyle=':')
            
            if v < len(vars)-2:
                ax[var].tick_params(axis='x', which='both', labelbottom=False)
            else:
                ax[var].set_xlim([dates_local[0], dates_local[-1]])
                ax[var].xaxis.set_major_formatter(mdates.DateFormatter('%b'))

            # add a legend
            if i == 0 and v == 1:
                ax[var].legend(loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=2)
            

    plt.show()


# -----------------------------------------------------------------------------

# # DIFFERENCE BETWEEN LOADING AND NO-LOADING TRANSPORTS (% of no-loading transport)
# for b,basin in enumerate(basins):
        
#     # Initialize figure
#     fig = plt.figure(figsize=(12,8))
#     ax = generate_axes(fig)

#     plt.suptitle(r'Loading $-$ No-loading'+'\n'+r'($\%$ change relative to max no-loading value)',fontsize=14)

#     # add map
#     # All Puget Sound
#     ax['map'].pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho), vmin=0, vmax=1.1, cmap='bone' )
#     # basins
#     if basin == 'Hood Canal':
#         mask = mask_hc
#     elif basin == 'South Sound':    
#         mask = mask_ss
#     elif basin == 'Whidbey Basin':
#         mask = mask_wb
#     elif basin == 'Main Basin':
#         mask = mask_mb
#     ax['map'].pcolormesh(plon, plat, np.where(mask == 0, np.nan, mask), vmin=0, vmax=1.5, cmap='bone' )
#     # fix aspect ratio
#     ax['map'].set_ylabel('Latitude', fontsize=12)
#     ax['map'].set_xlabel('Longitude', fontsize=12)
#     ax['map'].tick_params(axis='both', labelsize=12)
#     ax['map'].set_title('Basin', loc='left', fontsize=12, fontweight='bold')
#     pfun.dar(ax['map'])

#     # get boundaries
#     basin_bounds = basin_boundaries[basin]

#     # loop through boundaries
#     for i, boundary in enumerate(basin_bounds):

#     # --------------------------- get TEF terms ----------------------------------------
#         in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'cps' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (boundary + '.nc')
#         bulk_loading = xr.open_dataset(in_dir)
#         tef_df_loading, vn_list, vec_list = get_two_layer.get_two_layer(bulk_loading)

#         in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1noDIN_x11b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (boundary + '.nc')
#         bulk_NOloading = xr.open_dataset(in_dir) 
#         tef_df_NOloading, vn_list, vec_list = get_two_layer.get_two_layer(bulk_NOloading)

#         # determine which direction goes into the basin (positive) and which goes out (negative)
#         flip = 1 # default is positive goes into the basin
#         if basin == 'South Sound':
#             if boundary == 'ss':
#                 flip = -1 # positive is going from South Sound to Main Basin, so we need to flip
#         if basin == 'Hood Canal':
#             if boundary == 'hc':
#                 flip = -1 # positive is going from Hood Canal to Main Basin, so we need to flip
#         if basin == 'Whidbey Basin':
#             if boundary == 'wb':
#                 flip = 1 # positive is going from Main Basin to Whidbey Basin, so NO flip
#             elif boundary == 'dp':
#                 flip = -1 # positive is going from Whidbey to Straits, so we need to flip
#         if basin == 'Main Basin':
#             if boundary == 'ai':
#                 flip = 1 # positive is going from Straits into Main Basin, so NO flip
#             elif boundary == 'hc':
#                 flip = 1 # positive is going from Hood Canal into Main Basin, so NO flip
#             elif boundary == 'ss':
#                 flip = 1 # positive is going from South Sound into Main Basin, so NO flip
#             elif boundary == 'wb':
#                 flip = -1 # positive is going from Whidbey into Main Basin, so we need to flip

#         # get transports (positive goes into the basin)
#         Qin_loading    = flip * tef_df_loading['q_p'] # Qin [m3/s]
#         Qout_loading   = flip * tef_df_loading['q_m'] # Qout [m3/s]
#         Qin_NOloading  = flip * tef_df_NOloading['q_p'] # Qin [m3/s]
#         Qout_NOloading = flip * tef_df_NOloading['q_m'] # Qout [m3/s]

#         # plot basin boundaries
#         df = pd.read_pickle( Ldir['LOo'] / 'extract' / 'tef2' / 'sections_cas7_cps' / (boundary + '.p'))
#         ax['map'].plot(df['x'], df['y'], color=boundary_colors[i], linewidth=2, alpha=1)

#         # loop through variables and plot the transports
#         for v, var in enumerate(vars):
#             ax[var].text(0.05,0.85,var, fontsize=12, fontweight='bold',transform=ax[var].transAxes)

#             QinCin_loading     = tef_df_loading[var+'_p'].values * Qin_loading.values
#             QoutCout_loading   = tef_df_loading[var+'_m'].values * Qout_loading.values
#             QinCin_NOloading   = tef_df_NOloading[var+'_p'].values * Qin_NOloading.values
#             QoutCout_NOloading = tef_df_NOloading[var+'_m'].values * Qout_NOloading.values

#             # difference between loading and no-loading
#             diff = (QinCin_loading + QoutCout_loading) - (QinCin_NOloading + QoutCout_NOloading)
#             diff_norm = diff / np.nanmax(np.abs(QinCin_NOloading + QoutCout_NOloading)) * 100
#             ax[var].plot(dates_local_daily,zfun.lowpass(diff_norm,n=10),
#                     color=noloading_color[i],alpha=1,linewidth=2,label='Loading minus No-Loading')
#             ax[var].axhline(y=0, color='black', alpha=0.5, linewidth=1, linestyle=':')
            
#             if v < len(vars)-2:
#                 ax[var].tick_params(axis='x', which='both', labelbottom=False)
#             else:
#                 ax[var].set_xlim([dates_local[0], dates_local[-1]])
#                 ax[var].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
            

#     plt.show()

