"""
Plot transport of NPZD+O terms all basins
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import get_two_layer
from matplotlib.patches import Polygon
from matplotlib.transforms import Affine2D

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
boundaries = ['ai','dp','hc','ss','wb']
boundary_dict = {'ai':'Admiralty Inlet', 'dp':'Deception Pass', 'hc':'Hood Canal',
                  'ss':'South Sound', 'wb':'Whidbey Basin'}

# state variables
vars = ['oxygen', 'NO3', 'NH4', 'phytoplankton', 'zooplankton', 'SdetritusN', 'LdetritusN'] # 'salt'


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
##           Helper Functions to draw arrow             ##
##########################################################

def polygon_area(verts):
    x, y = verts[:, 0], verts[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

def polygon_centroid(verts):
    """
    Centroid of a non-self-intersecting polygon (shoelace formula).
    """
    x = verts[:, 0]
    y = verts[:, 1]
    x2 = np.roll(x, -1)
    y2 = np.roll(y, -1)

    cross = x * y2 - x2 * y
    A = 0.5 * np.sum(cross)
    Cx = (1 / (6 * A)) * np.sum((x + x2) * cross)
    Cy = (1 / (6 * A)) * np.sum((y + y2) * cross)
    return np.array([Cx, Cy])

def make_centered_unit_arrow():
    """
    Arrow points along +x. Geometry has rectangular shaft + triangular head.
    Then we recenter polygon so centroid is exactly at (0,0).
    """
    L_shaft = 0.25
    L_head  = 0.25
    W_shaft = 0.22
    W_head  = 0.46

    # Tail at x=0, tip at x=L_shaft+L_head before recentering
    verts = np.array([
        [0.0,             -W_shaft/2],   # tail lower
        [L_shaft,         -W_shaft/2],   # shaft lower-right
        [L_shaft,         -W_head/2],    # head base lower
        [L_shaft+L_head,   0.0],         # tip
        [L_shaft,          W_head/2],    # head base upper
        [L_shaft,          W_shaft/2],   # shaft upper-right
        [0.0,              W_shaft/2],   # tail upper
    ], dtype=float)

    # Shift so centroid is at origin: this makes (x,y) the true center anchor
    c = polygon_centroid(verts)
    verts -= c
    return verts

def add_flux_arrow(ax, x, y, angle_deg, flux, k_area=0.04,facecolor='tab:red', edgecolor='none', lw=0.6, alpha=0.6):
    """
    Draw one flux arrow centered exactly on (x,y).
    Area is proportional to flux: A = k_area * flux.
    """
    base = make_centered_unit_arrow()
    A0 = polygon_area(base)
    A_target = k_area * flux

    # Uniform scaling -> preserves angles (keeps right angles right)
    s = np.sqrt(A_target / A0)

    trans = (Affine2D()
             .scale(s, s)            # isotropic scale only
             .rotate_deg(angle_deg)  # rigid rotation
             .translate(x, y))       # center lands at (x,y)

    patch = Polygon(base, closed=True,facecolor=facecolor, edgecolor=edgecolor, lw=lw, alpha=alpha,transform=trans + ax.transData)
    ax.add_patch(patch)

##########################################################
##            Get all variables for analysis            ##
##########################################################

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

# Initialize figure
fig,ax = plt.subplots(1,len(vars),figsize=(18,8),sharey=True,sharex=True)

fig.suptitle('2017 mean fluxes',fontsize=14)

# add map
for axis in ax:
    # All Puget Sound
    axis.pcolormesh(plon, plat, np.where(mask_rho == 0, np.nan, mask_rho),
        vmin=0, vmax=1.2, cmap='bone' )
    # fix aspect ratio
    axis.set_xlabel('Longitude', fontsize=12)
    axis.tick_params(axis='both', labelsize=12)
    # axis.set_ylim([46.95,48.5])
    pfun.dar(axis)
ax[0].set_ylabel('Latitude', fontsize=12)

# loop through variables and plot the transports
for v, var in enumerate(vars):

    ai_norm = 0

    # loop through boundaries
    for i, boundary in enumerate(boundaries):

    # --------------------------- get TEF terms ----------------------------------------
        in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'cps' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (boundary + '.nc')
        bulk_loading = xr.open_dataset(in_dir)
        tef_df_loading, vn_list, vec_list = get_two_layer.get_two_layer(bulk_loading)

        in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1noDIN_x11b' / 'tef2' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (boundary + '.nc')
        bulk_NOloading = xr.open_dataset(in_dir) 
        tef_df_NOloading, vn_list, vec_list = get_two_layer.get_two_layer(bulk_NOloading)

        # # determine which direction goes into the basin (positive) and which goes out (negative)

        # get transports (positive goes into the basin)
        Qin_loading    = tef_df_loading['q_p'] # Qin [m3/s]
        Qout_loading   = tef_df_loading['q_m'] # Qout [m3/s]
        Qin_NOloading  = tef_df_NOloading['q_p'] # Qin [m3/s]
        Qout_NOloading = tef_df_NOloading['q_m'] # Qout [m3/s]

        QinCin_loading     = tef_df_loading[var+'_p'].values * Qin_loading.values
        QoutCout_loading   = tef_df_loading[var+'_m'].values * Qout_loading.values
        QinCin_NOloading   = tef_df_NOloading[var+'_p'].values * Qin_NOloading.values
        QoutCout_NOloading = tef_df_NOloading[var+'_m'].values * Qout_NOloading.values

        # calculate annual mean fluxes
        loading_mean     = np.nanmean(QinCin_loading+QoutCout_loading)
        noloading_mean   = np.nanmean(QinCin_NOloading+QoutCout_NOloading)
        diff = loading_mean - noloading_mean
        # arrow scale (normalize by biggest value of admiralty inlet no-loading run)
        if boundary == 'ai':
            ai_norm = np.max(np.abs(noloading_mean))
        scale = 1/ai_norm * 5e-2

        # get boundary location
        df = pd.read_pickle( Ldir['LOo'] / 'extract' / 'tef2' / 'sections_cas7_cps' / (boundary + '.p'))
        x_lon = df['x'][round(len(df['x'])/2)]
        y_lat = df['y'][round(len(df['y'])/2)]
        # plot basin boundaries
        # ax[0].plot(df['x'], df['y'], color='teal', linewidth=1, alpha=1)

        # determine direction and rotation of arrow
        # get rotation angle for arrow
        if boundary == 'ai':
            angle = -45
        elif boundary == 'dp':
            angle = 180
        elif boundary == 'hc':
            angle = 35
        elif boundary == 'ss':
            angle = 50
        elif boundary == 'wb':
            angle = 48
        # flip angle if the flux is negative
        if noloading_mean < 0:
            meanflip = 180
        else:
            meanflip = 0
        # get direction and color of change due to WWTP loading
        diffflip = meanflip
        if np.sign(loading_mean) ==  np.sign(noloading_mean) and np.abs(loading_mean) > np.abs(noloading_mean):
            # no flip
            diffflip += 0
            diff_color = 'red'
        else:
            # flip
            diffflip += 180
            diff_color = 'blue'

        # plot arrows
        # No-loading run
        add_flux_arrow(ax[v], x_lon, y_lat, angle + meanflip, np.abs(noloading_mean), k_area=scale, facecolor='gray', edgecolor='gray')
        # Loading minus no-loading
        # ax[1].scatter(x_lon, y_lat, s=np.abs(percentage_change)*50, facecolor=diff_color, alpha=0.6)
        add_flux_arrow(ax[v], x_lon, y_lat, angle + diffflip, np.abs(diff), k_area=scale, facecolor=diff_color)

        # set titles
        # ax[0].set_title('No-Loading', fontsize=12)
        # ax[1].set_title(r'Loading $-$ No-Loading', fontsize=12)
        ax[v].set_title(var + '\n' + str(round(ai_norm/1000,2)) + ' mol/s',fontsize=12)
        

    #     # total TEF (inflow + outflow)
    #     ax[var].plot(dates_local_daily,zfun.lowpass(QinCin_NOloading + QoutCout_NOloading,n=10),
    #             color='cadetblue',alpha=0.4,linewidth=3,label='No-Loading')
    #     ax[var].plot(dates_local_daily,zfun.lowpass(QinCin_loading + QoutCout_loading,n=10),
    #             color='teal',alpha=1,linewidth=1,linestyle='--',label='Loading')
    #     ax[var].axhline(y=0, color='black', alpha=0.5, linewidth=1, linestyle=':')
        
    #     if v < len(vars)-2:
    #         ax[var].tick_params(axis='x', which='both', labelbottom=False)
    #     else:
    #         ax[var].set_xlim([dates_local[0], dates_local[-1]])
    #         ax[var].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
            

    plt.subplots_adjust(left=0.05,right=0.99,wspace=0)
    plt.show()

