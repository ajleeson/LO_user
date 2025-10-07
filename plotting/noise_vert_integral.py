"""
Code to compare runs to evaluate noise in runs with
biogeochemical forcing differences.
"""

from subprocess import Popen as Po
from subprocess import PIPE as Pi
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean
from datetime import datetime
from mpl_toolkits.axes_grid1 import make_axes_locatable

from lo_tools import Lfun, zrfun, zfun
from lo_tools import plotting_functions as pfun

import sys

plt.close('all')

Ldir = Lfun.Lstart()

###################################################################
##                          User Inputs                          ##  
################################################################### 

# USER OPTIONS ----------------------------------------------------

d0 = '2012.10.07'
d1 = '2013.09.27'

list_type = 'weeklyaverage' #'weekly', 'daily', 'hourly ', 'allhours'


# # dstr = 'f2012.10.07'
# # dstr = 'f2013.04.04'
# dstr = 'f2013.12.27'

# runs
# gtx1 = 'cas7_t1_x11ab'
# gtx2 = 'cas7_t1noDIN_x11ab'
# his = 'ocean_his_0002.nc'

# gtagex of files to difference
Ldir_WWTP   = Lfun.Lstart(gridname='cas7', tag='t1', ex_name='x11ab')
Ldir_noWWTP = Lfun.Lstart(gridname='cas7', tag='t1noDIN', ex_name='x11ab')


###################################################################
##       Custom get fn list to do weekly from average files      ##  
################################################################### 

# format used for naming day folders
ds_fmt = '%Y.%m.%d'

def get_fn_list(list_type, Ldir, ds0, ds1, his_num=2):
    """
    INPUT:
    A function for getting lists of history files.
    List items are Path objects
    
    NEW 2023.10.05: for list_type = 'hourly', if you pass his_num = 1
    it will start with ocean_his_0001.nc on the first day instead of the default which
    is to start with ocean_his_0025.nc on the day before.

    NEW 2025.06.20: for list_type = 'hourly0'
    which will start with ocean_his_0001.nc on the first day instead of the default which
    is to start with ocean_his_0025.nc on the day before.
    This is identical to passing his_num = 1, but may be more convenient, especially
    as we move to "continuation" start_type, which always writes an 0001 file.
    """
    dt0 = datetime.strptime(ds0, ds_fmt)
    dt1 = datetime.strptime(ds1, ds_fmt)
    dir0 = Ldir['roms_out'] / Ldir['gtagex']
    if list_type == 'snapshot':
        # a single file name in a list
        his_string = ('0000' + str(his_num))[-4:]
        fn_list = [dir0 / ('f' + ds0) / ('ocean_his_' + his_string + '.nc')]
    elif list_type == 'hourly':
        # list of hourly files over a date range
        fn_list = Lfun.fn_list_utility(dt0,dt1,Ldir,his_num=his_num)
    elif list_type == 'hourly0':
        # list of hourly files over a date range, starting with 0001 of dt0.
        fn_list = Lfun.fn_list_utility(dt0,dt1,Ldir,his_num=1)
    elif list_type == 'daily':
        # list of history file 21 (Noon PST) over a date range
        fn_list = []
        date_list = Lfun.date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_his_0021.nc'
            fn_list.append(fn)
    elif list_type == 'lowpass':
        # list of lowpassed files (Noon PST) over a date range
        fn_list = []
        date_list = Lfun.date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'lowpassed.nc'
            fn_list.append(fn)
    elif list_type == 'average':
        # list of daily averaged files (Noon PST) over a date range
        fn_list = []
        date_list = Lfun.date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_avg_0001.nc'
            fn_list.append(fn)
    elif list_type == 'weekly':
        # like "daily" but at 7-day intervals
        fn_list = []
        date_list = Lfun.date_list_utility(dt0, dt1, daystep=7)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_his_0021.nc'
            fn_list.append(fn)
    elif list_type == 'weeklyaverage':
        # like "daily" but at 7-day intervals
        fn_list = []
        date_list = Lfun.date_list_utility(dt0, dt1, daystep=7)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_his_0002.nc'
            fn_list.append(fn)
    elif list_type == 'allhours':
        # a list of all the history files in a directory
        # (this is the only list_type that actually finds files)
        in_dir = dir0 / ('f' + ds0)
        fn_list = [ff for ff in in_dir.glob('ocean_his*nc')]
        fn_list.sort()

    return fn_list

# User Parker's LO_roms output on apogee for the with loading run (long hindcast)
Ldir_WWTP['roms_out'] = Ldir['roms_out5']

# get list of history files to plot
fn_list_WWTP   = get_fn_list(list_type, Ldir_WWTP, d0, d1)
fn_list_noWWTP = get_fn_list(list_type, Ldir_noWWTP, d0, d1)

# prepare output directory for results
outdir = Ldir['LOo'] / 'AL_custom_plots' / ('noise_test_loading_runs_' + list_type)
Lfun.make_dir(outdir)

###################################################################
##                      Loop through files                       ##  
################################################################### 

for i,fn_WWTP in enumerate(fn_list_WWTP):

    # fn1 = Ldir['roms_out'] / gtx1 / dstr / his
    # fn2 = Ldir['roms_out'] / gtx2 / dstr / his

    # get model output
    fn_noWWTP = fn_list_noWWTP[i]
    ds_model1 = xr.open_dataset(fn_WWTP)
    ds_model2 = xr.open_dataset(fn_noWWTP)

    fn1 = fn_WWTP
    fn2 = fn_noWWTP

    ds1 = xr.open_dataset(fn_WWTP)
    ds2 = xr.open_dataset(fn_noWWTP)

    G, S, T = zrfun.get_basic_info(fn1)
    zr, zw = zrfun.get_z(G['h'],ds1.zeta[0,:,:].to_numpy(),S)
    dz = np.diff(zw, axis=0)

    vn = 'TN'

    if vn == 'DIN':
        fld1 = ds1['NH4'][0,:,:,:].to_numpy() + ds1['NO3'][0,:,:,:].to_numpy()
        fld2 = ds2['NH4'][0,:,:,:].to_numpy() + ds2['NO3'][0,:,:,:].to_numpy()
    if vn == 'TN':
        vvn_list = ['NO3','NH4','phytoplankton','zooplankton','LdetritusN','SdetritusN']
        ii = 0
        for vvn in vvn_list:
            if ii == 0:
                fld1 = ds1[vvn][0,:,:,:].to_numpy()
                fld2 = ds2[vvn][0,:,:,:].to_numpy()
            else:
                fld1 += ds1[vvn][0,:,:,:].to_numpy()
                fld2 += ds2[vvn][0,:,:,:].to_numpy()
            ii += 1
    else:
        fld1 = ds1[vn][0,:,:,:].to_numpy()
        fld2 = ds2[vn][0,:,:,:].to_numpy()

    # vertical integrals
    c1 = (fld1 * dz).sum(axis=0)
    c2 = (fld2 * dz).sum(axis=0)

    # get dA to get total integral (total N in moles)
    # get grid cell area
    DX = (ds1.pm.values)**-1
    DY = (ds1.pn.values)**-1
    DA = DX*DY # get area in m^2

    # Get total integral
    c1 = c1 * DA
    c2 = c2 * DA

    # get delta between the runs
    dc = c1 - c2 

    #########################################
    ##           Aurora's plots            ##
    #########################################

    plt.close('all')
    fig, ax = plt.subplots(1,1,figsize = (9,9))

    dcmax = np.nanmax(dc)

    maxval = np.nanmax(dc)
    if maxval == 0:
        maxval = 1

    # crop colormap so we can see extremes
    cmap = cmocean.tools.crop_by_percent(cmocean.cm.balance, 15, which='both')
    cmap.set_bad('grey', alpha=1.0)

    # Map plot
    x,y = pfun.get_plon_plat(G['lon_rho'],G['lat_rho'])
    vv = 10#1e-1#1e-3
    cs = ax.pcolormesh(x,y,dc/maxval*100,cmap=cmap,vmin=-vv,vmax=vv)
    # cs = ax.pcolormesh(x,y,dc/np.nansum(dc)*100,cmap=cmap,vmin=-vv,vmax=vv)
    # cs = ax.pcolormesh(x,y,dc,cmap=cmap)#,vmin=-vv,vmax=vv)
    # format figure
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    # ax.set_title(r'Percent of $\Delta$' + ' Normalized Total Integrated %s [mmol]' % (vn))
    ax.set_title(r'$\%$ of max($\Delta$' + ' Integrated %s [mmol])' % (vn))

    # Put colorbar on the left
    cbar = fig.colorbar(cs, ax=ax, location='left', fraction=0.03, pad=0.04)
    cbar.outline.set_visible(False)


    # ==============================================================
    # CUMULATIVE PLOTS
    # ==============================================================

    # Compute cumulative values
    # dc_masked = np.nan_to_num(dc)
    lon = G['lon_rho'][0, :]
    lat = G['lat_rho'][:, 0]
    cum_lon = np.nansum(dc, axis=0)
    cum_lat = np.nansum(dc, axis=1)

    # Create a divider linked to the main axis
    divider = make_axes_locatable(ax)

    # get positive and negative TN integrals
    pos_lon_frac = cum_lon[cum_lon >= 0]
    pos_lon      = lon[cum_lon >= 0]
    neg_lon_frac = cum_lon[cum_lon < 0]
    neg_lon      = lon[cum_lon < 0]

    pos_lat_frac = cum_lat[cum_lat >= 0]
    pos_lat      = lat[cum_lat >= 0]
    neg_lat_frac = cum_lat[cum_lat < 0]
    neg_lat      = lat[cum_lat < 0]

    # Bottom axis, sharing x
    ax_bottom = divider.append_axes("bottom", size="15%", pad=0.1, sharex=ax)
    # ax_bottom.plot(lon, cum_lon/np.nansum(dc)*100, color='orchid', linewidth=2, alpha=0.5)
    ax_bottom.scatter(pos_lon, pos_lon_frac/np.nanmax(pos_lon_frac)*100, s=5, alpha=0.3,color='royalblue',
                label=r'Positive $\%$',zorder=5)
    ax_bottom.scatter(neg_lon, neg_lon_frac/np.nanmax(pos_lon_frac)*100, s=5, alpha=0.3,color='crimson',
                label=r'Negative $\%$',zorder=5)
    ax_bottom.set_ylim(1e-1,1e2)
    ax_bottom.set_yscale('log')
    ax_bottom.set_xlabel('Longitude')
    ax_bottom.set_ylabel(r'$\%$ of max$(\Delta$TN)'+'\nper Lon')
    ax_bottom.grid(True, color='gainsboro')
    plt.setp(ax_bottom.get_xticklabels(), rotation=45)

    # Right axis, sharing y
    ax_right = divider.append_axes("right", size="15%", pad=0.1, sharey=ax)
    # ax_right.plot(cum_lat/np.nansum(dc)*100, lat, color='orchid', linewidth=2, alpha=0.5)
    # ax_right.set_xlim(-0.1,5)
    ax_right.scatter(pos_lat_frac/np.nanmax(pos_lat_frac)*100, pos_lat, s=5, alpha=0.3,color='royalblue',
                label=r'Positive $\%$',zorder=5)
    ax_right.scatter(neg_lat_frac/np.nanmax(pos_lat_frac)*100, neg_lat, s=5, alpha=0.3,color='crimson',
                label=r'Negative $\%$',zorder=5)
    ax_right.set_xlim(1e-1,1e2)
    ax_right.set_xscale('log')
    ax_right.set_ylabel('Latitude')
    ax_right.set_xlabel(r'$\%$ of max$(\Delta$TN)'+'\nper Lat')
    ax_right.grid(True, color='gainsboro')
    ax_right.yaxis.set_label_position('right')
    plt.setp(ax_right.get_yticklabels(), visible=False)

    # Add date to top of figure
    plt.suptitle(str(ds1.ocean_time.values[0].astype('datetime64[D]')))

    # prepare a directory for results
    nouts = ('0000' + str(i))[-4:]
    outname = 'plot_' + nouts + '.png'
    outfile = outdir / outname
    print('Plotting ' + str(fn_WWTP))
    sys.stdout.flush()
    plt.savefig(outfile)
    plt.close()

# make movie
if len(fn_list_WWTP) > 1:
    cmd_list = ['ffmpeg','-r','3','-i', str(outdir)+'/plot_%04d.png', '-vcodec', 'libx264',
        '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir)+'/movie.mp4']
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    if len(stdout) > 0:
        print('\n'+stdout.decode())
    if len(stderr) > 0:
        print('\n'+stderr.decode())