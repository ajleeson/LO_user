'''
Movie comparing surface phytoplankton between two runs
The first panel shows the basecase surface phytoplankton
The second panel shows the difference between the basecase and the test condition
'''

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import numpy as np
import xarray as xr
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import pinfo
import os, sys

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
import pinfo
from importlib import reload
reload(pfun)
reload(pinfo)

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu

pth = Path(__file__).absolute().parent.parent.parent / 'LO_user' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun

Gr = gfun.gstart()

Ldir = Lfun.Lstart()

import matplotlib as mpl
# mpl.use('QtAgg')
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.close('all')

##############################################################
##                       USER INPUTS                        ##
##############################################################

list_type = 'weekly' #'snapshot', 'daily', 'hourly ', 'allhours', 'weekly'
his_num = 16

year = '2013'

# Provide information about models to compare
# basecase (N-less run)
noN_gtagex = 'cas7_t0noN_x4b'
# long hindcast
hindcast_gtagex = 'cas7_t0_x4b'
gtagexes = [noN_gtagex, hindcast_gtagex]

# Variables for plotting
vn = 'phytoplankton'

#'Puget Sound'
region = 'Puget Sound'

# Show WWTP locations?
WWTP_loc = True

d0 = '2013.01.01'
d1 = '2013.12.31'

##############################################################
##                     Pre-processing                       ##
##############################################################
# gtagex of files to difference
Ldir_hindcast   = Lfun.Lstart(gridname='cas7', tag='t0', ex_name='x4b')
Ldir_noN = Lfun.Lstart(gridname='cas7', tag='t0noN', ex_name='x4b')

HOME = Path.home()
try:
    HOSTNAME = os.environ['HOSTNAME']
except KeyError:
    HOSTNAME = 'BLANK'

# if running script from perigee, get output for hindcast from apogee
if (str(HOME) == '/home/auroral') & ('perigee' in HOSTNAME):
    Ldir_hindcast['roms_out'] = Path('/dat1/parker/LO_roms')

# get list of history files to plot
fn_list_hindcast   = Lfun.get_fn_list(list_type, Ldir_hindcast,
    d0, d1, his_num)
fn_list_noN = Lfun.get_fn_list(list_type, Ldir_noN,
    d0, d1, his_num)

# Read in WWTP locations
wwtp_loc = pd.read_csv('../../LO_data/grids/cas7/wwtp_info.csv')

# lon/lat limits (Puget Sound)
xmin = -123.2
xmax = -122.1
ymin = 46.93
ymax = 48.45

##############################################################
##                     WWTP locations                       ##
##############################################################

# helper function to convert Ecology name to LO name
def SSM2LO_name(rname):
    """
    Given a river name in LiveOcean, find corresponding river name in SSM
    """
    repeatrivs_fn = '../../LO_data/trapsD00/LiveOcean_SSM_rivers.xlsx'
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

##########################################################
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
    ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc')
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

    # Point Sources -------------------------------------------------------------
    # Prepare data for spatial summary plots

    # get flow, nitrate, and ammonium values
    fp_wwtps = '../../LO_output/pre/trapsP00/point_sources/lo_base/Data_historical/'
    flowdf_wwtps = pd.read_pickle(fp_wwtps+'CLIM_flow.p')    # m3/s
    no3df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NO3.p')      # mmol/m3
    nh4df_wwtps = pd.read_pickle(fp_wwtps+'CLIM_NH4.p')      # mmol/m3

    # calculate total DIN concentration in mg/L
    dindf_wwtps = (no3df_wwtps + nh4df_wwtps)/71.4    # mg/L

    # calculate daily loading timeseries in kg/d
    dailyloaddf_wwtps = 86.4*dindf_wwtps*flowdf_wwtps # kg/d = 86.4 * mg/L * m3/s

    # calculate average daily load over the year (kg/d)
    avgload_wwtps = dailyloaddf_wwtps.mean(axis=0).to_frame(name='avg-daily-load(kg/d)')

    # add row and col index for plotting on LiveOcean grid
    griddf0_wwtps = pd.read_csv('../../LO_data/grids/cas6/wwtp_info.csv')
    griddf_wwtps = griddf0_wwtps.set_index('rname') # use point source name as index
    avgload_wwtps = avgload_wwtps.join(griddf_wwtps['row_py']) # add row to avg load df (uses rname to index)
    avgload_wwtps = avgload_wwtps.join(griddf_wwtps['col_py']) # do the same for cols

    # get point source lat and lon
    lon_wwtps = [X[int(col)] for col in avgload_wwtps['col_py']]
    lat_wwtps = [Y[int(row)] for row in avgload_wwtps['row_py']]
    
    # define marker sizes (minimum size is 10 so dots don't get too small)
    sizes_wwtps = [max(0.3*load,30) for load in avgload_wwtps['avg-daily-load(kg/d)']]

##############################################################
##                        Plotting                          ##
##############################################################
outdir0 = Ldir['LOo'] / 'AL_custom_plots' / (hindcast_gtagex + '_MINUS_' + noN_gtagex) / 'surface_videos'
Lfun.make_dir(outdir0)

if len(fn_list_hindcast) > 1:
    # prepare a directory for results if making a movie
    outdir = outdir0 / (region + '_' + list_type+ '_'+vn)
    Lfun.make_dir(outdir, clean=True)

# Loop through history files to create movie
for his_file_ind,fn_hindcast in enumerate(fn_list_hindcast):
    fs = 10
    # plt.tight_layout()
    pfun.start_plot(fs=fs, figsize=(36,27))
    fig = plt.figure()
    gs = fig.add_gridspec(nrows=1, ncols=2, left=0.05, right=0.95, wspace=0.05, hspace=0.05)
    fn_noN = fn_list_noN[his_file_ind]
    ds_hindcast = xr.open_dataset(fn_hindcast)
    ds_noN = xr.open_dataset(fn_noN)

    lons = ds_hindcast.coords['lon_rho'].values
    lats = ds_hindcast.coords['lat_rho'].values
    px, py = pfun.get_plon_plat(lons,lats)

    # set axes range for different state variables
    if vn == 'NO3':
        vmin = 0
        vmax = 40
        cmap = cmocean.cm.matter
    elif vn == 'NH4':
        vmin = 0
        vmax = 6
        cmap = cmocean.cm.matter
    elif vn == 'phytoplankton':
        vmin = 0
        vmax = 30
        cmap = cmocean.cm.algae
    elif vn == 'oxygen':
        vmin = 0
        vmax = 10
        cmap = cmocean.cm.thermal#cmocean.cm.oxy
    elif vn == 'SdetritusN':
        vmin = 0
        vmax = 5
        cmap = cmocean.cm.matter
    elif vn == 'LdetritusN':
        vmin = 0
        vmax = 0.1
        cmap = cmocean.cm.matter

    # scale variable & get units
    scale =  pinfo.fac_dict[vn]
    units = pinfo.units_dict[vn]

    slev = -1

    for i,gtagex in enumerate(gtagexes):
    # Plot basecase map field
        if i == 0:
            ax = fig.add_subplot(1,2,1)
            v = ds_hindcast[vn][0,slev,:,:].values * scale
            cs = ax.pcolormesh(px,py,v, vmin=vmin, vmax=vmax, cmap=cmap)
            # add colorbar
            cbar = fig.colorbar(cs, location='left')
            cbar.ax.tick_params(labelsize=32)#,length=10, width=2)
            cbar.outline.set_visible(False)
            # format figure
            ax.set_xlim([xmin,xmax])
            ax.set_ylim([ymin,ymax])
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            ax.axis('off')
            # pfun.add_coast(ax)
            pfun.dar(ax)
            ax.set_title('(a) N-less run', fontsize=38)
            # save dataset for later use
            bc_ds = ds_hindcast
            bc = v

            # add 10 km bar
            lat0 = 46.94
            lon0 = -123.05
            lat1 = lat0
            lon1 = -122.91825
            distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
            x_dist_km = round(distances_m[0]/1000)
            ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=8)
            ax.text(lon0,lat0+0.01,'{} km'.format(x_dist_km),color='k',fontsize=28)

            # add puget sound map
            inset_map = plt.imread('puget_sound.png')
            imagebox = OffsetImage(inset_map)#, zoom = 0.15)
            ab = AnnotationBbox(imagebox, (-122.3, 47.05), frameon = False)
            ax.add_artist(ab)

        elif i == 1:
            # save test condition
            v = ds_noN[vn][0,slev,:,:].values * scale # we can take time = 0, because there is a snapshot of time (length of 1 in time dimension)
            c1_ds = ds_noN
            c1 = v

    # plot the pcolormesh difference 
    plt.subplots_adjust(wspace=0.01)
    ax = fig.add_subplot(1,2,2)
    diff = (c1 - bc)

    # fixed min and max difference so axis don't change frame to frame
    mindiff = -1
    maxdiff = 6

    # make sure the colorbar is always centered about zero
    cmap = cmocean.tools.crop(cmocean.cm.balance_r, mindiff, maxdiff, 0)
    cs = ax.pcolormesh(px,py,diff, vmin=mindiff, vmax=maxdiff, cmap=cmap)
    cbar = fig.colorbar(cs, location='right')
    cbar.ax.tick_params(labelsize=32)
    cbar.outline.set_visible(False)
    # format everything else
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_title('(b) Anomaly: hindcast minus N-less run', fontsize=38)
    ax.set_xlim([xmin,xmax])
    ax.set_ylim([ymin,ymax])
    ax.axis('off')
    # pfun.add_coast(ax)
    pfun.dar(ax)

    # add date
    sub1 = "/f"
    sub2 = "/ocean"
    # get date between substrings (because file name is very long)
    date_str = str(fn_list_hindcast[his_file_ind])
    idx1 = date_str.find(sub1)
    idx2 = date_str.find(sub2)
    date = date_str[idx1 + len(sub1): idx2]
    ax.text(0.1, 0.01, date,
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes, fontsize=30)

    # add wwtp locations
    if WWTP_loc == True:
        ax.scatter(lon_wwtps,lat_wwtps,color='none', edgecolors='k', linewidth=3, s=sizes_wwtps, label='WWTPs')
        leg_szs = [100, 1000, 10000]
        szs = [0.3*(leg_sz) for leg_sz in leg_szs]
        l0 = plt.scatter([],[], s=szs[0], color='none', edgecolors='k', linewidth=3)
        l1 = plt.scatter([],[], s=szs[1], color='none', edgecolors='k', linewidth=3)
        l2 = plt.scatter([],[], s=szs[2], color='none', edgecolors='k', linewidth=3)
        labels = ['< 100', '1,000', '10,000']
        legend = ax.legend([l0, l1, l2], labels, fontsize = 24, markerfirst=False,
            title='WWTP N loading \n'+r' (kg N d$^{-1}$)',loc='lower right', labelspacing=1, borderpad=0.8)
        plt.setp(legend.get_title(),fontsize=28)

                        
    # Add colormap title
    plt.suptitle( year + ' ' + list_type + ' surface ' + vn + ' ' + units,
                fontsize=44, fontweight='bold', y=0.95)

#-------------------------------------------------------

    if len(fn_list_hindcast) == 1:
        plt.savefig(outdir0 / (list_type + '_'
                + 'withWWTPs_minus_noWWTPs.png'))

    elif len(fn_list_hindcast) > 1:
        nouts = ('0000' + str(his_file_ind))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir / outname
        print('Plotting ' + str(fn_hindcast))
        sys.stdout.flush()
        plt.savefig(outfile)
        plt.close()

# make movie
if len(fn_list_hindcast) > 1:
    cmd_list = ['ffmpeg','-r','8','-i', str(outdir)+'/plot_%04d.png', '-vcodec', 'libx264',
        '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir)+'/movie.mp4']
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    if len(stdout) > 0:
        print('\n'+stdout.decode())
    if len(stderr) > 0:
        print('\n'+stderr.decode())

