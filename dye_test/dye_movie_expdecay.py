"""

Figures saved in LO_output/AL_custom_plots/dye_release/

"""

###################################################################
##                       import packages                         ##  
###################################################################      

from subprocess import Popen as Po
from subprocess import PIPE as Pi
import numpy as np
import xarray as xr
import cmocean
import matplotlib.pylab as plt

from lo_tools import Lfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
# from pathlib import Path
# pth = Path(__file__).absolute().parent.parent / 'LO' / 'pgrid'
# if str(pth) not in sys.path:
#     sys.path.append(str(pth))
# import gfun_utility as gfu
# import gfun

# Gr = gfun.gstart()
Ldir = Lfun.Lstart()

plt.close('all')

###################################################################
##                          User Inputs                          ##  
################################################################### 

# USER OPTIONS ----------------------------------------------------

d0= '2012.10.07'
d1 = '2012.10.07'

list_type = 'hourly' #'snapshot', 'daily', 'hourly ', 'allhours'
his_num = 26# 11
timestep_interval = 60

# ----------------------------------------------------------------

# gtagex of files to difference
Ldir_dye   = Lfun.Lstart(gridname='cas7', tag='exdye', ex_name='x11exdye')

# get list of history files to plot (and skip ocean_his_0025 from previous day)
fn_list   = Lfun.get_fn_list(list_type, Ldir_dye,
    d0, d1, his_num)[0:his_num:]

fn_list[0] = Ldir['roms_out'] / 'cas7_exdye_x11exdye' / 'f2012.10.07' / 'ocean_his_0001.nc'

vns = ['dye_01','dye_02']
date = '2012.10.07'

outdir0 = Ldir['LOo'] / 'AL_custom_plots' / 'expdecay_dye' / 'dye_01_02_sidebyside'
Lfun.make_dir(outdir0)


# PLOTTING
###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values

###################################################################
##                    Colormap differences                       ##  
################################################################### 

for i,fn in enumerate(fn_list):

    print(i)

    # get model output
    ds_model1 = xr.open_dataset(fn)

    h_dims = ('eta_rho', 'xi_rho')
    lon = ds_model1['lon_rho']
    lat = ds_model1['lat_rho']


    # set bounds
    vmin = 0
    vmax = 0.001

    # Get model1 data
    dye_01_val = ds_model1[vns[0]][0,0,:,:].values
    dye_02_val = ds_model1[vns[1]][0,0,:,:].values

    # Initialize figure
    fig = plt.figure(figsize=(16,8)) # 15,11 for Puget sound and 18,8 for Salish Sea
    plt.tight_layout()

    subplotnums = [121,122]
    stexts = vns
    values = [dye_01_val,dye_02_val]

    newcmap = cmocean.cm.matter

    # loop through all of the plots we need to make
    for j,stext in enumerate(stexts):

        # add water/land
        ax = fig.add_subplot(subplotnums[j])

        # plot values
        cs = ax.pcolormesh(lon, lat,values[j],vmin=vmin, vmax=vmax, cmap=newcmap)
        cbar = fig.colorbar(cs)
        cbar.ax.tick_params(labelsize=14)
        cbar.outline.set_visible(False)
        # format figure
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        ax.axis('off')
        ax.set_ylim([47.4,48])
        ax.set_xlim([-122.9,-122.1])
        # add West Point location
        wwtp_fn = Ldir['data'] / 'trapsD01' / 'processed_data'/ 'wwtp_data_wasielewski_etal_2024.nc'
        wwtp_data_ds = xr.open_dataset(wwtp_fn)
        WP_lat = wwtp_data_ds.lat.sel(source = wwtp_data_ds.name=='King County West Point WWTP').values
        WP_lon = wwtp_data_ds.lon.sel(source = wwtp_data_ds.name=='King County West Point WWTP').values
        ax.scatter(WP_lon,WP_lat,s=80, facecolors='none', edgecolors='deeppink')
        ax.scatter(WP_lon,WP_lat,s=80, facecolors='none', edgecolors='deeppink')
        # pfun.add_coast(ax, color='k')
        pfun.dar(ax)
        ax.set_title(stext + ' bottom concentration [kg/m3]', fontsize=16)
        
        if j == 1:
            ax.text(0.6, 0.1, 'Hour {}'.format(i), color='black', fontweight='bold', fontsize=12,
            transform=ax.transAxes)
        
    # prepare a directory for results
    nouts = ('0000' + str(i))[-4:]
    outname = 'plot_' + nouts + '.png'
    outfile = outdir0 / outname
    print('Plotting ' + str(fn))
    sys.stdout.flush()
    plt.savefig(outfile)
    plt.close()

# make movie
if len(fn_list) > 1:
    cmd_list = ['ffmpeg','-r','3','-i', str(outdir0)+'/plot_%04d.png', '-vcodec', 'libx264',
        '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir0)+'/movie.mp4']
    proc = Po(cmd_list, stdout=Pi, stderr=Pi)
    stdout, stderr = proc.communicate()
    if len(stdout) > 0:
        print('\n'+stdout.decode())
    if len(stderr) > 0:
        print('\n'+stderr.decode())