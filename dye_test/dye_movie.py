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
Ldir_dye   = Lfun.Lstart(gridname='cas7', tag='dyetest2', ex_name='x11bd')

# get list of history files to plot (and skip ocean_his_0025 from previous day)
fn_list   = Lfun.get_fn_list(list_type, Ldir_dye,
    d0, d1, his_num)[0:his_num:]

fn_list[0] = Ldir['roms_out'] / 'cas7_dyetest2_x11bd' / 'f2012.10.07' / 'ocean_his_0001.nc'

# # add additional files
# if his_num > 25:
#     for i in np.linspace(26,his_num,his_num-26+1):
#         oceanhis = 'ocean_his_00' +str(int(i)) + '.nc'
#         fn_list.append(Ldir['roms_out'] / 'cas7_ats00_debugx11ab' / 'f2012.10.07' / oceanhis)

vn = 'dye_01' # one variable at a time
date = '2012.10.07'

outdir0 = Ldir['LOo'] / 'AL_custom_plots' / 'dye_release_2' / vn
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

    # Get data, and get rid of ocean_time dim (because this is at a single time)
    v1 = ds_model1[vn].squeeze()

    if ('eta_rho' in v1.dims) and ('xi_rho' in v1.dims):
        h_dims = ('eta_rho', 'xi_rho')
        lon = ds_model1['lon_rho']
        lat = ds_model1['lat_rho']
    elif ('eta_u' in v1.dims) and ('xi_u' in v1.dims):
        h_dims = ('eta_u', 'xi_u')
        lon = ds_model1['lon_u']
        lat = ds_model1['lat_u']
    elif ('eta_v' in v1.dims) and ('xi_v' in v1.dims):
        h_dims = ('eta_v', 'xi_v')
        lon = ds_model1['lon_v']
        lat = ds_model1['lat_v']
    elif ('eta_psi' in v1.dims) and ('xi_psi' in v1.dims):
        h_dims = ('eta_psi', 'xi_psi')
        lon = ds_model1['lon_psi']
        lat = ds_model1['lat_psi']
    else:
        raise ValueError(f"Unknown grid type for variable '{vn}'.")

    # set bounds
    vmin = 0
    vmax = 0.0000001

    # Get model1 data
    surf_val = ds_model1[vn][0,-1,:,:].values
    bott_val = ds_model1[vn][0,0,:,:].values

    # Initialize figure
    fig = plt.figure(figsize=(16,8)) # 15,11 for Puget sound and 18,8 for Salish Sea
    plt.tight_layout()

    subplotnums = [121,122]
    stexts = ['Surface','Bottom']
    values = [surf_val,bott_val]

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
        ax.set_title(stext + ' ' + vn + ' concentration [kg/m3]', fontsize=16)
        
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