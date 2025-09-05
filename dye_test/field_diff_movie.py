"""
Plot difference movie between surface/bottom values of specified state variable.
Calculates difference between two different runs
(Written to compare long model1 to N-less run)

From ipython: run field_diff_movie
Figures saved in LO_output/AL_custom_plots/noise_WestPoint/[vn]_diff.mp3

"""

###################################################################
##                       import packages                         ##  
###################################################################      

import numpy as np
import xarray as xr
import matplotlib.pylab as plt
import csv
from scipy.optimize import curve_fit

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys

Ldir = Lfun.Lstart()

plt.close('all')

###################################################################
##                          User Inputs                          ##  
################################################################### 

# Which models to compare??

model1 = 'cas7_exdye2_x11exdye2'            # Dye 2 decays, Same inputs
# model1 = 'cas7_exdye2duplicate_x11exdye2'   # Duplicate of above
# model1 = 'cas7_twindye_x11twindye'          # Both dyes decay, Different inputs
# model1 = 'cas7_twindyeduplicate_x11twindye' # Duplicate of above

model2 = 'cas7_exdye2_x11exdye2'            # Dye 2 decays, Same inputs
# model2 = 'cas7_exdye2duplicate_x11exdye2'   # Duplicate of above
# model2 = 'cas7_twindye_x11twindye'          # Both dyes decay, Different inputs
# model2 = 'cas7_twindyeduplicate_x11twindye' # Duplicate of above

# Which variables to compare??
vn1 = 'dye_01'
vn2 = 'dye_02'

# USER OPTIONS ----------------------------------------------------

d0= '2012.10.07'
d1 = '2012.10.07'

list_type = 'hourly'
his_num = 25
timestep_interval = 60

# ----------------------------------------------------------------

# Some logic for naming things

# split gtagex
grid1, tag1, exec1 = model1.split('_', 2)
grid2, tag2, exec2 = model2.split('_', 2)

# ----------------------------------------------------------------

# gtagex of files to difference
Ldir_model1 = Lfun.Lstart(gridname=grid1, tag=tag1, ex_name=exec1)
Ldir_model2 = Lfun.Lstart(gridname=grid2, tag=tag2, ex_name=exec2)

# get list of history files to plot (and skip ocean_his_0025 from previous day)
fn_list_model1   = Lfun.get_fn_list(list_type, Ldir_model1,
    d0, d1, his_num)[0:his_num:]
fn_list_model2 = Lfun.get_fn_list(list_type, Ldir_model2,
    d0, d1, his_num)[0:his_num:]

fn_list_model1[0] = Ldir['roms_out'] / model1 / ('f' + d0) / 'ocean_his_0001.nc'
fn_list_model2[0] = Ldir['roms_out'] / model2 / ('f' + d0) / 'ocean_his_0001.nc'

# PLOTTING
outdir0 = Ldir['LOo'] / 'AL_custom_plots' / 'dye_noise_tests' / (model1 + '&' + model2)
Lfun.make_dir(outdir0)

if len(fn_list_model1) > 1:
    # prepare a directory for results if making a movie
    outdir = outdir0
    Lfun.make_dir(outdir / 'binary', clean=True)
    Lfun.make_dir(outdir / 'pcolormesh', clean=True)

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

# ###################################################################
# ##                   Binary differences movie                    ##  
# ################################################################### 

# for i,fn_model1 in enumerate(fn_list_model1):

#     # get model output
#     fn_model2 = fn_list_model2[i]
#     ds_model1 = xr.open_dataset(fn_model1)
#     ds_model2 = xr.open_dataset(fn_model2)

#     # Get data, and get rid of ocean_time dim (because this is at a single time)
#     v1 = ds_model1[vn1].squeeze()
#     v2 = ds_model2[vn2].squeeze()

#     # Identify vertical and horizontal dims
#     if 's_rho' in v1.dims:
#         vert_dim = 's_rho'
#     elif 's_w' in v1.dims:
#         vert_dim = 's_w'
#     else:
#         raise ValueError(f"No vertical dimension found for {vn1}")

#     if ('eta_rho' in v1.dims) and ('xi_rho' in v1.dims):
#         h_dims = ('eta_rho', 'xi_rho')
#         lon = ds_model1['lon_rho']
#         lat = ds_model1['lat_rho']
#     elif ('eta_u' in v1.dims) and ('xi_u' in v1.dims):
#         h_dims = ('eta_u', 'xi_u')
#         lon = ds_model1['lon_u']
#         lat = ds_model1['lat_u']
#     elif ('eta_v' in v1.dims) and ('xi_v' in v1.dims):
#         h_dims = ('eta_v', 'xi_v')
#         lon = ds_model1['lon_v']
#         lat = ds_model1['lat_v']
#     elif ('eta_psi' in v1.dims) and ('xi_psi' in v1.dims):
#         h_dims = ('eta_psi', 'xi_psi')
#         lon = ds_model1['lon_psi']
#         lat = ds_model1['lat_psi']
#     else:
#         raise ValueError(f"Unknown grid type for variable '{vn1}'.")

#     # Compute strict difference (True where different, False where equal)
#     diff_mask = (v1 != v2) | (v1.isnull() != v2.isnull())
#     diff_2d = diff_mask.any(dim=vert_dim)

#     # Mask out locations where both are NaN at all depths
#     both_nan = v1.isnull() & v2.isnull()
#     diff_2d = diff_2d.where(~both_nan.any(dim=vert_dim))

#     # Convert to numeric: 0 = same, 1 = diff, NaN = both missing
#     plot_data = diff_2d.astype(float)

#     # Set up colormap: black = 0, lightblue = 1, white = NaN
#     cmap = mcolors.ListedColormap(['paleturquoise', 'black'])
#     bounds = [-0.5, 0.5, 1.5]
#     norm = mcolors.BoundaryNorm(bounds, cmap.N)

#     # Plotting -------------------------------------------------- 

#     # Initialize figure
#     fig, ax = plt.subplots(1,1, figsize=(10, 8))

#     # plot
#     plt.pcolormesh(lon, lat, plot_data, cmap=cmap, norm=norm, shading='auto')

#     # add West Point location
#     wwtp_fn = Ldir['data'] / 'trapsD01' / 'processed_data'/ 'wwtp_data_wasielewski_etal_2024.nc'
#     wwtp_data_ds = xr.open_dataset(wwtp_fn)
#     WP_lat = wwtp_data_ds.lat.sel(source = wwtp_data_ds.name=='King County West Point WWTP').values
#     WP_lon = wwtp_data_ds.lon.sel(source = wwtp_data_ds.name=='King County West Point WWTP').values
#     ax.scatter(WP_lon,WP_lat,s=80, facecolors='none', edgecolors='deeppink')

#     # format figure
#     ax.text(0.04, 0.95, 'Differences', color='black', fontweight='bold', fontsize=12,
#             transform=ax.transAxes)
#     ax.text(0.04, 0.92, 'No Differences', color='turquoise', fontweight='bold', fontsize=12,
#             transform=ax.transAxes)
#     # plt.suptitle('Locations where dye_01 and dye_02 differs between runs at any s-level\nWithin same run',
#     #              fontsize=14,fontweight='bold')
#     plt.suptitle('Locations where {} and {} differs between runs at any s-level\n{} & {}'.format(vn1,vn2,model1,model2),
#                  fontsize=14,fontweight='bold')
#     plt.xlabel('Lon', fontsize=12)
#     plt.ylabel('Lat', fontsize=12)
#     pfun.dar(ax)

#     ax.text(0.7, 0.95, 'Hour {}'.format(i), color='black', fontweight='bold', fontsize=12,
#             transform=ax.transAxes)

#     # Salish Sea
#     ax.set_ylim([46.5,50])
#     ax.set_xlim([-126.5,-122])

#     # # West Point
#     # ax.set_ylim([47.4,48])
#     # ax.set_xlim([-122.9,-122.1])

#     plt.tight_layout()

#     # prepare a directory for results
#     nouts = ('0000' + str(i))[-4:]
#     outname = 'plot_' + nouts + '.png'
#     outfile = outdir /'binary' / outname
#     print('Plotting ' + str(fn_model1))
#     sys.stdout.flush()
#     plt.savefig(outfile)
#     plt.close()

# # make movie
# if len(fn_list_model1) > 1:
#     cmd_list = ['ffmpeg','-r','3','-i', str(outdir / 'binary')+'/plot_%04d.png', '-vcodec', 'libx264',
#         '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir / 'binary')+'/movie.mp4']
#     proc = Po(cmd_list, stdout=Pi, stderr=Pi)
#     stdout, stderr = proc.communicate()
#     if len(stdout) > 0:
#         print('\n'+stdout.decode())
#     if len(stderr) > 0:
#         print('\n'+stderr.decode())

# ###################################################################
# ##                    Colormap differences                       ##  
# ################################################################### 

# for i,fn_model1 in enumerate(fn_list_model1):

#     # get model output
#     fn_model2 = fn_list_model2[i]
#     ds_model1 = xr.open_dataset(fn_model1)
#     ds_model2 = xr.open_dataset(fn_model2)

#     # Get data, and get rid of ocean_time dim (because this is at a single time)
#     v1 = ds_model1[vn1].squeeze()
#     v2 = ds_model2[vn2].squeeze()

#     if ('eta_rho' in v1.dims) and ('xi_rho' in v1.dims):
#         h_dims = ('eta_rho', 'xi_rho')
#         lon = ds_model1['lon_rho']
#         lat = ds_model1['lat_rho']
#     elif ('eta_u' in v1.dims) and ('xi_u' in v1.dims):
#         h_dims = ('eta_u', 'xi_u')
#         lon = ds_model1['lon_u']
#         lat = ds_model1['lat_u']
#     elif ('eta_v' in v1.dims) and ('xi_v' in v1.dims):
#         h_dims = ('eta_v', 'xi_v')
#         lon = ds_model1['lon_v']
#         lat = ds_model1['lat_v']
#     elif ('eta_psi' in v1.dims) and ('xi_psi' in v1.dims):
#         h_dims = ('eta_psi', 'xi_psi')
#         lon = ds_model1['lon_psi']
#         lat = ds_model1['lat_psi']
#     else:
#         raise ValueError(f"Unknown grid type for variable '{vn1}'.")

#     # set bounds
#     vmin = -0.00001
#     vmax =  0.00001

#     # Get model1 data
#     surf_vn_model1 = ds_model1[vn1][0,-1,:,:].values
#     bott_vn_model1 = ds_model1[vn1][0,0,:,:].values
#     # Get model2 data
#     surf_vn_model2 = ds_model2[vn2][0,-1,:,:].values
#     bott_vn_model2 = ds_model2[vn2][0,0,:,:].values
#     # Get difference
#     surf_diff = (surf_vn_model1 - surf_vn_model2)
#     bott_diff = (bott_vn_model1 - bott_vn_model2)

#     # Initialize figure
#     fig = plt.figure(figsize=(16,8)) # 15,11 for Puget sound and 18,8 for Salish Sea
#     plt.tight_layout()

#     subplotnums = [121,122]
#     stexts = ['Surface','Bottom']
#     values = [surf_diff,bott_diff]

#     newcmap = cmocean.tools.crop_by_percent(cmocean.cm.balance_r, 20, which='both', N=None)

#     # loop through all of the plots we need to make
#     for j,stext in enumerate(stexts):

#         # add water/land
#         ax = fig.add_subplot(subplotnums[j])

#         # plot values
#         cs = ax.pcolormesh(lon, lat,values[j],vmin=vmin, vmax=vmax, cmap=newcmap)
#         cbar = fig.colorbar(cs)
#         cbar.ax.tick_params(labelsize=14)
#         cbar.outline.set_visible(False)
#         # format figure
#         # ax.set_xlim([xmin,xmax])
#         # ax.set_ylim([ymin,ymax])
#         ax.set_yticklabels([])
#         ax.set_xticklabels([])
#         ax.axis('off')
#         ax.set_ylim([47.4,48])
#         ax.set_xlim([-122.9,-122.1])
#         ax.scatter(WP_lon,WP_lat,s=80, facecolors='none', edgecolors='deeppink')
#         # pfun.add_coast(ax, color='k')
#         pfun.dar(ax)
#         ax.set_title('Dye difference at ' + stext + '[kg/m3]', fontsize=16)
#         fig.suptitle('{} minus {}\nBetween {} & {}'.format(vn1,vn2,model1,model2) + d0,
#                     fontsize=11, fontweight='bold')
        
#         if j == 1:
#             ax.text(0.6, 0.1, 'Hour {}'.format(i), color='black', fontweight='bold', fontsize=12,
#             transform=ax.transAxes)
        
#     # prepare a directory for results
#     nouts = ('0000' + str(i))[-4:]
#     outname = 'plot_' + nouts + '.png'
#     outfile = outdir /'pcolormesh' / outname
#     print('Plotting ' + str(fn_model1))
#     sys.stdout.flush()
#     plt.savefig(outfile)
#     plt.close()

# # make movie
# if len(fn_list_model1) > 1:
#     cmd_list = ['ffmpeg','-r','3','-i', str(outdir / 'pcolormesh')+'/plot_%04d.png', '-vcodec', 'libx264',
#         '-pix_fmt', 'yuv420p', '-crf', '25', str(outdir / 'pcolormesh')+'/movie.mp4']
#     proc = Po(cmd_list, stdout=Pi, stderr=Pi)
#     stdout, stderr = proc.communicate()
#     if len(stdout) > 0:
#         print('\n'+stdout.decode())
#     if len(stderr) > 0:
#         print('\n'+stderr.decode())

###################################################################
##                         Residual plot                         ##  
################################################################### 

wwtp_fn = Ldir['data'] / 'trapsD01' / 'processed_data'/ 'wwtp_data_wasielewski_etal_2024.nc'
wwtp_data_ds = xr.open_dataset(wwtp_fn)
WP_lat = wwtp_data_ds.lat.sel(source = wwtp_data_ds.name=='King County West Point WWTP').values
WP_lon = wwtp_data_ds.lon.sel(source = wwtp_data_ds.name=='King County West Point WWTP').values

# West Point position
lat0, lon0 = WP_lat, WP_lon

ds_model1_lasthour = xr.open_dataset(fn_list_model1[-1])
ds_model2_lasthour = xr.open_dataset(fn_list_model2[-1])

# Compute residual and squeeze ocean_time
res = ds_model1_lasthour[vn1].squeeze() - ds_model2_lasthour[vn2].squeeze()  # shape (s_rho, eta_rho, xi_rho)

# Load lat/lon on rho grid (shape (eta_rho, xi_rho))
lat = ds_model1_lasthour['lat_rho'].values
lon = ds_model1_lasthour['lon_rho'].values

# get residual data: (s_rho, eta_rho, xi_rho)
res_data = res.values  # (30, 1302, 663)

# Broadcast lat/lon to match res vertical dimension for flattening
# lat and lon are 2D (eta_rho, xi_rho)
# We want 3D (s_rho, eta_rho, xi_rho) by repeating lat/lon vertically
lat3d = np.broadcast_to(lat, res_data.shape)
lon3d = np.broadcast_to(lon, res_data.shape)

# Flatten all arrays to 1D for scatter plot
flat_res = res_data.flatten()
flat_lat = lat3d.flatten()
flat_lon = lon3d.flatten()

# Mask to keep only finite, nonzero residuals
valid = np.isfinite(flat_res) & (flat_res != 0)
flat_res = flat_res[valid]
flat_lat = flat_lat[valid]
flat_lon = flat_lon[valid]

# Get distance from West Point
x, y = zfun.ll2xy(flat_lon, flat_lat, lon0, lat0)
distance = np.sqrt(x**2 + y**2) / 1000  # convert to km

residuals = flat_res
# get positive and negative residual values, and take absolute value
pos_residuals = flat_res[flat_res > 0]
pos_distances = distance[flat_res > 0]
neg_residuals = flat_res[flat_res < 0] * -1
neg_distances = distance[flat_res < 0]

# # get mean residual value
# mean_res = np.nanmean(flat_res) * scale

# Plot
fig, axes = plt.subplots(2,1, figsize=(8, 8),sharex=True)
ax = axes.ravel()

# Residual
ax[0].scatter(distance, residuals, s=5, alpha=0.3,color='black',zorder=5)
# format figure
ax[0].set_ylabel('dye residual [kg/m3]',fontsize=12)
ax[0].grid(True,color='gainsboro')
# ax[0].set_ylim([-0.2,3.6])


# Absolute value and log transform
ax[1].scatter(pos_distances, pos_residuals, s=5, alpha=0.3,color='royalblue',
            label='Positive Residual',zorder=5)
ax[1].scatter(neg_distances, neg_residuals, s=5, alpha=0.3,color='crimson',
            label='Negative Residual',zorder=5)
# format figure
ax[1].set_yscale('log')
ax[1].set_xlim([0,25])
# ax[1].set_ylim([10e-12,10e1])
ax[1].set_xlabel('Distance from West Point [km]',fontsize=12)
ax[1].set_ylabel('dye residual [kg/m3]',fontsize=12)
ax[1].grid(True,color='gainsboro')
ax[1].legend(loc='upper right',fontsize=12)


plt.suptitle('{} - {} vs. distance from West Point\n({} - {})'.format(vn1,vn2,model1,model2),fontsize=12)
plt.tight_layout()
plt.savefig(outdir0/'residuals')
plt.close()

###################################################################
##                 Verify exponential decay                      ##  
###################################################################

# initialize empty list of dye ratios
dye_mass_1 = []
dye_mass_2 = []
seconds = np.linspace(0,86400,25)

# get thickness of dye in watercolumn at every lat/lon cell
# units are in m (thickness of hypoxic layer)
# get S for the whole grid
Sfp = Ldir['data'] / 'grids' / 'cas7' / 'S_COORDINATE_INFO.csv'
reader = csv.DictReader(open(Sfp))
S_dict = {}
for row in reader:
    S_dict[row['ITEMS']] = row['VALUES']
S = zrfun.get_S(S_dict)
# get cell thickness
h = grid_ds['h'].values # height of water column
z_rho, z_w = zrfun.get_z(h, 0*h, S) 
dzr = np.diff(z_w, axis=0) # vertical thickness of all cells [m] 

# get horizontal gridcell area
DX = (grid_ds.pm.values)**-1
DY = (grid_ds.pn.values)**-1
DA = DX*DY # get area in m2

# loop through days
for i,fn_model1 in enumerate(fn_list_model1):

    # get model output
    fn_model2 = fn_list_model2[i]
    ds_model1 = xr.open_dataset(fn_model1)
    ds_model2 = xr.open_dataset(fn_model2) 

    # Vertical dye thickness
    # Now get dye values at every grid cell
    dye_kgm3_model1 = ds_model1[vn1].values
    dye_kgm3_model2 = ds_model2[vn2].values
    # Multiple cell height array by dye concentration 
    dye_scell_thickness_1 = dzr * dye_kgm3_model1 # kg/m2
    dye_scell_thickness_2 = dzr * dye_kgm3_model2 # kg/m2
    # Sum along z to get thickness of dye layer
    dye_thick_1 = np.nansum(dye_scell_thickness_1,axis=1)
    dye_thick_2 = np.nansum(dye_scell_thickness_2,axis=1)

    # total dye mass
    dye_mass_temp_1 = np.sum(dye_thick_1 * DA, axis=(1, 2)) # kg
    dye_mass_temp_2 = np.sum(dye_thick_2 * DA, axis=(1, 2)) # kg
    
    # add to lists
    dye_mass_1.extend(dye_mass_temp_1)
    dye_mass_2.extend(dye_mass_temp_2)

# Plot actual mass of dye
fig, ax = plt.subplots(1,1, figsize=(8, 4))
ax.plot(seconds, dye_mass_1, 'o', markersize=5, linestyle='None',
        alpha=0.5,color='hotpink',label='{} total mass'.format(vn1))
ax.plot(seconds, dye_mass_2, 'o', markersize=5, linestyle='None',
        alpha=0.5,color='royalblue',label='{} total mass'.format(vn2))

# Plot expected mass of dye
Q = 4.07676581702 # m3/s from West Point
C = 10 # kg/m3
a = 1e-5 # 1/s
time = np.linspace(0,86500,1000)
expected_dye01 = Q*C*time
expected_dye02 = Q*C/a * (1-np.exp(-1*a*time))

ax.plot(seconds, expected_dye01, color='crimson',label='Expected {}'.format(vn1))
ax.plot(seconds, expected_dye02, color='navy',label='Expected {}'.format(vn2))


# # expected decay
# time = np.linspace(0,86500,1000)
# expected = np.exp(-1*1e-5*time)
# ax.plot(time,expected,color='black',linewidth=3,
#         label='Expected exponential decay')

# # Fit decay
# time = np.linspace(0,86500,1000)
# # Convert to arrays and automatically remove NaNs/Infs
# mask = np.isfinite(all_sec) & np.isfinite(dye_ratios)
# all_sec_arr = np.array(all_sec)[mask]
# dye_ratios_arr = np.array(dye_ratios)[mask]
# # Fit exponential decay with A=1
# def exp_decay_fixedA(t, k):
#     return np.exp(-k * t)
# k_fit, _ = curve_fit(exp_decay_fixedA, all_sec_arr, dye_ratios_arr, p0=[1e-5])
# # Fitted curve
# fit = np.exp(-k_fit * time)
# ax.plot(time,fit,color='royalblue',linewidth=3,
#         linestyle='--',label='Fitted exponential decay')

ax.set_ylim(0,3.5e6)
ax.set_xlim(0,86400)
ax.set_xlabel('Seconds',fontsize=12)
ax.set_ylabel('Dye Mass [kg]',fontsize=12)
ax.grid(True,color='gainsboro')
ax.legend(loc='best')

ax.set_title('Total mass of dye')

plt.savefig(outdir0/'exp_decay_check')
plt.close()



