"""
Plot difference between surface/bottom values of specified state variable.
Calculates difference between two different runs
(Written to compare long hindcast to N-less run)

From ipython: run model_field_diff
Figures saved in LO_output/AL_custom_plots/[vn]_diff.png

"""

###################################################################
##                       import packages                         ##  
###################################################################      

import numpy as np
import xarray as xr
import matplotlib.patches as patches
import matplotlib.pylab as plt
import cmocean
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

Gr = gfun.gstart()
Ldir = Lfun.Lstart()

###################################################################
##                          User Inputs                          ##  
################################################################### 

vn = 'oxygen' # u, v, w, oxygen, salt, temp
# date = '2013.01.01'
date = '2013.01.14'
# date = '2013.05.01'
# date = '2014.01.09'
# date = '2013.03.01'
# date = '2013.09.01'

###################################################################
##          load output folder, grid data, model output          ##  
################################################################### 

gtagex_ant = 'cas7_t0_x4b'
gtagex_noN = 'cas7_t0noN_x4b'

# where to put output figures
out_dir = Ldir['LOo'] / 'AL_custom_plots'
Lfun.make_dir(out_dir)

# Get grid data
G = zrfun.get_basic_info(Ldir['data'] / 'grids/cas7/grid.nc', only_G=True)
grid_ds = xr.open_dataset(Ldir['data'] / 'grids/cas7/grid.nc')
lon = grid_ds.lon_rho.values
lat = grid_ds.lat_rho.values
lon_u = grid_ds.lon_u.values
lat_u = grid_ds.lat_u.values
lon_v = grid_ds.lon_v.values
lat_v = grid_ds.lat_v.values
X = lon[0,:] # grid cell X values
Y = lat[:,0] # grid cell Y values

# layer, text, axis limits
slev_surf = -1
slev_bott = 0

# Puget Sound boundaries
PS_imin = 460
PS_imax = 642
PS_jmin = 600
PS_jmax = 900

# # Salish Sea boundaries
SS_imin = 300
SS_imax = 652
SS_jmin = 590
SS_jmax = 1150

# west boundary
west_imin = 0
west_imax = 21
west_jmin = 0
west_jmax = -1

# get model output
fp_ant = Ldir['roms_out'] / gtagex_ant / ('f'+date) / 'ocean_his_0025.nc'
fp_noN = Ldir['roms_out'] / gtagex_noN / ('f'+date) / 'ocean_his_0025.nc'
ds_ant = xr.open_dataset(fp_ant)
ds_noN = xr.open_dataset(fp_noN)

###################################################################
##                     Calculate differences                     ##  
################################################################### 

# get coordinates for pcolormesh
if vn == 'u':
    px, py = pfun.get_plon_plat(lon_u,lat_u)
elif vn == 'v':
    px, py = pfun.get_plon_plat(lon_v,lat_v)
else:
    px, py = pfun.get_plon_plat(lon,lat)

# sum nitrate and ammonium for DIN
if vn == 'DIN':
    vn_name = 'NO3'
    vn_no3 = 'NO3'
    vn_nh4 = 'NH4'
    # vmin = -5
    # vmax =  5
    vmin = -0.001
    vmax =  0.001
elif vn == 'u' or vn == 'v':
    vn_name = vn
    vmin = -0.00001#-0.01
    vmax =  0.00001#0.01
elif vn == 'oxygen':
    vn_name = vn
    vmin = -0.001
    vmax =  0.001
elif vn == 'salt':
    vn_name = vn
    vmin = -0.00001
    vmax =  0.00001
elif vn == 'temp':
    vn_name = vn
    vmin = -0.00001
    vmax =  0.00001
else:
    print('vmin and vmax not provided for '+ vn)

# scale variable
scale =  pinfo.fac_dict[vn_name]

# Get anthropogenic data
surf_vn_ant_fulldomain = ds_ant[vn][0,slev_surf,:,:].values
bott_vn_ant_fulldomain = ds_ant[vn][0,slev_bott,:,:].values
surf_vn_ant_PS = ds_ant[vn][0,slev_surf,PS_jmin:PS_jmax,PS_imin:PS_imax].values
bott_vn_ant_PS = ds_ant[vn][0,slev_bott,PS_jmin:PS_jmax,PS_imin:PS_imax].values
surf_vn_ant_SS = ds_ant[vn][0,slev_surf,SS_jmin:SS_jmax,SS_imin:SS_imax].values
bott_vn_ant_SS = ds_ant[vn][0,slev_bott,SS_jmin:SS_jmax,SS_imin:SS_imax].values
surf_vn_ant_west = ds_ant[vn][0,slev_surf,west_jmin:west_jmax,west_imin:west_imax].values
bott_vn_ant_west = ds_ant[vn][0,slev_bott,west_jmin:west_jmax,west_imin:west_imax].values

# Get noN data
surf_vn_noN_fulldomain = ds_noN[vn][0,slev_surf,:,:].values
bott_vn_noN_fulldomain = ds_noN[vn][0,slev_bott,:,:].values
surf_vn_noN_PS = ds_noN[vn][0,slev_surf,PS_jmin:PS_jmax,PS_imin:PS_imax].values
bott_vn_noN_PS = ds_noN[vn][0,slev_bott,PS_jmin:PS_jmax,PS_imin:PS_imax].values
surf_vn_noN_SS = ds_noN[vn][0,slev_surf,SS_jmin:SS_jmax,SS_imin:SS_imax].values
bott_vn_noN_SS = ds_noN[vn][0,slev_bott,SS_jmin:SS_jmax,SS_imin:SS_imax].values
surf_vn_noN_west = ds_noN[vn][0,slev_surf,west_jmin:west_jmax,west_imin:west_imax].values
bott_vn_noN_west = ds_noN[vn][0,slev_bott,west_jmin:west_jmax,west_imin:west_imax].values

# Get difference
surf_diff_fulldomain = (surf_vn_noN_fulldomain - surf_vn_ant_fulldomain) * scale
bott_diff_fulldomain = (bott_vn_noN_fulldomain - bott_vn_ant_fulldomain) * scale
surf_diff_PS = (surf_vn_noN_PS - surf_vn_ant_PS) * scale
bott_diff_PS = (bott_vn_noN_PS - bott_vn_ant_PS) * scale
surf_diff_SS = (surf_vn_noN_SS - surf_vn_ant_SS) * scale
bott_diff_SS = (bott_vn_noN_SS - bott_vn_ant_SS) * scale
surf_diff_west = (surf_vn_noN_west - surf_vn_ant_west) * scale
bott_diff_west = (bott_vn_noN_west - bott_vn_ant_west) * scale

###################################################################
##                  Plotting and saving figure                   ##  
################################################################### 

# define colors
PS_color = 'red'
SS_color = 'purple'
west_color = 'royalblue'

plt.close('all')

# Initialize figure
fig = plt.figure(figsize=(12, 7))
ax_map = plt.subplot(1,2,1)
ax_surf = plt.subplot(2,2,2)
ax_bott = plt.subplot(2,2,4)
axes = [ax_map, ax_surf, ax_bott]

# plot map
# pfun.add_coast(ax_map)
# difference
ax_map.set_title('(a) Surface difference',loc='left')
newcmap = cmocean.cm.balance_r
cs = ax_map.pcolormesh(px,py,surf_diff_fulldomain,vmin=vmin, vmax=vmax, cmap=newcmap)
cbar = fig.colorbar(cs, location='left')
# cbar.ax_map.tick_params(labelsize=14)
cbar.outline.set_visible(False)
# Puget Sound
PS_rect = patches.Rectangle((X[PS_imin], Y[PS_jmin]), X[PS_imax]-X[PS_imin], Y[PS_jmax]-Y[PS_jmin],
                          edgecolor=PS_color, facecolor='none', linewidth=2)
ax_map.add_patch(PS_rect)
# Salish Sea
SS_rect = patches.Rectangle((X[SS_imin], Y[SS_jmin]), X[SS_imax]-X[SS_imin], Y[SS_jmax]-Y[SS_jmin],
                          edgecolor=SS_color, facecolor='none', linewidth=2)
ax_map.add_patch(SS_rect)
# Western boundary
west_rect = patches.Rectangle((X[west_imin], Y[west_jmin]), X[west_imax]-X[west_imin], Y[west_jmax]-Y[west_jmin],
                          edgecolor=west_color, facecolor='none', linewidth=2)
ax_map.add_patch(west_rect)
pfun.dar(ax_map)
ax_map.set_xlim([X[0],X[-1]])
ax_map.set_ylim([Y[0],Y[-1]])

# get values to put in histogram (flatten to 1D array)
surf_diff_fulldomain_1d = surf_diff_fulldomain.flatten().tolist()
bott_diff_fulldomain_1d = bott_diff_fulldomain.flatten().tolist()
surf_diff_PS_1d = surf_diff_PS.flatten().tolist()
bott_diff_PS_1d = bott_diff_PS.flatten().tolist()
surf_diff_SS_1d = surf_diff_SS.flatten().tolist()
bott_diff_SS_1d = bott_diff_SS.flatten().tolist()
surf_diff_west_1d = surf_diff_west.flatten().tolist()
bott_diff_west_1d = bott_diff_west.flatten().tolist()
# remove nans and zero
surf_diff_fulldomain_1d = [x for x in surf_diff_fulldomain_1d if str(x) != 'nan' and x != 0]
bott_diff_fulldomain_1d = [x for x in bott_diff_fulldomain_1d if str(x) != 'nan' and x != 0]
surf_diff_PS_1d = [x for x in surf_diff_PS_1d if str(x) != 'nan' and x != 0]
bott_diff_PS_1d = [x for x in bott_diff_PS_1d if str(x) != 'nan' and x != 0]
surf_diff_SS_1d = [x for x in surf_diff_SS_1d if str(x) != 'nan' and x != 0]
bott_diff_SS_1d = [x for x in bott_diff_SS_1d if str(x) != 'nan' and x != 0]
surf_diff_west_1d = [x for x in surf_diff_west_1d if str(x) != 'nan' and x != 0]
bott_diff_west_1d = [x for x in bott_diff_west_1d if str(x) != 'nan' and x != 0]
# add to list
surf_diff_all = [surf_diff_fulldomain_1d,
                 surf_diff_PS_1d,
                 surf_diff_SS_1d,
                 surf_diff_west_1d]
bott_diff_all = [bott_diff_fulldomain_1d,
                 bott_diff_PS_1d,
                 bott_diff_SS_1d,
                 bott_diff_west_1d]
diff_names = ['Full domain',
            'Puget Sound',
            'Salish Sea',
            'Western boundary']

# print values
print('Full domain ========================')
print('    min:{}'.format(np.nanmin(surf_diff_fulldomain_1d)))
print('    max:{}'.format(np.nanmax(surf_diff_fulldomain_1d)))
print('    median:{}'.format(np.nanmedian(surf_diff_fulldomain_1d)))
print('    mean:{}'.format(np.nanmean(surf_diff_fulldomain_1d)))
print('Puget Sound ========================')
print('    min:{}'.format(np.nanmin(surf_diff_PS_1d)))
print('    max:{}'.format(np.nanmax(surf_diff_PS_1d)))
print('    median:{}'.format(np.nanmedian(surf_diff_PS_1d)))
print('    mean:{}'.format(np.nanmean(surf_diff_PS_1d)))
print('Salish Sea ========================')
print('    min:{}'.format(np.nanmin(surf_diff_SS_1d)))
print('    max:{}'.format(np.nanmax(surf_diff_SS_1d)))
print('    median:{}'.format(np.nanmedian(surf_diff_SS_1d)))
print('    mean:{}'.format(np.nanmean(surf_diff_SS_1d)))
print('Western boundary ========================')
print('    min:{}'.format(np.nanmin(surf_diff_west_1d)))
print('    max:{}'.format(np.nanmax(surf_diff_west_1d)))
print('    median:{}'.format(np.nanmedian(surf_diff_west_1d)))
print('    mean:{}'.format(np.nanmean(surf_diff_west_1d)))


bins = [x*0.02 - 0.14 for x in np.linspace(0,8,10)]

# plot surface difference histogram
ax_surf.set_title('(b) Surface difference',loc='left')
# ax_surf.hist(surf_diff_fulldomain_1d,bins=bins,alpha=0.5, color='royalblue')
ax_surf.axhline(y=0, xmin=-0.5, xmax=1.05,color='silver',linewidth=1,linestyle='--')
bplot = ax_surf.boxplot(surf_diff_all, patch_artist=True, labels=diff_names,
                showmeans=True, showfliers=False, boxprops={'color':'darkgray'}, meanprops=
                {'marker': 'o', 'markerfacecolor': 'navy', 'markersize': 7, 'markeredgecolor': 'none'},
                flierprops=
                {'marker': '.', 'markerfacecolor': 'black', 'markersize': 1, 'markeredgecolor': 'none'})
# format boxplot
for patch in bplot['boxes']:
    patch.set_facecolor('whitesmoke')
for element in ['whiskers', 'medians', 'caps']:
        plt.setp(bplot[element], color='darkgray')
# ax_surf.hist(surf_diff_PS_1d,bins=bins,alpha=0.5, color=PS_color)
# ax_surf.hist(surf_diff_west_1d,bins=bins,alpha=0.5, color=west_color)
# ax_surf.set_yscale('log')

# # plot bottom difference histogram
ax_bott.set_title('(c) Bottom difference', loc='left')
# # ax_bott.hist(bott_diff_fulldomain_1d,bins=bins,alpha=0.5, color='royalblue')
ax_bott.axhline(y=0, xmin=-0.5, xmax=1.05,color='silver',linewidth=1,linestyle='--')
bplot = ax_bott.boxplot(bott_diff_all, patch_artist=True, labels=diff_names,
                showmeans=True, showfliers=False, boxprops={'color':'darkgray'}, meanprops=
                {'marker': 'o', 'markerfacecolor': 'navy', 'markersize': 7, 'markeredgecolor': 'none'},
                flierprops=
                {'marker': '.', 'markerfacecolor': 'black', 'markersize': 1, 'markeredgecolor': 'none'})
# format boxplot
for patch in bplot['boxes']:
    patch.set_facecolor('whitesmoke')
for element in ['whiskers', 'medians', 'caps']:
        plt.setp(bplot[element], color='darkgray')
# ax_bott.hist(bott_diff_PS_1d,bins=bins,alpha=0.5, color=PS_color)
# ax_bott.hist(bott_diff_west_1d,bins=bins,alpha=0.5, color=west_color)
# ax_bott.set_yscale('log')
# ax_bott.set_xlabel('Natural - Anthropogenic [mg/L]')

plt.tight_layout
plt.suptitle(date + '\nNatural - Anthropogenic DO [mg/L]')

# subplotnums = [121,122]
# stexts = [stext_surf,stext_bott]
# values = [surf_diff_fulldomain,bott_diff_fulldomain]

# newcmap = cmocean.cm.balance_r

# # loop through all of the plots we need to make
# for i,stext in enumerate(stexts):

#     # add water/land
#     ax = fig.add_subplot(subplotnums[i])

#     # plot values
#     cs = ax.pcolormesh(px,py,values[i],vmin=vmin, vmax=vmax, cmap=newcmap)
#     cbar = fig.colorbar(cs)
#     cbar.ax.tick_params(labelsize=14)
#     cbar.outline.set_visible(False)
#     # format figure
#     # ax.set_xlim([xmin,xmax])
#     # ax.set_ylim([ymin,ymax])
#     ax.set_yticklabels([])
#     ax.set_xticklabels([])
#     ax.axis('off')
#     # pfun.add_coast(ax, color='k')
#     pfun.dar(ax)
#     ax.set_title(vn + ' difference at ' + stext + pinfo.units_dict[vn_name], fontsize=16)
#     fig.suptitle('Anthropogenic - Natural\n' + date + ' ocean_his_0025',
#                 fontsize=18, fontweight='bold')

#     # add 10 km bar
#     lat0 = 47
#     lon0 = -122.4
#     lat1 = lat0
#     lon1 = -122.27
#     distances_m = zfun.ll2xy(lon1,lat1,lon0,lat0)
#     x_dist_km = round(distances_m[0]/1000)
#     # ax.plot([lon0,lon1],[lat0,lat1],color='k',linewidth=6)
#     # ax.text((lon0+lon1)/2,lat0+0.02,'{} km'.format(x_dist_km),color='k',
#     #         horizontalalignment='center', fontsize=15)
    
#     # # add WWTP locations
#     # ax.scatter(wwtp_lon,wwtp_lat,s=30,alpha=0.5,
#     #         facecolors='none',edgecolors='k')

# # Generate plot
# plt.tight_layout
# plt.subplots_adjust(wspace=0.05)
# plt.savefig(out_dir / (date+'_'+vn+'_diff.png'))