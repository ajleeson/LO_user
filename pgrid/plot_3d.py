"""
Plot grid in 3D to look at
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='',
        type=str)
args = parser.parse_args()

import cmocean
from lo_tools import zfun,Lfun

import gfun
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
from lo_tools import plotting_functions as pfun
import gfun_plotting as gfp

import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pickle

testing = True
if testing:
    from importlib import reload
    reload(gfun)
    reload(gfp)

Ldir = Lfun.Lstart()
out_dir = Ldir['LOo'] / 'analytical_grids'
Lfun.make_dir(out_dir)

# select grid file
in_fn = gfun.select_file(Gr)

# load the default choices
dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))

# get river info if it exists
do_riv = False
ri_fn = Gr['gdir'] / 'roms_river_info.csv'
if ri_fn.is_file():
    rri_df = pd.read_csv(ri_fn, index_col='rname')
    do_riv = True

# load the data
ds = xr.open_dataset(in_fn)
z = -ds.h.values
mask_rho = ds.mask_rho.values

lon = ds.lon_rho.values
lat = ds.lat_rho.values

plon, plat = pfun.get_plon_plat(lon,lat)
pad = 0.05*(plat[-1,0]-plat[0,0])
ax_lims = (plon[0,0]-pad, plon[0,-1]+pad, plat[0,0]-pad, plat[-1,0]+pad)

# make a version of z with nans where masked
zm = z.copy()
zm[mask_rho == 0] = 0#np.nan

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(24,24))

# zoom in
# latmin = round(0.5*np.shape(plat)[1]) - round(0.4*np.shape(plat)[1])
# latmax = round(0.5*np.shape(plat)[1]) + round(0.4*np.shape(plat)[1])
# lonmin = round(0.5*np.shape(plon)[0]) - round(0.5*np.shape(plon)[0])
# lonmax = round(0.5*np.shape(plon)[0]) + round(0.4*np.shape(plon)[0])
latmin = round(0.5*np.shape(plat)[1]) - round(0.5*np.shape(plat)[1])
latmax = round(0.5*np.shape(plat)[1]) + round(0.5*np.shape(plat)[1])
lonmin = round(0.5*np.shape(plon)[0]) - round(0.5*np.shape(plon)[0])
lonmax = round(0.5*np.shape(plon)[0]) + round(0.5*np.shape(plon)[0])

x, y = zfun.ll2xy(lon, lat, -122.6, 48)

# bathymetry
fig = plt.figure(figsize=(8,15))
ax = fig.add_subplot(2,1,1, projection='3d')
cs = ax.plot_surface(lon, lat, zm, rstride=1, cstride=1,
                    edgecolor='none', alpha = 1, cmap=cmocean.cm.deep_r)
cset = ax.contourf(lon, lat, zm, zdir='z', offset=np.min(zm), cmap=cmocean.cm.deep_r)
# cs = ax.plot_surface(lon[lonmin:lonmax,latmin:latmax], lat[lonmin:lonmax,latmin:latmax], zm[lonmin:lonmax,latmin:latmax],
#  rstride=1, cstride=1, edgecolor='none', alpha = 1, cmap=cmocean.cm.deep_r)
# cset = ax.contourf(lon[lonmin:lonmax,latmin:latmax], lat[lonmin:lonmax,latmin:latmax], zm[lonmin:lonmax,latmin:latmax],
#  zdir='z', offset=np.min(zm[lonmin:lonmax,latmin:latmax]), cmap=cmocean.cm.deep_r)
# fig.colorbar(cs, ax=ax)
if dch['analytical'] == True:
    pass
else:
    pfun.add_coast(ax)
#pfun.dar(ax)
#ax.axis(ax_lims)
#ax.axis([-1,1,44, 47])
# ax.set_title(in_fn.name)
# make the panes transparent
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# make the grid lines transparent
ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
# label axis
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Depth [m]')
# ax.view_init(elev=30, azim=20)
ax.view_init(elev=40, azim=45)
#ax.text(.05, .95, Gr['gridname'], transform=ax.transAxes)
#ax.text(.95, .05, str(mask_rho.shape), ha='right', transform=ax.transAxes)

# plot cross-section
ax = fig.add_subplot(6,1,4)
ax.plot(y[:,0]/1000,zm[:,95])
ax.set_xlabel('distance from mouth [km]')
ax.set_ylabel('depth [m]')
# ax.set_xlim([47,48.5])

if do_riv:
    # gfp.add_river_tracks(Gr, ds, ax)
    print('Skip plotting river tracks')

if False:    
    # mask
    fig = plt.figure()
    ax = fig.add_subplot(111)
    tt = ['rho', 'u', 'v']
    sym = dict(zip(['rho', 'u', 'v'],['o','>','^']))
    c = dict(zip(['rho', 'u', 'v'],['b','orange','r']))
    for t in tt:
        x = ds['lon_'+t].values
        y = ds['lat_'+t].values
        m = ds['mask_'+t].values
        ax.plot(x, y, sym[t], c=c[t], alpha=.2, ms=3)
        ax.plot(x[m==1], y[m==1], sym[t], c=c[t], ms=3)
    
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(ax_lims)
    ax.set_title(in_fn.name)
    ax.text(.05, .95, Gr['gridname'], transform=ax.transAxes)
        
ds.close()
plt.savefig(out_dir / '{}_3d.png'.format(Gr['gridname']))
# plt.show()
