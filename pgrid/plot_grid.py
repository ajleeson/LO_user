"""
Plot grid to have a look at it. Accepts an optional command line argument
to look at a grid other than the one set in gfun.py.
"""
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pickle
import cmocean
from lo_tools import Lfun

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='',type=str) # e.g. cas6
parser.add_argument('-dmax', default=5, type=int) # max depth for colormap [m]
parser.add_argument('-small', default=False, type=Lfun.boolean_string) # True for laptop size
args = parser.parse_args()
zmin = -args.dmax


import gfun
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
from lo_tools import plotting_functions as pfun
import gfun_plotting as gfp

Ldir = Lfun.Lstart()
out_dir = Ldir['LOo'] / 'analytical_grids'
Lfun.make_dir(out_dir)

testing = True
if testing:
    from importlib import reload
    reload(gfun)
    reload(gfp)

# select grid file
in_fn = gfun.select_file(Gr)

plt.close('all')

# load the default choices
try:
    dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))
except FileNotFoundError:
    # you could fill this in by hand if you wanted
    dch = {'analytical': False} # hack to make cas6 work

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
zm[mask_rho == 0] = np.nan

# PLOTTING
plt.close('all')
if args.small:
    figsize = (8,8)
else:
    figsize = (12,12)
pfun.start_plot(figsize=figsize)

# bathymetry
fig = plt.figure()
ax = fig.add_subplot(111)
cs = ax.pcolormesh(plon, plat, zm, vmin=zmin, vmax=0, cmap=cmocean.cm.deep_r)#cmap='Spectral_r')
# cs = ax.pcolormesh(plon, plat, zm, vmin=-120, vmax=-100, cmap='Spectral_r')
ax.tick_params(axis='x', labelrotation = 30)
fig.colorbar(cs, ax=ax)
if dch['analytical'] == True:
    pass
else:
    pfun.add_coast(ax)
pfun.dar(ax)
ax.axis(ax_lims)
ax.set_title(in_fn.name)
ax.text(.05, .95, Gr['gridname'], transform=ax.transAxes)
ax.text(.95, .05, str(mask_rho.shape), ha='right', transform=ax.transAxes, bbox=pfun.bbox)
if do_riv:
    gfp.add_river_tracks(Gr, ds, ax)

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
plt.savefig(out_dir / '{}.png'.format(Gr['gridname']))
plt.close('all')
