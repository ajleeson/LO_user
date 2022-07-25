"""
Plot grid to have a look at it. Accepts an optional command line argument
to look at a grid other than the one set in gfun.py.
"""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='',
        type=str)
args = parser.parse_args()

import cmocean

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
zm[mask_rho == 0] = np.nan

# PLOTTING
plt.close('all')
pfun.start_plot(figsize=(12,12))

# bathymetry
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)
cs = ax.pcolormesh(plon, plat, zm, vmin=-150, vmax=0, cmap=plt.get_cmap(cmocean.cm.deep_r))
fig.colorbar(cs, ax=ax)
if dch['analytical'] == True:
    pass
else:
    pfun.add_coast(ax)
pfun.dar(ax)
ax.axis(ax_lims)
ax.set_title(in_fn.name)
ax.text(.05, .95, Gr['gridname'], transform=ax.transAxes)
ax.text(.95, .05, str(mask_rho.shape), ha='right', transform=ax.transAxes)
if do_riv:
    gfp.add_river_tracks(Gr, ds, ax)

# plot wwtps if they exist
do_wwtp = False
wwtp_fn = Gr['wwtp_dir'] / 'wwtp_loc_info.csv'
# read wwtp lat lon info
if wwtp_fn.is_file():
    do_wwtp = True
    wwtp_df = pd.read_csv(wwtp_fn)
    # print(wwtp_df)
if do_wwtp:
    # plot wwtp locations on grid
    ax.scatter(wwtp_df['lon'],wwtp_df['lat'], color='black', label='wwtps')
    # print labels
    for i,wwtp in enumerate(wwtp_df['dname']):
        wwtp_lon = wwtp_df['lon'][i]
        wwtp_lat = wwtp_df['lat'][i]+0.03
        ax.text(wwtp_lon, wwtp_lat, wwtp, fontsize=12, horizontalalignment='center')

# plot point sources linked to the wwtp if the point sources have been created
do_ps = False
ps_fn = Gr['gdir'] / 'roms_wwtp_info.csv'
# read point source location data
if ps_fn.is_file():
    do_ps = True
    ps_df = pd.read_csv(ps_fn)
if do_ps:
    # plot point source locations on grid
    X = lon[0,:]
    Y = lat[:,0]
    ps_lon = [X[int(ind)] for ind in ps_df['col_py']]
    ps_lat = [Y[int(ind)] for ind in ps_df['row_py']]
    ax.scatter(ps_lon,ps_lat, color='deeppink', marker='x', s=20, label='point source')
    for i,ps in enumerate(ps_df['wname']):
        ax.plot([wwtp_df['lon'][i], ps_lon[i]],
        [wwtp_df['lat'][i], ps_lat[i]],
        color='deeppink', linewidth=0.5)
        ax.legend(loc='best',fontsize=12)
     

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

plt.show()
