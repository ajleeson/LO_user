"""
Module of plotting functions.

Each function creates, and optionally saves, a plot of fields
from a ROMS history file.

INPUT: in_dict: a tuple with information to pass to the plot, such as:
- fn: text string with the full path name of the history file to plot
- fn_out: text string with full path of output file name
- auto_vlims: a boolean governing how color limits are set
- testing: a boolean for testing (e.g. shorter, faster particle tracking)
OUTPUT: either a screen image or a graphics file

"""
from matplotlib import markers
import numpy as np
import xarray as xr
import pickle
from datetime import datetime, timedelta
import pandas as pd
from cmocean import cm
import matplotlib.dates as mdates
import argparse
import math
import scipy.interpolate as interp
from scipy.optimize import curve_fit

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
if '_mac' in Ldir['lo_env']: # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
    
def P_basic(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'temp']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            # pfun.add_info(ax, in_dict['fn']) I commented this out so it is easier to see the point sources. Add back in later. --------------
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=False)

            # # plot wwtps if they exist
            # do_wwtp = False
            # wwtp_fn = Gr['wwtp_dir'] / 'wwtp_loc_info.csv'
            # # read wwtp lat lon info
            # if wwtp_fn.is_file():
            #     do_wwtp = True
            #     wwtp_df = pd.read_csv(wwtp_fn)
            #     # print(wwtp_df)
            # if do_wwtp:
            #     # plot wwtp locations on grid
            #     ax.scatter(wwtp_df['lon'],wwtp_df['lat'], color='black', label='wwtps')
            #     # print labels
            #     for i,wwtp in enumerate(wwtp_df['dname']):
            #         wwtp_lon = wwtp_df['lon'][i]
            #         wwtp_lat = wwtp_df['lat'][i]+0.05
            #         ax.text(wwtp_lon, wwtp_lat, wwtp, fontsize=14, horizontalalignment='center')

            # # plot point sources linked to the wwtp if the point sources have been created
            # do_ps = False
            # ps_fn = Ldir['data']/ 'grids'/ Gr['gridname'] / 'wwtp_info.csv'
            # # read point source location data
            # if ps_fn.is_file():
            #     do_ps = True
            #     ps_df = pd.read_csv(ps_fn)
            # if do_ps:
            #     # plot point source locations on grid
            #     lon = ds.lon_rho.values
            #     lat = ds.lat_rho.values
            #     X = lon[0,:]
            #     Y = lat[:,0]
            #     ps_lon = [X[int(ind)] for ind in ps_df['col_py']]
            #     ps_lat = [Y[int(ind)] for ind in ps_df['row_py']]
            #     ax.scatter(ps_lon,ps_lat, color='deeppink', marker='x', s=40, label='point sources')
            #     for i,ps in enumerate(ps_df['wname']):
            #         ax.plot([wwtp_df['lon'][i], ps_lon[i]],
            #         [wwtp_df['lat'][i], ps_lat[i]],
            #         color='deeppink', linewidth=1)
            #         ax.legend(loc='best',fontsize=14)

        elif ii == 2:
            ax.set_yticklabels([])
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_Fb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'Fb']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if vn == 'salt':
            cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                    cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                    vlims_fac=pinfo.range_dict[vn])
            ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
                    fontsize=1.2*fs)
        elif vn == 'Fb':
            # plot vertically integrateed buoyancy flux
            # calculate potential density
            import seawater as sw
            rho = sw.dens0(ds.salt.values.squeeze(), ds.temp.values.squeeze())
            # calculate vertically-integrated buoyancy flux
            Fb = -9.8 * np.sum(ds.AKs[0,1:-1,:,:].squeeze() * np.diff(rho, axis=0), axis=0).values
            Fb[ds.mask_rho.values==0] = np.nan
            plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
            cs = ax.pcolormesh(plon,plat, Fb, vmin=0, vmax=.5, cmap='YlOrRd')
            ax.set_title('Vertically Integrated Fb [W/m2]')
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        elif ii == 2:
            ax.set_yticklabels([])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_fancy(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'temp']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        if vn == 'salt':
            cmap = 'jet'
            vlims_fac = .5
        elif vn == 'temp':
            cmap = 'RdYlBu_r'
            vlims_fac = 1
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=cmap, fac=pinfo.fac_dict[vn], vlims_fac=vlims_fac)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            ax.set_yticklabels([])
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_dive_vort(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    fig = plt.figure()
    
    # create fields
    u = ds.u[0,-1,:,:].values
    v = ds.v[0,-1,:,:].values
    dx = 1/ds.pm.values
    dy = 1/ds.pn.values
    # dive is on the trimmed rho grid
    dive = np.diff(u[1:-1,:], axis=1)/dx[1:-1,1:-1] + np.diff(v[:,1:-1],axis=0)/dy[1:-1,1:-1]
    # vort is on the psi grid (plot with lon_rho, lat_rho)
    vort = np.diff(v,axis=1)/dx[1:,1:] - np.diff(u,axis=0)/dy[1:,1:]
    
    # set color limits
    vv = 2*np.nanstd(vort)
    
    # PLOT CODE
    if in_dict['auto_vlims']:
        pinfo.vlims_dict['vort'] = (-vv, vv)
        pinfo.vlims_dict['dive'] = (-vv, vv)
        
    vmin = pinfo.vlims_dict['vort'][0]
    vmax = pinfo.vlims_dict['vort'][1]
    
    for ii in [1,2]:
        ax = fig.add_subplot(1, 2, ii)
        cmap = 'RdYlBu_r'
        if ii == 1:
            plon, plat = pfun.get_plon_plat(ds.lon_rho[1:-1,1:-1].values, ds.lat_rho[1:-1,1:-1].values)
            cs = plt.pcolormesh(plon, plat, dive, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Divergence $[s^{-1}]$', fontsize=1.2*fs)
        elif ii == 2:
            cs = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, vort, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Vorticity $[s^{-1}]$', fontsize=1.2*fs)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            pass
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_dive_vort2(in_dict):
    # same as dive_vort but focused on a specific region
    
    # JdF:
    aa = [-125, -122.3, 47.8, 48.8]
    
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    # aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 6
    pfun.start_plot(fs=fs, figsize=(10,10))
    fig = plt.figure()
    
    # create fields
    u = ds.u[0,-1,:,:].values
    v = ds.v[0,-1,:,:].values
    dx = 1/ds.pm.values
    dy = 1/ds.pn.values
    # dive is on the trimmed rho grid
    dive = np.diff(u[1:-1,:], axis=1)/dx[1:-1,1:-1] + np.diff(v[:,1:-1],axis=0)/dy[1:-1,1:-1]
    # vort is on the psi grid (plot with lon_rho, lat_rho)
    vort = np.diff(v,axis=1)/dx[1:,1:] - np.diff(u,axis=0)/dy[1:,1:]
    
    # set color limits
    vv = 4*np.nanstd(vort)
    
    # PLOT CODE
    if in_dict['auto_vlims']:
        pinfo.vlims_dict['vort'] = (-vv, vv)
        pinfo.vlims_dict['dive'] = (-vv, vv)
        
    vmin = pinfo.vlims_dict['vort'][0]
    vmax = pinfo.vlims_dict['vort'][1]
    
    for ii in [1,2]:
        ax = fig.add_subplot(2, 1, ii)
        cmap = 'RdYlBu_r'
        if ii == 1:
            plon, plat = pfun.get_plon_plat(ds.lon_rho[1:-1,1:-1].values, ds.lat_rho[1:-1,1:-1].values)
            cs = plt.pcolormesh(plon, plat, dive, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Divergence $[s^{-1}]$', fontsize=1.2*fs)
        elif ii == 2:
            cs = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, vort, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Vorticity $[s^{-1}]$', fontsize=1.2*fs)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_ylabel('Latitude')
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            #pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            ax.set_xlabel('Longitude')
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_ri(in_dict):
    """
    Simplified Richardson number
    """
    # START
    fs = 10
    pfun.start_plot(fs=fs, figsize=(20,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PLOT CODE
    xrho = ds['lon_rho'][0,:].values
    yrho = ds['lat_rho'][:,0].values

    # define box
    aa = [-123.25, -122.1, 47, 48.75]
    ix0 = zfun.find_nearest_ind(xrho, aa[0])
    ix1 = zfun.find_nearest_ind(xrho, aa[1])
    iy0 = zfun.find_nearest_ind(yrho, aa[2])
    iy1 = zfun.find_nearest_ind(yrho, aa[3])

    h = ds.h[iy0:iy1, ix0:ix1].values
    rho_bot = ds.rho[0, 0, iy0:iy1, ix0:ix1].values
    rho_top = ds.rho[0, -1, iy0:iy1, ix0:ix1].values
    drho = rho_bot - rho_top
    u = ds.ubar[0, iy0:iy1, ix0-1:ix1].values
    v = ds.vbar[0, iy0-1:iy1, ix0:ix1].values
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    uu = (u[:, 1:] + u[:, :-1])/2
    vv = (v[1:, :] + v[:-1, :])/2
    spd2 = uu**2 + vv**2
    spd2[np.isnan(drho)] = np.nan
    spd2[spd2 < .001] = .001 # avoid divide by zero errors

    # approximate Richardson number
    rho0 = ds.rho0.values
    g = 9.8
    Ri = g * drho * h / (rho0 * spd2)

    # psi_grid coordinates
    x, y = np.meshgrid(ds.lon_u.values[0,ix0-1:ix1], ds.lat_v.values[iy0-1:iy1,0])

    # PLOTTING
    plt.close('all')
    pfun.start_plot(fs=10, figsize=(18,10))
    fig = plt.figure()

    xt = [-123.2, -122.2]
    yt = [47, 47.5, 48, 48.5]

    ax = fig.add_subplot(131)
    cs = ax.pcolormesh(x, y, drho, vmin=0, vmax=5, cmap=cm.dense)
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    ax.set_title(r'$\Delta\rho\ [kg\ m^{-3}]$')
    ax.set_xticks(xt)
    ax.set_yticks(yt)

    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(x, y, np.sqrt(spd2), vmin=0, vmax=2, cmap=cm.speed)
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    ax.set_title(r'Speed $[m\ s^{-1}]$')
    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.set_yticklabels([])

    ax = fig.add_subplot(133)
    cs = ax.pcolormesh(x, y, 4*Ri, vmin=0, vmax = 2, cmap='RdYlBu')
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    ax.set_title(r'$4 x Ri$')
    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.set_yticklabels([])
        
    fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_Chl_DO(in_dict):
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn_list = ['phytoplankton', 'oxygen']
    fs = 14
    ii = 1
    for vn in vn_list:
        if vn == 'phytoplankton':
            slev = -1
            stext = 'Surface'
        elif vn == 'oxygen':
            slev = 0
            stext = 'Bottom'
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                vlims_fac=pinfo.range_dict[vn], do_mask_edges=True)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('%s %s %s' % (stext, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        pfun.add_bathy_contours(ax, ds, txt=True)
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        ii += 1
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_DO_WA_shelf(in_dict):
    # Focus on bottom DO on the WA shelf
    aa = [-126.1, -123.7, 45.8, 48.8]
    xtl = [-126, -125, -124]
    ytl = [46, 47, 48]
    
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(7,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'oxygen'
    slev = 0
    stext = 'Bottom'
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    ax = fig.add_subplot(111)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
            vlims_fac=pinfo.range_dict[vn], do_mask_edges=True)
    fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('%s %s %s' % (stext, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
    ax.set_xlabel('Longitude')
    pfun.add_bathy_contours(ax, ds, txt=False)
    ax.set_ylabel('Latitude')
    ax.set_xticks(xtl)
    ax.set_yticks(ytl)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    
    pfun.add_windstress_flower(ax, ds, t_scl=0.5, t_leglen=0.1, center=(.85,.65), fs=12)
    # ADD MEAN WINDSTRESS VECTOR
    # t_scl: scale windstress vector (smaller to get longer arrows)
    # t_leglen: # Pa for wind stress vector legend
        
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_ths(in_dict):
    # Plot property-property plots, like theta vs. s
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(10,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    # make a potential density field
    import seawater as sw
    s0 = 25; s1 = 35
    th0 = 0; th1 = 20
    SS, TH = np.meshgrid(np.linspace(s0, s1, 50), np.linspace(th0, th1, 50))
    SIG = sw.dens0(SS, TH) - 1000
    S = zrfun.get_basic_info(in_dict['fn'], only_S=True)
    h = ds['h'].values
    z = zrfun.get_z(h, 0*h, S, only_rho=True)
    s = ds['salt'].values.squeeze()
    th = ds['temp'].values.squeeze()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Theta (deg C)')
    ax.contour(SS, TH, SIG, 20)
    nsub = 500
    alpha = .1
    mask = z > -10
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.r', alpha=alpha)
    mask = (z < -10) & (z > -200)
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.g', alpha=alpha)
    mask = z < -200
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.b', alpha=alpha)
    ax.set_xlim(s0, s1)
    ax.set_ylim(th0, th1)
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_debug(in_dict):
    # Focused on debugging
    vn_list = ['u', 'v', 'zeta']
    do_wetdry = False
    
    # START
    fs = 10
    pfun.start_plot(fs=fs, figsize=(8*len(vn_list),10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    ii = 1
    for vn in vn_list:
        if 'lon_rho' in ds[vn].coords:
            tag = 'rho'
        if 'lon_u' in ds[vn].coords:
            tag = 'u'
        if 'lon_v' in ds[vn].coords:
            tag = 'v'
        x = ds['lon_'+tag].values
        y = ds['lat_'+tag].values
        px, py = pfun.get_plon_plat(x,y)
        if vn in ['u', 'v']:
            v = ds[vn][0,-1,:,:].values
            vmin = -2
            vmax = 2
            cmap='hsv_r'
        elif vn == 'zeta':
            v = ds[vn][0,:,:].values
            h = ds.h.values
            mr = ds.mask_rho.values
            v[mr==0] = np.nan
            h[mr==0] = np.nan
            v = v + h
            vn = 'depth'
            vmin = 2
            vmax = 4
            cmap='RdYlGn'
        else:
            v = ds[vn][0, -1,:,:].values
        ax = fig.add_subplot(1, len(vn_list), ii)
        ax.set_xticks([])
        ax.set_yticks([])
        cs = ax.pcolormesh(px, py, v, cmap=cmap, vmin=vmin, vmax=vmax)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'], his_num=True)
        vmax, vjmax, vimax, vmin, vjmin, vimin = pfun.maxmin(v)
        ax.plot(x[vjmax,vimax], y[vjmax,vimax],'*y', mec='k', markersize=15)
        ax.plot(x[vjmin,vimin], y[vjmin,vimin],'oy', mec='k', markersize=10)
        ax.set_title(('%s ((*)max=%0.1f, (o)min=%0.1f)' % (vn, vmax, vmin)))
        ii += 1

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_layer(in_dict):
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn_list = ['oxygen', 'temp']
    z_level = -250
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, z_level)
        v_scaled = pinfo.fac_dict[vn]*laym
        vlims = pinfo.vlims_dict[vn]
        if len(vlims) == 0:
            vlims = pfun.auto_lims(v_scaled)
            pinfo.vlims_dict[vn] = vlims
        cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], v_scaled[1:-1,1:-1],
                           vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
        cb = fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        ax.set_title('%s %s on Z = %d (m)' % (pinfo.tstr_dict[vn], pinfo.units_dict[vn], z_level))
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'])
            ax.set_ylabel('Latitude')
            pfun.add_windstress_flower(ax, ds)
        if ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], zlev=z_level)
        ii += 1
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_bpress(in_dict):
    """
    Specialized plot related to bottom pressure anomalies.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    sta_dict = {
        'CE01':(-124.095, 44.6598), # Oregon Inshore (25 m)
        'CE02':(-124.304, 44.6393), # Oregon Shelf (80 m)
        'CE04':(-124.956, 44.3811), # Oregon Offshore (588 m)
        'PN01A':(-125.3983, 44.5096), # Slope Base (2905 m)
        }
    vn_list = ['salt', 'temp']
    z_level = -300
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, z_level)
        v_scaled = pinfo.fac_dict[vn]*laym
        vlims = pinfo.vlims_dict[vn]
        if len(vlims) == 0:
            if vn == 'salt':
                vlims = pfun.auto_lims(v_scaled, vlims_fac=0.3)
            elif vn == 'temp':
                vlims = pfun.auto_lims(v_scaled, vlims_fac=2)
            else:
                vlims = pfun.auto_lims(v_scaled)
            pinfo.vlims_dict[vn] = vlims
        cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], v_scaled[1:-1,1:-1],
                           vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
        cb = fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        ax.set_title('%s %s on Z = %d (m)' % (pinfo.tstr_dict[vn], pinfo.units_dict[vn], z_level))
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'])
            ax.set_ylabel('Latitude')
            pfun.add_windstress_flower(ax, ds)
        if ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], zlev=z_level, v_scl=5, v_leglen=0.3)
            for sta in sta_dict.keys():
                ax.plot(sta_dict[sta][0], sta_dict[sta][1], 'o', mfc='y', mec='k', ms=10)
        ii += 1
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_sect_al0(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -30

        # x = np.linspace(lon.min(), lon.max(), 500)
        # y = 45 * np.ones(x.shape)

        # ----------------------------------------------------------------------------
        lon_list = [-2, 0, 1, 2]
        x_res_list = [2500, 500, 500, 2500]
        lat_list = [43, 44.9, 45.1, 47]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list, lat_list, y_res_list)
        lonvals, latvals = np.meshgrid(Lon_vec, Lat_vec)

        x_intermed0, y_intermed = zfun.ll2xy(lonvals, latvals, 0, 45)
        #y = y[int(len(y)/2)::]

        x_intermed = (y_intermed/80000)*np.exp(-y_intermed/2e4)*np.sin(y_intermed/2e3)

        x_withinrange = np.abs(x_intermed[:,0]) <= 2

        x = x_intermed[:,0] * x_withinrange
        y = latvals[:,0]

        # only upper half
        x = x[int(len(x)/2):int(len(x)/1.3)]
        y = y[int(len(y)/2):int(len(y)/1.3)]

        # ----------------------------------------------------------------------------

    # or read in a section (or list of sections)
    else:
        tracks_path = Ldir['data'] / 'section_lines'
        tracks = ['Line_jdf_v0.p', 'Line_ps_main_v0.p']
        zdeep = -300
        xx = np.array([])
        yy = np.array([])
        for track in tracks:
            track_fn = tracks_path / track
            # get the track to interpolate onto
            pdict = pickle.load(open(track_fn, 'rb'))
            xx = np.concatenate((xx,pdict['lon_poly']))
            yy = np.concatenate((yy,pdict['lat_poly']))
        for ii in range(len(xx)-1):
            x0 = xx[ii]
            x1 = xx[ii+1]
            y0 = yy[ii]
            y1 = yy[ii+1]
            nn = 20
            if ii == 0:
                x = np.linspace(x0, x1, nn)
                y = np.linspace(y0,y1, nn)
            else:
                x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
                y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    #cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
     #       cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section


    # pfun.add_coast(ax)
    # aaf = [-125.5, -122.1, 46.8, 50.3] # focus domain
    # ax.axis(aaf)
    # pfun.dar(ax)

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    #pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    ratio = ((hgt*2.5/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    #ax.axis(pfun.get_aa(ds))
    #pfun.dar(ax)

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    # ax.set_xticks([-2, 0, 2])
    # ax.set_yticks([44, 45, 46, 47])


    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_sect_alpe(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -30

        y = np.linspace(1.04*lat.min(), 0.975*lat.max(), 500)
        x = np.zeros(y.shape)

    # or read in a section (or list of sections)
    else:
        tracks_path = Ldir['data'] / 'section_lines'
        tracks = ['Line_jdf_v0.p', 'Line_ps_main_v0.p']
        zdeep = -300
        xx = np.array([])
        yy = np.array([])
        for track in tracks:
            track_fn = tracks_path / track
            # get the track to interpolate onto
            pdict = pickle.load(open(track_fn, 'rb'))
            xx = np.concatenate((xx,pdict['lon_poly']))
            yy = np.concatenate((yy,pdict['lat_poly']))
        for ii in range(len(xx)-1):
            x0 = xx[ii]
            x1 = xx[ii+1]
            y0 = yy[ii]
            y1 = yy[ii+1]
            nn = 20
            if ii == 0:
                x = np.linspace(x0, x1, nn)
                y = np.linspace(y0,y1, nn)
            else:
                x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
                y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    #cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
     #       cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section


    # pfun.add_coast(ax)
    # aaf = [-125.5, -122.1, 46.8, 50.3] # focus domain
    # ax.axis(aaf)
    # pfun.dar(ax)

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    #pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    ratio = ((hgt*2.5/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    #ax.axis(pfun.get_aa(ds))
    #pfun.dar(ax)

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    # ax.set_xticks([-2, 0, 2])
    # ax.set_yticks([44, 45, 46, 47])


    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_alpe2(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'temp'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -20

        x = np.linspace(-1.4, 1.4, 500)
        y = 45.07 * np.ones(x.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    #cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
     #       cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    #pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    ratio = ((hgt*1.02/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    #ax.axis(pfun.get_aa(ds))
    #pfun.dar(ax)

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    # ax.set_xticks([-2, 0, 2])
    # ax.set_yticks([44, 45, 46, 47])


    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.5)

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_upw_v_vel(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'v'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_v']
        lat = G['lat_v']
        zdeep = -100

        y = np.linspace(44.8, 44, 500)
        x = np.zeros(y.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(2, 2, (1,2))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
            # cmap=cm.balance, fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*0.2/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)

    # section
    ax = fig.add_subplot(2, 2, (3,4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.5)

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_upw_u_vel(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'u'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_u']
        lat = G['lat_u']
        zdeep = -100

        y = np.linspace(44.8, 44, 500)
        x = np.zeros(y.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(2, 2, (1,2))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
            # cmap=cm.balance, fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*0.2/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)

    # section
    ax = fig.add_subplot(2, 2, (3,4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.5)

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_upw_spline(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """

    # Function to get interface depth from numerical model results
    def get_interface(v3,dist,sf):
        # get rho depth values
        zrho = [arr[0] for arr in v3['zrf']]
        # initialize array to save values
        zeta_depths = []

        # loop through for every distance away from the shelf
        for i in range(len(dist)):

            # apply cubic spline interpolation to salinity
            sal_spline = interp.CubicSpline(zrho, sf[:,i],
             bc_type='not-a-knot', extrapolate=None)

            # check if the watercolumn is more or less uniform in salinity (within 2%)
            if math.isclose(sf[0,i], sf[-1,i], rel_tol=0.02, abs_tol=0.0):
                # if so, then there is no interface
                zeta_depths = zeta_depths + [np.nan]

            else:
                # calculate depth of greatest salinity gradient
                ds_dz_spline = sal_spline.derivative()
                ds_dz_vals = ds_dz_spline(range(-100,0))

                # calculate depth of greatest ds/dz slope
                zeta_index = np.where(ds_dz_vals == np.min(ds_dz_vals))
                # create array from 0 to -100, and index into it
                depths = np.linspace(-100,-1,100)
                zeta_depth_curr = depths[zeta_index[0]]
                zeta_depths = zeta_depths + [zeta_depth_curr[0]]

        # Remove nans from data
        valid = ~(np.isnan(dist) | np.isnan(zeta_depths))
        p0 = [1,-0.5,-24.5] # initial guess
        # fit an exponential curve to data
        popt, pcov = curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, dist[valid], np.array(zeta_depths)[valid], p0)
        a = popt[0]
        b = popt[1]
        c = popt[2]
        # return the exponential fit of our interface based on the calculated interface from our spline.
        zeta_depths = a * np.exp(b * dist) + c

        return zeta_depths #zeta_depth_spline

    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -100

        y = np.linspace(44.8, 44, 500)
        x = np.zeros(y.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(2, 2, (1,2))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*0.2/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)

    # section
    ax = fig.add_subplot(2, 2, (3,4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # plot interface (numerical model)
    zeta_depth = get_interface(v3,dist,sf)
    ax.plot(dist,zeta_depth,color = 'cyan',linewidth = 2, label='interface depth (numerical)')

    # calculate interface (analytical model)
    kms = np.linspace(0,90,91) # km
    H = 24.5 # m
    tau = 0.1 # N m-2
    rho = 1023 # kg m-3
    F = tau/(rho*H) # m s-2
    f = 1e-4 # s-1
    g = 9.8 # m s-1
    gprime = g*(2.5/rho) # m s-1
    lambda_ = np.sqrt(gprime*H)/f # m
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    curr_time = T['dt']
    start_time = datetime(2020, 1, 1, 0, 0)
    t = (curr_time - start_time).total_seconds()
    # calculate zeta (includes offset from the starting interface depth)
    zeta_ana = [-H + F*H/(f*lambda_)*t*np.exp(-1*(km*1000)/lambda_) for km in kms]
    # plot interface (analytical model)
    ax.plot(kms,zeta_ana,color = 'xkcd:bubblegum pink',linewidth = 2,
     linestyle = '--', label='interface depth (analytical)')

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.5)

    # add legend
    ax.legend(fancybox=True, loc='lower left', framealpha = 0)

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_upw_avgdiff(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """

    # Function to get interface depth from numerical model results
    def get_interface(v3,dist,sf):
        # get rho depth values
        zrho = [arr[0] for arr in v3['zrf']]
        # initialize array to save values
        zeta_depths = []

        # loop through for every distance away from the shelf
        for i in range(len(dist)):

            # apply cubic spline interpolation to salinity
            sal_spline = interp.CubicSpline(zrho, sf[:,i],
             bc_type='not-a-knot', extrapolate=None)

            # initialize array of differences
            avg_diff = []

            # check if the watercolumn is more or less uniform in salinity (within 2%)
            if math.isclose(sf[0,i], sf[-1,i], rel_tol=0.02, abs_tol=0.0):
                # if so, then there is no interface
                zeta_depths = zeta_depths + [np.nan]
                avg_diff = avg_diff + [0]

            else:
                # calculate depth of greatest salinity gradient
                sal_vals = sal_spline(range(-100,0))
                # else, loop through depth in 1 m increments
                for j in range(-100,0):
                    # calculate magnitude of difference between average salinity above and below a depth
                    avg_bottom = np.mean(sal_vals[0:j])
                    avg_top = np.mean(sal_vals[j::])
                    diff = abs(avg_bottom - avg_top)
                    avg_diff = avg_diff + [diff]
                
                # calculate depth of greatest ds/dz slope
                zeta_index = avg_diff.index(np.nanmax(avg_diff))
                # create array from 0 to -100, and index into it
                depths = np.linspace(-100,-1,100)
                zeta_depth_curr = depths[zeta_index]
                zeta_depths = zeta_depths + [zeta_depth_curr]
        
        # # fit cubic spline to interface so we get a smooth line
        # zeta_depth_spline = interp.CubicSpline(dist, zeta_depths)
        # #  bc_type='not-a-knot', extrapolate=None)

        return zeta_depths #zeta_depth_spline

    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -100

        y = np.linspace(44.8, 44, 500)
        x = np.zeros(y.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(2, 2, (1,2))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*0.2/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)

    # section
    ax = fig.add_subplot(2, 2, (3,4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # plot interface (numerical model)
    zeta_depth = get_interface(v3,dist,sf)
    ax.plot(dist,zeta_depth,color = 'cyan',linewidth = 2, label=r'$\zeta$ numerical')

    # calculate interface (analytical model)
    kms = np.linspace(0,90,91) # km
    H = 24.5 # m
    tau = 0.1 # N m-2
    rho = 1023 # kg m-3
    F = tau/(rho*H) # m s-2
    f = 1e-4 # s-1
    g = 9.8 # m s-1
    gprime = g*(2.5/rho) # m s-1
    lambda_ = np.sqrt(gprime*H)/f # m
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    curr_time = T['dt']
    start_time = datetime(2020, 1, 1, 0, 0)
    t = (curr_time - start_time).total_seconds()
    # calculate zeta (includes offset from the starting interface depth)
    zeta_ana = [-H + F*H/(f*lambda_)*t*np.exp(-1*(km*1000)/lambda_) for km in kms]
    # plot interface (analytical model)
    ax.plot(kms,zeta_ana,color = 'xkcd:bubblegum pink',linewidth = 2,
     linestyle = '--', label=r'$\zeta$ analytical')

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.5)

    # add legend
    ax.legend(fancybox=True, loc='lower left', framealpha = 0)

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_upw_grad(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """

    # Function to get interface depth from numerical model results
    def get_interface(v3,dist,sf):
        # get rho depth values
        zrho = [arr[0] for arr in v3['zrf']]
        # initialize array to save values
        zeta_depths = []

        # calculate gradient of salinity array
        grad = np.gradient(sf)
        # sum gradient in x and y direction
        grad_sum = grad[0] + 10*grad[1]

        # loop through for every distance away from the shelf
        for i in range(len(dist)):

            # check if the watercolumn is more or less uniform in salinity (within 2%)
            if math.isclose(sf[0,i], sf[-1,i], rel_tol=0.02, abs_tol=0.0):
                # if so, then there is no interface
                zeta_depths = zeta_depths + [np.nan]

            else:
                # calculate depth of greatest salinity gradient in the water column
                grad_col = grad_sum[:,i]
                zeta_index = np.where(grad_col == np.min(grad_col))
                zeta_depth_curr = zrho[zeta_index[0][0]]
                zeta_depths = zeta_depths + [zeta_depth_curr]

        # Remove nans from data
        valid = ~(np.isnan(dist) | np.isnan(zeta_depths))
        p0 = [1,-0.5,-24.5] # initial guess
        # fit an exponential curve to data
        popt, pcov = curve_fit(lambda t, a, b, c: a * np.exp(b * t) + c, dist[valid], np.array(zeta_depths)[valid], p0)
        a = popt[0]
        b = popt[1]
        c = popt[2]
        # return the exponential fit of our interface based on the calculated interface from our spline.
        zeta_depths = a * np.exp(b * dist) + c

        return zeta_depths #zeta_depth_spline

    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -100

        y = np.linspace(44.8, 44, 500)
        x = np.zeros(y.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(2, 2, (1,2))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*0.2/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)

    # section
    ax = fig.add_subplot(2, 2, (3,4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # plot interface (numerical model)
    zeta_depth = get_interface(v3,dist,sf)
    ax.plot(dist,zeta_depth,color = 'cyan',linewidth = 2, label='interface depth (numerical)')

    # calculate interface (analytical model)
    kms = np.linspace(0,90,91) # km
    H = 24.5 # m
    tau = 0.1 # N m-2
    rho = 1023 # kg m-3
    F = tau/(rho*H) # m s-2
    f = 1e-4 # s-1
    g = 9.8 # m s-1
    gprime = g*(2.5/rho) # m s-1
    lambda_ = np.sqrt(gprime*H)/f # m
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    curr_time = T['dt']
    start_time = datetime(2020, 1, 1, 0, 0)
    t = (curr_time - start_time).total_seconds()
    # calculate zeta (includes offset from the starting interface depth)
    zeta_ana = [-H + F*H/(f*lambda_)*t*np.exp(-1*(km*1000)/lambda_) for km in kms]
    # plot interface (analytical model)
    ax.plot(kms,zeta_ana,color = 'xkcd:bubblegum pink',linewidth = 2,
     linestyle = '--', label='interface depth (analytical)')

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.5)

    # add legend
    ax.legend(fancybox=True, loc='lower left', framealpha = 0)

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_upw_contour(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    """

    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(9,13))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -100

        y = np.linspace(44.8, 44, 500)
        x = np.zeros(y.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(3, 2, (1,2))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*0.2/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)

    # section
    ax = fig.add_subplot(3, 2, (3,4))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs1 = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs1, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # create contour line at salinity = 31.5 to get interface (numerical model)
    cs2 = ax.contour(v3['distf'], v3['zrf'], sf, levels=[31.5], alpha=0)
    # get values from contour line
    points = cs2.collections[0].get_paths()[0]
    verts = points.vertices
    kms_num = verts[:,0]
    zeta_num = verts[:,1]
    # plot interface (analytical model)
    ax.plot(kms_num,zeta_num,color = 'cyan',linewidth = 2,
     linestyle = '-', label='interface depth (numerical)')

    # calculate interface (analytical model)
    kms = np.linspace(0,90,91) # km
    H = 24.5 # m
    tau = 0.1 # N m-2
    rho = 1023 # kg m-3
    F = tau/(rho*H) # m s-2
    f = 1e-4 # s-1
    g = 9.8 # m s-1
    gprime = g*(2.5/rho) # m s-1
    lambda_ = np.sqrt(gprime*H)/f # m
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    curr_time = T['dt']
    start_time = datetime(2020, 1, 1, 0, 0)
    t = (curr_time - start_time).total_seconds()
    # calculate zeta (includes offset from the starting interface depth)
    zeta_ana = [-H + F*H/(f*lambda_)*t*np.exp(-1*(km*1000)/lambda_) for km in kms]
    zeta_ana_to_compare = [-H + F*H/(f*lambda_)*t*np.exp(-1*(km*1000)/lambda_) for km in kms_num]
    # plot interface (analytical model)
    ax.plot(kms,zeta_ana,color = 'xkcd:bubblegum pink',linewidth = 2,
     linestyle = '--', label='interface depth (analytical)')

    # add legend
    ax.legend(fancybox=True, loc='lower left', framealpha = 0)

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.5)

    # -------------------------------------------------------------------
    # comparison of interface depth
    ax = fig.add_subplot(3, 2, (5,6))
    num_min_ana = zeta_num - zeta_ana_to_compare
    ax.plot(kms_num,num_min_ana)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(-8,H)
    ax.set_title(r'$\zeta_{num}-\zeta_{ana}$')
    ax.set_ylabel('Difference (m)')
    ax.set_xlabel('Distance (km)')
    ax.axhline(y=0,color='k',linestyle=':')

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.3)

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_alpe_withTides(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -30

        y = np.linspace(1.04*lat.min(), 0.975*lat.max(), 500)
        x = np.zeros(y.shape)


    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(2, 3, (1,4))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*2.5/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    # hardcoded latitude of mooring location from superplot
    lat = 45.4

    # get zeta and latitude
    zetas = v2['zeta']
    lats = v2['lat']
    
    # difference from lats to desired lat
    lats_diff = abs(lats - lat*np.ones(np.shape(lats)))

    # lat_index
    lat_index = np.where(lats_diff == np.min(lats_diff))

    # Location of zeta
    zeta_lat = lat_index
    

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    ax.plot(x[idist0],lat*np.ones(np.shape(y[idist0])),marker='*',markersize=15,color='hotpink')

    # section -------------------------------------------------------------------
    ax = fig.add_subplot(2, 3, (2,3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.plot(dist[zeta_lat],zetas[zeta_lat],marker='*',markersize=25,color='hotpink')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # Get time
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt_local = pfun.get_dt_local(T['dt'])

    # Tidal Amplitude -------------------------------------------------------------

    ax = fig.add_subplot(2, 3, (5,6))
    # ax.plot(dt_local,zetas[zeta_lat],marker='*',markersize=8,color='hotpink')

    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)

    # # get the day
    tm = T['dt'] # datetime
    TM = datetime(tm.year, tm.month, tm.day, tm.hour)

    # Tides
    alpha = 1
    # convert zeta timestamp to local time
    local_t_tides = [pfun.get_dt_local(x) for x in fdf['Timestamp']]

    ax.plot(local_t_tides, fdf['Tide Height (m)'].values, '-k',
        lw=1.5, alpha=alpha)

    # time marker
    ax.plot(fdf.loc[TM, 'Timestamp'], fdf.loc[TM, 'Tide Height (m)'],
        marker='*', color='hotpink', markersize=20)
    ax.set_ylim(-2,2)
    ax.set_title('Sea Surface Height (m)')
    ax.set_ylabel(r'$\zeta$ (m)')
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%D"))
    ax.tick_params('x', labelrotation=45)
    fig.tight_layout()

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_alpe2_withTides(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -30

        # CAUTION: the multipliers will change the start and end points of the section. Sometimes, I've accidentally
        # set it such that the multiplier times the min is larger than the multiplier times the max,
        # so the section line looked flipped. 
        y = np.linspace(44.8, 45.8, 500)
        x = np.zeros(y.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(2, 3, (1,4))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*1.0/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    # hardcoded latitude of mooring location from superplot
    lat = 45.3

    # get zeta and latitude
    zetas = v2['zeta']
    lats = v2['lat']
    
    # difference from lats to desired lat
    lats_diff = abs(lats - lat*np.ones(np.shape(lats)))

    # lat_index
    lat_index = np.where(lats_diff == np.min(lats_diff))

    # Location of zeta
    zeta_lat = lat_index
    

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    ax.plot(x[idist0],lat*np.ones(np.shape(y[idist0])),marker='*',markersize=15,color='hotpink')

    # section -------------------------------------------------------------------
    ax = fig.add_subplot(2, 3, (2,3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.plot(dist[zeta_lat],zetas[zeta_lat],marker='*',markersize=25,color='hotpink')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # Get time
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt_local = pfun.get_dt_local(T['dt'])

    # Tidal Amplitude -------------------------------------------------------------

    ax = fig.add_subplot(2, 3, (5,6))
    # ax.plot(dt_local,zetas[zeta_lat],marker='*',markersize=8,color='hotpink')

    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)

    # # get the day
    tm = T['dt'] # datetime
    TM = datetime(tm.year, tm.month, tm.day, tm.hour)

    # Tides
    alpha = 1
    # convert zeta timestamp to local time
    local_t_tides = [pfun.get_dt_local(x) for x in fdf['Timestamp']]

    ax.plot(local_t_tides, fdf['Tide Height (m)'].values, '-k',
        lw=1.5, alpha=alpha)

    # time marker
    ax.plot(fdf.loc[TM, 'Timestamp'], fdf.loc[TM, 'Tide Height (m)'],
        marker='*', color='hotpink', markersize=20)
    ax.set_ylim(-2,2)
    ax.set_title('Sea Surface Height (m)')
    ax.set_ylabel(r'$\zeta$ (m)')
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%D"))
    ax.tick_params('x', labelrotation=45)
    #set aspect ratio
    ax.set_aspect(5/2)
    fig.tight_layout()

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_soundspeed(in_dict):
    """
    Soundspeed section plot
    """
    import gsw

    ds = xr.open_dataset(in_dict['fn'])
    # create track by hand
    x = np.linspace(-124.85,-124.2, 100) # shelf only
    #x = np.linspace(-126,-124.2, 100) # shows SOFAR channel
    y = 47 * np.ones(x.shape)
    v2, v3, dist, idist0 = pfun.get_section(ds, 'salt', x, y, in_dict)
    s = v3['sectvarf']
    v2, v3, dist, idist0 = pfun.get_section(ds, 'temp', x, y, in_dict)
    th = v3['sectvarf']

    X = v3['distf']
    Z = v3['zrf']
    # adjust Z so surface is at 0
    Z = Z - Z[-1,:]

    p = gsw.p_from_z(Z, 47)
    SA = gsw.SA_from_SP(s, p, -125, 47)
    CT = gsw.CT_from_pt(SA, th)
    spd = gsw.sound_speed(SA, CT, p)

    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(16,9))
    fig, axes = plt.subplots(nrows=3, ncols=2)

    ax = axes[0,0]
    cs = ax.pcolormesh(X, Z, SA, cmap='jet')
    fig.colorbar(cs, ax=ax)
    ax.text(.95, .05, 'Absolute Salinity', transform=ax.transAxes, ha='right')

    ax = axes[1,0]
    cs = ax.pcolormesh(X, Z, CT, cmap='jet')
    fig.colorbar(cs, ax=ax)
    ax.text(.95, .05, 'Conservative Temperature', transform=ax.transAxes, ha='right')

    ax = axes[2,0]
    cs = ax.pcolormesh(X, Z, spd, cmap='jet')
    fig.colorbar(cs, ax=ax)
    ax.text(.95, .05, 'Soundspeed [m/s]', transform=ax.transAxes, ha='right')

    ax = axes[0,1]
    ax.plot(SA,Z, alpha=.2)
    ax.text(.05, .05, 'Absolute Salinity', transform=ax.transAxes, ha='left')

    ax = axes[1,1]
    ax.plot(CT,Z, alpha=.2)
    ax.text(.95, .05, 'Conservative Temperature', transform=ax.transAxes, ha='right')

    ax = axes[2,1]
    ax.plot(spd,Z, alpha=.2)
    ax.text(.95, .05, 'Soundspeed [m/s]', transform=ax.transAxes, ha='right')

    fig.suptitle(str(in_dict['fn']))

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    
def P_splash(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  Eventually I could automate making this new every day.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from PyCO2SYS import CO2SYS
    import seawater as sw
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from warnings import filterwarnings
    filterwarnings('ignore') # skip a warning from PyCO2SYS

    Ldir = Lfun.Lstart()

    do_topo = True

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x = ds.lon_psi.values
    y = ds.lat_psi.values
    th = ds['temp'][0,-1,1:-1,1:-1].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 3
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    def get_arag(ds, fn, aa, nlev):
        G = zrfun.get_basic_info(fn, only_G=True)
        # find indices that encompass region aa
        i0 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[0]) - 1
        i1 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[1]) + 2
        j0 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[2]) - 1
        j1 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[3]) + 2
        px = G['lon_psi'][j0:j1-1, i0:i1-1]
        py = G['lat_psi'][j0:j1-1, i0:i1-1]
        lat = G['lat_rho'][j0:j1,i0:i1] # used in sw.pres
        # first extract needed fields and save in v_dict
        v_dict = {}
        vn_in_list = ['temp', 'salt' , 'rho', 'alkalinity', 'TIC']
        for cvn in vn_in_list:
            L = ds[cvn][0,nlev,j0:j1,i0:i1].values
            v_dict[cvn] = L
        # ------------- the CO2SYS steps -------------------------
        # create pressure
        Ld = G['h'][j0:j1,i0:i1]
        Lpres = sw.pres(Ld, lat)
        # get in situ temperature from potential temperature
        Ltemp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, Lpres)
        # convert from umol/L to umol/kg using in situ dentity
        Lalkalinity = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
        Lalkalinity[Lalkalinity < 100] = np.nan
        LTIC = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
        LTIC[LTIC < 100] = np.nan
        CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
            Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
        # PH = CO2dict['pHout']
        # PH = zfun.fillit(PH.reshape((v_dict['salt'].shape)))
        ARAG = CO2dict['OmegaARout']
        ARAG = ARAG.reshape((v_dict['salt'].shape))
        ARAG = ARAG[1:-1, 1:-1]
        return px, py, ARAG

    # LARGE MAP
    ax = fig.add_subplot(121)
    cmap = 'RdYlBu_r'
    cs = ax.pcolormesh(x,y,th, cmap=cmap, vmin=11, vmax=20)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmap = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmap, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.2,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.98,.95,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
         ha='right', va='top', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Willapa and Grays Harbor
    aa = [-124.6, -123.65, 46, 47.2]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='g', alpha=1, linewidth=2, inset=0)

    # SMALL. MAP
    ax = fig.add_subplot(122)
    px, py, ARAG = get_arag(ds, fn, aa, nlev)
    cs = ax.pcolormesh(px,py,ARAG, cmap='coolwarm_r', vmin=0, vmax=3)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-124.5, -124])
    ax.set_yticks([46, 47])
    ax.set_xlabel('Longitude')
    ax.text(.98,.99,'Bottom water\nAragonite\nSaturation\nState',
         ha='right', va='top', weight='bold', transform=ax.transAxes)

    fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_splash2(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  This one is focused on the Salish Sea.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x = ds.lon_psi.values
    y = ds.lat_psi.values
    th = ds['temp'][0,-1,1:-1,1:-1].values
    ox = pinfo.fac_dict['oxygen'] * ds['oxygen'][0,0,1:-1,1:-1].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 3
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    # LARGE MAP
    ax = fig.add_subplot(121)
    cmap = 'RdYlBu_r'
    cs = ax.pcolormesh(x,y,th, cmap=cmap, vmin=11, vmax=20)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmap = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmap, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.5,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.03,.45,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
         ha='left', va='bottom', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Salish Sea
    aa = [-125.3, -122.1, 46.8, 50.3]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='g', alpha=1, linewidth=2, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    ax.text(-123.072,46.7866,'Washington', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-122.996,44.5788,'Oregon', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    
    ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
        style='italic',ha='center',va='center',rotation=-10)
    

    # SMALL MAP
    from cmocean import cm
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,ox, cmap=cm.oxy, vmin=0, vmax=10)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-125, -124, -123])
    ax.set_yticks([47, 48, 49, 50])
    ax.set_xlabel('Longitude')
    ax.text(.84,.95,'Salish Sea\n\nBottom Oxygen\n$[mg\ L^{-1}]$',
         ha='right', va='top', weight='bold', transform=ax.transAxes)
         
    # add labels
    ax.text(-122.8,49.335,'Fraser\nRiver',size=fs2,
        style='italic',ha='center',va='center',rotation=0)
    ax.text(-123.7,49.2528,'Strait of Georgia',size=fs2,
        style='italic',ha='center',va='center',rotation=-30)
    ax.text(-123.5,48.28,'Strait of Juan de Fuca',size=fs2,
        style='italic',ha='center',va='center',rotation=0,
        color='w')
    ax.text(-123.3,47.6143,'Puget\nSound',size=fs2,
        style='italic',ha='center',va='center',rotation=+55)
    ax.text(-122.3,48.48,'Skagit\nRiver',size=fs3,
        style='italic',ha='center',va='center',
        bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    ax.text(-123.173,48.44,'Haro\nStrait',size=fs3,
        style='italic',ha='center',va='center',
        color='w')

    fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_splash3(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  This one is focused on the Puget Sound.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(15,12))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x,y = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
    th = ds['temp'][0,-1,:,:].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 3
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    # LARGE MAP
    ax = fig.add_subplot(121)
    cs = ax.pcolormesh(x,y,th, cmap='RdYlBu_r', vmin=11, vmax=20)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmap = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap='gist_earth', shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.5,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.03,.45,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
         ha='left', va='bottom', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Puget Sound and the San Juans
    aa = [-123.4, -122, 47, 48.8]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='m', alpha=1, linewidth=2, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    ax.text(-123.072,46.7866,'Washington', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-122.996,44.5788,'Oregon', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    
    ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
        style='italic',ha='center',va='center',rotation=-10)
    

    # SMALL MAP
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,th, cmap='RdYlBu_r', vmin=11, vmax=20)
    # Inset colorbar
    # cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    # fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-123, -122.5, -122])
    ax.set_yticks([47, 48])
    ax.set_xlabel('Longitude')
    ax.text(.03,.5,'Puget Sound &\nSan Juans',
         ha='left', va='center', weight='bold', transform=ax.transAxes)
             

    #fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_splash4(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  This one is focused on the Puget Sound Salinity.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(15,12))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True
    cmap = 'Spectral_r'
    vlims = (25, 33) # full map
    vlims2 = (25, 33) # PS map

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x,y = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
    salt = ds['salt'][0,-1,:,:].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 1
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    # LARGE MAP
    ax = fig.add_subplot(121)
    cs = ax.pcolormesh(x,y,salt, cmap=cmap, vmin=vlims[0], vmax=vlims[1])
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmapt = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmapt, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.5,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.03,.45,'Surface Salinity $[g\ kg^{-1}]$\n'+tstr,
         ha='left', va='bottom', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Puget Sound and the San Juans
    aa = [-123.4, -122, 47, 48.8]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='m', alpha=1, linewidth=2, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    ax.text(-123.072,46.7866,'Washington', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-122.996,44.5788,'Oregon', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    
    ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
        style='italic',ha='center',va='center',rotation=-10)
    

    # SMALL MAP
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,salt, cmap=cmap, vmin=vlims2[0], vmax=vlims2[1])
    # Inset colorbar
    # cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    # fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmapt = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmapt, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-123, -122.5, -122])
    ax.set_yticks([47, 48])
    ax.set_xlabel('Longitude')
    ax.text(.03,.5,'Puget Sound &\nSan Juans',
         ha='left', va='center', weight='bold', transform=ax.transAxes)
             

    #fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_salt(in_dict):
    # Plot salinity maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'salt'
    vlims = (28.5, 33) # full map
    vlims2 = (22, 32) # PS map
    vlims3 = (29, 32) # PS section
    cmap = 'Spectral_r'

    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'].values
    lat = ds['lat_psi'].values
    v =ds[vn][0, -1, 1:-1, 1:-1].values
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nSalinity\n'
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')
    
    ax.text(.99,.97,'S range\n'+ str(vlims), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
        va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(0,1.5)
    ax.set_axis_off()

    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-5,20)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_oxygen(in_dict):
    # Plot bottom oxygen maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'oxygen'
    vlims = (0, 10) # full map
    vlims2 = (0, 10) # PS map
    vlims3 = (0, 10) # PS section
    from cmocean import cm
    cmap = cm.oxy

    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['Line_HC_thalweg_long.p']
    zdeep = -250
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'].values
    lat = ds['lat_psi'].values
    v =ds[vn][0, 0, 1:-1, 1:-1].values
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nBottom Oxygen\n'
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    cb.ax.tick_params(labelsize=.85*fs)
    ax.text(1, .85, r'$[mg\ L^{-1}]$', transform=ax.transAxes, fontsize=fs, ha='right')
    # ax.text(.99,.97,'S range\n'+ str(vlims), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    # ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    # ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
    #     va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nHood Canal', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(0,1.5)
    ax.set_axis_off()

    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-5,20)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_chl(in_dict):
    # Plot phytoplankton maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'phytoplankton'
    vlims = (0, 25) # full map
    vlims2 = (0, 25) # PS map
    vlims3 = (0, 25) # PS section
    cmap = 'Spectral_r'
    
    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'].values
    lat = ds['lat_psi'].values
    v =ds[vn][0, -1, 1:-1, 1:-1].values
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nPhytoplankton\n'+pinfo.units_dict[vn]+'\n'
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.99,.97,'range\n'+ str(vlims), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    # ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    # ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
    #     va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(0,1.5)
    ax.set_axis_off()

    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-5,20)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_tidal_avg_vel_alpev40d(in_dict):
    # Plot tidally averaged velocity profile from mooring data.

    # get gtagex
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # Define dates
    d_str = '2020.02.01_2020.02.10'

    # get gridname
    gridname = 'alpe'

    # mooring
    fn = Ldir['LOo'] / 'extract' / gtagex / 'moor' / gridname / ('superplot_' + d_str + '.nc')
    moor = xr.open_dataset(fn)

    v = moor['v']
    z_w = moor['z_w']
    z_rho = moor['z_rho']

    year = '2020'
    dates = pd.date_range(start='2/1/'+year, end='2/11/'+year, freq= '1H')[12::24]
    dates_local = [pfun.get_dt_local(x) for x in dates]

    # Calculate first difference of z_w
    dz = np.diff(z_w)

    # Calculate transport velocity
    vdz = v * dz

    # time average using 24-hour hanning window filter, and subsample per day (At noon)
    vdz_filtered = zfun.lowpass(np.array(vdz), f='hanning', n=24, nanpad=True)[12::24]
    dz_filtered = zfun.lowpass(np.array(dz), f='hanning', n=24, nanpad=True)[12::24]

    print('shape(vdz_filtered) = {}'.format(np.shape(vdz_filtered)))
    print('shape(dz_filtered) = {}'.format(np.shape(dz_filtered)))

    # divide average tranpsort velocity by average dz
    v_filtered = vdz_filtered/dz_filtered
    print('shape(v_filtered) = {}'.format(np.shape(v_filtered)))

    # calculate average depth
    avgz = zfun.lowpass(np.array(z_rho), f='hanning', n=24, nanpad=True)[12::24]
    print('shape(avgz) = {}'.format(np.shape(avgz)))

    # plot
    pfun.start_plot(figsize=(15,7))
    fig = plt.figure()
    ax = fig.add_subplot(1, 3, (1,2))

    n = np.shape(v_filtered)[0]

    # Make color gradient
    colors = cm.haline(np.linspace(0.1, 0.9, n))

    for i, c in zip(range(n),colors):
        ax.plot(v_filtered[i,:],avgz[i,:], label = dates_local[i].strftime('%b-%d'), color = c)

    ax.legend(loc = 'best')
    ax.axvline(0, color='k', linestyle='--')
    ax.set_ylabel('Z (m)')
    ax.set_title(r'Daily Average Velocity Profile ($m \ s^{-1}$)')

    # Calculate tidally-averaged velocity profile
    v_7day = v[72::,:] # get one week of half spring and half neap tide
    z_7day = z_rho[72::,:] 
    depth_7day = z_w[72::,:]
    n = np.shape(v_7day)[0]
    # average for whole day
    v_7dayavg = zfun.lowpass(np.array(v_7day), f='hanning', n=n, nanpad=True)[int(n/2)]
    z_7dayavg = zfun.lowpass(np.array(z_7day), f='hanning', n=n, nanpad=True)[int(n/2)]

    # Calculate average depth
    depth_avg = zfun.lowpass(np.array(depth_7day), f='hanning', n=n, nanpad=True)[int(n/2)]
    depth_avg = -1*depth_avg[0]

    # Add comparison of model and analytical velocity profile from CEWA 570
    ax = fig.add_subplot(1, 3, 3)

    # INITIALIZE CONSTANTS
    s1 = 0         # [ppt] river salinity
    s2 = 30        # [ppt] PS salinity
    del_s = s2-s1  # [ppt]
    Q = 1000       # [m^3/s] river discharge
    beta = 0.77e-3 # [ppt^-1] 
    H = depth_avg     # [m] estuary depth
    print(H)
    B = 5767       # [m] estuary width
    L = 10e3       # [m] estuary length
    nu_e = 1e-2    # [m^2/s] eddy viscosity
    g = 9.8        # [m/s^2] gravity

    # x3 vector
    x3 = np.linspace(-H,0,10000)

    # Define a temporary wind velocity scale u star
    u_star = 0 # 1e-2        # [m/s]

    # Calculate component velocity scales
    u_R = Q/(B*H)
    u_W = np.square(u_star)*H/nu_e
    u_E = g*beta*del_s*np.power(H,3)/(48*nu_e*L)

    # Calculate component velocity profiles
    u_river = -1*u_R*(3/2 - (3/2)*np.square(x3/H))
    u_wind = -1*u_W*((3/4)*np.square(x3/H) + (x3/H) + 1/4)
    u_density = -1*u_E*(1 - 9*np.square(x3/H) - 8*np.power(x3/H,3))

    # Calculate tidally-averaged velocity profile
    u_1 = u_river + u_wind + u_density

    ax.plot(u_1,x3,label = 'analytical')
    ax.plot(v_7dayavg,z_7dayavg, label = 'numerical')
    ax.axvline(0, color='k', linestyle='--')
    ax.legend(loc = 'best')

    # FINISH
    moor.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_tidal_avg_sal_alpev40d(in_dict):

    
    # Plot tidally averaged velocity profile from mooring data.

    # get gtagex
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # Define dates
    d_str = '2020.02.01_2020.02.10'

    # get gridname
    gridname = 'alpe'

    # mooring
    fn = Ldir['LOo'] / 'extract' / gtagex / 'moor' / gridname / ('superplot_' + d_str + '.nc')
    moor = xr.open_dataset(fn)

    z_rho = moor['z_rho']
    salt = moor['salt']

    year = '2020'
    dates = pd.date_range(start='2/1/'+year, end='2/11/'+year, freq= '1H')[12::24]
    dates_local = [pfun.get_dt_local(x) for x in dates]


    # time average using 24-hour hanning window filter, and subsample per day (At noon)
    sal_filtered = zfun.lowpass(np.array(salt), f='hanning', n=24, nanpad=True)[12::24]
    print('shape(sal_filtered) = {}'.format(np.shape(sal_filtered)))

    # calculate average depth
    avgz = zfun.lowpass(np.array(z_rho), f='hanning', n=24, nanpad=True)[12::24]
    print('shape(avgz) = {}'.format(np.shape(avgz)))

    # plot
    fig, ax = plt.subplots(1,1,figsize = (7,5))

    n = np.shape(sal_filtered)[0]

    # Make color gradient
    colors = cm.haline(np.linspace(0.1, 0.9, n))

    for i, c in zip(range(n),colors):
        plt.plot(sal_filtered[i,:],avgz[i,:], label = dates_local[i].strftime('%b-%d'), color = c)

    plt.legend(loc = 'best')
    plt.ylabel('Z (m)')
    plt.title(r'Daily Average Salinity Profile ($g \ kg^{-1}$)')


    # FINISH
    moor.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_tidal_avg_vel_alpe2v2mon(in_dict):
    # Plot tidally averaged velocity profile from mooring data.

    # get gtagex
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # Define dates
    d_str = '2020.02.01_2020.03.01'

    # get gridname
    gridname = 'alpe'

    # mooring
    fn = Ldir['LOo'] / 'extract' / gtagex / 'moor' / gridname / ('superplot_' + d_str + '.nc')
    moor = xr.open_dataset(fn)

    v = moor['v']
    z_w = moor['z_w']
    z_rho = moor['z_rho']

    year = '2020'
    dates = pd.date_range(start='2/1/'+year, end='3/2/'+year, freq= '1H')[36:-24:72]
    dates_local = [pfun.get_dt_local(x) for x in dates]

    # Calculate first difference of z_w
    dz = np.diff(z_w)

    # Calculate transport velocity
    vdz = v * dz

    # time average using 24-hour hanning window filter, and subsample per day (At noon)
    vdz_filtered = zfun.lowpass(np.array(vdz), f='godin', nanpad=True)[36:-24:72]
    dz_filtered = zfun.lowpass(np.array(dz), f='godin', nanpad=True)[36:-24:72]

    print('shape(vdz_filtered) = {}'.format(np.shape(vdz_filtered)))
    print('shape(dz_filtered) = {}'.format(np.shape(dz_filtered)))

    # divide average tranpsort velocity by average dz
    v_filtered = vdz_filtered/dz_filtered
    print('shape(v_filtered) = {}'.format(np.shape(v_filtered)))

    # calculate average depth
    avgz = zfun.lowpass(np.array(z_rho), f='godin', nanpad=True)[36:-24:72]
    print('shape(avgz) = {}'.format(np.shape(avgz)))

    # plot
    pfun.start_plot(figsize=(15,7))
    fig = plt.figure()
    ax = fig.add_subplot(1, 3, (1,2))

    n = np.shape(v_filtered)[0]

    # Make color gradient
    colors = cm.haline(np.linspace(0.1, 0.9, n))

    print('n={}'.format(n))
    print('num dates={}'.format(np.shape(dates_local)))

    for i, c in zip(range(n),colors):
        ax.plot(v_filtered[i,:],avgz[i,:], label = dates_local[i].strftime('%b-%d'), color = c)

    ax.legend(loc = 'best')
    ax.axvline(0, color='k', linestyle='--')
    ax.set_ylabel('Z (m)')
    ax.set_title(r'Daily Average Velocity Profile ($m \ s^{-1}$)')

    # Calculate tidally-averaged velocity profile ------------------------------------------------------------------
    v_7day = v[72::,:] # get one week of half spring and half neap tide
    z_7day = z_rho[72::,:] 
    depth_7day = z_w[72::,:]
    n = np.shape(v_7day)[0]
    # average for whole day
    v_7dayavg = zfun.lowpass(np.array(v_7day), f='hanning', n=n, nanpad=True)[int(n/2)]
    z_7dayavg = zfun.lowpass(np.array(z_7day), f='hanning', n=n, nanpad=True)[int(n/2)]

    # Calculate average depth
    depth_avg = zfun.lowpass(np.array(depth_7day), f='hanning', n=n, nanpad=True)[int(n/2)]
    depth_avg = -1*depth_avg[0]

    # Add comparison of model and analytical velocity profile from CEWA 570
    ax = fig.add_subplot(1, 3, 3)

    # INITIALIZE CONSTANTS
    s1 = 0         # [ppt] river salinity
    s2 = 30        # [ppt] PS salinity
    del_s = s2-s1  # [ppt]
    Q = 1000       # [m^3/s] river discharge
    beta = 0.77e-3 # [ppt^-1] 
    H = depth_avg     # [m] estuary depth
    print(H)
    B = 5767       # [m] estuary width
    L = 10e3       # [m] estuary length
    nu_e = 1e-2    # [m^2/s] eddy viscosity
    g = 9.8        # [m/s^2] gravity

    # x3 vector
    x3 = np.linspace(-H,0,10000)

    # Define a temporary wind velocity scale u star
    u_star = 0 # 1e-2        # [m/s]

    # Calculate component velocity scales
    u_R = Q/(B*H)
    u_W = np.square(u_star)*H/nu_e
    u_E = g*beta*del_s*np.power(H,3)/(48*nu_e*L)

    # Calculate component velocity profiles
    u_river = -1*u_R*(3/2 - (3/2)*np.square(x3/H))
    u_wind = -1*u_W*((3/4)*np.square(x3/H) + (x3/H) + 1/4)
    u_density = -1*u_E*(1 - 9*np.square(x3/H) - 8*np.power(x3/H,3))

    # Calculate tidally-averaged velocity profile
    u_1 = u_river + u_wind + u_density

    ax.plot(u_1,x3,label = 'analytical')
    ax.plot(v_7dayavg,z_7dayavg, label = 'numerical')
    ax.axvline(0, color='k', linestyle='--')
    ax.legend(loc = 'best')

    # FINISH
    moor.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_tidal_avg_sal_alpe2v2mon(in_dict):

    
    # Plot tidally averaged velocity profile from mooring data.

    # get gtagex
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # Define dates
    d_str = '2020.02.01_2020.03.01'

    # get gridname
    gridname = 'alpe'

    # mooring
    fn = Ldir['LOo'] / 'extract' / gtagex / 'moor' / gridname / ('superplot_' + d_str + '.nc')
    moor = xr.open_dataset(fn)

    z_rho = moor['z_rho']
    salt = moor['salt']

    year = '2020'
    dates = pd.date_range(start='2/1/'+year, end='3/02/'+year, freq= '1H')[36:-24:72]
    dates_local = [pfun.get_dt_local(x) for x in dates]


    # time average using godin filter, and subsample per day (At noon)
    sal_filtered = zfun.lowpass(np.array(salt), f='godin', nanpad=True)[36:-24:72]
    print('shape(sal_filtered) = {}'.format(np.shape(sal_filtered)))

    # calculate average depth
    avgz = zfun.lowpass(np.array(z_rho), f='godin', nanpad=True)[36:-24:72]
    print('shape(avgz) = {}'.format(np.shape(avgz)))

    # plot
    fig, ax = plt.subplots(1,1,figsize = (7,5))

    n = np.shape(sal_filtered)[0]

    # Make color gradient
    colors = cm.haline(np.linspace(0.1, 0.9, n))

    print('n={}'.format(n))
    print('num dates={}'.format(np.shape(dates_local)))

    for i, c in zip(range(n),colors):
        plt.plot(sal_filtered[i,:],avgz[i,:], label = dates_local[i].strftime('%b-%d'), color = c)

    plt.legend(loc = 'best')
    plt.ylabel('Z (m)')
    plt.title(r'Daily Average Salinity Profile ($g \ kg^{-1}$)')


    # FINISH
    moor.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_upw_v_vel(in_dict):
    """
    This plots compares the surface velocity along a section of the 
    simulation to an analytical solution for upwelling
    
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(9,12))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'v'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_v']
        lat = G['lat_v']
        zdeep = -100

        y = np.linspace(44.8, 44, 500)
        x = np.zeros(y.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(3, 2, (1,2))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
            # cmap=cm.balance, fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
    fig.colorbar(cs, ax=ax)

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*0.2/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)

    # plot cross-shore velocity (numerical model)
    ax = fig.add_subplot(3, 2, (3,4))
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(-0.1, 0.1)
    v_num = sf[-1,:] # surface layer
    ax.plot(dist,v_num, color = 'cyan', linewidth = 2, label='cross-shore velocity (numerical,surface)')
    # for i in range(1,15):
    #     ax.plot(dist,sf[-1*i,:], alpha=0.3)

    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # calculate cross-shore velocity (analytical model)
    H = 24.5 # m
    tau = 0.1 # N m-2
    rho = 1023 # kg m-3
    F = tau/(rho*H) # m s-2
    f = 1e-4 # s-1
    g = 9.8 # m s-1
    gprime = g*(2.5/rho) # m s-1
    lambda_ = np.sqrt(gprime*H)/f # m
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    curr_time = T['dt']
    start_time = datetime(2020, 1, 1, 0, 0)
    t = (curr_time - start_time).total_seconds()
    # cross-shore velocity
    v_ana = [-1*(F/f)*(1-np.exp(-1*(km*1000)/lambda_)) for km in dist] # negative because southward
    # plot velocity (analytical model)
    ax.plot(dist,v_ana,color = 'xkcd:bubblegum pink',linewidth = 2,
     linestyle = '--', label='cross-shore velocity (analytical)')

    # line at zero
    ax.axhline(y=0,color='k',linestyle=':')

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.3)

    # add legend
    ax.legend(fancybox=True, loc='upper right', framealpha = 0)

    # -------------------------------------------------------------------
    # comparison of velocity
    ax = fig.add_subplot(3, 2, (5,6))
    num_min_ana = v_num - v_ana
    ax.plot(dist,num_min_ana)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(-0.07,0.07)
    ax.set_title(r'$v_{num}-v_{ana}$')
    ax.set_ylabel(r'Difference (m $s^{-1}$)')
    ax.set_xlabel('Distance (km)')
    ax.axhline(y=0,color='k',linestyle=':')

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.3)

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_upw_u_vel(in_dict):
    """
    This plots compares the surface velocity along a section of the 
    simulation to an analytical solution for upwelling
    
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(9,12))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'u'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if True:
        lon = G['lon_u']
        lat = G['lat_u']
        zdeep = -100

        y = np.linspace(44.8, 44, 500)
        x = np.zeros(y.shape)

    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING

    # map with section line
    ax = fig.add_subplot(3, 2, (1,2))
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
            # cmap=cm.balance, fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
    fig.colorbar(cs, ax=ax)

    # -----------------------------------------------------------
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    ratio = ((hgt*0.2/AR)/(hgt))

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    # -----------------------------------------------------------

    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    print('xlims: {},{}'.format(x.min(),x.max()))
    print('ylims: {},{}'.format(y.min(),y.max()))
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)

    # plot cross-shore velocity (numerical model)
    ax = fig.add_subplot(3, 2, (3,4))
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(0, 1)
    # u_num = sf[-15:-1,:] # top 15 layers
    # u_num = np.mean(u_num, axis = 0)
    u_num = sf[-1,:] # surface layer
    ax.plot(dist,u_num, color = 'cyan', linewidth = 2, label='along-shore velocity (numerical,surface)')
    # for i in range(1,15):
    #     ax.plot(dist,sf[-1*i,:], alpha=0.3)

    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()

    # calculate cross-shore velocity (analytical model)
    H = 24.5 # m
    tau = 0.1 # N m-2
    rho = 1023 # kg m-3
    F = tau/(rho*H) # m s-2
    f = 1e-4 # s-1
    g = 9.8 # m s-1
    gprime = g*(2.5/rho) # m s-1
    lambda_ = np.sqrt(gprime*H)/f # m
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    curr_time = T['dt']
    start_time = datetime(2020, 1, 1, 0, 0)
    t = (curr_time - start_time).total_seconds()
    # cross-shore velocity
    u_ana = [F*t*np.exp(-1*(km*1000)/lambda_) for km in dist] # negative because southward
    # plot velocity (analytical model)
    ax.plot(dist,u_ana,color = 'xkcd:bubblegum pink',linewidth = 2,
     linestyle = '--', label='along-shore velocity (analytical)')

    # line at zero
    ax.axhline(y=0,color='k',linestyle=':')

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.3)

    # add legend
    ax.legend(fancybox=True, loc='upper right', framealpha = 0)

    # -------------------------------------------------------------------
    # comparison of velocity
    ax = fig.add_subplot(3, 2, (5,6))
    num_min_ana = u_num - u_ana
    ax.plot(dist,num_min_ana)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(-0.8,0.8)
    ax.set_title(r'$u_{num}-u_{ana}$')
    ax.set_ylabel(r'Difference (m $s^{-1}$)')
    ax.set_xlabel('Distance (km)')
    ax.axhline(y=0,color='k',linestyle=':')

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*0.3)

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()