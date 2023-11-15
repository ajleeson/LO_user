"""
Code to plot the processed ORCA mooring data, combining all stations
on a single plot, and showing seasonal profiles of one property per plot.
"""

import pandas as pd
import xarray as xr
import numpy as np
import gsw
import cmocean
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from warnings import filterwarnings
filterwarnings('ignore') # skip warning about all-nan layers

from lo_tools import plotting_functions as pfun
from lo_tools import Lfun, zrfun
Ldir = Lfun.Lstart()
from lo_tools import obs_functions as obfun

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-test','--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
testing = args.testing

year = '2017'
source = 'orca'
otype = 'moor' # introducing a new "otype" beyond ctd and bottle
gtx = 'cas7_trapsV00_meV00'
# out_dir = Ldir['LOo'] / 'obs' / source / otype
out_dir_model = Ldir['LOo'] / 'extract' / gtx / otype / source
out_dir_orca = Ldir['LOo'] / 'obs' / source / otype

plot_out_dir = Ldir['LOo'] / 'obs' / source / otype / 'figs' # general archive of figures
Lfun.make_dir(plot_out_dir)

sn_list = ['CI','PW','NB','DB','HP','TW']
axnum_dict = dict(zip(sn_list,[1,2,4,5,7,8]))

# get the grid data
ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
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

##########################################################
##                  Plot model output                   ##
##########################################################

axlim_dict = {'salt':(24,32),'temp':(5,18),'oxygen':(0,400),'SIG0':(18,24)}

if testing:
    vn_list = ['SA']
else:
    vn_list = ['salt','temp','oxygen', 'SIG0']
    
plt.close('all')
pfun.start_plot(figsize=(16,8))

for vn in vn_list:
    print( vn)
    fig = plt.figure()

    if vn == 'salt':
        var = 'SA'
    elif vn == 'temp':
        var = 'CT'
    elif vn == 'oxygen':
        var = 'DO (uM)'
    elif vn == 'SIG0':
        var = vn

    jj = 1
    for sn in sn_list:
        out_fn = out_dir_model / (sn + '_2017.01.01_2017.12.31.nc')
        ds = xr.open_dataset(out_fn)

        lat = ds.lat_rho.values
        lon = ds.lon_rho.values

        z = ds.z_rho.values
        t = ds.ocean_time.to_index()
    
        if vn != 'SIG0':
            fld = np.transpose(ds[vn].values) # pack as (z,time) for plotting

        # get pressure, salt, and temp
        p = [gsw.p_from_z(depth,lat) for depth in z.transpose()]
        temp = ds['temp'].transpose() # potential temperature
        salt = ds['salt'].transpose() # practical salinity

        # calculate absolute salinity, conservative temperature, and density
        salt_abs = gsw.conversions.SA_from_SP(salt, p, lon, lat)
        temp_cons = gsw.conversions.CT_from_pt(salt_abs, temp)
        if vn == 'salt':
            fld = salt_abs
        elif vn == 'temp':
            fld = temp_cons
        elif vn == 'SIG0':
            fld = gsw.sigma0(salt_abs,temp_cons)
    
        # time-mean profiles, by season
        ax = fig.add_subplot(3,3,axnum_dict[sn])
        # get season info dicts
        c_dict = obfun.season_c_dict
        daylist = obfun.season_daylist
        for ii in range(len(daylist)-1): # loop over seasons 0-3
            tmask = (t.dayofyear>=daylist[ii]) & (t.dayofyear<daylist[ii+1])
            this_fld = fld[:,tmask]
            # mask out depths where poor time coverage makes jumpy averages
            nans = np.isnan(this_fld).sum(axis=1)
            nans_norm = nans/np.min(nans)
            fld_mean = np.nanmean(this_fld,axis=1)
            fld_mean[nans_norm > 1.25] = np.nan
            # plot line for this season
            ax.plot(fld_mean,np.mean(z,0),'-',color=c_dict[ii],lw=4, alpha = 0.3)

##########################################################
##                  Plot observations                   ##
##########################################################

        print(' - '+sn)
        out_fn = out_dir_orca / (sn + '_daily.nc')
        ds_orca = xr.open_dataset(out_fn)
        # crop to year
        ds_orca = ds_orca.sel(time=slice(year+'-01-01', year+'-12-31'))

        z = ds_orca.z.values
        t = ds_orca.time.to_index()
    
        fld = np.transpose(ds_orca[var].values) # pack as (z,time) for plotting

        # map
        if jj == 1:
            axm = fig.add_subplot(133)
            # pfun.add_coast(axm)
            lonmin = -123.3
            lonmax = -122
            latmin = 47
            latmax = 48.3
            axm.axis([lonmin,lonmax,latmin,latmax])
            # pfun.dar(axm)
            axm.set_title('ORCA Station Location')
            axm.set_xlabel('Longitude')
            axm.set_ylabel('Latitude')
            newcmap = cmocean.cm.ice
            newcmap.set_bad('#EEEEEE',1.) # background color
            plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-10, vmax=0, cmap=newcmap)
            pfun.dar(axm)
            pfun.add_coast(axm,color='gray')
            for border in ['top','right','bottom','left']:
                axm.spines[border].set_visible(False)
            plt.xticks(rotation=30)
            plt.yticks(rotation=30)
        axm.plot(ds_orca.attrs['lon'],ds_orca.attrs['lat'],'or',ms=10,alpha=.6)
        axm.text(ds_orca.attrs['lon'],ds_orca.attrs['lat'],sn+': '+ds_orca.attrs['Station Name'],
            bbox=pfun.bbox,rotation=20)
        

        # get season info dicts
        c_dict = obfun.season_c_dict
        daylist = obfun.season_daylist
        for ii in range(len(daylist)-1): # loop over seasons 0-3
            tmask = (t.dayofyear>=daylist[ii]) & (t.dayofyear<daylist[ii+1])
            this_fld = fld[:,tmask]
            # mask out depths where poor time coverage makes jumpy averages
            nans = np.isnan(this_fld).sum(axis=1)
            nans_norm = nans/np.min(nans)
            fld_mean = np.nanmean(this_fld,axis=1)
            fld_mean[nans_norm > 1.25] = np.nan
            # plot line for this season
            ax.plot(fld_mean,z,'--',color=c_dict[ii],lw=2)
        
##########################################################
##                     Format figure                    ##
##########################################################

            # name of season
            if jj == 1:
                ax.text(.95,.35 - ii*.1,obfun.season_name_dict[ii],color=c_dict[ii],
                    fontweight='bold',ha='right',transform=ax.transAxes)
            # model or obs
            if jj == 2:
                ax.text(.5,.2,'Model: light bold solid',color='gray', fontsize = 10,
                    fontweight='bold',ha='center',transform=ax.transAxes)
                ax.text(.5,.1,'Observations: dark thin dashed',color='k', fontsize = 10,
                    ha='center',transform=ax.transAxes)
        ax.text(.05,.15,sn,color='k',
            fontweight='bold',transform=ax.transAxes)
        ax.set_ylim(-120,0)
        ax.set_xlim(axlim_dict[vn][0], axlim_dict[vn][1])
        ax.grid(True,color='w',linewidth=2)
        ax.set_facecolor('#EEEEEE')
        for border in ['top','right','bottom','left']:
            ax.spines[border].set_visible(False)
        plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.90, wspace=0.2, hspace=0.1)
        plt.tight_layout
        if jj in [5,6]:
            ax.set_xlabel('Seasonal Mean\n' + ds_orca[var].attrs['long_name']+' ['+ds_orca[var].attrs['units']+']')
        else:
            ax.set_xticklabels([])
        if jj in [1,3,5]:
            ax.set_ylabel('Z [m]')
        else:
            ax.set_yticklabels([])

        plt.suptitle(year + ' ' + vn)
        
        jj += 1
    
    if testing:
        plt.show()
    else:
        out_name = 'model_All_sta_seasonal_' + vn.replace(' ','_') + '.png'
        plt.savefig(plot_out_dir / out_name)
        plt.close()

