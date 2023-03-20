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
from matplotlib.patches import Rectangle
from matplotlib import pyplot as plt

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
# import pinfo
from importlib import reload
reload(pfun)
# reload(pinfo)


fn = 'results/roms_his_og.nc'

# START
ds = xr.open_dataset(fn)
# print(list(ds.keys()))

# find aspect ratio of the map
# aa = pfun.get_aa(ds)
# AR is the aspect ratio of the map: Vertical/Horizontal
# AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
fs = 14
hgt = 10
pfun.start_plot(fs=fs, figsize=(12,8))#, figsize=(int(hgt*2.5/AR),int(hgt)))
fig = plt.figure()
# PLOT CODE
vn_list = ['w', 'temp']
# print(ds['salt'])
ii = 1
for vn in vn_list:
    # if in_dict['auto_vlims']:
    #     pinfo.vlims_dict[vn] = ()
    ax = fig.add_subplot(1, len(vn_list), ii)
    x = ds['xi_rho'].values
    y = ds['eta_rho'].values
    v = ds[vn][0,-1,:,:].values
    # px, py = pfun.get_plon_plat(x,y)
    cs = ax.pcolormesh(x, y, v)#, cmap=cmap, vmin=vmin, vmax=vmax)
    # cs = pfun.add_map_field(ax, ds, vn)#, pinfo.vlims_dict,
            # cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
    fig.colorbar(cs)
    # pfun.add_coast(ax)
    # ax.axis(pfun.get_aa(ds))
    # Puget Sound:
    # ax.set(xlim=(-123.5, -122), ylim=(46.8, 49))
    plt.locator_params(axis='x', nbins=3)
    pfun.dar(ax)
    ax.set_title('Surface {}'.format(vn))
    # ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
    ax.set_xlabel('Longitude')
    if ii == 1:
        ax.set_ylabel('Latitude')
        # pfun.add_info(ax, in_dict['fn']) I commented this out so it is easier to see the point sources. Add back in later. --------------
        #pfun.add_windstress_flower(ax, ds)
        # pfun.add_bathy_contours(ax, ds, txt=False)

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
        # pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
    ii += 1
#fig.tight_layout()
# FINISH
ds.close()
pfun.end_plot()
# if len(str(in_dict['fn_out'])) > 0:
#     plt.savefig(in_dict['fn_out'])
#     plt.close()
# else:
plt.show()
