"""
Create map of terminal inlets
"""

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cmocean
import pandas as pd
from matplotlib.patches import Rectangle
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun

Ldir = Lfun.Lstart()

year = '2017'
# set up dates
startdate = year + '.01.01'
enddate = year + '.12.31'

jobname = 'twentyoneinlets'
# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')
# Get stations:
sta_dict = job_lists.get_sta_dict(jobname)

# colors
shallow = 'steelblue'
hypoxic = 'hotpink'
oxygenated = 'dimgrey'

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../../../LO_data/grids/cas7/grid.nc')
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

# Create map
plt.close('all')
fig = plt.figure(figsize=(7.7,8.7))
ax = fig.add_subplot(1,1,1)
lon_low = -123.6#42
lon_high =  -122
lat_low = 46.93
lat_high = 48.46
ax.add_patch(Rectangle((lon_low, lat_low), lon_high-lon_low,lat_high-lat_low, facecolor='#EEEEEE'))
ax.tick_params(axis='both', labelsize=12)
plt.xticks(rotation=30)
# plt.pcolormesh(plon, plat, zm, linewidth=0.5, vmin=-1.2, vmax=0, cmap=plt.get_cmap('Greys'))
plt.pcolormesh(plon, plat, zm, vmin=-8, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))

# get terminal inlet locations
for stn,station in enumerate(sta_dict): # stations: 
    if station == 'lynchcove2':
        continue
    # get segment information
    seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
    seg_df = pd.read_pickle(seg_name)
    ji_list = seg_df[station+'_p']['ji_list']
    jj = [x[0] for x in ji_list]
    ii = [x[1] for x in ji_list]
    # set everything to nan that is not in the inlet
    # first make everything a nan
    inlet_loc = np.full(zm.shape, np.nan) 
    # set values of 1 for everything that is in the inlet
    inlet_loc[jj,ii] = 20
    # add inlet locations
    # plt.pcolormesh(plon, plat, inlet_loc, linewidth=0.5, vmin=0, vmax=65, cmap=plt.get_cmap('coolwarm'))
    plt.pcolormesh(plon, plat, inlet_loc, linewidth=0.5, vmin=0, vmax=35, cmap=plt.get_cmap(cmocean.cm.ice))

# oxygenated inlets
for stn,station in enumerate(['sinclair','quartermaster','dyes','crescent','carr','elliot','commencement']):#sta_dict): # stations: 
    # get segment information
    seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
    seg_df = pd.read_pickle(seg_name)
    ji_list = seg_df[station+'_p']['ji_list']
    jj = [x[0] for x in ji_list]
    ii = [x[1] for x in ji_list]
    # set everything to nan that is not in the inlet
    # first make everything a nan
    inlet_loc = np.full(zm.shape, np.nan) 
    # set values of 1 for everything that is in the inlet
    inlet_loc[jj,ii] = 70
    # add inlet locations
    # plt.pcolormesh(plon, plat, inlet_loc, linewidth=0.5, vmin=0, vmax=65, cmap=plt.get_cmap('coolwarm'))
    plt.pcolormesh(plon, plat, inlet_loc, linewidth=0.5, vmin=0, vmax=100, cmap=plt.get_cmap('Greys'))

# Hypoxic inlets
for stn,station in enumerate(['lynchcove','holmes','dabob','penn','portsusan','case']): # stations: 
    # get segment information
    seg_name = Ldir['LOo'] / 'extract' / 'tef2' / 'seg_info_dict_cas7_c21_traps00.p'
    seg_df = pd.read_pickle(seg_name)
    ji_list = seg_df[station+'_p']['ji_list']
    jj = [x[0] for x in ji_list]
    ii = [x[1] for x in ji_list]
    # set everything to nan that is not in the inlet
    # first make everything a nan
    inlet_loc = np.full(zm.shape, np.nan) 
    # set values of 1 for everything that is in the inlet
    inlet_loc[jj,ii] = 40
    # add inlet locations
    plt.pcolormesh(plon, plat, inlet_loc, linewidth=0.5, vmin=0, vmax=100, cmap=plt.get_cmap('spring'))

# format
# ax.axes.xaxis.set_visible(False)
# ax.axes.yaxis.set_visible(False)
# ax.set_xlim(-123.42, -122) # Puget Sound
# ax.set_ylim(46.93, 48.46) # Puget Sound
ax.set_xlim(lon_low, lon_high) # Puget Sound
ax.set_ylim(lat_low, lat_high) # Puget Sound

pfun.dar(ax)

# remove border
for border in ['top','right','bottom','left']:
    ax.spines[border].set_visible(False)

# Add inlet labels
for sta in sta_dict:
    if sta == 'lynchcove2':
        continue
    sta_lon = sta_dict[sta][0]
    sta_lat = sta_dict[sta][1]
    lon_off = 0
    lat_off = 0.04
    ha = 'center'
    color = shallow
    if sta in ['sinclair','quartermaster','dyes','crescent','carr','elliot','commencement']:
        color = oxygenated
    if sta in ['lynchcove','holmes','dabob','penn','portsusan','case']:
        color = hypoxic
    # move text to make things more visible
    if sta_lon > -122.61 or sta in ['henderson','budd']:
        lon_off = 0.03
        lat_off = 0
        ha = 'left'
    if sta_lon < -123 or sta in ['oak','sinclair','case']:
        lon_off = -0.03
        lat_off = 0
        ha = 'right'
    if sta in ['carr']:
        lat_off = 0.1
    if sta in ['commenecement','elliot','quartermaster']:
        lon_off = 0.05
    if sta in ['budd']:
        lon_off = -0.05
        lat_off = -0.12
    if sta in ['eld']:
        lon_off = -0.1
        lat_off = -0.08
    if sta in ['totten']:
        lon_off = -0.07
    if sta in ['henderson']:
        lat_off = -0.02
    if sta in ['lynchcove']:
        lat_off = 0.06
    if sta in ['portsusan']:
        lon_off = 0.05
    if sta in ['holmes']:
        lon_off = -0.05
        lat_off = -0.05
    if sta in ['penn','dabob']:
        ha = 'right'
        lon_off = -0.05
        lat_off = -0.016
    ax.text(sta_lon+lon_off,sta_lat+lat_off,sta,va='center',ha=ha,color=color,fontsize = 10, fontweight='bold')

# add labels
ax.add_patch(Rectangle((-123.58, 48.12), 0.7,0.3, facecolor='white',alpha=0.7,edgecolor='gray'))
ax.text(-123.54,48.37,'Shallow inlets\n(mean depth < 10 m)',va='center',ha='left',
        color=shallow,fontsize = 11, fontweight='bold')
ax.text(-123.54,48.27,'Hypoxic deep inlets\n(mean depth > 10 m)',va='center',ha='left',
        color=hypoxic,fontsize = 11, fontweight='bold')
ax.text(-123.54,48.17,'Oxygenated deep inlets\n(mean depth > 10 m)',va='center',ha='left',
        color=oxygenated,fontsize = 11, fontweight='bold')

##########################################################
##                 Plot station locations               ##
##########################################################


# # dictionary of ecology stations
# ecol_stn = {
#     'sinclair': 'SIN001',
#     'elliot': 'ELB015',
#     'lynchcove': 'HCB007',
#     'commencement': 'CMB003',
#     'hammersley': 'OAK004',
#     # 'totten': 'TOT002',
#     'budd': 'BUD005'
# }

# # Plot cast locations
# for sta in sta_dict:
#     sta_lon = sta_dict[sta][0]
#     sta_lat = sta_dict[sta][1]
#     lon_off = 0
#     lat_off = 0.04
#     ha = 'center'
#     ecol = ''
#     # format depending on whether there is an ecology station or not
#     if sta in ecol_stn:
#          color = 'mediumorchid'
#          ecol = '\n('+ ecol_stn[sta] +')'
#          ax.plot(sta_lon,sta_lat,linestyle='none',marker='o',markersize=6,
#                  color='plum',markeredgecolor='k',markeredgewidth=0.8)
#     else:
#         color = 'k'

#     # move text to make things more visible
#     if sta_lon > -122.61 or sta in ['henderson','budd']:
#         lon_off = 0.03
#         lat_off = 0
#         ha = 'left'
#     if sta_lon < -123 or sta in ['oak','sinclair','case']:
#         lon_off = -0.03
#         lat_off = 0
#         ha = 'right'
#     if sta in ['carr']:
#         lat_off = 0.06
#     if sta in ['commenecement','elliot','quartermaster']:
#         lon_off = 0.05
#     if sta in ['budd']:
#         lon_off = -0.05
#         lat_off = -0.12
#     if sta in ['eld']:
#         lon_off = -0.1
#         lat_off = -0.08
#     if sta in ['totten']:
#         lon_off = -0.07
#     if sta in ['henderson']:
#         lat_off = -0.02
#     if sta in ['lynchcove']:
#         lat_off = 0.06
#     if sta in ['portsusan']:
#         lon_off = 0.05
#     if sta in ['holmes']:
#         lon_off = -0.05
#         lat_off = -0.05
#     if sta in ['penn','dabob']:
#         ha = 'right'
#         lon_off = -0.05
#         lat_off = -0.016
#     ax.text(sta_lon+lon_off,sta_lat+lat_off,sta+ecol,va='center',ha=ha,color=color,fontsize = 11)

# # add labels
# ax.add_patch(Rectangle((-123.39, 48.27), 0.6,0.15, facecolor='white',alpha=0.7,edgecolor='gray'))
# ax.text(-123.37,48.37,'Ecology monitoring\nstation (ID)',va='center',ha='left',fontweight='bold',color='mediumorchid',fontsize = 12)
# ax.text(-123.37,48.3,'No observations',va='center',ha='left',fontweight='bold',color='k',fontsize = 12)

plt.title('Inlet Locations',fontsize = 14)
plt.tight_layout()


# where to put output figures
out_dir = Ldir['LOo'] / 'pugetsound_DO' / ('DO_budget_'+startdate+'_'+enddate) / '2layer_figures'
Lfun.make_dir(out_dir)
plt.subplots_adjust(left=0.05, top=0.95, bottom=0.1, right=0.95)
plt.savefig('terminal_inlet_map.png')#,transparent=True)
