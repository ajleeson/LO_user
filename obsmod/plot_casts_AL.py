"""
Code to plot obs and mod casts at a given station, typically from ecology because
Those are monthly time series at a named location.
"""

import sys
import pandas as pd
import numpy as np
import pickle
from matplotlib.patches import Rectangle
import xarray as xr
import cmocean
import matplotlib.pyplot as plt
from lo_tools import plotting_functions as pfun
from lo_tools import Lfun, zfun, zrfun
Ldir = Lfun.Lstart()

year = '2017'
in_dir = Ldir['parent'] / 'LO_output' / 'obsmod'

# choices
sta_name = 'HCB004'
vn = 'DO (uM)'
# vn = 'SA'
# vn = 'CT'
# vn = 'Chl (mg m-3)'

"""
HCB003 is around Hoodsport
HCB004 is near Alderbrook
HCB007 is closer to the head of Lynch Cove
"""

# specify input (created by process_multi_bottle.py and process_multi_ctd.py)
otype = 'ctd'#'bottle'
in_fn = in_dir / ('multi_' + otype + '_' + year + '.p')
df_dict = pickle.load(open(in_fn, 'rb'))

source = 'ecology'
for gtx in df_dict.keys():
    df_dict[gtx] = df_dict[gtx].loc[df_dict[gtx].source==source,:]

# get station lat/lons
stn_df = df_dict['obs'].groupby('name', as_index=False).first()
stn_names = stn_df['name'].values
stn_lat = stn_df['lat'].values
stn_lon = stn_df['lon'].values
    
for gtx in df_dict.keys():
    df_dict[gtx] = df_dict[gtx].loc[df_dict[gtx].name==sta_name,:]

# dicts to help with plotting
lim_dict = {'SA':(15,36),'CT':(0,20),'DO (uM)':(0,500),
    'NO3 (uM)':(0,50),'NH4 (uM)':(0,10),'DIN (uM)':(0,50),
    'DIC (uM)':(1500,2500),'TA (uM)':(1500,2500),'Chl (mg m-3)':(0,20)}

c_list = ['k']#['darkturquoise','k']
c_dict = {'obs':'mediumorchid'}
ii = 0
for gtx in df_dict.keys():
    print(gtx)
    if gtx == 'obs' or gtx == 'cas7_trapsV00_meV00_AugVFCinis' or gtx == 'cas6_traps2_x2b':
        pass
    else:
        c_dict[gtx] = c_list[ii]
        ii += 1
    
    
# plotting
plt.close('all')

cid_list = df_dict['obs'].cid.unique()
cid_list.sort()


################################################################ 
##                 Plot all cast locations                    ##
################################################################ 

# Puget Sound only
lat_low = 47
lat_high = 48.15 #48.5
lon_low = -123.3
lon_high = -122.1

# # Hood Canal only
# lat_low = 47.3
# lat_high = 47.9
# lon_low = -123.3
# lon_high = -122.6

# get the grid data
ds = xr.open_dataset('../../LO_data/grids/cas7/grid.nc')
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

pfun.start_plot(figsize=(6,9))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.subplots_adjust(wspace=0, hspace=0.1)

# Plot map
ax.add_patch(Rectangle((lon_low, lat_low), lon_high-lon_low,lat_high-lat_low, facecolor='#EEEEEE'))
plt.pcolormesh(plon, plat, zm, vmin=-5, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
pfun.add_coast(ax, color='gray')
ax.axis([lon_low,lon_high,lat_low,lat_high])
pfun.dar(ax)
ax.set_xlabel('')
ax.set_ylabel('')
plt.xticks(rotation=30)
for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)

# Plot cast locations
ax.plot(stn_lon,stn_lat,linestyle='none',marker='o',color='k')
for i,stn in enumerate(stn_names):
    ax.text(stn_lon[i],stn_lat[i]+0.01,stn,va='bottom',ha='center',fontweight='bold', fontsize = 8)

plt.title('Cast Locations',fontsize = 16)
fig.tight_layout

plt.show()

################################################################ 
##                     Plot cast profiles                     ##
################################################################ 

# pfun.start_plot(figsize=(20,12))
pfun.start_plot(figsize=(15,8))
fig = plt.figure()

plt.subplots_adjust(wspace=0, hspace=0.1)

labels = ['(a) January', '(b) February', '(c) March', '(d) April', '(e) May', '(f) June',
          '(g) July', '(h) August', '(i) September', '(j) October', '(k) November', '(l) December']

runnames = ['Observations','Updated Model']#['Observations','Current LiveOcean', 'Updated Model']

# styles = ['-','-','--']

widths = ['3','1']#['3','3','1']

zbot = 0
ax_dict = dict()
for i,cid in enumerate(cid_list):
    ii = i + 1
    # mo = df_dict[gtx].loc[df_dict[gtx].cid==cid,'time'].to_list()[0].month
    ax = fig.add_subplot(2,6,ii)
    j = 0
    for gtx in df_dict.keys():
        if gtx == 'cas7_trapsV00_meV00_AugVFCinis' or gtx == 'cas6_traps2_x2b':
            pass
        else: 
            x = df_dict[gtx].loc[df_dict[gtx].cid==cid,vn].to_numpy()
            y = df_dict[gtx].loc[df_dict[gtx].cid==cid,'z'].to_numpy()
            zbot = np.min((zbot,np.min(y)))
            ax.plot(x,y,linestyle='-',c=c_dict[gtx],linewidth=widths[j],label=runnames[j]) # label = gtx
            ax.text(.1,.05,labels[i],transform=ax.transAxes,bbox=pfun.bbox)
            ax.set_xlim(lim_dict[vn])
            j += 1
    ax_dict[ii] = ax
    ax.set_xlim(lim_dict[vn])
    if ii not in [1,7]:
        ax.set_yticklabels([])
    else:
        ax.set_ylabel('Z [m]')
    if ii < 7:
        ax.set_xticklabels([])
    plt.tick_params(axis='x',rotation=30)
    
    
for ii in ax_dict.keys():
    ax = ax_dict[ii]
    ax.set_ylim(zbot,0)
    # format figure
    ax.set_facecolor('#EEEEEE')
    ax.grid(True,color='w',linewidth=2)
    for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
    
fig.tight_layout()
fig.subplots_adjust(top=0.86, bottom=0.1)

handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc=(0.25, 0.88), ncol=3,
           frameon = False, labelcolor = 'linecolor', prop=dict(weight='bold'))
    
fig.suptitle('Station: %s, %s, %s' % (sta_name,vn,year),
             fontweight='bold', fontsize=16)
plt.show()
