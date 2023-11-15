"""
This focuses on property-property plots and obs-mod plots.

It specializes on model-model-obs comparisons, then with a bunch of
choices for filtering the data based on source, season, and depth.

Hence it is primarily a tool for model development: is one version
different of better than another?
"""
import sys
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
from lo_tools import plotting_functions as pfun
from matplotlib.patches import Rectangle
import xarray as xr
import cmocean
from lo_tools import Lfun, zfun, zrfun
Ldir = Lfun.Lstart()

testing = False

#############################
##       User Inputs       ##
#############################

month = 'Dec'
vn = 'DO (uM)' #'DIN (uM)'

#############################

# define fontsizes
fs_header = 24
fs_label = 20

year = '2017'
in_dir = Ldir['parent'] / 'LO_output' / 'obsmod'

plt.close('all')

# Salish Sea
lat_low = 46.9
lat_high = 49.7
lon_low = -124.5
lon_high = -122.1

# # Puget Sound only
# lat_low = 46.9
# lat_high = 48.2
# lon_low = -123.3
# lon_high = -122.1

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

# specify input (created by process_multi_bottle.py and process_multi_ctd.py)
otype = 'bottle'
in_fn = in_dir / ('multi_' + otype + '_' + year + '.p')
df0_dict = pickle.load(open(in_fn, 'rb'))

# where to put output figures
out_dir = Ldir['parent'] / 'LO_output' / 'obsmod_plots'
Lfun.make_dir(out_dir)

# add DIN field
for gtx in df0_dict.keys():
    df0_dict[gtx]['DIN (uM)'] = df0_dict[gtx]['NO3 (uM)'] + df0_dict[gtx]['NH4 (uM)']      
        
source = 'all'
time_range = 'all'
            
df_dict = df0_dict.copy()

# ===== FILTERS ======================================================
f_str = otype + ' ' + year + '\n\n' # a string to put for info on the map
ff_str = otype + '_' + year # a string for the output .png file name

# limit which sources to use
if source == 'all':
    # use df_dict as-is
    f_str += 'Source = all\n'
    ff_str += '_all'
        
# ====================================================================

# Plotting

if vn == 'DIN (uM)':
    vmin = -25
    vmax = 25
elif vn == 'DO (uM)':
    vmin = -155
    vmax = 155

pfun.start_plot(figsize=(8,12), fs=fs_label)

gtx_list = ['cas7_trapsV00_meV00']

alpha = 0.3
fig = plt.figure()

dti = pd.DatetimeIndex(df_dict[gtx].time)
# specify month
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
month_num = months.index(month) + 1
# get data from the correct month
mask = (dti.month==month_num)
obs_df = df_dict['obs'].loc[mask,:]
gtx_df = df_dict[gtx].loc[mask,:]
# sort by depth-- deepest to shallowest-- smallest to largest (for plotting later)
obs_df = obs_df.sort_values(by=['z'])
gtx_df = gtx_df.sort_values(by=['z'])


# Format station map
ax = fig.add_subplot(1,1,1)
ax.set_xticks([])
ax.set_yticks([])
ax.add_patch(Rectangle((lon_low, lat_low), lon_high-lon_low,lat_high-lat_low, facecolor='#FFFFFF'))
pfun.add_coast(ax, color='gray')
ax.axis([lon_low,lon_high,lat_low,lat_high])
pfun.dar(ax)
ax.set_xlabel('')
ax.set_ylabel('')
ax.text(.05,0.05,'bottle 2017',va='bottom',transform=ax.transAxes,fontweight='bold', fontsize = fs_label)
for border in ['top','right','bottom','left']:
        ax.spines[border].set_visible(False)
plt.xticks(fontsize=fs_label,rotation=30)
plt.yticks(fontsize=fs_label)

# Get station locations
x_obs_lon = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
                (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), 'lon'].to_numpy()
y_obs_lat = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
                (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), 'lat'].to_numpy()

# Get model data
obs_val = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
            (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), vn].to_numpy()
mod_val = gtx_df.loc[(gtx_df['lat']>lat_low) & (gtx_df['lat']<lat_high) &
                (gtx_df['lon']>lon_low) & (gtx_df['lon']<lon_high), vn].to_numpy()
# model minus observations (used for colorbar)
mod_minus_obs = mod_val - obs_val

# Get depths
depths = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
            (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), 'z'].to_numpy()
depths = depths * -5
sizes = [max(depth,50) for depth in depths]

# Plot data
newcmap = cmocean.cm.balance# cmocean.tools.crop(cmocean.cm.balance, vmin=vmin, vmax=vmax, pivot = 0)
map = plt.scatter(x_obs_lon,y_obs_lat, c=mod_minus_obs, cmap = newcmap, s = depths)
cbar = plt.colorbar(map)
plt.clim(vmin,vmax) 
cbar.ax.tick_params(labelsize=fs_label)

# Create legend
leg_szs = [10, 50, 150]
# szs = leg_szs
szs = [5*(leg_sz) for leg_sz in leg_szs]
# szs = [10*np.sqrt(leg_sz) for leg_sz in leg_szs]
l0 = plt.scatter([],[], s=szs[0], color='grey', edgecolors='none', alpha=0.5)
l1 = plt.scatter([],[], s=szs[1], color='grey', edgecolors='none', alpha=0.5)
l2 = plt.scatter([],[], s=szs[2], color='grey', edgecolors='none', alpha=0.5)
labels = ['< 10 m', '50 m', '150 m', '']
legend = ax.legend([l0, l1, l2], labels, fontsize = fs_label, markerfirst=False,
    title=r'Depth $\propto$ area', loc='best', labelspacing=1.5, borderpad=0.7)
plt.setp(legend.get_title(),fontsize=fs_label)


# format and save figure

fig.tight_layout()
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.93, wspace=0.05, hspace=0.2)
plt.suptitle(month + ' ' + vn + ': Model minus Observation', fontsize = fs_header)

print('Plotting ' + ff_str)
sys.stdout.flush()

plt.savefig(out_dir / (ff_str + '_bullseye_map' '.png'))
# plt.close('all')

    
