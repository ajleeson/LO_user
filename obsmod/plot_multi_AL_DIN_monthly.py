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

# define fontsizes
fs_header = 15
fs_label = 14

year = '2017'
in_dir = Ldir['parent'] / 'LO_output' / 'obsmod'

plt.close('all')

# Salish Sea
lat_low = 46.9
lat_high = 50.16
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
for otype in ['bottle']:#, 'ctd']:
    in_fn = in_dir / ('multi_' + otype + '_' + year + '.p')
    df0_dict = pickle.load(open(in_fn, 'rb'))

    # where to put output figures
    out_dir = Ldir['parent'] / 'LO_output' / 'obsmod_plots'
    Lfun.make_dir(out_dir)

    if otype == 'bottle':
        # add DIN field
        for gtx in df0_dict.keys():
            if gtx == 'cas6_v0_live':
                df0_dict[gtx]['DIN (uM)'] = df0_dict[gtx]['NO3 (uM)']
            else:
                df0_dict[gtx]['DIN (uM)'] = df0_dict[gtx]['NO3 (uM)'] + df0_dict[gtx]['NH4 (uM)']

    # loop over a variety of choices

    if otype == 'bottle':
        if True:
            source_list = ['all']
        else:
            source_list = ['nceiCoastal', 'nceiSalish', 'dfo1', 'ecology']
            #source_list = ['nceiSalish']
        
    elif otype == 'ctd':
        if True:
            source_list = ['all']
        else:
            source_list = ['dfo1', 'ecology']
            
    if True:
        time_range_list = ['all']
    else:
        time_range_list = ['spring','summer']
        
    if False:
        depth_range_list = ['all']
    else:
        depth_range_list = ['shallow','deep']
        
        
    for source in source_list:
        for depth_range in depth_range_list:
            for time_range in time_range_list:
            
                df_dict = df0_dict.copy()

                # ===== FILTERS ======================================================
                f_str = otype + ' ' + year + '\n\n' # a string to put for info on the map
                ff_str = otype + '_' + year # a string for the output .png file name

                # limit which sources to use
                if source == 'all':
                    # use df_dict as-is
                    f_str += 'Source = all\n'
                    ff_str += '_all'
                else:
                    # use just one source
                    f_str += 'Source = ' + source + '\n'
                    ff_str += '_' + source
                    for gtx in df_dict.keys():
                        df_dict[gtx] = df_dict[gtx].loc[df_dict[gtx].source==source,:]

                # depth range
                if depth_range == 'all':
                    pass
                elif depth_range == 'shallow':
                    # shallow water
                    zz = -30
                    f_str += 'Z above ' + str(zz) + ' [m]\n'
                    ff_str += '_shallow'
                    for gtx in df_dict.keys():
                        df_dict[gtx] = df_dict[gtx].loc[df_dict[gtx].z >= zz,:]
                elif depth_range == 'deep':
                    # deep water
                    zz = -30
                    f_str += 'Z below ' + str(zz) + ' [m]\n'
                    ff_str += '_deep'
                    for gtx in df_dict.keys():
                        df_dict[gtx] = df_dict[gtx].loc[df_dict[gtx].z <= zz,:]
        
                # time range
                if time_range == 'all':
                    pass
                elif time_range == 'spring':
                    # specific months
                    f_str += 'Months = [4,5,6]\n'
                    ff_str += '_spring'
                    for gtx in df_dict.keys():
                        dti = pd.DatetimeIndex(df_dict[gtx].time)
                        mask = (dti.month==4) | (dti.month==5) | (dti.month==6)
                        df_dict[gtx] = df_dict[gtx].loc[mask,:]
                elif time_range == 'summer':
                    # specific months
                    f_str += 'Months = [7,8,9]\n'
                    ff_str += '_summer'
                    for gtx in df_dict.keys():
                        dti = pd.DatetimeIndex(df_dict[gtx].time)
                        mask = (dti.month==7) | (dti.month==8) | (dti.month==9)
                        df_dict[gtx] = df_dict[gtx].loc[mask,:]
                # ====================================================================

                # Plotting

                fs = 12
                pfun.start_plot(figsize=(20,12), fs=fs)

                gtx_list = ['cas7_trapsV00_meV00']

                alpha = 0.3
                fig = plt.figure()
                # plt.subplots_adjust(wspace=0.1, hspace=0.1)

                if otype == 'bottle':
                    month_list = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
                    jj_list = [1,2,3,4,6,7,8,9,11,12,13,14] # indices for the data plots

                lim_dict = {'DIN (uM)':(0,50)}

                for ii,month in enumerate(month_list):
                    jj = jj_list[ii]
                    mon_num = ii + 1
                    if otype == 'bottle':
                        ax = fig.add_subplot(3,5,jj)
                    elif otype == 'ctd':
                        ax = fig.add_subplot(2,3,jj)
                    vn = 'DIN (uM)'
                    
                    # Get and plot model vs. observational data in the region
                    gtx = gtx_list[0]

                    # get data from the correct month
                    dti = pd.DatetimeIndex(df_dict[gtx].time)
                    mask = (dti.month==mon_num)
                    obs_df = df_dict['obs'].loc[mask,:]#df_dict['obs']
                    gtx_df = df_dict[gtx].loc[mask,:]#df_dict[gtx]
                    x = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
                                (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), vn].to_numpy()
                    y = gtx_df.loc[(gtx_df['lat']>lat_low) & (gtx_df['lat']<lat_high) &
                                    (gtx_df['lon']>lon_low) & (gtx_df['lon']<lon_high), vn].to_numpy()
                    ax.plot(x,y,marker='.',ls='',color='darkturquoise', alpha=alpha, markersize = 20)
    
                    # Add statistics to plot
                    if (not np.isnan(x).all()) and (not np.isnan(y).all()) and (len(x) > 0) and (len(y) > 0):
                        bias = np.nanmean(y-x)
                        rmse = np.sqrt(np.nanmean((y-x)**2))
                        ax.text(.9,0.05,'bias=%0.1f, rmse=%0.1f' % (bias,rmse),c='k',
                            transform=ax.transAxes, ha='right', fontweight='bold', bbox=pfun.bbox,
                            fontsize=fs_label,style='italic')

                    # Set tick labels
                    if otype == 'bottle':
                        if jj in [11,12,13,14]:
                            ax.set_xlabel('Observed', fontsize = fs_header)
                        else:
                            ax.tick_params(axis='x', colors='white')
                        if jj in [1,6,11]:
                            ax.set_ylabel('Modeled', fontsize = fs_header)
                        else:
                            ax.tick_params(axis='y', colors='white')
                    plt.xticks(fontsize=fs_label)
                    plt.yticks(fontsize=fs_label)
        
                    # Add label and 1-1 line
                    var = vn
                    ax.text(.07,.87,month,transform=ax.transAxes, fontweight='bold', fontsize = fs_header)
                    ax.axis([lim_dict[vn][0], lim_dict[vn][1], lim_dict[vn][0], lim_dict[vn][1]])
                    ax.plot([lim_dict[vn][0], lim_dict[vn][1]], [lim_dict[vn][0], lim_dict[vn][1]],'-k')
                    ax.grid(True,color='w',linewidth=2)

                    # format figure
                    plt.gca().set_aspect('equal', adjustable='box')
                    ax.set_facecolor('#EEEEEE')
                    for border in ['top','right','bottom','left']:
                        ax.spines[border].set_visible(False)
    
                # station map
                if otype == 'bottle':
                    ax = fig.add_subplot(1,5,5)
                # Plot stations locations
                obs_df = df_dict['obs']
                x_obs = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
                                (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), 'lon'].to_numpy()
                y_obs = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
                                (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), 'lat'].to_numpy()
                ax.plot(x_obs,y_obs,marker='.',ls='',color='k', alpha=alpha)
                ax.set_xticks([])
                ax.set_yticks([])


                ax.add_patch(Rectangle((lon_low, lat_low), lon_high-lon_low,lat_high-lat_low, facecolor='#EEEEEE'))
                plt.pcolormesh(plon, plat, zm, vmin=-10, vmax=0, cmap=plt.get_cmap(cmocean.cm.ice))
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


                # format and save figure

                fig.tight_layout()
                plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.93, wspace=0.05, hspace=0.2)

                plt.suptitle(depth_range + ' ' + vn, fontsize = 24)
                
                print('Plotting ' + ff_str)
                sys.stdout.flush()
                
                if testing:
                    plt.show()
                else:
                    plt.savefig(out_dir / (ff_str + '_DIN_monthly' '.png'))
                    # plt.close('all')

    
