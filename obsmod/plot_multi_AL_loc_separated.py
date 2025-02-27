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

# get the grid data
ds = xr.open_dataset('../../LO_data/grids/cas6/grid.nc')
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
                pfun.start_plot(figsize=(17,12), fs=fs)

                gtx_list = ['cas7_trapsV00_meV00']
                label_list = ['Main Basin', 'South Sound', 'Hood Canal']
                c_dict = dict(zip(label_list,['navy', 'orchid', 'coral']))#['navy','olivedrab','orchid']))
                t_dict = dict(zip(label_list,[.05,0.15,0.25]))#,0.25])) # vertical position of stats text

                alpha = 0.3
                fig = plt.figure()
                plt.subplots_adjust(wspace=0, hspace=0)

                if otype == 'bottle':
                    vn_list = ['SA','CT','DO (uM)','NO3 (uM)','NH4 (uM)','DIN (uM)',
                        'DIC (uM)', 'TA (uM)', 'Chl (mg m-3)']
                    jj_list = [1,2,3,5,6,7,9,10,11] # indices for the data plots
                elif otype == 'ctd':
                    vn_list = ['SA','CT','DO (uM)','Chl (mg m-3)']
                    jj_list = [1,2,4,5] # indices for the data plots

                lim_dict = {'SA':(12,36),'CT':(0,24),'DO (uM)':(0,600),
                    'NO3 (uM)':(0,50),'NH4 (uM)':(0,10),'DIN (uM)':(0,50),
                    'DIC (uM)':(1500,2500),'TA (uM)':(1500,2500),'Chl (mg m-3)':(0,24)}

                for ii in range(len(vn_list)):
                    jj = jj_list[ii]
                    if otype == 'bottle':
                        ax = fig.add_subplot(3,4,jj)
                    elif otype == 'ctd':
                        ax = fig.add_subplot(2,3,jj)
                    vn = vn_list[ii]
                    # x = df_dict['obs'][vn].to_numpy()
                    obs_df = df_dict['obs']
                    
                    # Loop through the different regions ----------------------------------------------
                    for label in label_list:

                        if label == 'Strait of Georgia':
                            lat_low = 48.8
                            lat_high = 50.16
                            lon_low = -124.9
                            lon_high = -122.5

                        elif label == 'Main Basin':
                            lat_low = 47.3
                            lat_high = 48
                            lon_low = -122.6
                            lon_high = -122

                        elif label == 'South Sound':
                            lat_low = 47
                            lat_high = 47.3
                            lon_low = -123.3
                            lon_high = -122

                        elif label == 'Hood Canal':
                            lat_low = 47.3
                            lat_high = 47.9
                            lon_low = -123.3
                            lon_high = -122.6

                        # y = df_dict[gtx][vn].to_numpy()
                        gtx = gtx_list[0]
                        gtx_df = df_dict[gtx]
                        x = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
                                    (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), vn].to_numpy()
                        y = gtx_df.loc[(gtx_df['lat']>lat_low) & (gtx_df['lat']<lat_high) &
                                        (gtx_df['lon']>lon_low) & (gtx_df['lon']<lon_high), vn].to_numpy()
                        ax.plot(x,y,marker='.',ls='',color=c_dict[label], alpha=alpha)
        
                        if (not np.isnan(x).all()) and (not np.isnan(y).all()) and (len(x) > 0) and (len(y) > 0):
                            bias = np.nanmean(y-x)
                            rmse = np.sqrt(np.nanmean((y-x)**2))
                            ax.text(.9,t_dict[label],'bias=%0.1f, rmse=%0.1f' % (bias,rmse),c=c_dict[label],
                                transform=ax.transAxes, ha='right', fontweight='bold', bbox=pfun.bbox,
                                fontsize=fs_label,style='italic')

                    if otype == 'bottle':
                        if jj in [9,10,11]:
                            ax.set_xlabel('Observed', fontsize = fs_header)
                        if jj in [1,5,9]:
                            ax.set_ylabel('Modeled', fontsize = fs_header)
                    elif otype == 'ctd':
                        if jj in [4,5]:
                            ax.set_xlabel('Observed', fontsize = fs_header)
                        if jj in [1,4]:
                            ax.set_ylabel('Modeled', fontsize = fs_header)
                    plt.xticks(fontsize=fs_label)
                    plt.yticks(fontsize=fs_label)
        
                    # add labels to identify the model runs with the colors
                    if jj == 1:
                        yy = 0
                        for label in c_dict.keys():
                            ax.text(.05, .85 - 0.08*yy, label, c=c_dict[label], transform=ax.transAxes,
                                fontweight='bold', ha='left', fontsize = fs_label)
                            yy += 1
            
                    if vn == 'SA':
                        var = 'Salinity'
                    elif vn =='CT':
                        var = 'Cons. Temp'
                    else:
                        var = vn
                    ax.text(.05,.93,var,transform=ax.transAxes, fontweight='bold', fontsize = fs_header)
                    ax.axis([lim_dict[vn][0], lim_dict[vn][1], lim_dict[vn][0], lim_dict[vn][1]])
                    ax.plot([lim_dict[vn][0], lim_dict[vn][1]], [lim_dict[vn][0], lim_dict[vn][1]],'-k')
                    ax.grid(True,color='w',linewidth=2)

                    # format figure
                    # ax.set_box_aspect(1)
                    plt.gca().set_aspect('equal', adjustable='box')
                    ax.set_facecolor('#EEEEEE')
                    for border in ['top','right','bottom','left']:
                        ax.spines[border].set_visible(False)
    
                # station map
                if otype == 'bottle':
                    ax = fig.add_subplot(1,4,4)
                elif otype == 'ctd':
                    ax = fig.add_subplot(1,3,3)

                # df_dict['obs'].plot(x='lon',y='lat',style='.', color='deeppink',legend=False, ax=ax)

                # Plot stations in different regions
                for label in label_list:
                        if label == 'Strait of Georgia':
                            lat_low = 48.8
                            lat_high = 50.16
                            lon_low = -124.9
                            lon_high = -122.5

                        elif label == 'Main Basin':
                            lat_low = 47.3
                            lat_high = 48
                            lon_low = -122.6
                            lon_high = -122

                        elif label == 'South Sound':
                            lat_low = 47
                            lat_high = 47.3
                            lon_low = -123.3
                            lon_high = -122

                        elif label == 'Hood Canal':
                            lat_low = 47.3
                            lat_high = 47.9
                            lon_low = -123.3
                            lon_high = -122.6

                        obs_df = df_dict['obs']
                        x_obs = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
                                        (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), 'lon'].to_numpy()
                        y_obs = obs_df.loc[(obs_df['lat']>lat_low) & (obs_df['lat']<lat_high) &
                                        (obs_df['lon']>lon_low) & (obs_df['lon']<lon_high), 'lat'].to_numpy()
                        ax.plot(x_obs,y_obs,marker='.',ls='',color=c_dict[label], alpha=alpha)

                # # Salish Sea
                # lat_low = 46.9
                # lat_high = 50.16
                # lon_low = -124.5
                # lon_high = -122.1

                # Puget Sound only
                lat_low = 46.9
                lat_high = 48.2
                lon_low = -123.3
                lon_high = -122.1

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

                fig.tight_layout()
                
                print('Plotting ' + ff_str)
                sys.stdout.flush()
                
                if testing:
                    plt.show()
                else:
                    plt.savefig(out_dir / (ff_str + '_diff_locs' '.png'))
                    # plt.close('all')

    
