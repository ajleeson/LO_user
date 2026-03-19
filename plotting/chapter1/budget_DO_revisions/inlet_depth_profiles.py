"""
Generate salinity depth profiles
to determine interface depth of 13 inlets
"""

from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from scipy.stats import pearsonr
from datetime import datetime, timedelta
import pandas as pd
import cmocean
import matplotlib.pylab as plt
import pinfo
import seawater as sw
import get_two_layer

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

##########################################################
##                    Define inputs                     ##
##########################################################

gtagex = 'cas7_t1_x11ab'
jobname = 'twentyoneinlets'
startdate = '2017.01.01'
enddate = '2017.12.31'

# season
start = '06-01'
end = '08-31'
year = '2017'

##########################################################
##              Get stations and gtagexes               ##
##########################################################

# parse gtagex
gridname, tag, ex_name = gtagex.split('_')
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)

# find job lists from the extract moor
job_lists = Lfun.module_from_file('job_lists', Ldir['LOu'] / 'extract' / 'moor' / 'job_lists.py')

# Get mooring stations:
sta_dict = job_lists.get_sta_dict(jobname)
del sta_dict['lynchcove2']

# # thirteen inlets in different basins
# inlets = 

# inlets sorted by depth
inlets = ['oak','henderson','hammersley','similk','totten','killsut','budd',
          'eld','sinclair','quartermaster','lynchcove','penn','crescent','dyes',
          'case','holmes','elliot','carr','portsusan','commencement','dabob']

# remove inlets with mean depth < 10 m:
inlets = ['sinclair','quartermaster','lynchcove','penn','crescent','dyes',
          'case','holmes','elliot','carr','portsusan','commencement','dabob']

# ##########################################################
# ##                Plot density profiles                 ##
# ##########################################################

# # fig,axes = plt.subplots(3,7, figsize=(15,8.5), sharex=True)
# fig,axes = plt.subplots(3,5, figsize=(15,8.5))
# ax = axes.ravel()
    
# # plot density profiles
# for i,station in enumerate(inlets): #enumerate(['elliot']):
#     # download .nc files
#     fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
#     ds = xr.open_dataset(fn)
#     # # crop to season
#     # ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

#     # # print(station + ' = {}'.format(ds['h'].values))

#     # loop through and plot each day of the year
#     for day in range(len(ds.ocean_time)):
#         depth = ds['z_rho'][day,:]
#         salt = ds['salt'][day,:]
#         temp = ds['temp'][day,:]
#         rho = sw.dens0(salt, temp) # potential density
#         ax[i].plot(rho,depth,alpha=0.1,color='silver')

#     # # calculate average depth profile and standard deviation
#     # avg_depth = np.nanmean(ds['z_rho'],axis=0)
#     # avg_salt = np.nanmean(ds['salt'],axis=0)
#     # avg_temp = np.nanmean(ds['temp'],axis=0)
#     # avg_rho = sw.dens0(avg_salt, avg_temp)

#     # ax[i].plot(avg_rho,avg_depth,linewidth=2,color='navy')

#     # crop to season
#     seasons = {'JFM':['01-01','03-30'],
#                'AMJ':['04-01','06-30'],
#                'JAS':['07-01','09-30'],
#                'OND':['10-01','12-31']}
#     seasons = {'JAS':['07-01','09-30']}
#     colors = ['black']#['royalblue','hotpink','black','darkorange']
#     for s,season in enumerate(seasons):
#         ds_season = ds.sel(ocean_time=slice(np.datetime64(year+'-'+seasons[season][0]),
#                                      np.datetime64(year+'-'+seasons[season][1])))
#         # calculate average depth profile and standard deviation
#         avg_depth = np.nanmean(ds_season['z_rho'],axis=0)
#         avg_salt = np.nanmean(ds_season['salt'],axis=0)
#         avg_temp = np.nanmean(ds_season['temp'],axis=0)
#         avg_rho = sw.dens0(avg_salt, avg_temp)
#         if season == 'JAS':
#             alpha=1
#             linewidth=3
#             # calculate and plot interface location (inflection point)
#             # determine drho/dz
#             drho_dz = np.diff(avg_rho)/np.diff(avg_depth)
#             minindx = drho_dz.tolist().index(min(drho_dz.tolist()))
#             interface_depth = np.nanmean(avg_depth[minindx:minindx+2])
#         else:
#             alpha=0.5
#             linewidth=2
#         ax[i].plot(avg_rho,avg_depth,linewidth=linewidth,color=colors[s],alpha=alpha)
#         ax[i].axhline(interface_depth,0,1030,color='crimson',linewidth=1)
#         ax[i].text(np.nanmedian(avg_rho) - 1.3*(np.nanmax(avg_rho)-np.nanmin(avg_rho)),
#                    interface_depth+0.04*np.nanmin(avg_depth),
#                    str(round(interface_depth,1)) + ' m',color='crimson',ha='left',va='top',
#                    fontweight='bold')

#     ax[i].set_title(station,fontweight='bold')
#     # ax[i].axhline(-5,0,1030,color='crimson',linewidth=1)

#     ax[0].text(1018,-11,'Daily averages',color='silver',ha='left')
#     ax[0].text(1018,-12.1,'Jul-Sep average',color='black',fontweight='bold',ha='left')

#     if i in [0,5,10]:
#         ax[i].set_ylabel('z (m)')
#     if i >= 10:
#         ax[i].set_xlabel('potential density [kg/m3]')

#     ax[i].axhline(-5,0,1030,color='royalblue',linewidth=1,linestyle='--')


# plt.tight_layout()
# plt.show()


##########################################################
##             Plot density & DO profiles               ##
##########################################################

# fig,axes = plt.subplots(3,7, figsize=(15,8.5), sharex=True)
fig,axes = plt.subplots(3,5, figsize=(15,8.5))
ax = axes.ravel()
    
# plot density profiles
for i,station in enumerate(inlets): #enumerate(['elliot']):
    # download .nc files
    fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)
    # # crop to season
    # ds = ds.sel(ocean_time=slice(np.datetime64(year+'-'+start),np.datetime64(year+'-'+end)))

    # # print(station + ' = {}'.format(ds['h'].values))

    # loop through and plot each day of the year
    for day in range(len(ds.ocean_time)):
        depth = ds['z_rho'][day,:]
        salt = ds['salt'][day,:]
        temp = ds['temp'][day,:]
        rho = sw.dens0(salt, temp) # potential density
        ax[i].plot(salt,depth,alpha=0.1,color='silver')

    # # calculate average depth profile and standard deviation
    # avg_depth = np.nanmean(ds['z_rho'],axis=0)
    # avg_salt = np.nanmean(ds['salt'],axis=0)
    # avg_temp = np.nanmean(ds['temp'],axis=0)
    # avg_rho = sw.dens0(avg_salt, avg_temp)

    # ax[i].plot(avg_rho,avg_depth,linewidth=2,color='navy')

    # crop to season
    seasons = {'JFM':['01-01','03-30'],
               'AMJ':['04-01','06-30'],
               'JAS':['07-01','09-30'],
               'OND':['10-01','12-31']}
    seasons = {'JAS':['07-01','09-30']}
    colors = ['black']#['royalblue','hotpink','black','darkorange']
    for s,season in enumerate(seasons):
        ds_season = ds.sel(ocean_time=slice(np.datetime64(year+'-'+seasons[season][0]),
                                     np.datetime64(year+'-'+seasons[season][1])))
        # calculate average depth profile and standard deviation
        avg_depth = np.nanmean(ds_season['z_rho'],axis=0)
        avg_salt = np.nanmean(ds_season['salt'],axis=0)
        avg_temp = np.nanmean(ds_season['temp'],axis=0)
        avg_rho = sw.dens0(avg_salt, avg_temp)
        if season == 'JAS':
            alpha=1
            linewidth=3
            # # calculate and plot interface location (inflection point)
            # # determine drho/dz
            # drho_dz = np.diff(avg_rho)/np.diff(avg_depth)
            # minindx = drho_dz.tolist().index(min(drho_dz.tolist()))
            # interface_depth = np.nanmean(avg_depth[minindx:minindx+2])

            # calculate and plot interface location (inflection point)
            # determine dsalt/dz
            dsalt_dz = np.diff(avg_salt)/np.diff(avg_depth)
            minindx = dsalt_dz.tolist().index(min(dsalt_dz.tolist()))
            interface_depth = np.nanmean(avg_depth[minindx:minindx+2])

            # Add DO
            ax2 = ax[i].twiny()
            avg_DO = np.nanmean(ds_season['oxygen'],axis=0)* 32/1000 # [mg/L]
            ax2.plot(avg_DO,avg_depth,linewidth=1,color='mediumturquoise',alpha=alpha)
            ax2.tick_params(axis='x', colors='mediumturquoise')

        else:
            alpha=0.5
            linewidth=2
        ax[i].plot(avg_salt,avg_depth,linewidth=linewidth,color=colors[s],alpha=alpha)
        ax[i].axhline(interface_depth,0,40,color='crimson',linewidth=1)
        # label calculated interface depths
        ax[i].text(0.3,0.2,r'max d$\rho$/dz = '+ str(round(interface_depth,1)) + ' m',
                   color='crimson',ha='left',transform=ax[i].transAxes, fontweight='bold')
        
        # get 1/3 depth
        # create dictionaries with interface depths
        interface_dict = dict()
        with open('interface_depths.csv', 'r') as f:
            for line in f:
                inlet, interface_depth = line.strip().split(',')
                interface_dict[inlet] = interface_depth # in meters. NaN means that it is one-layer
        z_interface = float(interface_dict[station])
        ax[i].text(0.25,0.1,r'1/3 depth = '+ str(z_interface) + ' m', fontweight='bold',
                   color='royalblue',ha='left',transform=ax[i].transAxes)
        ax[i].axhline(z_interface,0,40,color='royalblue',linewidth=1,linestyle='--')

    ax[i].set_title(station,fontweight='bold')
    # ax[i].axhline(-5,0,1030,color='crimson',linewidth=1)

    ax[0].text(24,-2,'Daily averages',color='grey',ha='left',va='top')
    ax[0].text(24,-3.1,'Jul-Sep average',color='black',fontweight='bold',ha='left',va='top')
    ax[0].text(24,-4.2,'Jul-Sep avg DO [mg/L]',color='teal',fontweight='bold',ha='left',va='top')

    if i in [0,5,10]:
        ax[i].set_ylabel('z (m)')
    if i >= 10:
        ax[i].set_xlabel('Salinity [g/kg]')

    


plt.tight_layout()
plt.show()


##########################################################
##              Explore interface depths                ##
##########################################################

# get dividing salinity based on TEF

TEF_salt_dict = {}

for station in inlets:

    # --------------------------- get TEF exchange flow terms ----------------------------------------
    in_dir = Ldir['LOo'] / 'extract' / 'cas7_t1_x11b' / 'tef2' / 'c21' / ('bulk_'+year+'.01.01_'+year+'.12.31') / (station + '.nc')
    bulk = xr.open_dataset(in_dir)
    tef_df, vn_list, vec_list = get_two_layer.get_two_layer(bulk)
    Q_p = tef_df['q_p'] # Qin [m3/s]
    Q_m = tef_df['q_m'] # Qout [m3/s]
    S_p = tef_df['salt_p'] # Sin [g/kg]
    S_m = tef_df['salt_m'] # Sout [g/kg]

    # calculate annual average dividing salinity
    TEF_div_salt = np.nanmean([np.nanmean(S_m),np.nanmean(S_p)])

    # add dividing salinity to dictionary 
    TEF_salt_dict[station] = TEF_div_salt

    # print('===========\n'+station)
    # print('Sin = {}'.format(np.nanmean(S_p)))
    # print('Sout = {}'.format(np.nanmean(S_m)))
    # print('midpoint = {}'.format(np.nanmean([np.nanmean(S_m),np.nanmean(S_p)])))


##########################################################
##             Plot density & DO profiles               ##
##########################################################

# fig,axes = plt.subplots(3,7, figsize=(15,8.5), sharex=True)
fig,axes = plt.subplots(3,5, figsize=(15,8.5))
ax = axes.ravel()
    
# plot density profiles
for i,station in enumerate(inlets): #enumerate(['elliot']):
    # download .nc files
    fn = '../../../../LO_output/extract/' + gtagex + '/moor/' + jobname + '/' + station + '_' + startdate + '_' + enddate + '.nc'
    ds = xr.open_dataset(fn)

    # loop through and plot each day of the year
    for day in range(len(ds.ocean_time)):
        depth = ds['z_rho'][day,:]
        salt = ds['salt'][day,:]
        temp = ds['temp'][day,:]
        rho = sw.dens0(salt, temp) # potential density
        ax[i].plot(salt,depth,alpha=0.1,color='silver')

    # crop to season
    seasons = {'JFM':['01-01','03-30'],
               'AMJ':['04-01','06-30'],
               'JAS':['07-01','09-30'],
               'OND':['10-01','12-31']}
    seasons = {'JAS':['07-01','09-30']}
    colors = ['black']#['royalblue','hotpink','black','darkorange']
    for s,season in enumerate(seasons):
        ds_season = ds.sel(ocean_time=slice(np.datetime64(year+'-'+seasons[season][0]),
                                     np.datetime64(year+'-'+seasons[season][1])))
        # calculate average depth profile and standard deviation
        avg_depth = np.nanmean(ds_season['z_rho'],axis=0)
        avg_salt = np.nanmean(ds_season['salt'],axis=0)
        avg_temp = np.nanmean(ds_season['temp'],axis=0)
        avg_rho = sw.dens0(avg_salt, avg_temp)
        if season == 'JAS':
            alpha=1
            linewidth=3

            # Add DO
            ax2 = ax[i].twiny()
            avg_DO = np.nanmean(ds_season['oxygen'],axis=0)* 32/1000 # [mg/L]
            ax2.plot(avg_DO,avg_depth,linewidth=1,color='teal',alpha=alpha)
            ax2.tick_params(axis='x', colors='teal')

        else:
            alpha=0.5
            linewidth=2
        ax[i].plot(avg_salt,avg_depth,linewidth=linewidth,color=colors[s],alpha=alpha)


    ax[i].set_title(station,fontweight='bold')

    # add dividing salinity based on TEF
    interface_salt = TEF_salt_dict[station]
    index_of_interface = min(range(len(avg_salt)), key=lambda i: abs(avg_salt[i]-interface_salt))
    interface_depth = avg_depth[index_of_interface]
    ax[i].axhline(interface_depth,0,40,color='crimson',linewidth=1)
    ax[i].text(0.05,0.5,r'from TEF = '+ str(np.round(interface_depth,1)) + ' m', fontweight='bold',
                   color='crimson',ha='left',transform=ax[i].transAxes)

    # label lines
    ax[0].text(23.2,-9,'Daily averages',color='grey',ha='left',va='top')
    ax[0].text(23.2,-10.1,'Jul-Sep average',color='black',fontweight='bold',ha='left',va='top')
    ax[0].text(23.2,-11.2,'Jul-Sep avg DO [mg/L]',color='teal',fontweight='bold',ha='left',va='top')

    if i in [0,5,10]:
        ax[i].set_ylabel('z (m)')
    if i >= 10:
        ax[i].set_xlabel('Salinity [g/kg]')

    


plt.tight_layout()
plt.show()
