"""
Make climatologies for pre-existing LO rivers.
Biogeochemisty variables only.

Based on Ecology's timeseries, using data stored in 
LO_data/traps/all_nonpoint_source_data.nc

To run, from ipython:
run make_climatology_LOrivbio.py
"""

plotting = True

#################################################################################
#                              Import packages                                  #
#################################################################################

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import datetime
import matplotlib.dates as mdates
import datetime
    

#################################################################################
#                     Get data and set up dataframes                            #
#################################################################################

# define year range to create climatologies
year0 = 1999
year1 = 2017

# location to save file
clim_dir = Ldir['LOo'] / 'pre' / 'traps' / 'LO_rivbio' /'Data_historical'

# get flow and loading data
triv_fn = Ldir['data'] / 'traps' / 'all_nonpoint_source_data.nc'
ecology_data_ds = xr.open_dataset(triv_fn)

# get triv names
trivnames_all = ecology_data_ds['name'].values

# Remove pre-existing LO rivers
# read overlapping rivers
repeatrivs_fn = Ldir['data'] / 'traps' / 'LiveOcean_SSM_rivers.xlsx'
repeatrivs_df = pd.read_excel(repeatrivs_fn)
SSM_repeats = repeatrivs_df['SSM_rname'].values
# remove nans
SSM_repeats = [x for x in SSM_repeats if str(x) != 'nan']
# remove repeat river names from list of river names
trivnames = [river for river in trivnames_all if river in SSM_repeats]

# # just test 5 trivs for now -------------------------------------------------
# trivnames = trivnames[28:33]

# separate rivers with weird biogeochem
weird_biogeochem = ['Alberni Inlet', 'Chehalis R', 'Gold River',
                    'Willapa R', 'Columbia R', 'Comox']
# get list of rivers with more realistic biogeochem
trivnames_realisticBGC = [river for river in trivnames
                          if river not in weird_biogeochem]

# initialize dataframes for all trivs
NO3_clim_df  = pd.DataFrame()
NH4_clim_df  = pd.DataFrame()
TIC_clim_df  = pd.DataFrame()
Talk_clim_df = pd.DataFrame()
DO_clim_df   = pd.DataFrame()

# variable names
vns = ['DO(mg/L)','NO3(mmol/m3)', 'NH4(mmol/m3)','TIC(mmol/m3)','Talk(meq/m3)']
letters = ['(a)','(b)','(c)','(d)','(e)']

# create one-year date range for plotting
yrday = pd.date_range(start ='1/1/2020', end ='12/31/2020', freq ='D')

#################################################################################
#                          Calculate climatologies                              #
#################################################################################

print('Calculating climatologies...')

# loop through all nonpoint sources
for i,rname in enumerate(trivnames_realisticBGC):

    # turn dataset information for this triv into a dataframe
    # so it's easier to manipulate
    d = {'Date': ecology_data_ds.date.values,
         'NO3(mmol/m3)':ecology_data_ds.NO3[ecology_data_ds.name==rname,:].values[0],
         'NH4(mmol/m3)':ecology_data_ds.NH4[ecology_data_ds.name==rname,:].values[0],
         'TIC(mmol/m3)':ecology_data_ds.TIC[ecology_data_ds.name==rname,:].values[0],
         'Talk(meq/m3)': ecology_data_ds.Talk[ecology_data_ds.name==rname,:].values[0],
         'DO(mmol/m3)': ecology_data_ds.DO[ecology_data_ds.name==rname,:].values[0]}
    triv_df = pd.DataFrame(data=d)
    # replace all zeros with nans, so zeros don't bias data
    triv_df = triv_df.replace(0, np.nan)

    # add day of year column
    triv_df['day_of_year'] = triv_df.apply(lambda row: row.Date.dayofyear, axis = 1)
    # add year column
    triv_df['year'] = pd.DatetimeIndex(triv_df['Date']).year

    # calculate averages
    # (compress 1999-2017 timeseries to single year, with an average for each day)
    triv_avgs_df = triv_df.groupby('day_of_year').mean().reset_index()
    # calculate standard deviation
    triv_sds_df = triv_df.groupby('day_of_year').std(ddof=0).reset_index()

    # replace all nans with zeros, so I'm no longer injecting nans
    triv_avgs_df = triv_avgs_df.replace(np.nan,0)
    triv_sds_df = triv_sds_df.replace(np.nan,0)

#################################################################################
#                            Create climatologies                               #
#################################################################################

    # Add data to climatology dataframes
    NO3_clim_df  = pd.concat([NO3_clim_df, pd.Series(triv_avgs_df['NO3(mmol/m3)'], name=rname)], axis = 1)          # [mmol/m3]
    NH4_clim_df  = pd.concat([NH4_clim_df, pd.Series(triv_avgs_df['NH4(mmol/m3)'], name=rname)], axis = 1)          # [mmol/m3]
    TIC_clim_df  = pd.concat([TIC_clim_df, pd.Series(triv_avgs_df['TIC(mmol/m3)'], name=rname)], axis = 1)          # [mmol/m3]
    Talk_clim_df = pd.concat([Talk_clim_df, pd.Series(triv_avgs_df['Talk(meq/m3)'], name=rname)], axis = 1)         # [meq/m3]
    DO_clim_df   = pd.concat([DO_clim_df, pd.Series(triv_avgs_df['DO(mmol/m3)'], name=rname)], axis = 1)            # [mmol/m3]

# Calculate summary statistics for all trivs
clim_avgs = pd.DataFrame()
clim_max = pd.DataFrame()
clim_min = pd.DataFrame()
clim_sds = pd.DataFrame()
# list climatology dfs
clim_df_list = [DO_clim_df,NO3_clim_df,NH4_clim_df,TIC_clim_df,Talk_clim_df]

# rename variables
vns = ['DO(mmol/m3)','NO3(mmol/m3)','NH4(mmol/m3)','TIC(mmol/m3)','Talk(meq/m3)']
for i,vn in enumerate(vns):
    # average values of all nonpoint sources
    clim_avgs[vn] = clim_df_list[i].mean(axis=1)
    # max climatology values
    clim_max[vn] = clim_df_list[i].max(axis=1)
    # min climatology values
    clim_min[vn] = clim_df_list[i].min(axis=1)
    # standard deviation of all nonpoint sources
    clim_sds[vn] = clim_df_list[i].std(axis=1)


#################################################################################
#                           Plot summary statistics                             #
#################################################################################

# Plot Summary Statistics
fig, axes = plt.subplots(3,2, figsize=(16, 8), sharex=True)
ax = axes.ravel()
for j,vn in enumerate(vns):
    i = j+1

    # convert DO from mmol/m3 to mg/L for plotting
    if vn == 'DO(mmol/m3)':
        scale = 1/31.26 
        var = 'DO(mg/L)'
    else:
        scale = 1
        var = vn

    # label subplot
    ax[i].text(0.05,0.85,letters[j]+' '+var,transform=ax[i].transAxes,fontsize=14)
    # Plot average
    ax[i].plot(yrday,clim_avgs[vn].values*scale, label='Average of all Sources', color='mediumpurple', linewidth=1.5)
    # Plot error shading
    upper_bound = [min(clim_avgs[vn].values[ii]+clim_sds[vn].values[ii],clim_max[vn].values[ii])*scale for ii in range(366)] # don't go higher than max value
    lower_bound = [max(clim_avgs[vn].values[ii]-clim_sds[vn].values[ii],clim_min[vn].values[ii])*scale for ii in range(366)] # don't go lower than min value
    ax[i].fill_between(yrday,upper_bound,lower_bound,label='One SD',color='mediumpurple',alpha=0.2,edgecolor='none')
    # Plot max
    ax[i].plot(yrday,clim_max[vn].values*scale, label='Max Value', color='firebrick', linestyle='--', linewidth=1)
    # Plot min
    ax[i].plot(yrday,clim_min[vn].values*scale, label='Min Value', color='cornflowerblue', linestyle='--', linewidth=1)
    # fontsize of tick labels
    ax[i].tick_params(axis='both', which='major', labelsize=12)
    ax[i].tick_params(axis='x', which='major', rotation=30)
    ax[i].set_xlim([datetime.date(2020, 1, 1), datetime.date(2020, 12, 31)])
    ax[i].set_ylim([0,1.3*max(clim_max[vn].values)])
    if i < 5:
        ax[i].set_ylim([0,1.3*max(clim_max[vn].values)*scale])
    # create legend
    if i ==5:
        ax[i].set_ylim([0,1.3*max(clim_max['TIC(mmol/m3)'].values)])
        handles, labels = ax[5].get_legend_handles_labels()
        ax[0].legend(handles, labels, loc='center', ncol = 2,fontsize=14)
        ax[0].axis('off')
    # Define the date format
    if i >= 4:
        date_form = mdates.DateFormatter("%b")
        ax[i].xaxis.set_major_formatter(date_form)
# plot title is name of source
plt.suptitle('Pre-existing LO River Climatology Summary (n={})'.format(len(trivnames_realisticBGC)),fontsize=18)
# Save figure
figname = 'LOrivbio_river_summary.png'
save_path = clim_dir / figname
fig.savefig(save_path)
# plt.close('all')
# plt.show()

#################################################################################
#                          Plot all river climatologies                         #
#################################################################################

if plotting == True:
    print('Climatologies done\n')
    print('Plotting...')

    for i,rname in enumerate(trivnames_realisticBGC):

        print('{}/{}: {}'.format(i+1,len(trivnames_realisticBGC),rname))

        # turn dataset information for this triv into a dataframe
        # so it's easier to manipulate
        d = {'Date': ecology_data_ds.date.values,
            'NO3(mmol/m3)':ecology_data_ds.NO3[ecology_data_ds.name==rname,:].values[0],
            'NH4(mmol/m3)':ecology_data_ds.NH4[ecology_data_ds.name==rname,:].values[0],
            'TIC(mmol/m3)':ecology_data_ds.TIC[ecology_data_ds.name==rname,:].values[0],
            'Talk(meq/m3)': ecology_data_ds.Talk[ecology_data_ds.name==rname,:].values[0],
            'DO(mmol/m3)': ecology_data_ds.DO[ecology_data_ds.name==rname,:].values[0]}
        triv_df = pd.DataFrame(data=d)
        # replace all zeros with nans, so zeros don't bias data
        triv_df = triv_df.replace(0, np.nan)

        # add day of year column
        triv_df['day_of_year'] = triv_df.apply(lambda row: row.Date.dayofyear, axis = 1)
        # add year column
        triv_df['year'] = pd.DatetimeIndex(triv_df['Date']).year

        # Plotting
        fig, axes = plt.subplots(3,2, figsize=(16, 8), sharex=True)
        ax = axes.ravel()
        for j,var in enumerate(vns):

            # convert DO from mmol/m3 to mg/L for plotting
            if var == 'DO(mmol/m3)':
                scale = 1/31.26 
                vn = 'DO(mg/L)'
            else:
                scale = 1
                vn = var

            i = j+1
            # label subplot
            ax[i].set_title(vn,fontsize=14)
            # Plot individual years
            for yr in range(1999,2017):
                triv_yr_df = triv_df.loc[triv_df['year'] == yr]
                values_to_plot = triv_yr_df[var].values*scale
                values_to_plot = values_to_plot.tolist()
                # skip leap years
                if np.mod(yr,4) != 0:
                    # pad Feb 29th with nan
                    values_to_plot = values_to_plot[0:60] + [np.nan] + values_to_plot[60::]
                if yr == 2017:
                    yrday_17 = pd.date_range(start ='1/1/2020', end ='8/02/2020', freq ='D') # don't have full 2017 dataset
                    ax[i].plot(yrday_17,values_to_plot,alpha=0.5, label=yr, linewidth=1)
                else:
                    ax[i].plot(yrday,values_to_plot,alpha=0.5, label=yr, linewidth=1)
            # Plot climatology
            ax[i].plot(yrday,clim_df_list[j][rname].values*scale, label='climatology', color='black', linewidth=1.5)
            # fontsize of tick labels
            ax[i].tick_params(axis='both', which='major', labelsize=12)
            ax[i].tick_params(axis='x', which='major', rotation=30)
            ax[i].set_xlim([datetime.date(2020, 1, 1), datetime.date(2020, 12, 31)])
            # create legend
            if i ==5:
                handles, labels = ax[5].get_legend_handles_labels()
                ax[0].legend(handles, labels, loc='center', ncol = 4,fontsize=14)
                ax[0].axis('off')
            # Define the date format
            if i >= 4:
                date_form = mdates.DateFormatter("%b")
                ax[i].xaxis.set_major_formatter(date_form)
        # plot title is name of source
        plt.suptitle(rname,fontsize=18)
        # Save figure
        figname = rname + '.png'
        save_path = clim_dir / 'climatology_plots' / figname
        fig.savefig(save_path)
        plt.close('all')


#################################################################################
#                             Save climatologies                                #
#################################################################################

# check for missing values:
if pd.isnull(NO3_clim_df).sum().sum() != 0:
    print('Warning, there are missing nitrate values!')
if pd.isnull(NH4_clim_df).sum().sum() != 0:
    print('Warning, there are missing ammonium values!')
if pd.isnull(TIC_clim_df).sum().sum() != 0:
    print('Warning, there are missing TIC values!')
if pd.isnull(Talk_clim_df).sum().sum() != 0:
    print('Warning, there are missing alkalinity values!')
if pd.isnull(DO_clim_df).sum().sum() != 0:
    print('Warning, there are missing oxygen values!')

# save results
NO3_clim_df.to_pickle(clim_dir / ('CLIM_NO3_' + str(year0) + '_' + str(year1) + '.p'))
NH4_clim_df.to_pickle(clim_dir / ('CLIM_NH4_' + str(year0) + '_' + str(year1) + '.p'))
TIC_clim_df.to_pickle(clim_dir / ('CLIM_TIC_' + str(year0) + '_' + str(year1) + '.p'))
Talk_clim_df.to_pickle(clim_dir / ('CLIM_Talk_' + str(year0) + '_' + str(year1) + '.p'))
DO_clim_df.to_pickle(clim_dir / ('CLIM_DO_' + str(year0) + '_' + str(year1) + '.p'))

print('Done')

# # Print statements for testing
# for i,clim in enumerate(clim_df_list):
#     print(vns[i]+'=================================================')
#     print(clim)