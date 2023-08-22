"""
Plot timeseries for all WWTPs
So I can identify open and close dates for each of them

To create individual climatology figures, run from ipython with:
run wwtp_flow_timeseries.py -test True
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import xarray as xr
import matplotlib.pyplot as plt
import argparse
import matplotlib.dates as mdates
    

#################################################################################
#                     Get data and set up dataframes                            #
#################################################################################

# read arguments
parser = argparse.ArgumentParser()
# -test True will output plots
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

# location to save file
clim_dir = Ldir['LOo'] / 'pre' / 'traps' / 'point_sources' /'Data_historical'
Lfun.make_dir(clim_dir)

# get flow and loading data
wwtp_fn = Ldir['data'] / 'traps' / 'all_point_source_data.nc'
ecology_data_ds = xr.open_dataset(wwtp_fn)

# get wwtp names and wwtp ids
wwtpnames = ecology_data_ds['name'].values


#################################################################################
#              Plot all individual point source climatologies                   #
#################################################################################

if args.testing == True:
    print('Plotting...')

    # generate directory to save files
    fig_dir = clim_dir / 'wwtp_flow_timeseries'
    Lfun.make_dir(fig_dir)

    for i,wname in enumerate(wwtpnames):

        print('{}/{}: {}'.format(i+1,len(wwtpnames),wname))
        dates = ecology_data_ds.date.values
        flows = ecology_data_ds.flow[ecology_data_ds.name == wname].values[0]

        # Plot flow timeseries for each source
        fig,ax = plt.subplots(figsize=(10,5))
        ax.plot(dates,flows)
        # format figure
        ax.grid(True)
        ax.xaxis.set_major_locator(mdates.YearLocator(base=1))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        plt.xticks(rotation=45)
        ax.set_xlim([ecology_data_ds.date.values[0],ecology_data_ds.date.values[-1]])
        ax.set_ylim([0,1.3*max(flows)])
        ax.set_ylabel('Flow [m3/s]')

        # plot title is name of source
        plt.suptitle(wname,fontsize=18)

        # Save figure
        figname = wname + '.png'
        save_path = fig_dir / figname
        fig.savefig(save_path)
        plt.close('all')

    print('Done')
    