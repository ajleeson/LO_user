#%%
"""
Plot fields in one or more history files.

Examples:

Plot a single figure to the screen with default arguments:
run pan_plot
(it will prompt for the plot type)

Here is an example with explicit flags:
run pan_plot -gtx cas6_v0_live -ro 0 -0 2019.07.04 -lt snapshot -pt P_Chl_DO -avl False

When using the default -avl True for multiple plots (e.g. when making a movie)
the color limits will all be set to match those set by auto_lims() from the first plot.

Use -avl False to have color limits all set to match those set by pinfo.vlims_dict.

Using -test True will create the object "ds", an xarray Dataset of the most recent history file.

"""
import os, sys
import argparse
from datetime import datetime, timedelta
from subprocess import Popen as Po
from subprocess import PIPE as Pi

from lo_tools import Lfun
import roms_plots
from importlib import reload
reload(roms_plots)

########################
# helper functions

# format used for naming day folders
ds_fmt = '%Y.%m.%d'

def get_fn_list(list_type, Ldir, ds0, ds1, his_num=2):
    """
    INPUT:
    A function for getting lists of history files.
    List items are Path objects
    
    NEW 2023.10.05: for list_type = 'hourly', if you pass his_num = 1
    it will start with ocean_his_0001.nc on the first day instead of the default which
    is to start with ocean_his_0025.nc on the day before.

    NEW 2025.06.20: for list_type = 'hourly0'
    which will start with ocean_his_0001.nc on the first day instead of the default which
    is to start with ocean_his_0025.nc on the day before.
    This is identical to passing his_num = 1, but may be more convenient, especially
    as we move to "continuation" start_type, which always writes an 0001 file.
    """
    dt0 = datetime.strptime(ds0, ds_fmt)
    dt1 = datetime.strptime(ds1, ds_fmt)
    dir0 = Ldir['roms_out'] / Ldir['gtagex']
    if list_type == 'snapshot':
        # a single file name in a list
        his_string = ('0000' + str(his_num))[-4:]
        fn_list = [dir0 / ('f' + ds0) / ('ocean_his_' + his_string + '.nc')]
    elif list_type == 'hourly':
        # list of hourly files over a date range
        fn_list = fn_list_utility(dt0,dt1,Ldir,his_num=his_num)
    elif list_type == 'hourly0':
        # list of hourly files over a date range, starting with 0001 of dt0.
        fn_list = fn_list_utility(dt0,dt1,Ldir,his_num=1)
    elif list_type == 'daily':
        # list of history file 21 (Noon PST) over a date range
        fn_list = []
        date_list = date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_his_0021.nc'
            fn_list.append(fn)
    elif list_type == 'lowpass':
        # list of lowpassed files (Noon PST) over a date range
        fn_list = []
        date_list = date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'lowpassed.nc'
            fn_list.append(fn)
    elif list_type == 'average':
        # list of daily averaged files (Noon PST) over a date range
        fn_list = []
        date_list = date_list_utility(dt0, dt1)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_avg_0001.nc'
            fn_list.append(fn)
    elif list_type == 'weekly':
        # like "daily" but at 7-day intervals
        fn_list = []
        date_list = date_list_utility(dt0, dt1, daystep=7)
        for dl in date_list:
            f_string = 'f' + dl
            fn = dir0 / f_string / 'ocean_his_0001.nc'
            fn_list.append(fn)
    elif list_type == 'allhours':
        # a list of all the history files in a directory
        # (this is the only list_type that actually finds files)
        in_dir = dir0 / ('f' + ds0)
        fn_list = [ff for ff in in_dir.glob('ocean_his*nc')]
        fn_list.sort()

    return fn_list

def date_list_utility(dt0, dt1, daystep=1):
    """
    INPUT: start and end datetimes
    OUTPUT: list of LiveOcean formatted dates
    """
    date_list = []
    dt = dt0
    while dt <= dt1:
        date_list.append(dt.strftime(ds_fmt))
        dt = dt + timedelta(days=daystep)
    return date_list

def fn_list_utility(dt0, dt1, Ldir, hourmax=24, his_num=2):
    """
    INPUT: start and end datetimes
    OUTPUT: list of all history files expected to span the dates
    - list items are Path objects
    """
    dir0 = Ldir['roms_out'] / Ldir['gtagex']
    fn_list = []
    date_list = date_list_utility(dt0, dt1)
    if his_num == 1:
        # New scheme 2023.10.05 to work with new or continuation start_type,
        # by assuming we want to start with ocean_his_0001.nc of dt0
        fn_list.append(dir0 / ('f'+dt0.strftime(ds_fmt)) / 'ocean_his_0001.nc')
        for dl in date_list:
            f_string = 'f' + dl
            hourmin = 1
            for nhis in range(hourmin+1, hourmax+2):
                nhiss = ('0000' + str(nhis))[-4:]
                fn = dir0 / f_string / ('ocean_his_' + nhiss + '.nc')
                fn_list.append(fn)
    else:
        # For any other value of his_num we assume this is a perfect start_type
        # and so there is no ocean_his_0001.nc on any day and we start with
        # ocean_his_0025.nc of the day before.
        dt00 = (dt0 - timedelta(days=1))
        fn_list.append(dir0 / ('f'+dt00.strftime(ds_fmt)) / 'ocean_his_0025.nc')
        for dl in date_list:
            f_string = 'f' + dl
            hourmin = 1
            for nhis in range(hourmin+1, hourmax+2):
                nhiss = ('0000' + str(nhis))[-4:]
                fn = dir0 / f_string / ('ocean_his_' + nhiss + '.nc')
                fn_list.append(fn)
    return fn_list

###########################



#%%

parser = argparse.ArgumentParser()

# which run to use
parser.add_argument('-gtx', '--gtagex', default='cas6_v0_live', type=str)
parser.add_argument('-ro', '--roms_out_num', default=0, type=int)
# 2 = Ldir['roms_out2'], etc.

# select time period and frequency
parser.add_argument('-0', '--ds0', default='2019.07.04', type=str)
parser.add_argument('-1', '--ds1', type=str)
parser.add_argument('-lt', '--list_type', default='snapshot', type=str)
# snapshot, hourly, or daily

# arguments that allow you to bypass the interactive choices
parser.add_argument('-hn', '--his_num', default=1, type=int)
parser.add_argument('-pt', '--plot_type', type=str)

# arguments that influence other behavior
#  e.g. make a movie, override auto color limits
parser.add_argument('-mov', '--make_movie', default=False, type=Lfun.boolean_string)
parser.add_argument('-avl', '--auto_vlims', default=True, type=Lfun.boolean_string)
parser.add_argument('-save', '--save_plot', default=False, type=Lfun.boolean_string)
parser.add_argument('-test', '--testing', default=False, type=Lfun.boolean_string)

# do things with the arguments
args = parser.parse_args()
argsd = args.__dict__
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# set second date string if omitted (need to have Ldir['ds0'])
if Ldir['ds1'] == None:
    Ldir['ds1'] = Ldir['ds0']
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num']) + '_blowup']
    print(Ldir['roms_out'])

# choose the type of list to make
if Ldir['list_type'] == None:
    print(' pan_plot '.center(60,'='))
    print('\n%s\n' % '** Choose List type (return for snapshot) **')
    lt_list = ['snapshot', 'daily', 'hourly ', 'allhours']
    Nlt = len(lt_list)
    lt_dict = dict(zip(range(Nlt), lt_list))
    for nlt in range(Nlt):
        print(str(nlt) + ': ' + lt_list[nlt])
    my_nlt = input('-- Input number -- ')
    if len(my_nlt)==0:
        Ldir['list_type'] = 'snapshot'
    else:
        Ldir['list_type'] = lt_dict[int(my_nlt)]
else:
    pass

# choose the type of plot to make
if Ldir['plot_type'] == None:
    print('\n%s\n' % '** Choose Plot type (return for P_basic) **')
    pt_list_raw = dir(roms_plots)
    pt_list = []
    for pt in pt_list_raw:
        if pt[:2] == 'P_':
            pt_list.append(pt)
    Npt = len(pt_list)
    pt_dict = dict(zip(range(Npt), pt_list))
    for npt in range(Npt):
        print(str(npt) + ': ' + pt_list[npt])
    my_npt = input('-- Input number -- ')
    if len(my_npt)==0:
        Ldir['plot_type'] = 'P_basic'
    else:
        Ldir['plot_type'] = pt_dict[int(my_npt)]
else:
    pass
    
whichplot = getattr(roms_plots, Ldir['plot_type'])

in_dict = dict()
in_dict['auto_vlims'] = Ldir['auto_vlims']
in_dict['testing'] = Ldir['testing']

# get list of history files to plot
fn_list = Lfun.get_fn_list(Ldir['list_type'], Ldir,
    Ldir['ds0'], Ldir['ds1'], his_num=Ldir['his_num'])
    
if (Ldir['list_type'] == 'allhours') and Ldir['testing']:
    fn_list = fn_list[:4]

# PLOTTING
outdir0 = Ldir['LOo'] / 'plots'
Lfun.make_dir(outdir0)
if '_mac' in Ldir['lo_env']: # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
plt.close('all')

# #---------------------------------------------------------------------
# fn_list = fn_list[0:-5]

if len(fn_list) == 1:
    # plot a single image to screen
    fn = fn_list[0]
    in_dict['fn'] = fn
    if Ldir['save_plot'] == True:
        in_dict['fn_out'] = outdir0 / (Ldir['list_type'] + '_'
            + Ldir['plot_type'] + '_' + Ldir['gtagex'] + '.png')
    else:
        in_dict['fn_out'] = ''
    whichplot(in_dict)
    
elif len(fn_list) > 1:
    # fn_list = fn_list[0:20] # I previously used this line to shorted hourly movie to only 20 hours, due to blow up
    # prepare a directory for results
    outdir = outdir0 / (Ldir['list_type'] + '_' + Ldir['plot_type'] + '_' + Ldir['gtagex'])
    Lfun.make_dir(outdir, clean=True)
    # plot to a folder of files
    jj = 0
    for fn in fn_list[1::]: # changed from fn_list so that we don't try to make video with hour 25 from previous day
        nouts = ('0000' + str(jj))[-4:]
        outname = 'plot_' + nouts + '.png'
        outfile = outdir / outname
        print('Plotting ' + str(fn))
        sys.stdout.flush()
        in_dict['fn'] = fn
        in_dict['fn_out'] = outfile
        whichplot(in_dict)
        # after the first plot we no longer change vlims
        in_dict['auto_vlims'] = False
        jj += 1
    # and make a movie
    if Ldir['make_movie']:
        cmd_list = ['ffmpeg','-r','4','-i', str(outdir)+'/plot_%04d.png', '-vcodec', 'libx264',
            '-pix_fmt', 'yuv420p', '-crf', '21', str(outdir)+'/movie.mp4']
        proc = Po(cmd_list, stdout=Pi, stderr=Pi)
        stdout, stderr = proc.communicate()
        if len(stdout) > 0:
            print('\n'+stdout.decode())
        if len(stderr) > 0:
            print('\n'+stderr.decode())
            
if Ldir['testing']:
    import xarray as xr
    ds = xr.open_dataset(fn)
