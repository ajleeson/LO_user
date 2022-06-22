# %%
from matplotlib import markers
import numpy as np
import pickle
from datetime import datetime, timedelta
import pandas as pd
import xarray as xr
from cmocean import cm
import netCDF4 as nc
import matplotlib.pyplot as plt
import pytz

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'lo_tools'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
# import gfun_utility as gfu
# import gfun

# from lo_tools import Lfun, zfun, zrfun
# from lo_tools import plotting_functions as pfun
# import pinfo
# from importlib import reload
# reload(pfun)
# reload(pinfo)

# Ldir = Lfun.Lstart()
# if '_mac' in Ldir['lo_env']: # mac version
#     pass
# else: # remote linux version
#     import matplotlib as mpl
#     mpl.use('Agg')
# import matplotlib.pyplot as plt


def get_dt_local(dt, tzl='US/Pacific'):
    # take a model datetime (assumed to be UTC) and return local datetime
    tz_utc = pytz.timezone('UTC')
    tz_local = pytz.timezone(tzl)
    dt_utc = dt.replace(tzinfo=tz_utc)
    dt_local = dt_utc.astimezone(tz_local)
    return dt_local

# %%

gtagex = 'alpe_v40d_uu1k'
gridname = 'alpe'
year_str = '2020'
d_str = '2020.02.01_2020.02.10'

loo = Path(__file__).absolute().parent.parent.parent / 'LO_output'

# get forcing fields
ffn = loo / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
fdf = pd.read_pickle(ffn)


# %% SALINITY PROFILE
# mooring
fn = loo / 'extract' / gtagex / 'moor' / gridname / ('superplot_' + d_str + '.nc')
moor = xr.open_dataset(fn)

zeta = moor['zeta']
salt = moor['salt'].transpose()
z_rho = moor['z_rho'].transpose()

year = '2020'
dates = pd.date_range(start='2/1/'+year, end='2/11/'+year, freq= '1H')
dates_local = [get_dt_local(x) for x in dates]

n = 241
fig, ax = plt.subplots(1,1,figsize = (12,4))
cs = ax.pcolormesh(dates_local,z_rho[:,0:n],salt[:,0:n])
fig.colorbar(cs)
plt.ylabel('Z (m)')
plt.title(r'Salinity Depth Profile ($g \ kg^{-1}$)')
#ax.xaxis.set_major_formatter(mdates.DateFormatter("%D"))
ax.tick_params('x', labelrotation=45)

# %% VELOCITY  PROFILE

# mooring
fn = loo / 'extract' / gtagex / 'moor' / gridname / ('superplot_' + d_str + '.nc')
moor = xr.open_dataset(fn)

v = moor['v']
z_w = moor['z_w']
z_rho = moor['z_rho']

year = '2020'
dates = pd.date_range(start='2/1/'+year, end='2/11/'+year, freq= '1H')
dates_local = [get_dt_local(x) for x in dates]

# Calculate first difference of z_w
zw_fd = np.diff(z_w)

# Calculate transport velocity
vdz = v * zw_fd

# time average with filter
# use a godin like the superplot rms calculation

# calculate average depth
avgz = np.mean(z_rho,0)

n = 241
fig, ax = plt.subplots(1,1,figsize = (6,6))
plt.ylabel('Z (m)')
plt.title(r'Velcotiy Profile ($m \ s^{-1}$)')


# %%
