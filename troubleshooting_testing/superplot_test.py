# %%
from matplotlib import markers
import numpy as np
import pickle
from datetime import datetime, timedelta
import pandas as pd
import xarray as xr
from cmocean import cm
import netCDF4 as nc

import sys
from pathlib import Path
# pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'lo_tools'
# if str(pth) not in sys.path:
#     sys.path.append(str(pth))
# import gfun_utility as gfu
# import gfun

# from lo_tools import Lfun, zfun, zrfun
# from lo_tools import plotting_functions as pfun
# import pinfo
# from importlib import reload
# reload(pfun)
# reload(pinfo)

Ldir = Lfun.Lstart()
if '_mac' in Ldir['lo_env']: # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt

# %%

gtagex = 'alpe_v40d_uu1k'
gridname = 'alpe'
year_str = '2020'
d_str = '2020.02.01_2020.02.10'

loo = Path(__file__).absolute().parent.parent.parent / 'LO_output'

# get forcing fields
ffn = loo / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
fdf = pd.read_pickle(ffn)

# mooring
fn = loo / 'extract' / gtagex / 'moor' / gridname / ('superplot_' + d_str + '.nc')
moor = nc.Dataset(fn)