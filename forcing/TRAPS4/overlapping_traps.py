"""
To test from ipython:

run overlapping_traps.py -g cas6 -r backfill -s continuation -d 2021.01.01 -f TRAPS2
"""

from pathlib import Path
import sys, os
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from lo_tools import forcing_argfun2 as ffun

Ldir = ffun.intro() # this handles all the argument passing

import xarray as xr
from lo_tools import Lfun, zrfun
import numpy as np
import pandas as pd
import rivfun
import trapsfun

# get the list of pre-existing rivers and indices for this grid
lo_fn = Ldir['grid'] / 'river_info.csv'
lo_df = pd.read_csv(lo_fn, index_col='rname')

# get the list of tiny rivers and indices for this grid
tr_fn = Ldir['grid'] / 'triv_info.csv'
tr_df = pd.read_csv(tr_fn, index_col='rname')

# get the list of point sources and indices for this grid
ps_fn = Ldir['grid'] / 'wwtp_info.csv'
ps_df = pd.read_csv(ps_fn, index_col='rname')

# merge the dataframes
traps_df = tr_df.append(ps_df)
all_df = traps_df.append(lo_df)

# get list of sources that discharge to the same cell (overlapping)
duplicate_df = all_df[all_df.duplicated(['row_py','col_py'], keep=False) == True]
print(duplicate_df)