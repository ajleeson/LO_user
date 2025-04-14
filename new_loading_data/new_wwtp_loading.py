"""
Script to plot and compare new WWTP data and my older climatologies.
"""

#################################################################################
#                              Import packages                                  #
#################################################################################

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import datetime
import matplotlib.dates as mdates
import datetime
import traps_helper
import os
from pathlib import Path
    

#################################################################################
#                          Get old WWTP climatology data                        #
#################################################################################

# get flow
fn_oldWWTP_flow = Ldir['LOo'] / 'pre' / 'trapsP00' / 'point_sources' / 'lo_base' / 'Data_historical' / 'CLIM_flow.p'
df_oldWWTP_flow = pd.read_pickle(fn_oldWWTP_flow)

print(df_oldWWTP_flow)