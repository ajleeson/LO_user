"""
Calculate upper limit of BOD effect from WWTPs
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import get_two_layer
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from matplotlib import colormaps
from matplotlib.colors import ListedColormap
from scipy.stats import pearsonr
from scipy.linalg import lstsq
from matplotlib.ticker import AutoMinorLocator

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

########################
# set inputs

# max average value from wwtps (five day demand)
bod = 30 # mg/L/5day
bod_1day = bod/5 # mg/L 

# max wwtp flowrate (west point)
Q = 0.1 # 5 # m3/s

########################
# calcualate DO consumption

# Get flowrate in L/s
Q_Ls = Q * 1000 # m3/s * 1000 L/m3

# convert to L/d
Q_Ld = Q_Ls * 60 * 60 * 24 # L/s * 60 s/min * 60 min/hr * 24 hr/d

# calculate DO consumption rate
DO_consumption_rate = bod_1day * Q_Ld # mg/L * L/d = mg/d

# divide by volume of nominal bottom layer
vol_m3 = 500 * 500 * 0.9 # m3
# convert to L
vol_L = vol_m3 * 1000 # m3 * 1000 L/m3

# determine DO consumed in one bottom grid cell
DO_consumption_bottom_cell_mgLd = DO_consumption_rate / vol_L # mg/d / L = mg/L/d

print(DO_consumption_bottom_cell_mgLd)
