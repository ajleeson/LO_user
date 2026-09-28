"""
Make pie charts comparie NH4 and NO3 percentage in WWTPs compared to all of Puget Sound

"""

# import things
from subprocess import Popen as Po
from subprocess import PIPE as Pi
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
import xarray as xr
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator
import matplotlib.patches as patches
from matplotlib.offsetbox import (OffsetImage, AnnotationBbox)
import matplotlib.image as image
import pandas as pd
import cmocean
import matplotlib.pylab as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patheffects as PathEffects
import pinfo

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun

plt.close('all')

#########################
# Data from change_in_N.py

# Puget Sound No-loading NO3 and NH4 average concentrations
PS_noload_NO3_avg_mmol_m3 = 26.53 # [mmol/m3]
PS_noload_NH4_avg_mmol_m3 =  0.78 # [mmol/m3]

#########################
# Data from change_in_N.py

WWWTP_NO3_avg_mmol_m3 =  440.64068171910645 # [mmol/m3]
WWWTP_NH4_avg_mmol_m3 = 1858.4869055732954 # [mmol/m3]

#########################
# plot pie charts

# labels = 'NO3', 'NH4'
# sizes_pugetsound = [PS_noload_NO3_avg_mmol_m3, PS_noload_NH4_avg_mmol_m3]
# sizes_wwtps = [WWWTP_NO3_avg_mmol_m3, WWWTP_NH4_avg_mmol_m3]
# colors = ['purple','crimson']

# fig, ax = plt.subplots(1,2)
# ax[0].pie(sizes_pugetsound, labels=labels, colors=colors,
#         textprops={'size': 'large'}, autopct='%1.1f%%')
# ax[1].pie(sizes_wwtps, labels=labels, colors=colors,
#         textprops={'size': 'large'}, autopct='%1.1f%%')


labels = 'NO3', 'NH4'
sizes_pugetsound = [PS_noload_NO3_avg_mmol_m3, PS_noload_NH4_avg_mmol_m3]
sizes_wwtps = [WWWTP_NO3_avg_mmol_m3, WWWTP_NH4_avg_mmol_m3]
colors = ['purple','crimson']

fig, ax = plt.subplots(1,2)
wedges_pugetsound, texts_pugetsound, autotexts_pugetsound = ax[0].pie(sizes_pugetsound, labels=None, colors=colors,
        textprops={'color':'white', 'fontsize':24, 'fontweight':'bold'}, autopct='%1.1f%%')
wedges_wwtps, texts_wwtps, autotexts_wwtps = ax[1].pie(sizes_wwtps, labels=None, colors=colors,
        textprops={'color':'white', 'fontsize':24, 'fontweight':'bold'}, autopct='%1.1f%%')

ax[0].set_title('Puget Sound\n(No-loading)', fontsize=24, fontweight='bold')
ax[1].set_title('WWTPs', fontsize=24, fontweight='bold')

# ax[1].legend(wedges_pugetsound, labels, loc='upper left', fontsize=14)
plt.tight_layout(rect=[0, 0.1, 1, 1])