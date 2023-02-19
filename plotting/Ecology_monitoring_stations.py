"""
Code to plot Ecology monitoring stations
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as col
import pandas as pd
import math
import xarray as xr
import cmocean
from lo_tools import plotting_functions as pfun


# Read monitoring site locations
fp = '../../LO_user/downloaded_data/EIM_monitoring_locations.csv'
df = pd.read_csv(fp)

lat = df['Calculated_Latitude_Decimal_Degrees_NAD83HARN']
lon = df['Calculated_Longitude_Decimal_Degrees_NAD83HARN']

fig = plt.figure(figsize=(6,8))
ax = fig.add_subplot(111)
# pfun.add_coast(ax,color='black')
pfun.dar(ax)
ax.scatter(lon,lat,s=50, edgecolors='k', color='deeppink' )
ax.set_xlim(-123.5,-122)
ax.set_ylim(47,48.9)

plt.show()