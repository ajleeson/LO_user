
'''
Script to show where Salish Sea is located on the globe
'''

import cartopy
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import xarray as xr
import cmocean
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()


########################################################
##                Salish Sea on Globe                ##
########################################################

plt.close('all')

# pick where to center globe
# ax = plt.axes(projection=ccrs.Orthographic(-69, 30))
ax = plt.axes(projection=ccrs.Orthographic(-100, 30))


# color land and water
ax.add_feature(cartopy.feature.OCEAN, zorder=0, color='#B1DEE2')
ax.add_feature(cartopy.feature.LAND, zorder=0, color='white')

# set firnat to globe, and add axis
ax.set_global()
ax.gridlines(linewidth=0.2,color='#389198')

# add rectangle around Salish Sea
ax.plot([-124.9854945365928, -124.9854945365928, -122, -122, -124.9854945365928],
        [46.81655190352876, 50.132729215720914, 50.132729215720914, 46.81655190352876, 46.81655190352876],
         color='black', linewidth=3, transform=ccrs.Geodetic())

plt.show()