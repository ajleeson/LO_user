"""
Code to draw plot regional bathymetry
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import xarray as xr
import cmocean
from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()


background = 'white'

# Get LiveOcean grid info --------------------------------------------------

# get the grid data
ds = xr.open_dataset('../../../LO_data/grids/cas7/grid.nc')
z = -ds.h.values
mask_rho = np.transpose(ds.mask_rho.values)
lon = ds.lon_rho.values
lat = ds.lat_rho.values
X = lon[0,:] # grid cell X values
Y = lat[:,0] # grid cell Y values
plon, plat = pfun.get_plon_plat(lon,lat)
# make a version of z with nans where masked
zm = z.copy()
zm[np.transpose(mask_rho) == 0] = np.nan

# Create bathymetry plot --------------------------------------------------------------

# Initialize figure
plt.close('all')
fs = 10
pfun.start_plot(fs=fs, figsize=(7,8))
fig = plt.figure()
plt.subplots_adjust(wspace=0, hspace=0)
# create colormap

zm[np.transpose(mask_rho) != 0] = -1
newcmap = plt.get_cmap('Greys_r')
newcmap.set_bad(background,1.)

ax0 = fig.add_subplot(1,1,1)
cs = ax0.pcolormesh(plon, plat, zm, vmin=-5, vmax=0, cmap=newcmap)

# TEST BED LOCATIONS ------------------------------------------
testbed_color = 'deeppink'
ax0.scatter(-123.457488,48.129378,      color=testbed_color,s=30,edgecolor='None') # Port Angeles (Ebb Carbon)
ax0.scatter(-123.045000,48.078611,      color=testbed_color,s=30,edgecolor='None') # Sequim Bay (Ebb Carbon)
ax0.scatter(-123.0121,48.5453,          color=testbed_color,s=30,edgecolor='None') # FHLOO
ax0.scatter(-122.342908,47.584536,      color=testbed_color,s=30,edgecolor='None') # Terminal 30 near Port of Seattle
ax0.scatter(-122.43131423,47.60079806,  color=testbed_color,s=30,edgecolor='None') # South KC WWTP (~same size as West Point, and near Port of Seattle)
ax0.text(-123.0121,48.5453       +0.05,'1',     color=testbed_color, size=11,fontweight='bold', ha='center',va='bottom') # FHLOO
ax0.text(-123.457488,48.129378   +0.05,'2',     color=testbed_color, size=11,fontweight='bold', ha='center',va='bottom') # Port Angeles
ax0.text(-123.045000,48.078611   +0.05,'3',     color=testbed_color, size=11,fontweight='bold', ha='center',va='bottom') # Sequim Bay
ax0.text(-122.43131423 -0.08,47.60079806-0.05,'4', color=testbed_color, size=11,fontweight='bold', ha='right',va='top') # South King WWTP
ax0.text(-122.342908   +0.08,47.584536  -0.05,'5', color=testbed_color, size=11,fontweight='bold', ha='left', va='top') # Terminal 30 (Port of Seattle)
# label locations
ax0.text(-126.5,49.6,'Testbed locations', #48.4
         color=testbed_color, fontweight='bold',size=11, ha='left', va='top')
ax0.text(-126.4,49.4,'1. FHLOO\n2. Port Angeles PNNL\n3. Seqium Bay PNNL\n4. South King WWTP\n5. Port of Seattle (Terminal 30)',
         color=testbed_color, size=11, ha='left', va='top')

ax0.scatter(-123.457488,48.129378,      color='None',s=100,edgecolor=testbed_color) # Port Angeles (Ebb Carbon)
ax0.scatter(-123.045000,48.078611,      color='None',s=100,edgecolor=testbed_color) # Sequim Bay (Ebb Carbon)
ax0.scatter(-123.0121,48.5453,          color='None',s=100,edgecolor=testbed_color) # FHLOO
ax0.scatter(-122.342908,47.584536,      color='None',s=100,edgecolor=testbed_color) # Terminal 30 near Port of Seattle
ax0.scatter(-122.43131423,47.60079806,  color='None',s=100,edgecolor=testbed_color) # South KC WWTP (~same size as West Point, and near Port of Seattle)

# NATURAL LOCATIONS -------------------------------------------
natural_color = 'royalblue'
ax0.scatter(-125.2220,50.1160,          color=natural_color,s=30,edgecolor='None') # Quadra Island Field Station at Hyacinthe Bay (good based on dye, bad based on surf mixed layer)
ax0.scatter(-123.20643409,49.22495359,  color=natural_color,s=30,edgecolor='None') # Iona WWTP
ax0.scatter(-123.06369538,47.35532783,  color=natural_color,s=30,edgecolor='None') # Alderbrook Resort & Spa
ax0.scatter(-123.946536,46.199841,      color=natural_color,s=30,edgecolor='None') # Point Adams Research Station, Northwest Fisheries Science Center
ax0.text(-125.2220,50.1160         +0.05,'6',color=natural_color, size=11,fontweight='bold', ha='center',va='bottom') # Quadra Island Field Station at Hyacinthe Bay (good based on dye, bad based on surf mixed layer)
ax0.text(-123.20643409,49.22495359 +0.05,'7',color=natural_color, size=11,fontweight='bold', ha='center',va='bottom') # Iona WWTP (Southern Strait of Georgia; bad based on dye, good based on surf mixed layer)
ax0.text(-124.202353,48.344473     +0.05,'8',color=natural_color, size=11,fontweight='bold', ha='center',va='bottom') # Middle of Strait of Juan de Fuca (mixed? maybe bad?)
ax0.text(-123.06369538,47.35532783 +0.05,'9',color=natural_color, size=11,fontweight='bold', ha='center',va='bottom') # Alderbrook Resort & Spa (Hood Canal; good location)
ax0.text(-123.946536,46.199841     +0.05,'10',color=natural_color, size=11,fontweight='bold', ha='center',va='bottom') # Point Adams Research Station; Columbia River (good location)

# label locations
ax0.text(-126.5,48.6,'Natural locations', #47.4
         color=natural_color, fontweight='bold',size=11, ha='left', va='top')
ax0.text(-126.4,48.4,'6. Northern SoG (Hakai Station)\n7. Southern SoG (Iona WWTP)\n8. SJdF\n9. Hood Canal (Alderbrook WWTP)\n10. Columbia River (Point Adams)',
         color=natural_color, size=11, ha='left', va='top')


ax0.scatter(-125.2220,50.1160,          color='None',s=100,edgecolor=natural_color) # Quadra Island Field Station at Hyacinthe Bay (good based on dye, bad based on surf mixed layer)
ax0.scatter(-123.20643409,49.22495359,  color='None',s=100,edgecolor=natural_color) # Iona WWTP (Southern Strait of Georgia; bad based on dye, good based on surf mixed layer)
ax0.scatter(-123.06369538,47.35532783,  color='None',s=100,edgecolor=natural_color) # Alderbrook Resort & Spa (Hood Canal; good location)
ax0.scatter(-123.946536,46.199841,      color='None',s=100,edgecolor=natural_color) # Point Adams Research Station; Columbia River (good location)
ax0.scatter(-124.202353,48.344473,      color='None',s=100,edgecolor=natural_color) # Middle of Strait of Juan de Fuca (mixed? maybe bad?)


# ADDITIONAL LOCATIONS NOT IN PROPOSED PLAN
additional_color = 'black'
ax0.scatter(-122.42066969,47.78164567,marker='x',color=additional_color,s=40) # Brightwater WWTP (too far south of Whidbey)
ax0.scatter(-122.45016076,47.66183568,marker='x',color=additional_color,s=40) # West Point WWTP  (too far from Port of Seattle)
ax0.scatter(-122.90528374,47.05846051,marker='x',color=additional_color,s=40) # LOTT WWTP 
ax0.scatter(-123.450444,48.650500,    marker='x',color=additional_color,s=40) # IOS (Institute of Ocean Sciences)
ax0.scatter(-122.684473,48.508155,    marker='x',color=additional_color,s=40) # Shannon Point Marine Center
ax0.scatter(-122.324606,47.348850,    marker='x',color=additional_color,s=40) # MaST Center Aquarium

# ax0.scatter(-123.40467062,48.12833167,    marker='x',color=additional_color,s=40) # Port Angeles STP


# format figure
pfun.dar(ax0)
# Set axis limits
ax0.set_xlim(-126.6,-122)#X[i2]) # Salish Sea
ax0.set_ylim(46,50.5) # Salish Sea


plt.xticks(rotation=30, color='gray', size=12)
plt.yticks(color='gray', size=12)



plt.show()
