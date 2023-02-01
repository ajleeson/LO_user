"""
Plot as-run river time series.

"""
from lo_tools import Lfun
from lo_tools import plotting_functions as pfun

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

Ldir = Lfun.Lstart(gridname='cas6', tag='v00')

# load extraction (an xarray Dataset)
fn = Ldir['LOo'] / 'pre' / 'river' / 'cas6' / 'Data_roms' / 'extraction_2021.01.01_2021.02.19.nc'
x = xr.load_dataset(fn)

# # get climatology
# clm_fn = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag'] / 'Data_historical' / 'CLIM_flow_1980_2020.p'
# dfc = pd.read_pickle(clm_fn)

# # add the climatology, for practice
# x['transport_clim'] = 0*x.transport
# x['yearday'] = (('time'), x.time.to_index().dayofyear.to_numpy())
# ydvec = x.yearday.values

# # add the climatology to the xarray dataset (maybe use groupby instead?)
# for rn in list(x.riv.values):
#     if rn in dfc.columns:
#         this_riv = dfc[rn] # a Series
#         this_riv_clim = 0 * ydvec
#         for ii in range(1,367):
#             this_riv_clim[ydvec==ii] = this_riv[ii]
#         x.transport_clim.loc[:,rn] = this_riv_clim
#     else:
#         print('Missing ' + rn)
        
# plotting
plt.close('all')
# pfun.start_plot()
# fig = plt.figure()

# # for plotting time series we are better off using pandas
# df = pd.DataFrame(index=x.time.values)
# ii = 1
# for rn in ['Oak Harbor Lagoon','Whidbey east']:
#     ax = fig.add_subplot(2,1,ii)
#     df.loc[:,'Q'] = x.transport.sel(riv=rn).values
#     # df.loc[:,'Qclim'] = x.transport_clim.sel(riv=rn).values
#     # if ii == 1:
#     #     leg = True
#     # else:
#     #     leg = False
#     df.plot(ax=ax, grid=True)#, legend=leg)
#     ax.set_ylim(bottom=0)
#     if ii == 1:
#         ax.set_ylim(top=0.1)
#     ax.set_xlim(x.time.values[0], x.time.values[-1])
#     ax.text(.7,.9,rn.title(), transform=ax.transAxes)
#     plt.title(r'Flowrate ($m^3 \ s^{-1}$)')
#     if ii == 1:
#         ax.set_xticklabels([])
#     ii += 1

# for plotting time series we are better off using pandas
df = pd.DataFrame(index=x.time.values)
fig, ax = plt.subplots(2,1,figsize = (10,6), sharex = True)
for i,rn in enumerate(['Oak Harbor Lagoon','Whidbey east']):
    Q = x.transport.sel(riv=rn).values
    ax[i].plot(x.time.values,Q)
    ax[i].set_ylim(bottom=0)
    if i == 0:
        ax[i].set_ylim(top=0.1)
        ax[i].text(0.05,0.2,'WWTP',transform=ax[i].transAxes)
    ax[i].set_xlim(x.time.values[0], x.time.values[-1])
    ax[i].set_title(rn, fontsize =14)
    if i == 1:
        ax[i].text(0.05,0.2,'Tiny River',transform=ax[i].transAxes)
        ax[i].set_xticks(x.time.values[::3])
        ax[i].set_xticklabels(ax[i].get_xticks(), rotation = 45)
        date_form = mdates.DateFormatter("%m/%d")
        ax[i].xaxis.set_major_formatter(date_form)
    plt.suptitle(r'Flowrate ($m^3 \ s^{-1}$)', fontsize = 16)


plt.show()
pfun.end_plot()
