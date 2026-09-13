"""
Compare budget terms between the different interface depth methods,
including the one-layer budget
"""

import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pylab as plt
import get_two_layer
import matplotlib.dates as mdates
from matplotlib.gridspec import GridSpec
from scipy import odr, stats
from matplotlib import colormaps
from matplotlib.colors import ListedColormap
from scipy.stats import pearsonr

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

Ldir = Lfun.Lstart()

plt.close('all')

###############################################################
##                   d/dt(DO) vs. DOdeep                     ##
###############################################################


# initialize figure
fig,axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
ax = axes.ravel()

interface_types = ['drdz','tef','og','halocline','oxycline','onelayer']


# loop through depth interface types
for t,type in enumerate(interface_types):
    fn = 'inlet_budgets_decline_period_mgL_day_' + type + '.csv'
    inlet_budget_df = pd.read_csv(fn)

    if type == 'onelayer':
        DO_type = 'SepOctInletDO[mg/L]'
        DO_type_err = 'SepOctInletDO_err[mg/L]'
    else:
        DO_type = 'SepOctDeepDO[mg/L]'
        DO_type_err = 'SepOctDeepDO_err[mg/L]'

    # plot points
    ax[t].errorbar(inlet_budget_df[DO_type][:-1],inlet_budget_df['d/dt(DO)'][:-1],
                    xerr=inlet_budget_df[DO_type_err][:-1],
                    yerr=inlet_budget_df['d/dt(DO)'+'_err'][:-1],
                    fmt='o',color='black')
    
    # plot points
    cmap_temp = colormaps['gnuplot2_r'].resampled(256)
    cmap_depth = ListedColormap(cmap_temp(np.linspace(0.1, 0.8, 256)))# get range of colormap
    cs = ax[t].scatter(inlet_budget_df[DO_type][:-1],inlet_budget_df['d/dt(DO)'][:-1],
                    s=50, zorder=5,c=inlet_budget_df['MeanDepth[m]'][:-1], cmap=cmap_depth, vmin=0, vmax=110)

    # create colorbarlegend
    if t == 0:
        cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
        cbar = fig.colorbar(cs, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(r'Inlet mean depth [m]', fontsize=12)
        cbar.outline.set_visible(False)

    # add zero line
    ax[t].axhline(0,0,8, color='gray',linestyle=':')

    ax[t].set_xlim([0,8])

    ax[t].tick_params(axis='both', labelsize=12)

    ax[t].set_xlim([0,8])


    if t in [0,3]:
        ax[t].set_ylabel('Decline period d/dt(DO)\n' + r'[mg L$^{-1}$ d$^{-1}$]',fontsize=12)
    if t >= 3 :
        ax[t].set_xlabel(r'Sep-Oct DO$_{deep}$ [mg L$^{-1}$]',fontsize=12)

    # label y-axis
    criteria = ''
    if type == 'og':
        criteria = 'Uniform 1/3'
    elif type == 'tef':
        criteria = 'TEF'
    elif type == 'drdz':
        criteria = r'd$\rho$/dz'
    elif type == 'halocline':
        criteria = 'Halocline'
    elif type == 'oxycline':
        criteria = 'Oxycline'
    elif type == 'onelayer':
        criteria = 'One-layer'
    ax[t].text(0.1,0.8,criteria + '\ninterface', ha='left',fontsize=12,
               transform=ax[t].transAxes, fontweight='bold')
    
    # calculate correlation
    r,p = pearsonr(inlet_budget_df[DO_type][:-1],inlet_budget_df['d/dt(DO)'][:-1])
    ax[t].text(0.03,0.16,'R = {}\np = {}'.format(round(r,2),round(p,2)),
                 transform=ax[t].transAxes, zorder=6, va='top')

plt.subplots_adjust(right=0.9)
# plt.tight_layout()
plt.show()


###############################################################
##         d/dt(DO) correlation between interfaces           ##
###############################################################

lims = [-0.06,0]

# initialize figure
fig,axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
ax = axes.ravel()

fig.suptitle('Different Interface Depth vs. One-Layer\nDecline period d/dt(DO) [mg/L/day]', fontsize=14)

interface_types = ['drdz','tef','og','halocline','oxycline','onelayer']
onelayer_budget_df = pd.read_csv('inlet_budgets_decline_period_mgL_day_onelayer.csv')


# loop through depth interface types
for t,type in enumerate(interface_types):
    fn = 'inlet_budgets_decline_period_mgL_day_' + type + '.csv'
    inlet_budget_df = pd.read_csv(fn)

    if type == 'onelayer':
        DO_type = 'SepOctInletDO[mg/L]'
        DO_type_err = 'SepOctInletDO_err[mg/L]'
    else:
        DO_type = 'SepOctDeepDO[mg/L]'
        DO_type_err = 'SepOctDeepDO_err[mg/L]'

    # # plot points
    # ax[t].errorbar(onelayer_budget_df['d/dt(DO)'][:-1],inlet_budget_df['d/dt(DO)'][:-1],
    #                 xerr=inlet_budget_df['d/dt(DO)'+'_err'][:-1],
    #                 yerr=onelayer_budget_df['d/dt(DO)'+'_err'][:-1],
    #                 fmt='o',color='black')
    
    # plot points
    cmap_temp = colormaps['rainbow_r'].resampled(256)
    cmap_DO = ListedColormap(cmap_temp(np.linspace(0, 0.9, 256)))# get range of colormap
    cs = ax[t].scatter(onelayer_budget_df['d/dt(DO)'][:-1],inlet_budget_df['d/dt(DO)'][:-1],
                    s=50, zorder=5,c=inlet_budget_df[DO_type][:-1], cmap=cmap_DO, vmin=0, vmax=8)
     # create colorbarlegend
    if t == 0:
        cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
        cbar = fig.colorbar(cs, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label('Sep & Oct deep DO [mg/L]', fontsize=12)
        cbar.outline.set_visible(False)

    # add 1-1 line
    x = lims
    ax[t].plot(x,x, color='gray',linestyle=':')

    ax[t].set_xlim(lims)
    ax[t].set_ylim(lims)

    ax[t].tick_params(axis='both', labelsize=12)
    
    # calculate correlation
    r,p = pearsonr(onelayer_budget_df['d/dt(DO)'][:-1],inlet_budget_df['d/dt(DO)'][:-1])
    ax[t].text(0.7,0.16,'R = {}\np = {}'.format(round(r,2),round(p,4)),
                 transform=ax[t].transAxes, zorder=6, va='top', ha='left')
    
    # label y-axis
    criteria = ''
    if type == 'og':
        criteria = 'Uniform 1/3'
    elif type == 'tef':
        criteria = 'TEF'
    elif type == 'drdz':
        criteria = r'd$\rho$/dz'
    elif type == 'halocline':
        criteria = 'Halocline'
    elif type == 'oxycline':
        criteria = 'Oxycline'
    elif type == 'onelayer':
        criteria = 'One-layer'
    ax[t].text(0.1,0.8,criteria + '\ninterface', ha='left',fontsize=12,
               transform=ax[t].transAxes, fontweight='bold')


###############################################################
##           DOdeep correlation between interfaces           ##
###############################################################

lims = [0,8]

# initialize figure
fig,axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
ax = axes.ravel()

fig.suptitle('Different Interface Depth vs. One-Layer\nSep & Oct Deep DO [mg/L]', fontsize=14)

interface_types = ['drdz','tef','og','halocline','oxycline','onelayer']
onelayer_budget_df = pd.read_csv('inlet_budgets_decline_period_mgL_day_onelayer.csv')


# loop through depth interface types
for t,type in enumerate(interface_types):
    fn = 'inlet_budgets_decline_period_mgL_day_' + type + '.csv'
    inlet_budget_df = pd.read_csv(fn)

    if type == 'onelayer':
        DO_type = 'SepOctInletDO[mg/L]'
        DO_type_err = 'SepOctInletDO_err[mg/L]'
    else:
        DO_type = 'SepOctDeepDO[mg/L]'
        DO_type_err = 'SepOctDeepDO_err[mg/L]'

    # # plot points
    # ax[t].errorbar(onelayer_budget_df['d/dt(DO)'][:-1],inlet_budget_df['d/dt(DO)'][:-1],
    #                 xerr=inlet_budget_df['d/dt(DO)'+'_err'][:-1],
    #                 yerr=onelayer_budget_df['d/dt(DO)'+'_err'][:-1],
    #                 fmt='o',color='black')
    
    # plot points
    cmap_temp = colormaps['rainbow_r'].resampled(256)
    cmap_DO = ListedColormap(cmap_temp(np.linspace(0, 0.9, 256)))# get range of colormap
    cs = ax[t].scatter(onelayer_budget_df['SepOctInletDO[mg/L]'][:-1],inlet_budget_df[DO_type][:-1],
                    s=50, zorder=5,c=inlet_budget_df['d/dt(DO)'][:-1], cmap=cmap_DO, vmin=-0.06, vmax=0)
     # create colorbarlegend
    if t == 0:
        cbar_ax = fig.add_axes([0.92, 0.1, 0.02, 0.8])
        cbar = fig.colorbar(cs, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label('Decrease period d/dt(DO) [mg/L/day]', fontsize=12)
        cbar.outline.set_visible(False)

    # add 1-1 line
    x = lims
    ax[t].plot(x,x, color='gray',linestyle=':')

    ax[t].set_xlim(lims)
    ax[t].set_ylim(lims)

    ax[t].tick_params(axis='both', labelsize=12)
    
    # calculate correlation
    r,p = pearsonr(onelayer_budget_df['SepOctInletDO[mg/L]'][:-1],inlet_budget_df[DO_type][:-1])
    ax[t].text(0.7,0.16,'R = {}\np = {}'.format(round(r,2),round(p,4)),
                 transform=ax[t].transAxes, zorder=6, va='top', ha='left')
    
    # label y-axis
    criteria = ''
    if type == 'og':
        criteria = 'Uniform 1/3'
    elif type == 'tef':
        criteria = 'TEF'
    elif type == 'drdz':
        criteria = r'd$\rho$/dz'
    elif type == 'halocline':
        criteria = 'Halocline'
    elif type == 'oxycline':
        criteria = 'Oxycline'
    elif type == 'onelayer':
        criteria = 'One-layer'
    ax[t].text(0.1,0.8,criteria + '\ninterface', ha='left',fontsize=12,
               transform=ax[t].transAxes, fontweight='bold')

