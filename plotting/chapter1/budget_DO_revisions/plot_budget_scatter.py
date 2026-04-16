"""
budget scatter plots using all 5 interface depth criteria
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


# ##########################################################
# ##            Get all variables for analysis            ##
# ##########################################################


# # initialize figure
# fig = plt.figure(figsize=(10, 5.5))
# gs = GridSpec(2, 12, figure=fig, height_ratios=[3, 4])#, wspace=2.4)
# axes = []  # list to store axes
# for i in range(4):
#     if i == 0:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3])
#     else:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3], sharex=axes[0])
#     axes.append(ax)
# for i in range(3):
#     ax = fig.add_subplot(gs[1, i*4:(i+1)*4], sharex=axes[0])
#     axes.append(ax)
# axes = [ax for ax in axes]

# # plot scatter points
# vars = ['ExchangeFlow','VerticalTransport',
#        'Photosynthesis','Consumption', 'd/dt(DO)',
#        'PhysicalResupply', 'NetEcosystemMetabolism']
# # colors = ['#0D4B91','#99C5F7','#8F0445','#FCC2DD','black','#488DDB','#F069A8']
# letters = ['(a) Exchange Flow','(b) Vertical','(c) Photosynthesis','(d) Consumption',
#            '(e) d/dt(DO)','(f) Physical Resupply', '(g) Net Ecosystem\nMetabolism']
# ylims = [ [-0.5,4], [-4.5,1], [-0.05,0.4], [-0.15,0.05],
#          [-0.2,0.1], [-0.5,0.2], [-0.06,0.25] ]

# # interface_types = ['onelaye','drdz','tef','og','halocline','oxycline']
# # inlet_colors = ['black','#62B6CB','#A8C256','#96031A','#957FEF','#F9627D']

# interface_types = ['drdz','tef','og','halocline','oxycline']
# inlet_colors = ['#62B6CB','#A8C256','#96031A','#957FEF','#F9627D']

# tef_budget_df = pd.read_csv('inlet_budgets_decline_period_mgL_day_tef.csv')

# # loop through depth interface types
# for t,type in enumerate(interface_types):
#     fn = 'inlet_budgets_decline_period_mgL_day_' + type + '.csv'
#     inlet_budget_df = pd.read_csv(fn)

#     for i,var in enumerate(vars):

#         varname = var

#         if type == 'onelaye':
#             if var == 'ExchangeFlow':
#                 varname = 'QinDOin'
#             elif var == 'VerticalTransport':
#                 varname = 'QoutDOout'

#         # # add error bars
#         # axes[i].errorbar(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df[var][:-1],
#         #                 xerr=inlet_budget_df['SepOctDeepDO_err[mg/L]'][:-1],
#         #                 yerr=inlet_budget_df[var+'_err'][:-1],
#         #                 fmt='o',color='black')
#         # plot points

#         if type == 'onelaye':
#             DOconc = inlet_budget_df['SepOctInletDO[mg/L]'][:-1]
#             # DOconc = inlet_budget_df['DOin-DOinlet[mg/L]'][:-1]
#         else:
#             DOconc = inlet_budget_df['SepOctDeepDO[mg/L]'][:-1]
#             # DOconc = inlet_budget_df['DOin-DOdeep[mg/L]'][:-1]

#         if i == 4:
#             axes[i].scatter(DOconc,inlet_budget_df[varname][:-1],
#                         color=inlet_colors[t], s=30, edgecolor='white',linewidth=0.5,zorder=5, label=type)
#         else:
#             axes[i].scatter(DOconc,inlet_budget_df[varname][:-1],
#                         color=inlet_colors[t], s=30, edgecolor='white',linewidth=0.5,zorder=5)
#         # # add mean line
#         # axes[i].axhline(inlet_budget_df[var].values[-1],0,8, color=colors[i])
#         # minval = inlet_budget_df[var].values[-1] - inlet_budget_df[var+'_err'].values[-1]
#         # maxval = inlet_budget_df[var].values[-1] + inlet_budget_df[var+'_err'].values[-1]
#         # axes[i].fill_between([0,9], [minval,minval], [maxval,maxval],
#         #                      color=colors[i],alpha=0.3)
#         # add zero line

#         if t == 0:
#             axes[i].axhline(0,0,8, color='gray',linestyle=':')
#             # axes[i].axvline(0,0,8, color='gray',linestyle=':')
#             # format panel
#             axes[i].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
#                         transform=axes[i].transAxes, zorder=6, va='top')
#             # axes[i].set_xlim([0,8])
#             # axes[i].set_ylim(ylims[i])

#             axes[i].tick_params(axis='both', labelsize=12)

#             if i in [0,4]:
#                 axes[i].set_ylabel('Decline period rates\n' + r'[mg L$^{-1}$ d$^{-1}$]',fontsize=12)
#             if i >=4 :
#                 axes[i].set_xlabel(r'Sep-Oct DO$_{deep}$ [mg L$^{-1}$]',fontsize=12)
#         if t == 4 and i == 4:
#             axes[i].legend(fontsize=12)

# plt.tight_layout()
# plt.show()

# # ----------------------- test new scatter plot (budget terms vs mean depth) -------------------------------

# # initialize figure
# fig = plt.figure(figsize=(10, 5.5))
# gs = GridSpec(2, 12, figure=fig, height_ratios=[3, 4])#, wspace=2.4)
# axes = []  # list to store axes
# for i in range(4):
#     if i == 0:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3])
#     else:
#         ax = fig.add_subplot(gs[0, i*3:(i+1)*3], sharex=axes[0])
#     axes.append(ax)
# for i in range(3):
#     ax = fig.add_subplot(gs[1, i*4:(i+1)*4], sharex=axes[0])
#     axes.append(ax)
# axes = [ax for ax in axes]

# # plot scatter points
# vars = ['ExchangeFlow','VerticalTransport',
#        'Photosynthesis','Consumption', 'd/dt(DO)',
#        'PhysicalResupply', 'NetEcosystemMetabolism']
# colors = ['#0D4B91','#99C5F7','#8F0445','#FCC2DD','black','#488DDB','#F069A8']
# letters = ['(a) Exchange Flow','(b) Vertical','(c) Photosynthesis','(d) Consumption',
#            '(e) d/dt(DO)','(f) Physical Resupply', '(g) Net Ecosystem\nMetabolism']
# # ylims = [ [-0.5,4], [-4.5,1], [-0.05,0.4], [-0.15,0.05],
# #          [-0.2,0.1], [-0.5,0.2], [-0.06,0.25] ]


# # loop through depth interface types
# for t,type in enumerate(interface_types):
#     fn = 'inlet_budgets_decline_period_mgL_day_' + type + '.csv'
#     inlet_budget_df = pd.read_csv(fn)

#     for i,var in enumerate(vars):

#         varname = var

#         if type == 'onelaye':
#             if var == 'ExchangeFlow':
#                 varname = 'QinDOin'
#             elif var == 'VerticalTransport':
#                 varname = 'QoutDOout'

#         # add error bars
#         # axes[i].errorbar(inlet_budget_df['MeanDepth[m]'][:-1],inlet_budget_df[var][:-1],
#         #                 yerr=inlet_budget_df[var+'_err'][:-1],
#         #                 fmt='o',color='black')
#         # plot points
#         if i == 4:
#             axes[i].scatter(tef_budget_df['MeanDepth[m]'][:-1],inlet_budget_df[varname][:-1], color=inlet_colors[t],
#                         s=30, edgecolor='white',linewidth=0.5,zorder=5, label=type)
#         else:
#             axes[i].scatter(tef_budget_df['MeanDepth[m]'][:-1],inlet_budget_df[varname][:-1], color=inlet_colors[t],
#                         s=30, edgecolor='white',linewidth=0.5,zorder=5)
#         # # add mean line
#         # axes[i].axhline(inlet_budget_df[var].values[-1],0,110, color=colors[i])
#         # minval = inlet_budget_df[var].values[-1] - inlet_budget_df[var+'_err'].values[-1]
#         # maxval = inlet_budget_df[var].values[-1] + inlet_budget_df[var+'_err'].values[-1]
#         # axes[i].fill_between([0,110], [minval,minval], [maxval,maxval],
#         #                     color=colors[i],alpha=0.3)
#         # add zero line
#         axes[i].axhline(0,0,110, color='gray',linestyle=':')
#         # format panel
#         axes[i].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
#                     transform=axes[i].transAxes, zorder=6, va='top')
#         # axes[i].set_xlim([0,110])
#         # axes[i].set_ylim(ylims[i])

#         axes[i].tick_params(axis='both', labelsize=12)

#         if t == 0:
#             axes[i].axhline(0,0,8, color='gray',linestyle=':')
#             # format panel
#             axes[i].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
#                         transform=axes[i].transAxes, zorder=6, va='top')
#             # axes[i].set_xlim([0,8])
#             # axes[i].set_ylim(ylims[i])

#             axes[i].tick_params(axis='both', labelsize=12)

#             if i in [0,4]:
#                 axes[i].set_ylabel('Decline period rates\n' + r'[mg L$^{-1}$ d$^{-1}$]',fontsize=12)
#             if i >=4 :
#                 axes[i].set_xlabel('Inlet mean depth [m]',fontsize=12)
#         if t == 5 and i == 4:
#             axes[i].legend(fontsize=12)

# plt.tight_layout()
# plt.show()

# # ----------------------- just d/dt(DO) -------------------------------

# # initialize figure
# fig, ax = plt.subplots(1,1,figsize=(7,5))

# # loop through depth interface types
# for t,type in enumerate(interface_types):
#     fn = 'inlet_budgets_decline_period_mgL_day_' + type + '.csv'
#     inlet_budget_df = pd.read_csv(fn)


#     # plot scatter points
#     var =  'd/dt(DO)'

#     if type == 'onelaye':
#             DOconc = inlet_budget_df['SepOctInletDO[mg/L]'][:-1]
#             DOconc_err = inlet_budget_df['SepOctInletDO_err[mg/L]'][:-1]
#             # DOconc = inlet_budget_df['DOin-DOinlet[mg/L]'][:-1]
#     else:
#         DOconc = inlet_budget_df['SepOctDeepDO[mg/L]'][:-1]
#         DOconc_err = inlet_budget_df['SepOctDeepDO_err[mg/L]'][:-1]
#         # DOconc = inlet_budget_df['DOin-DOdeep[mg/L]'][:-1]

#     # # add error bars
#     # ax.errorbar(DOconc,inlet_budget_df[var][:-1],
#     #                 xerr=inlet_budget_df['SepOctDeepDO_err[mg/L]'][:-1],
#     #                 yerr=inlet_budget_df[var+'_err'][:-1],
#     #                 fmt='o',color='black',elinewidth=0.5)
#     # plot points
#     ax.scatter(DOconc,inlet_budget_df[var][:-1], color=inlet_colors[t],
#                     s=30, edgecolor='None',linewidth=0.5,zorder=5)
#     # # add mean line
#     # ax.axhline(inlet_budget_df[var].values[-1],0,8, color='teal')
#     # minval = inlet_budget_df[var].values[-1] - inlet_budget_df[var+'_err'].values[-1]
#     # maxval = inlet_budget_df[var].values[-1] + inlet_budget_df[var+'_err'].values[-1]
#     # ax.fill_between([0,9], [minval,minval], [maxval,maxval],
#     #                         color='lightseagreen',alpha=0.3)


#     # 2. Set up the ODR model (Linear: y = m*x + b)
#     def linear_func(p, x):
#         m, b = p
#         return m * x + b

#     model = odr.Model(linear_func)
#     data = odr.RealData(DOconc, inlet_budget_df[var][:-1], sx=DOconc_err, sy=inlet_budget_df[var+'_err'][:-1])
#     my_odr = odr.ODR(data, model, beta0=[1., 0.])

#     # 3. Run regression and extract results
#     output = my_odr.run()
#     slope = output.beta[0]
#     slope_std_err = output.sd_beta[0]

#     # 4. Calculate Significance (t-test)
#     # Degrees of freedom = n - 2
#     df = len(DOconc) - 2
#     t_stat = slope / slope_std_err
#     p_value = 2 * (1 - stats.t.cdf(np.abs(t_stat), df))

#     print('============================\n{}'.format(type))
#     print(f"Slope: {slope:.4f}")
#     print(f"P-value: {p_value:.4f}")
#     print("Significant" if p_value < 0.05 else "Not Significant")

#     # # add linear fit
#     # input_array = np.array([DOconc, [1]*len(DOconc)]).T
#     # B,a,b,c = lstsq(input_array,inlet_budget_df[var][:-1])
#     # slope = B[0]
#     # intercept = B[1]
#     # x=[0,8]
#     # y=[intercept+slope*xval for xval in x]
#     # ax.plot(x,y,color=inlet_colors[t])
#     # # calculate r^2 and p value
#     # r,p = pearsonr(DOconc,inlet_budget_df[var][:-1])
#     # print('\n===============================================')
#     # print('{}: d/dt(DO) dependence on DO_deep'.format(type))
#     # print('   r = {}'.format(r))
#     # print('   R^2 = {}'.format(r**2))
#     # print('   p = {}'.format(p))
#     # print('===============================================\n')


#     # add zero line
#     ax.axhline(0,0,8, color='gray',linestyle=':')

#     ax.set_xlim([0,8])

#     ax.tick_params(axis='both', labelsize=12)
#     ax.set_ylabel(r'Decline period d/dt(DO) [mg L$^{-1}$ d$^{-1}$]',fontsize=12)
#     ax.set_xlabel(r'Hypoxic season DO$_{deep}$ [mg L$^{-1}$]',fontsize=12)

# plt.tight_layout()
# plt.show()


# # ----------------------- just DOin-DOdeep vs. depth -------------------------------

# # initialize figure
# fig, ax = plt.subplots(1,1,figsize=(7,5))

# # loop through depth interface types
# for t,type in enumerate(interface_types):
#     fn = 'inlet_budgets_decline_period_mgL_day_' + type + '.csv'
#     inlet_budget_df = pd.read_csv(fn)

#     DOconc = inlet_budget_df['DOin-DOdeep[mg/L]'][:-1]

#     # # add error bars
#     # ax.errorbar(DOconc,inlet_budget_df[var][:-1],
#     #                 xerr=inlet_budget_df['SepOctDeepDO_err[mg/L]'][:-1],
#     #                 yerr=inlet_budget_df[var+'_err'][:-1],
#     #                 fmt='o',color='black',elinewidth=0.5)
#     # plot points
#     ax.scatter(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],DOconc, color=inlet_colors[t],
#                     s=30, edgecolor='None',linewidth=0.5,zorder=5)
#     # # add mean line
#     # ax.axhline(inlet_budget_df[var].values[-1],0,8, color='teal')
#     # minval = inlet_budget_df[var].values[-1] - inlet_budget_df[var+'_err'].values[-1]
#     # maxval = inlet_budget_df[var].values[-1] + inlet_budget_df[var+'_err'].values[-1]
#     # ax.fill_between([0,9], [minval,minval], [maxval,maxval],
#     #                         color='lightseagreen',alpha=0.3)


#     # # add linear fit
#     # input_array = np.array([DOconc, [1]*len(DOconc)]).T
#     # B,a,b,c = lstsq(input_array,inlet_budget_df[var][:-1])
#     # slope = B[0]
#     # intercept = B[1]
#     # x=[0,8]
#     # y=[intercept+slope*xval for xval in x]
#     # ax.plot(x,y,color=inlet_colors[t])
#     # # calculate r^2 and p value
#     # r,p = pearsonr(DOconc,inlet_budget_df[var][:-1])
#     # print('\n===============================================')
#     # print('{}: d/dt(DO) dependence on DO_deep'.format(type))
#     # print('   r = {}'.format(r))
#     # print('   R^2 = {}'.format(r**2))
#     # print('   p = {}'.format(p))
#     # print('===============================================\n')


#     # add zero line
#     ax.axhline(0,0,8, color='gray',linestyle=':')

#     # ax.set_xlim([0,8])

#     ax.tick_params(axis='both', labelsize=12)
#     ax.set_ylabel(r'Decline period DO$_{in}$ - DO$_{deep}$ [mg L$^{-1}$]',fontsize=12)
#     # ax.set_xlabel(r'Mean Depth [m]',fontsize=12)
#     ax.set_xlabel(r'Hypoxic season DO$_{deep}$ [mg L$^{-1}$]',fontsize=12)

# plt.tight_layout()
# plt.show()


############################

# initialize figure
fig,axes = plt.subplots(2,3, figsize=(12,8), sharex=True, sharey=True)
ax = axes.ravel()

interface_types = ['drdz','tef','og','halocline','oxycline']


# loop through depth interface types
for t,type in enumerate(interface_types):
    fn = 'inlet_budgets_decline_period_mgL_day_' + type + '.csv'
    inlet_budget_df = pd.read_csv(fn)

    # plot points
    ax[t].errorbar(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df['d/dt(DO)'][:-1],
                    xerr=inlet_budget_df['SepOctDeepDO_err[mg/L]'][:-1],
                    yerr=inlet_budget_df['d/dt(DO)'+'_err'][:-1],
                    fmt='o',color='black')
    
    # plot points
    cmap_temp = colormaps['gnuplot2_r'].resampled(256)
    cmap_depth = ListedColormap(cmap_temp(np.linspace(0.1, 0.8, 256)))# get range of colormap
    cs = ax[t].scatter(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df['d/dt(DO)'][:-1],
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
    # # format panel
    # axes[t].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
    #              transform=axes[i].transAxes, zorder=6, va='top')
    ax[t].set_xlim([0,8])
    # axes[t].set_ylim(ylims[i])

    ax[t].tick_params(axis='both', labelsize=12)

    # if t == 0:
    #     axes[i].axhline(0,0,8, color='gray',linestyle=':')
    #     # format panel
    #     axes[i].text(0.03,0.96,letters[i],fontsize=11,fontweight='bold',
    #                 transform=axes[i].transAxes, zorder=6, va='top')
    ax[t].set_xlim([0,8])
    # ax[t].set_ylim(ylims[i])

    #     axes[i].tick_params(axis='both', labelsize=12)

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
    ax[t].text(0.1,0.8,criteria + '\ninterface', ha='left',fontsize=12,
               transform=ax[t].transAxes, fontweight='bold')
    
    # calculate correlation
    r,p = pearsonr(inlet_budget_df['SepOctDeepDO[mg/L]'][:-1],inlet_budget_df['d/dt(DO)'][:-1])
    ax[t].text(0.03,0.16,'R = {}\np = {}'.format(round(r,2),round(p,2)),
                 transform=ax[t].transAxes, zorder=6, va='top')

plt.subplots_adjust(right=0.9)
# plt.tight_layout()
plt.show()