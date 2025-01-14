"""
plot monthly means
"""
import numpy as np
import matplotlib.pylab as plt
from matplotlib import colormaps
from matplotlib.colors import ListedColormap

def plot_monthly_means(MONTHLYmean_DOdeep,
                        MONTHLYmean_DOin,
                        MONTHLYmean_Tflush,
                        MONTHLYmean_perchyp,
                        df_MONTHLYmean_DOdeep,
                        df_MONTHLYmean_DOin,
                        df_MONTHLYmean_Tflush):

    # initialize figure
    fig, ax = plt.subplots(2,2,figsize = (10,9))
    ax = ax.ravel()

    # define colormap for Tflush and DO concentration
    cmap_temp = colormaps['cubehelix_r'].resampled(256)
    cmap_tflush = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap
    cmap_oxy = ListedColormap(cmap_temp(np.linspace(0.2, 1, 256)))# get range of colormap


    #########################################################
    ##     Panel (a): DOdeep vs. DOin colored by Tflush    ##
    #########################################################

    # format figure
    ax[0].set_title('(a) All thirteen inlets', size=14, loc='left', fontweight='bold')
    ax[0].tick_params(axis='x', labelrotation=30)
    ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax[0].tick_params(axis='both', labelsize=12)
    ax[0].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    ax[0].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
    # plot
    ax[0].scatter(MONTHLYmean_DOin,MONTHLYmean_DOdeep,s=60, zorder=5, c='k', alpha=0.5)
    ax[0].plot([0,12],[0,12],color='gray')
    ax[0].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
    cs = ax[0].scatter(MONTHLYmean_DOin,MONTHLYmean_DOdeep,s=60, zorder=5, c=MONTHLYmean_Tflush, cmap=cmap_oxy)
    # create colorbarlegend
    cbar = fig.colorbar(cs)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
    cbar.outline.set_visible(False)
    ax[0].set_xlim([0,10])
    ax[0].set_ylim([0,10])

    #########################################################
    ##    Panel (b): DOin vs. Tflush colored by hyp vol    ##
    #########################################################

    # format figure
    ax[1].set_title('(b) All thirteen inlets', size=14, loc='left', fontweight='bold')
    ax[1].tick_params(axis='x', labelrotation=30)
    ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax[1].tick_params(axis='both', labelsize=12)
    ax[1].set_xlabel(r'Monthly mean T$_{flush}$ [days]', fontsize=12)
    ax[1].set_ylabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    ax[1].set_ylim([0,10])
    ax[1].set_xlim([0,85])
    # plot
    cmap_hyp = colormaps['gist_heat_r']
    cs_DO = ax[1].scatter(MONTHLYmean_Tflush,MONTHLYmean_DOin,s=60,zorder=5,edgecolor='gray',c=MONTHLYmean_perchyp,cmap=cmap_hyp)
    # create colorbarlegend
    cbar = fig.colorbar(cs_DO)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel('Monthly mean % hypoxic volume', rotation=90, fontsize=12)
    cbar.outline.set_visible(False)

    #########################################################
    ##               Panel (c): Crescent Bay               ##
    #########################################################

    # format figure 
    ax[2].set_title('(c) Crescent Bay', size=14, loc='left', fontweight='bold')
    ax[2].tick_params(axis='x', labelrotation=30)
    ax[2].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax[2].tick_params(axis='both', labelsize=12)
    ax[2].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    ax[2].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
    # plot
    ax[2].plot([0,11],[0,11],color='dimgray')
    ax[2].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
    # all inlets
    ax[2].scatter(MONTHLYmean_DOin,MONTHLYmean_DOdeep,s=60, zorder=5, color='gray',alpha=0.5, edgecolor='none')
    # highlight Crescent Bay
    cs = ax[2].scatter(df_MONTHLYmean_DOin['crescent'],df_MONTHLYmean_DOdeep['crescent'],marker='s',
                        s=150, zorder=6, c=df_MONTHLYmean_Tflush['crescent'], edgecolor='black',cmap=cmap_tflush,
                        linewidth=2, vmin=0, vmax=40)
    # create colorbarlegend
    cbar = fig.colorbar(cs)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
    cbar.outline.set_visible(False)
    ax[2].set_xlim([0,11])
    ax[2].set_ylim([0,11])

    #########################################################
    ##               Panel (d): Lynch Cove                 ##
    #########################################################

    # format figure
    ax[3].set_title('(d) Lynch Cove', size=14, loc='left', fontweight='bold')
    ax[3].tick_params(axis='x', labelrotation=30)
    ax[3].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax[3].tick_params(axis='both', labelsize=12)
    ax[3].set_xlabel(r'Monthly mean DO$_{in}$ [mg/L]', fontsize=12)
    ax[3].set_ylabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=12)
    # plot
    ax[3].plot([0,11],[0,11],color='dimgray')
    ax[3].text(0.9,0.9,'unity',rotation=45,va='center',ha='center',backgroundcolor='white',zorder=4, fontsize=10)
    # all inlets
    ax[3].scatter(MONTHLYmean_DOin,MONTHLYmean_DOdeep,s=60, zorder=5, color='gray',alpha=0.5, edgecolor='none')
    # highlight Lynch Cove
    cs = ax[3].scatter(df_MONTHLYmean_DOin['lynchcove'],df_MONTHLYmean_DOdeep['lynchcove'],marker='s',
                        s=150, zorder=6, c=df_MONTHLYmean_Tflush['lynchcove'], edgecolor='black',cmap=cmap_tflush,
                        linewidth=2, vmin=0, vmax=40)
    # create colorbarlegend
    cbar = fig.colorbar(cs)
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_ylabel(r'Monthly mean T$_{flush}$ [days]', rotation=90, fontsize=12)
    cbar.outline.set_visible(False)
    ax[3].set_xlim([0,11])
    ax[3].set_ylim([0,11])


    plt.tight_layout()
    plt.show()

    return