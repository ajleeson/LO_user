"""
plot monthly mean DOdeep vs. % hypoxic volume
and monthly mean DOdeep time series
"""
import matplotlib.pylab as plt
import matplotlib.dates as mdates
import helper_functions

def dodeep_hypvol_timeseries(MONTHLYmean_DOdeep,
                            MONTHLYmean_perchyp,
                            DOconcen_dict,
                            dates_local_daily,
                            dates_local_hrly,
                            inlets,minday,maxday):
    

    # initialize figure
    fig, ax = plt.subplots(1,2,figsize = (12,5),gridspec_kw={'width_ratios': [1, 1.5]})

    # Deep DO vs. % hypoxic volume
    ax[0].set_title('(a) All thirteen inlets', size=14, loc='left', fontweight='bold')
    ax[0].tick_params(axis='x', labelrotation=30)
    ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax[0].tick_params(axis='both', labelsize=12)
    ax[0].set_xlabel(r'Monthly mean DO$_{deep}$ [mg/L]', fontsize=14)
    ax[0].set_ylabel('Monthly mean % hypoxic volume', fontsize=14)
    # plot
    ax[0].scatter(MONTHLYmean_DOdeep,MONTHLYmean_perchyp,alpha=0.3,s=80,zorder=5,color='navy')
    ax[0].set_xlim([0,10])
    ax[0].set_ylim([0,100])

    # Deep DO timeseries
    ax[1].set_title('(b) All thirteen inlets', size=14, loc='left', fontweight='bold')
    # format grid
    ax[1].tick_params(axis='x', labelrotation=30)
    loc = mdates.MonthLocator(interval=1)
    ax[1].xaxis.set_major_locator(loc)
    ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax[1].tick_params(axis='both', labelsize=12)
    # add drawdown period
    ax[1].axvline(dates_local_daily[minday],0,12,color='grey')
    ax[1].axvline(dates_local_daily[maxday],0,12,color='grey')
    # loop through inlets
    for i,inlet in enumerate(inlets):
        # get average deep layer DO
        deep_lay_DO_alltime = DOconcen_dict[inlet]['Deep Layer DO']
        # 30-day hanning window
        deep_lay_DO_alltime = helper_functions.lowpass(deep_lay_DO_alltime.values,n=30)
        # plot
        ax[1].plot(dates_local_daily,deep_lay_DO_alltime,linewidth=1,color='navy',alpha=0.5)

    # format labels
    ax[1].set_xlim([dates_local_hrly[0],dates_local_hrly[-2]])
    ax[1].set_ylim([0,10])
    ax[1].set_ylabel(r'DO$_{deep}$ [mg/L]',fontsize=14)
    plt.tight_layout()
    plt.show()
    
    return