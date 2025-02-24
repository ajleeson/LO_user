"""
plot d/dt(DO) timeseries
"""
import matplotlib.pylab as plt
import matplotlib.dates as mdates
import helper_functions
import numpy as np

def d_dt_DO_timeseries(DOconcen_dict,
                        dates_local_daily,
                        dates_local_hrly,
                        inlets,minday,maxday,
                        dimensions_dict,deeplay_dict,
                        shallowlay_dict,):
    

    # initialize figure
    fig, ax = plt.subplots(1,2,figsize = (10,5), sharey=True)


    # Deep d/dt(DO) timeseries
    nwin = 30
    ax[0].set_title('Deep Layer\n{}-day Hanning Window'.format(nwin), size=14, loc='left', fontweight='bold')
    # format grid
    ax[0].tick_params(axis='x', labelrotation=30)
    loc = mdates.MonthLocator(interval=1)
    ax[0].xaxis.set_major_locator(loc)
    ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    # ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax[0].tick_params(axis='both', labelsize=12)
    # add zero line
    ax[0].axhline(0,0,100,color='k',linestyle=':')
    # add drawdown period
    ax[0].axvline(dates_local_daily[minday],0,12,color='grey')
    ax[0].axvline(dates_local_daily[maxday],0,12,color='grey')
    # loop through inlets
    for i,inlet in enumerate(inlets):
        # get average deep layer DO
        deep_lay_DO_alltime = DOconcen_dict[inlet]['Deep Layer DO']
        # # 30-day hanning window
        # deep_lay_DO_alltime = helper_functions.lowpass(deep_lay_DO_alltime.values,n=30)
        # # plot
        # ax[1].plot(dates_local_daily,deep_lay_DO_alltime,linewidth=1,color='navy',alpha=0.5)

        # calculate slope and plot slope
        d_dt_DO = np.diff(deep_lay_DO_alltime.values)
        # filter
        d_dt_DO = helper_functions.lowpass(d_dt_DO,n=nwin)
        # plot
        ax[0].plot(dates_local_daily[1::],d_dt_DO,linewidth=1, color='navy', alpha=0.5)

    # format labels
    ax[0].set_xlim([dates_local_hrly[0],dates_local_hrly[-2]])
    ax[0].set_ylim([-0.1,0.2])
    ax[0].set_ylabel('d/dt(DO) [mg/L per day]',fontsize=14)
    # plt.tight_layout()
    # plt.show()
    # ax[0].legend(loc='best',ncols=2,fontsize=10)




    ########################################################################

    ax[1].set_title('Entire inlet\n{}-day Hanning Window'.format(nwin), size=14, loc='left', fontweight='bold')
    # format grid
    ax[1].tick_params(axis='x', labelrotation=30)
    loc = mdates.MonthLocator(interval=1)
    ax[1].xaxis.set_major_locator(loc)
    ax[1].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    # ax[1].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax[1].tick_params(axis='both', labelsize=12)
    # add zero line
    ax[1].axhline(0,0,100,color='k',linestyle=':')
    # add drawdown period
    ax[1].axvline(dates_local_daily[minday],0,12,color='grey')
    ax[1].axvline(dates_local_daily[maxday],0,12,color='grey')
    # loop through inlets
    for i,inlet in enumerate(inlets):
        
        d_dt_DO = deeplay_dict[inlet]['d/dt(DO)'] + shallowlay_dict[inlet]['d/dt(DO)']
        total_volume = dimensions_dict[inlet]['Inlet volume'][0]
        storage_daily = d_dt_DO / total_volume # kmol O2 /s /m3
        
        # convert to mg/L per day
        storage_daily = storage_daily.values * 1000 * 32 * 60 * 60 * 24

        # filter
        d_dt_DO = helper_functions.lowpass(storage_daily,n=nwin)
        # plot
        ax[1].plot(dates_local_daily,d_dt_DO,linewidth=1,color='navy', alpha=0.5)

    # format labels
    ax[1].set_xlim([dates_local_hrly[0],dates_local_hrly[-2]])
    ax[1].set_ylim([-0.1,0.2])
    # ax[1].set_ylabel('d/dt(DO) [mg/L per day]',fontsize=14)
    plt.tight_layout()
    plt.show()
    # plt.legend(loc='best',ncols=4)
    
    return