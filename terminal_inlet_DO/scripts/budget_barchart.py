"""
Plots example budget for Lynch Cove
and plots barcharts of budget terms for hypoxic and oxygenated inlets
during the drawdown period (mid-Jul through mid-Aug).

Also conducts Welch's t-test to test whether biological drawdown rate
or net decrease rates are different between hypoxic and oxygenated inlets.
"""
import numpy as np
import matplotlib.pylab as plt
import matplotlib.dates as mdates
from scipy.stats import ttest_ind
import helper_functions

def budget_barchart(inlets,shallowlay_dict,deeplay_dict,
                    dates_local_hrly,dates_local_daily,hyp_inlets,
                    minday,maxday,kmolm3sec_to_mgLday): 

    # initialize figure
    fig, ax = plt.subplots(4,1,figsize=(9.1,9.5))

    ##########################################################
    ##   Panel (a): Lynch Cove example budget time series   ##
    ##########################################################

    # lynch cove example budget
    inlet = 'lynchcove'

    # format figure
    ax[0].text(0.02, 0.88,'(a) Lynch Cove',fontsize=12, fontweight='bold',transform=ax[0].transAxes,)
    ax[0].set_xlim([dates_local_hrly[0],dates_local_hrly[-25]])
    ax[0].set_ylabel('DO transport ' + r'[kmol O$_2$ s$^{-1}$]',size=10)
    ax[0].grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
    ax[0].tick_params(axis='x', labelrotation=30, labelsize=10)
    ax[0].tick_params(axis='y', labelsize=10)
    loc = mdates.MonthLocator(interval=1)
    ax[0].xaxis.set_major_locator(loc)
    ax[0].xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax[0].set_ylim([-0.23,0.18])

    # plot deep budget time series
    nwin = 10 # hanning window length
    ax[0].plot(dates_local_daily,helper_functions.lowpass(deeplay_dict[inlet]['d/dt(DO)'].values,n=nwin),color='k',
                linewidth=1,label=r'$\frac{d}{dt}\int_V$DO dV') #,alpha=0.6)
    ax[0].plot(dates_local_daily,helper_functions.lowpass(deeplay_dict[inlet]['Vertical Transport'].values + shallowlay_dict[inlet]['Vertical Transport'].values,n=nwin),
                color='darkorange', linewidth=1,label='Error')
    ax[0].plot(dates_local_daily,helper_functions.lowpass(deeplay_dict[inlet]['TEF Exchange Flow'].values,n=nwin),color='#0D4B91',
            linewidth=2,label='Exchange Flow')
    ax[0].plot(dates_local_daily,helper_functions.lowpass(deeplay_dict[inlet]['Vertical Transport'].values,n=nwin),color='#99C5F7',
            linewidth=2,label='Vertical')
    ax[0].plot(dates_local_daily,helper_functions.lowpass(deeplay_dict[inlet]['Photosynthesis'].values,n=nwin),color='#8F0445',
                linewidth=2, label='Photosynthesis')
    ax[0].plot(dates_local_daily,helper_functions.lowpass(deeplay_dict[inlet]['Bio Consumption'].values,n=nwin),color='#FCC2DD',
                linewidth=2,label='Consumption')
    ax[0].legend(loc='lower right',ncol=6, fontsize=9, handletextpad=0.15)

    # add drawdown period
    # mid Jul - mid Aug
    ax[0].axvline(dates_local_daily[minday],0,12,color='grey')
    ax[0].axvline(dates_local_daily[maxday],0,12,color='grey')

    ##########################################################
    ##    Panel (b): Lynch Cove example budget bar chart    ##
    ##########################################################
    #  
    ax[1].tick_params(axis='x', labelrotation=30)
    ax[1].axhline(y=0, xmin=-0.5, xmax=1.05,color='silver',linewidth=1,linestyle='--')
    ax[1].tick_params(axis='y', labelsize=10)
    ax[1].set_xticklabels([])
    ax[1].set_ylabel('mg/L per day',fontsize=10)
    ax[1].text(0.02, 0.88,'(b) Lynch Cove',fontsize=12, fontweight='bold',transform=ax[1].transAxes,)
    ax[1].set_xlim([-0.5,1.05])
    ax[1].set_ylim([-0.25,0.25])
    # set bar width
    width = 0.2
    # create bar chart
    for attribute, measurement in deeplay_dict[inlet].items():
        # skip variables we are not interested in
        if attribute in ['WWTPs',
                        'Exchange Flow & Vertical',
                        'Photosynthesis & Consumption',
                        'Volume',
                        'Qin m3/s']:
            continue
        # assign colors
        if attribute == 'TEF Exchange Flow':
            color = '#0D4B91'
            label = 'Exchange Flow'
            pos = 0.2
        if attribute == 'Vertical Transport':
            color = '#99C5F7'
            label = 'Vertical'
            pos = 0.4
        if attribute == 'Photosynthesis':
            color = '#8F0445'
            label = attribute
            pos = 0.7
        if attribute == 'Bio Consumption':
            color = '#FCC2DD'
            label = 'Consumption'
            pos = 0.9
        if attribute == 'd/dt(DO)':
            color = 'black'
            label = r'$\frac{d}{dt}$DO (net decrease)'
            pos = -0.2

        # calculate time average
        time_avg = np.nanmean(measurement[minday:maxday])
        # get volume average
        avg = time_avg/(np.nanmean(deeplay_dict[inlet]['Volume'][minday:maxday])) # kmol O2 /s /m3
        # convert to mg/L per day
        avg = avg * kmolm3sec_to_mgLday

        # plot bars
        ax[1].bar(pos, avg, width, zorder=5, align='center', edgecolor=color,color=color, label=label)
        if avg < 0:
            wiggle = 0.05
        if avg > 0:
            wiggle = -0.05
        ax[1].text(pos, wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                color='black',fontsize=10)

        ax[1].legend(bbox_to_anchor=(0.5, -0.3), loc='lower center', fontsize=9, ncol=5, handletextpad=0.15)

    ##########################################################
    ##   Panel (c): distinct physical and biological terms  ## 
    ##########################################################

    # set bar width
    width = 0.1

    # counters for spacing bars
    multiplier_deep1 = 0
    multiplier_deep2 = 0

    # format grid
    for axis in [ax[2],ax[3]]:
        axis.tick_params(axis='x', labelrotation=30)
        axis.axhline(y=0, xmin=-0.5, xmax=1.05,color='silver',linewidth=1,linestyle='--')
        axis.tick_params(axis='y', labelsize=10)
        axis.set_xticklabels([])
        axis.set_ylabel('mg/L per day',fontsize=10)
        axis.set_xlim([-0.5,1.05])
    ax[2].set_ylim([-2.5,2.5])
    ax[3].set_ylim([-0.35,0.25])

    # create a new dictionary of results
    oxy_dict = {}
    hyp_dict = {}
    oxy_dict_std = {}
    hyp_dict_std = {}

    for inlet in inlets:
        for attribute, measurement in deeplay_dict[inlet].items():
            # skip variables we are not interested in
            if attribute in ['WWTPs',
                            'Exchange Flow & Vertical',
                            'Photosynthesis & Consumption',
                            'Volume',
                            'Qin m3/s']:
                continue
            # calculate time average normalized by volume
            avg = np.nanmean(measurement[minday:maxday]/(deeplay_dict[inlet]['Volume'][minday:maxday])) # kmol O2 /s /m3
            # convert to mg/L per day
            avg = avg * kmolm3sec_to_mgLday

            # save values in dictionary
            if inlet in hyp_inlets:
                if attribute in hyp_dict.keys():
                    hyp_dict[attribute].append(avg)
                else:
                    hyp_dict[attribute] = [avg]
            else:
                if attribute in oxy_dict.keys():
                    oxy_dict[attribute].append(avg)
                else:
                    oxy_dict[attribute] = [avg]

    # t-test
    print('\n=====================Welch\'s t-test=======================\n')
    print('p < 0.05 means we reject null hypothesis\n\n')
    for attribute in oxy_dict:
        if attribute == 'd/dt(DO)':
            print('Null hypothesis: d/dt(DO) of hypoxic and oxygenated inlets is the same')
            print('    alpha = 0.05')
            a = oxy_dict[attribute]
            b = hyp_dict[attribute]
            ttest = ttest_ind(a, b, axis=0, equal_var=False)
            print('    '+str(ttest))
        else:
            continue
    print('\n')
    for i,dict in enumerate([oxy_dict,hyp_dict]):
    # average all oxygenated and hypoxic inlet rate values
        if i ==0:
            shift = 0
        else:
            shift = 0.1
        for attribute, measurement in dict.items():
            # choose color
            if attribute == 'TEF Exchange Flow':
                color = '#0D4B91'
                label = 'Exchange Flow'
                pos = 0.15
            if attribute == 'Vertical Transport':
                color = '#99C5F7'
                label = 'Vertical'
                pos = 0.35
            if attribute == 'Photosynthesis':
                color = '#8F0445'
                label = attribute
                pos = 0.65
            if attribute == 'Bio Consumption':
                color = '#FCC2DD'
                label = 'Consumption'
                pos = 0.85
            if attribute == 'd/dt(DO)':
                color = 'black'
                label = r'$\frac{d}{dt}$DO (net decrease)'
                pos = -0.25

            # calculate average
            avg = np.nanmean(measurement)

            if avg < 0:
                wiggle = 0.65
                ha = 'center'
            if avg > 0:
                wiggle = -0.7
                ha = 'center'
            # plot
            offset = width * multiplier_deep1
            if attribute == 'Storage':
                hatchcolor = 'white'
            else:
                hatchcolor = 'white'
            if i == 0:
                rects = ax[2].bar(pos + shift, avg, width, zorder=5, align='center', edgecolor=hatchcolor,color=color, hatch='xx')
                rects = ax[2].bar(pos + shift, avg, width, zorder=5, align='center', edgecolor=color,color='none')
                if attribute == 'd/dt(DO)':
                    ax[2].text(pos + shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
                            color=color,fontsize=10, fontweight='bold', rotation=45)
                else:
                    ax[2].text(pos + shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
                            color='gray',fontsize=10, rotation=45)
                multiplier_deep1 += 2
            elif i == 1:
                offset = width * multiplier_deep2
                rects = ax[2].bar(pos + shift, avg, width, zorder=5, align='center', edgecolor=color,color=color, label=label)
                if attribute == 'd/dt(DO)':
                    ax[2].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
                            color=color,fontsize=10, fontweight='bold', rotation=45)
                else:
                    ax[2].text(pos + shift, 0+wiggle, str(round(avg,3)),horizontalalignment=ha,verticalalignment='center',
                            color='gray',fontsize=10, rotation=45)
                multiplier_deep2 += 2
            

    ##########################################################
    ##   Panel (d): combined physical and biological terms  ## 
    ##########################################################

    # set bar width
    width = 0.2

    # counters for spacing bars
    multiplier_deep1 = 0
    multiplier_deep2 = 0

    # create a new dictionary of results
    oxy_dict = {}
    hyp_dict = {}

    for inlet in inlets:
        for attribute, measurement in deeplay_dict[inlet].items():
            # skip variables we are not interested in
            if attribute in ['TEF Exchange Flow',
                            'WWTPs',
                            'Vertical Transport',
                            'Photosynthesis',
                            'Bio Consumption',
                            'Volume',
                            'Qin m3/s']:
                continue
                # calculate time average normalized by volume
            avg = np.nanmean(measurement[minday:maxday]/(deeplay_dict[inlet]['Volume'][minday:maxday])) # kmol O2 /s /m3
            # convert to mg/L per day
            avg = avg * kmolm3sec_to_mgLday

            # save values in dictionary
            if inlet in ['penn','case','holmes','portsusan','lynchcove','dabob']:
                if attribute in hyp_dict.keys():
                    hyp_dict[attribute].append(avg)
                else:
                    hyp_dict[attribute] = [avg]
            else:
                if attribute in oxy_dict.keys():
                    oxy_dict[attribute].append(avg)
                else:
                    oxy_dict[attribute] = [avg]

    # t-test
    for attribute in oxy_dict:
        if attribute == 'Photosynthesis & Consumption':
            print('Null hypothesis: biological drawdown of hypoxic and oxygenated inlets is the same')
            print('    alpha = 0.05')
            a = oxy_dict[attribute]
            b = hyp_dict[attribute]
            ttest = ttest_ind(a, b, axis=0, equal_var=False)
            print('    '+str(ttest))
        else:
            continue

    for i,dict in enumerate([oxy_dict,hyp_dict]):
    # average all oxygenated and hypoxic inlet rate values
        if i ==0:
            shift = 0
        else:
            shift = 0.2
        for attribute, measurement in dict.items():
            # choose color
            if attribute == 'Exchange Flow & Vertical':
                color = '#488DDB'
                label = 'Exchange Flow & Vertical Transport'
                pos = 0.2
            if attribute == 'Photosynthesis & Consumption':
                color = '#F069A8'
                label = attribute + ' (drawdown)'
                pos = 0.7
            if attribute == 'd/dt(DO)':
                color = 'black'
                label = r'$\frac{d}{dt}$DO'
                pos = -0.3
        
            # calculate average
            avg = np.nanmean(measurement)

            if avg < 0:
                wiggle = 0.04
            if avg > 0:
                wiggle = -0.04
            # plot
            offset = width * multiplier_deep1
            hatchcolor = 'white'
            if i == 0:
                rects = ax[3].bar(pos + shift, avg, width, zorder=5, edgecolor=hatchcolor,color=color, hatch='xx')
                rects = ax[3].bar(pos + shift, avg, width, zorder=5, edgecolor=color,color='none')
                if attribute == 'd/dt(DO)':
                    ax[3].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color=color,fontsize=10, fontweight='bold')
                else:
                    ax[3].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color='gray',fontsize=10)
                multiplier_deep1 += 2
            elif i == 1:
                offset = width * multiplier_deep2
                rects = ax[3].bar(pos + shift, avg, width, zorder=5, edgecolor=color,color=color,label=label)
                if attribute == 'd/dt(DO)':
                    ax[3].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color=color,fontsize=10, fontweight='bold')
                else:
                    ax[3].text(pos+shift, 0+wiggle, str(round(avg,3)),horizontalalignment='center',verticalalignment='center',
                            color='gray',fontsize=10)
                multiplier_deep2 += 2

    ax[3].legend(loc='lower center', fontsize=9, ncol=3, handletextpad=0.15)

    ax[2].text(0.3, 0.12, 'HATCHED: oxygenated (n={})\nSOLID: hypoxic (n={})      '.format(len(oxy_dict['d/dt(DO)']),len(hyp_dict['d/dt(DO)'])),
            color='black', verticalalignment='bottom', horizontalalignment='right',zorder=6,
            transform=ax[2].transAxes, fontsize=9, fontweight='bold')

    ax[2].text(0.02, 0.88,'(c) All thirteen inlets',fontsize=12, fontweight='bold',transform=ax[2].transAxes,)
    ax[3].text(0.02, 0.88,'(d) All thirteen inlets',fontsize=12, fontweight='bold',transform=ax[3].transAxes,)

    plt.subplots_adjust(left=0.1, top=0.95, bottom=0.05, right=0.9, hspace=0.3)
    plt.show()

    return
