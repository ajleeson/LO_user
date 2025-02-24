"""
Plots boxplots of net decrease of oxygen 
in all thirteen terminal inlets
during the drawdown period (mid-Jul through mid-Aug).

Also conducts Welch's ANOVA test to test whether all inlets
have the same mean net decrease of oxgyen during the
drawdown period.
"""
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
import pingouin as pg

def net_decrease_boxplots(dimensions_dict,deeplay_dict,
                        shallowlay_dict,minday,maxday):
    

    print('\n=====================Welch\'s ANOVA (1layer)=======================\n')

    print('\nNull hypothesis: all inlets have same net decrease rate')
    print('p < 0.05 means we reject null hypothesis\n\n')

    # initialize figure
    fig, ax = plt.subplots(1,1,figsize = (10,5.5))

    # format figure
    ax.tick_params(axis='x', labelrotation=30)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_ylabel(r'$d/dt$DO [mg/L per day]', fontsize=12)
    ax.set_ylim([-0.55,0.4])

    # initialize lists to store net decrease rates
    storage_all = []
    storage_mean = []

    # mean_depth_sorted = dict(sorted(dimensions_dict.items(), key=lambda item: item[1]))
    # print(mean_depth_sorted)

    stations_sorted = ['sinclair','quartermaster','dyes',
                    'crescent','penn','case',
                    'lynchcove','carr','holmes',
                    'portsusan','elliot','commencement',
                    'dabob']

    for i,station in enumerate(stations_sorted):
        
        # get daily net decrease rate
        d_dt_DO = deeplay_dict[station]['d/dt(DO)'][minday:maxday] + shallowlay_dict[station]['d/dt(DO)'][minday:maxday]
        total_volume = dimensions_dict[station]['Inlet volume'][0] #deeplay_dict[station]['Volume'][minday:maxday] + shallowlay_dict[station]['Volume'][minday:maxday]
        storage_daily = d_dt_DO / total_volume # kmol O2 /s /m3
        
        # convert to mg/L per day
        storage_daily = storage_daily.values * 1000 * 32 * 60 * 60 * 24

        # add to array
        storage_all.append(list(storage_daily))
        storage_mean.append(np.nanmean(storage_daily))


    # create boxplot
    ax.axhline(y=0, xmin=-0.5, xmax=1.05,color='silver',linewidth=1,linestyle='--')
    bplot = plt.boxplot(storage_all, patch_artist=True, labels=stations_sorted,
                showmeans=True, showfliers=False, boxprops={'color':'darkgray'}, meanprops=
                {'marker': 'o', 'markerfacecolor': 'navy', 'markersize': 7, 'markeredgecolor': 'none'})

    # align x-axis labels
    plt.xticks(ha='right')

    # add mean depth
    for i,station in enumerate(stations_sorted):
        interval = 0.077
        ax.text(i*interval+interval/2, 0.03, str(round(dimensions_dict[station]['Mean depth'][0])) ,color='black',
                            horizontalalignment='center',transform=ax.transAxes, fontsize=10)
    ax.text(interval/3, 0.08, 'Mean depths [m]:' ,color='black', fontweight='bold',
                            horizontalalignment='left',transform=ax.transAxes, fontsize=10)

    # format boxplot
    for patch in bplot['boxes']:
        patch.set_facecolor('whitesmoke')
    for element in ['whiskers', 'medians', 'caps']:
            plt.setp(bplot[element], color='darkgray')

    # add mean of all inlets
    ax.axhline(y=np.nanmean(storage_mean), xmin=-0.5, xmax=1.05,color='deepskyblue',linewidth=2,
            label='mean: {} mg/L per day'.format(round(np.nanmean(storage_mean),3)))
    ax.legend(loc='best',frameon=False,fontsize=12)

    # condudct anova test
    flat_storage = [x for xs in storage_all for x in xs]
    repeats = maxday - minday
    df = pd.DataFrame({'storage':flat_storage,
                    'inlet': np.repeat(stations_sorted, repeats=repeats)})
    pingu = pg.welch_anova(data=df, dv='storage', between='inlet')
    print('Welch\'s ANOVA test of all thirteen inlets')
    print('    p = {:.2e}'.format((pingu['p-unc'].values[0])))
    print('\n')

    # remove Dabob Bay from analysis and repeat ANOVA
    df_nodabob = df[df['inlet'] != 'dabob']
    pingu_nodabob = pg.welch_anova(data=df_nodabob, dv='storage', between='inlet')
    print('OMITTING DABOB BAY: Welch\'s ANOVA test of all other twelve inlets')
    print('    p = {}'.format(round(pingu_nodabob['p-unc'].values[0],3)))
    print('\n')

    plt.title('d/dt(DO) in entire inlet (one-layer)')
    plt.tight_layout()
    plt.show()
    
    return