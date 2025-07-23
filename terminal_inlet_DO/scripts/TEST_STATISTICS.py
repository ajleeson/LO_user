"""
Plots example budget for Lynch Cove
and plots barcharts of budget terms for hypoxic and oxygenated inlets
during the drawdown period (mid-Jul through mid-Aug).

Also conducts Welch's t-test to test whether biological drawdown rate
or net decrease rates are different between hypoxic and oxygenated inlets.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf
from scipy.stats import ttest_ind
from scipy.stats import shapiro
from scipy.stats import bartlett
from scipy.stats import kruskal
import pingouin as pg
import helper_functions

def test_statistics(inlets,shallowlay_dict,deeplay_dict,
                    dates_local_hrly,dates_local_daily,hyp_inlets,
                    minday,maxday,kmolm3sec_to_mgLday): 
    
    print('\n===========================================================')
    print('====================Testing Statistics=====================')
    print('===========================================================\n')

    # some random inlet
    inlet = 'lynchcove'

    for attribute, measurement in deeplay_dict[inlet].items():
        # skip variables we are not interested in
        if attribute in ['WWTPs',
                        'Exchange Flow & Vertical',
                        'Photosynthesis & Consumption',
                        'Volume',
                        'Qin m3/s']:
            continue

        # calculate time average
        time_avg = np.nanmean(measurement[minday:maxday])
        # get volume average
        avg = time_avg/(np.nanmean(deeplay_dict[inlet]['Volume'][minday:maxday])) # kmol O2 /s /m3
        # convert to mg/L per day
        avg = avg * kmolm3sec_to_mgLday

    # ##########################################################
    # ##   Panel (c): distinct physical and biological terms  ## 
    # ##########################################################


    # create a new dictionary of results
    oxy_dict = {}
    hyp_dict = {}
    
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
    print('Looking at d/dt(DO) for oxygenated and hypoxic inlets...')
    for attribute in oxy_dict:
        if attribute == 'd/dt(DO)':
            
            # print('    alpha = 0.05')
            a = oxy_dict[attribute]
            b = hyp_dict[attribute]
            print('    First checking that mean d/dt(DO) rates of oxygenated and hypoxic')
            print('    inlet groups are normally distributed')
            print('      Shapiro-Wilk test (p < 0.05 means data are NOT normally distributed)')
            stat,shapiro_test_oxy_p = shapiro(a)
            stat,shapiro_test_hyp_p = shapiro(b)
            print('        p = {} for oxygenated inlets'.format(round(shapiro_test_oxy_p,3)))
            print('        p = {} for hypoxic inlets'.format(round(shapiro_test_hyp_p,3)))
            print('        => Data are normally distributed\n')

            # Perform Bartlett's test
            statistic, p_value = bartlett(a, b)

            print(f"Bartlett's test statistic: {statistic}")
            print(f"P-value: {p_value}")

            # Interpret the results
            alpha = 0.05
            if p_value < alpha:
                print("Reject the null hypothesis: Variances are significantly different.")
            else:
                print("Fail to reject the null hypothesis: Variances are not significantly different.")

            print('\nNull hypothesis: d/dt(DO) of hypoxic and oxygenated inlets is the same')
            print('p < 0.05 means we reject null hypothesis')
            ttest,p_value = ttest_ind(a, b, axis=0, equal_var=False)
            print('    p = {}'.format(round(p_value,3)))
        else:
            continue

     # t-test
    print('\n=====================Welch\'s t-test=======================\n')
    print('Looking at biology and physics for oxygenated and hypoxic inlets...')

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


    for attribute in oxy_dict:
        if attribute == 'Photosynthesis & Consumption':
            
            # print('    alpha = 0.05')
            a = oxy_dict[attribute]
            b = hyp_dict[attribute]
            print('    First checking that mean drawdown rates of oxygenated and hypoxic')
            print('    inlet groups are normally distributed')
            print('      Shapiro-Wilk test (p < 0.05 means data are NOT normally distributed)')
            stat,shapiro_test_oxy_p = shapiro(a)
            stat,shapiro_test_hyp_p = shapiro(b)
            print('        p = {} for oxygenated inlets'.format(round(shapiro_test_oxy_p,3)))
            print('        p = {} for hypoxic inlets'.format(round(shapiro_test_hyp_p,3)))
            print('        => Data are normally distributed\n')

            # Perform Bartlett's test
            statistic, p_value = bartlett(a, b)

            print(f"Bartlett's test statistic: {statistic}")
            print(f"P-value: {p_value}")

            # Interpret the results
            alpha = 0.05
            if p_value < alpha:
                print("Reject the null hypothesis: Variances are significantly different.")
            else:
                print("Fail to reject the null hypothesis: Variances are not significantly different.")

            print('\nNull hypothesis: drawdown of hypoxic and oxygenated inlets is the same')
            print('p < 0.05 means we reject null hypothesis')
            ttest,p_value = ttest_ind(a, b, axis=0, equal_var=False)
            print('    p = {}'.format(round(p_value,3)))
        else:
            continue

    print('--------------------------------')

    for attribute in oxy_dict:
        print(attribute)
        if attribute == 'Exchange Flow & Vertical':
            
            # print('    alpha = 0.05')
            a = oxy_dict[attribute]
            b = hyp_dict[attribute]
            print('    First checking that mean physical rates of oxygenated and hypoxic')
            print('    inlet groups are normally distributed')
            print('      Shapiro-Wilk test (p < 0.05 means data are NOT normally distributed)')
            stat,shapiro_test_oxy_p = shapiro(a)
            stat,shapiro_test_hyp_p = shapiro(b)
            print('        p = {} for oxygenated inlets'.format(round(shapiro_test_oxy_p,3)))
            print('        p = {} for hypoxic inlets'.format(round(shapiro_test_hyp_p,3)))
            print('        => Data are normally distributed\n')

            # Perform Bartlett's test
            statistic, p_value = bartlett(a, b)

            print(f"Bartlett's test statistic: {statistic}")
            print(f"P-value: {p_value}")

            # Interpret the results
            alpha = 0.05
            if p_value < alpha:
                print("Reject the null hypothesis: Variances are significantly different.")
            else:
                print("Fail to reject the null hypothesis: Variances are not significantly different.")

            print('\nNull hypothesis: physical of hypoxic and oxygenated inlets is the same')
            print('p < 0.05 means we reject null hypothesis')
            ttest,p_value = ttest_ind(a, b, axis=0, equal_var=False)
            print('    p = {}'.format(round(p_value,3)))
        else:
            continue

    # print('\n=====================Welch\'s ANOVA=======================\n')

    # print('\nNull hypothesis: all inlets have same net decrease rate')
    # print('p < 0.05 means we reject null hypothesis\n\n')

    # # initialize lists to store net decrease rates
    # storage_all = []
    # storage_mean = []

    # # mean_depth_sorted = dict(sorted(dimensions_dict.items(), key=lambda item: item[1]))
    # # print(mean_depth_sorted)

    # stations_sorted = ['sinclair','quartermaster','dyes',
    #                 'crescent','penn','case',
    #                 'lynchcove','carr','holmes',
    #                 'portsusan','elliot','commencement',
    #                 'dabob']

    # for i,station in enumerate(stations_sorted):
        
    #     # get daily net decrease rate
    #     storage_daily =  deeplay_dict[station]['d/dt(DO)'][minday:maxday]/(deeplay_dict[station]['Volume'][minday:maxday]) # kmol O2 /s /m3
        
    #     # convert to mg/L per day
    #     storage_daily = storage_daily.values * 1000 * 32 * 60 * 60 * 24

    #     # # subsample to weekly values
    #     # storage_weekly = storage_daily[::7]

    #     # # subsample to every three days
    #     # storage_daily = storage_daily[::3]

    #     # # Aggregate to weekly means
    #     # # Assume storage_daily is your 1D array with 61 daily values
    #     # n = len(storage_daily)
    #     # block_size = 3
    #     # # Truncate to full weeks only (optional, drops leftover days)
    #     # n_blocks = n // block_size
    #     # trimmed = storage_daily[:n_blocks * block_size]
    #     # # Reshape into (weeks, days_per_week) and compute means
    #     # storage_daily = trimmed.reshape(n_blocks, block_size).mean(axis=1)

    #     # add to array
    #     storage_all.append(list(storage_daily))
    #     storage_mean.append(np.nanmean(storage_daily))
    #     # storage_all.append(list(storage_weekly))
    #     # storage_mean.append(np.nanmean(storage_weekly))

    # # loop through all inlets adn check for normality
    # print('    First checking that daily d/dt(DO) rates of inlets are normally distributed')
    # print('      Shapiro-Wilk test (p < 0.05 means data are NOT normally distributed)')
    # for i,inlet_storage in enumerate(storage_all):
    #     stat,shapiro_test_daily = shapiro(inlet_storage)
    #     print('        p = {} for {}'.format(round(shapiro_test_daily,3),stations_sorted[i]))
    # print('        => Not all inlets have normally distributed daily d/dt(DO) rates!\n')

    # print('Kruskal-Wallis H-test of all thirteen inlets (weekly subsampling to reduce autocorrelation)')
    # stat,kruskal_p = kruskal(storage_all[0],
    #                          storage_all[1],
    #                          storage_all[2],
    #                          storage_all[3],
    #                          storage_all[4],
    #                          storage_all[5],
    #                          storage_all[6],
    #                          storage_all[7],
    #                          storage_all[8],
    #                          storage_all[9],
    #                          storage_all[10],
    #                          storage_all[11],
    #                          storage_all[12])
    # print('    p = {}'.format(kruskal_p))
    # print('\n')
    

    # print('---------------------------')
    # print('Checking autocorrelation of weekly data')
    
    # print(len(storage_all))
    # print(len(storage_all[0]))
    
    # fig, axes = plt.subplots(3,5,figsize=(10,8),sharex=True,sharey=True)
    # fig.supylabel('Autocorrelation')
    # fig.supxlabel('Lag (days)')
    # fig.suptitle('Autocorrelation of daily d/dt(DO) during drawdown period')
    # ax = axes.ravel()
    # for i,inlet_storage in enumerate(storage_all):
    #     # fig, ax = plt.subplots(1,1,figsize = (6,6))
    #     # pd.plotting.lag_plot(pd.Series(inlet_storage), lag=1)
    #     plot_acf(np.asarray(inlet_storage), title=stations_sorted[i], ax=ax[i])
    #     # ax.set_title(stations_sorted[i])


    # print('---------------------------')
    # print('Permutation test')

    # # # condudct anova test
    # # flat_storage = [x for xs in storage_all for x in xs]
    # # repeats = maxday - minday
    # # repeats = len(storage_all[0]) # for subsampled
    # # df = pd.DataFrame({'storage':flat_storage,
    # #                    'inlet': np.repeat(stations_sorted, repeats=repeats)})
    
    # # # Compute observed F-statistic
    # # obs = pg.anova(data=df, dv='storage', between='inlet', detailed=True)
    # # obs_F = obs.loc[0, 'F']
    # # n_perm = 10000
    # # perm_Fs = np.zeros(n_perm)
    # # for i in range(n_perm):
    # #     shuffled_df = df.copy()
    # #     shuffled_df['inlet'] = np.random.permutation(shuffled_df['inlet'])  # row-wise shuffle
    # #     result = pg.anova(data=shuffled_df, dv='storage', between='inlet', detailed=True)
    # #     perm_Fs[i] = result.loc[0, 'F']
    # # # Two-sided permutation p-value
    # # p_perm = np.mean(perm_Fs >= obs_F)
    # # print(f'Observed F: {obs_F:.3f}')
    # # print(f'Permutation-based p-value: {p_perm:.4f}')
    

    # # pingu = pg.welch_anova(data=df, dv='storage', between='inlet')
    # # print('Welch\'s ANOVA test of all thirteen inlets')
    # # print('    p = {:.2e}'.format((pingu['p-unc'].values[0])))
    # # print('\n')

    # # # remove Dabob Bay from analysis and repeat ANOVA
    # # df_nodabob = df[df['inlet'] != 'dabob']
    # # pingu_nodabob = pg.welch_anova(data=df_nodabob, dv='storage', between='inlet')
    # # print('OMITTING DABOB BAY: Welch\'s ANOVA test of all other twelve inlets')
    # # print('    p = {}'.format(round(pingu_nodabob['p-unc'].values[0],3)))
    # # print('\n')


    # return
