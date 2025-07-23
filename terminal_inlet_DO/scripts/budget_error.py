"""
calculate and print error of budget
expressed as a % of QinDOin and biological consumption
"""
import numpy as np

def budget_error(inlets,shallowlay_dict,deeplay_dict,
                 dimensions_dict,kmolm3sec_to_mgLday):

    print('\n=============================================================')
    print('========================Budget Error=========================')
    print('=============================================================\n')

    # initialize lists
    error_QinDOin_ann_avg = []
    error_consumption_ann_avg = []

    for inlet in inlets:

        # calculate budget error (mg/L per day)
        error_TEF = (shallowlay_dict[inlet]['Vertical Transport']+deeplay_dict[inlet]['Vertical Transport'])/ (
            dimensions_dict[inlet]['Inlet volume'].values) * kmolm3sec_to_mgLday
        inlet_error_ann_avg = np.nanmean(error_TEF)

        # calculate QinDOin (mg/L per day)
        QinDOin = (deeplay_dict[inlet]['TEF Exchange Flow'].values/dimensions_dict[inlet]['Inlet volume'].values) * kmolm3sec_to_mgLday
        inlet_QinDOin_ann_avg = np.nanmean(QinDOin)

        # calculate biological consumption in deep layer (mg/L per day)
        consumption = (deeplay_dict[inlet]['Bio Consumption'].values/dimensions_dict[inlet]['Inlet volume'].values) * kmolm3sec_to_mgLday
        inlet_consumption_ann_avg = np.nanmean(consumption)

        # add values to list
        error_QinDOin_ann_avg.append(inlet_error_ann_avg/inlet_QinDOin_ann_avg)
        error_consumption_ann_avg.append(inlet_error_ann_avg/inlet_consumption_ann_avg)
        
    # calculate bulk statistics
    error_QinDOin = np.abs(np.nanmean(error_QinDOin_ann_avg)) * 100
    error_consumption = np.abs(np.nanmean(error_consumption_ann_avg)) * 100

    # print bulk statistics
    print('(annual mean error)/(annual mean QinDOin) [expressed as percentage]')
    print('    {}%'.format(round(error_QinDOin,2)))
    print('\n')
    print('(annual mean error)/(annual mean deep consumption) [expressed as percentage]')
    print('    {}%'.format(round(error_consumption,2)))

    return