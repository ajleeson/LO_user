"""
calculate and print correlation between
mean biological consumption and QinDOin
for 13 terminal inlets during the drawdown period
(mid-Jul through mid-Aug)
"""
import numpy as np
from scipy.stats import pearsonr

def correl(inlets,deeplay_dict,minday,maxday,kmolm3sec_to_mgLday):
    
    print('\n=============================================================')
    print('=========Correlation between QinDOin and Consumption=========')
    print('=============================================================\n')

    print('UPDATE THIS TO INSTEAD BE DOin - DOdeep vs. TFLUSH CORRELATION,\n SINCE I REFERENCE THAT NUMBER IN THE MANUSCRIPT')

    exchange = np.array([])
    consumption = np.array([])

    for i,inlet in enumerate(inlets):
        
        # get QinDOin during drawdown period
        # calculate time average normalized by volume
        tef =  np.nanmean(deeplay_dict[inlet]['TEF Exchange Flow'][minday:maxday]/
                        (deeplay_dict[inlet]['Volume'][minday:maxday])) # kmol O2 /s /m3
        # convert to mg/L per day
        tef = tef * kmolm3sec_to_mgLday

        # get biological consumption during drawdown period
        # calculate time average normalized by volume
        cons =  np.nanmean(deeplay_dict[inlet]['Bio Consumption'][minday:maxday]/
                        (deeplay_dict[inlet]['Volume'][minday:maxday])) # kmol O2 /s /m3
        # convert to mg/L per day
        cons = cons * kmolm3sec_to_mgLday

        # add to array
        exchange = np.concatenate((exchange,[tef]))
        consumption = np.concatenate((consumption,[cons]))

    # calculate correlation coefficient (Pearson)
    r,p = pearsonr(exchange,consumption)
    print('correlation during drawdown period:')
    print('    r = {}'.format(round(r,3)))
    print('    p = {}'.format(round(p,3)))

    return