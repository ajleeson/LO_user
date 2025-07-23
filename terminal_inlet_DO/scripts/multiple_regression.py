"""
Calculate multiple linear regression
of DOdeep dependece on DOin and Tflush
"""
import numpy as np
from scipy.linalg import lstsq
from scipy.stats import pearsonr

def multiple_regression(MONTHLYmean_DOdeep,
                        MONTHLYmean_DOin,
                        MONTHLYmean_Tflush):
    

    print('\n=============================================================')
    print('==================Multiple Linear Regression=================')
    print('=============================================================\n')

    # create array of predictors
    input_array = np.array([MONTHLYmean_DOin, MONTHLYmean_Tflush, [1]*len(MONTHLYmean_DOin)]).T


    B,a,b,c = lstsq(input_array,MONTHLYmean_DOdeep)
    slope_DOin = B[0]
    slope_Tflush = B[1]
    intercept = B[2]


    # calculate r^2 and p value
    r,p = pearsonr(MONTHLYmean_DOin,MONTHLYmean_DOdeep)
    print('DO_deep dependence on DO_in')
    print('   r = {}'.format(round(r,3)))
    print('   R^2 = {}'.format(round((r**2),3)))
    print('   p = {:.2e}'.format(p))

    # calculate r^2 and p value
    r,p = pearsonr(MONTHLYmean_Tflush,MONTHLYmean_DOin-MONTHLYmean_DOdeep)
    print('\n(DO_in - DO_deep) dependence on T_flush')
    print('   r = {}'.format(round(r,3)))
    print('   R^2 = {}'.format(round((r**2),3)))
    print('   p = {:.2e}'.format(p))

    print('\nMean deep layer DO [mg/L] = {}*DOin + {}*Tflush + {}\n'.format(
        round(slope_DOin,2),round(slope_Tflush,2),round(intercept,2)))

    # calculate r^2 and p value
    predicted_DOdeep = slope_DOin * MONTHLYmean_DOin + slope_Tflush * MONTHLYmean_Tflush + intercept
    r,p = pearsonr(MONTHLYmean_DOdeep,predicted_DOdeep)
    print('DO_deep dependence on DO_in and T_flush')
    print('   r = {}'.format(round(r,3)))
    print('   R^2 = {}'.format(round((r**2),3)))
    print('   p = {:.2e}'.format(p))
    
    return