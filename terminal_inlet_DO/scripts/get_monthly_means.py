"""
calculate and return monthly mean
DOdeep, DOin, and Tflush 
for all terminal inlets

outputs:
    MONTHLYmean_XXXX are arrays of monthly mean values
    for all inlets, compressed into a single array

    df_MONTHLY_mean_XXX are dataframes, where each column
    is an individual inlet. All columns contain monthly
    mean values corresponding to the inlet (ie., 12 rows)
"""
import numpy as np
import pandas as pd

def get_monthly_means(deeplay_dict,DOconcen_dict,
                      dimensions_dict,inlets):

    # values for looping
    months = ['Jan','Feb','Mar',
            'Apr','May','Jun',
            'Jul','Aug','Sep',
            'Oct','Nov','Dec']
    intervals = len(months)

    # initialize arrays to store monthly mean values for all inlets
    MONTHLYmean_DOdeep = np.zeros(len(inlets)*intervals)
    MONTHLYmean_DOin = np.zeros(len(inlets)*intervals)
    MONTHLYmean_Tflush = np.zeros(len(inlets)*intervals)
    MONTHLYmean_perchyp = np.zeros(len(inlets)*intervals)

    # initialize dataframes to store monthly mean values for individual inlets
    df_MONTHLYmean_DOdeep = pd.DataFrame(columns=inlets)
    df_MONTHLYmean_DOin = pd.DataFrame(columns=inlets)
    df_MONTHLYmean_Tflush = pd.DataFrame(columns=inlets)
    df_MONTHLYmean_perchyp = pd.DataFrame(columns=inlets)

    # calculate monthly mean DOin, DOdeep, and Tflush in each inlet
    for i,inlet in enumerate(inlets):

        # initialize temporary array of monthly means for one inlet
        inlet_mean_DOdeep = []
        inlet_mean_DOin = []
        inlet_mean_Tflush = []
        inlet_mean_perchyp = []
        
        # Note that data extends from Jan 02 through Dec 31
        # So index = 0 corresponds to Jan 02
        # The data indices have been adjusts to align with the
        # start and end of every month, considering our date range.
        for month_index,month in enumerate(months):
            if month == 'Jan':
                MONTHminday = 0 
                MONTHmaxday = 30 
            elif month == 'Feb': 
                MONTHminday = 30
                MONTHmaxday = 58
            elif month == 'Mar':
                MONTHminday = 58
                MONTHmaxday = 89
            elif month == 'Apr':
                MONTHminday = 89
                MONTHmaxday = 119
            elif month == 'May':
                MONTHminday = 119
                MONTHmaxday = 150
            elif month == 'Jun':
                MONTHminday = 150
                MONTHmaxday = 180
            elif month == 'Jul':
                MONTHminday = 180
                MONTHmaxday = 211
            elif month == 'Aug':
                MONTHminday = 211
                MONTHmaxday = 242
            elif month == 'Sep': 
                MONTHminday = 242
                MONTHmaxday = 272
            elif month == 'Oct':
                MONTHminday = 272
                MONTHmaxday = 303
            elif month == 'Nov':
                MONTHminday = 303
                MONTHmaxday = 332
            elif month == 'Dec':
                MONTHminday = 332
                MONTHmaxday = 363

            # calculate monthly means
            mean_DOdeep = np.nanmean(DOconcen_dict[inlet]['Deep Layer DO'][MONTHminday:MONTHmaxday]) # mg/L
            mean_DOin = np.nanmean(DOconcen_dict[inlet]['DOin'][MONTHminday:MONTHmaxday]) # mg/L
            mean_Tflush = np.nanmean(dimensions_dict[inlet]['Inlet volume'][0]/deeplay_dict[inlet]['Qin m3/s'][MONTHminday:MONTHmaxday]) / (60*60*24) # days
            mean_perc_hyp_vol = np.nanmean(DOconcen_dict[inlet]['percent hypoxic volume'][MONTHminday:MONTHmaxday]) # percent

            # save values in arrays for all inlets
            MONTHLYmean_DOdeep[i*intervals+month_index] =  mean_DOdeep
            MONTHLYmean_DOin[i*intervals+month_index] = mean_DOin
            MONTHLYmean_Tflush[i*intervals+month_index] = mean_Tflush
            MONTHLYmean_perchyp[i*intervals+month_index] = mean_perc_hyp_vol

            # save values in temporary array for individual inlet
            inlet_mean_DOdeep.append(mean_DOdeep)
            inlet_mean_DOin.append(mean_DOin)
            inlet_mean_Tflush.append(mean_Tflush)
            inlet_mean_perchyp.append(mean_perc_hyp_vol)

        # save values in dataframes for individual inlets
        df_MONTHLYmean_DOdeep[inlet] = inlet_mean_DOdeep
        df_MONTHLYmean_DOin[inlet] = inlet_mean_DOin
        df_MONTHLYmean_Tflush[inlet] = inlet_mean_Tflush
        df_MONTHLYmean_perchyp[inlet] = inlet_mean_perchyp
    
    return [MONTHLYmean_DOdeep,
            MONTHLYmean_DOin,
            MONTHLYmean_Tflush,
            MONTHLYmean_perchyp,
            df_MONTHLYmean_DOdeep,
            df_MONTHLYmean_DOin,
            df_MONTHLYmean_Tflush,
            df_MONTHLYmean_perchyp]