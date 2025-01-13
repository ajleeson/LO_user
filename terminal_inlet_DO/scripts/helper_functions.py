"""
This module contains utility functions forfiltering.

Parker MacCready
"""

import numpy as np

def lowpass(data, f='hanning', n=40, nanpad=True):
    """
    A replacement for almost all previous filter code.
    f = 'hanning' (default) or 'godin'
    
    Input: ND numpy array, any number of dimensions, with time on axis 0.
    
    Output: Array of the same size, filtered with Hanning window of length n,
        or the Godin filter (hourly data only) padded with nan's.
    """
    if n == 1:
        return data
    else:
        if f == 'hanning':
            filt = hanning_shape(n=n)
        elif f == 'godin':
            filt = godin_shape()
        else:
            print('ERROR in filt_general(): unsupported filter ' + f)
            filt = np.nan
        npad = np.floor(len(filt)/2).astype(int)
        sh = data.shape
        df = data.flatten('F')
        dfs = np.convolve(df, filt, mode = 'same')
        smooth = dfs.reshape(sh, order='F')
        # note that the indexing below defaults to being on axis 0,
        # and correctly broadcasts without having to mention the other axes
        if nanpad:
            smooth[:npad] = np.nan
            smooth[-npad:] = np.nan
        else:
            smooth[:npad] = data[:npad]
            smooth[-npad:] = data[-npad:]
        return smooth
    
def godin_shape():
    """
    Based on matlab code of 4/8/2013  Parker MacCready
    Returns a 71 element numpy array that is the weights
    for the Godin 24-24-25 tildal averaging filter. This is the shape given in
    Emery and Thomson (1997) Eqn. (5.10.37)
    ** use ONLY with hourly data! **
    """
    k = np.arange(12)
    filt = np.NaN * np.ones(71)
    filt[35:47] = (0.5/(24*24*25))*(1200-(12-k)*(13-k)-(12+k)*(13+k))
    k = np.arange(12,36)
    filt[47:71] = (0.5/(24*24*25))*(36-k)*(37-k)
    filt[:35] = filt[:35:-1]
    return filt

def hanning_shape(n=40):
    """
    Returns a Hanning window of the specified length.
    """
    ff = np.cos(np.linspace(-np.pi,np.pi,n+2))[1:-1]
    filt = (1 + ff)/2
    filt = filt / filt.sum()
    return filt