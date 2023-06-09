import numpy as np

def intensity_correlation(ii, t1):
    """
    calculate the one-time intensity correlation for a single time-point t1
    
    :param ii: time-ensemble intensity distribution, I(t_1,k)
    :param t1: time to calculate correlation
    :returns: one-time correlation
    """

    a = np.ones(ii.shape[0])

    for itr in range(ii.shape[0]):

        a[itr] = np.correlate(ii[t1,:],ii[itr,:])

    return a/np.max(a)


def hanbury_brown_twiss(ii):
    """
    calculate the two-time correlation function of an intensity array at a single point in space
    
    :param ii: time-ensemble intensity distribution, I(t_1,k)
    :returns cc: two-time correlation
    """
    cc = np.zeros([ii.shape[0], ii.shape[0]])

    for itr in range(ii.shape[0]):

            cc[itr, :] = intensity_correlation(ii,itr)


    for itr in range(cc.shape[0]):

        cc[:,itr] = cc[itr,:]

    return cc

