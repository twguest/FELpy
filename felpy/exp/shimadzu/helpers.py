# -*- coding: utf-8 -*-

        
def shimadzu_test_data(nx, ny, npulses, ntrains, weight = 0.1,
                       processed = False):
    """
    generate data matching the structure of 'treated' shimadzu test data
    
    :param nx: number of horizontal pixels
    :param ny: number of vertical pixels
    :param npulses: number of pulses per train
    :param ntrains: number of trains
    :param weight: weighting of random modulation (1 leads to unccorelated fields)
     
    :returns data: numpy array of shape [nx, ny, npulses, ntrains]
    """
    
    data = np.ones([nx, ny, npulses, ntrains])*complex_gaussian_2d(nx,ny)
    data += np.random.rand(nx, ny npulses, ntrains)*weight
    
    return data
    
    
    

