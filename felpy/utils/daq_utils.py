# -*- coding: utf-8 -*-

"""
utilities to help with daq output at the EuXFEL. 
to be used in conjunction with extra/karabo_data
"""
from karabo_data import RunDirectory

def load_data(run, proposal, exp, VERBOSE = False):
    
    ddir = "/gpfs/exfel/exp/SPB/{}/{}/raw/{}/".format(exp, proposal, run)
    
    print("Loading data from: {}".format(ddir))
    
    data = RunDirectory(ddir)
    
    if VERBOSE: 
        data.info()

    return data

def shimadzu_reshape(arr, VERBOSE = False):
    """
    a verbose method of a describing theshape of the intensity data 
    """
    
    arr = arr.T
    
    if VERBOSE:
        print("# x pixels: {}".format(arr.shape[0]))
        print("# y pixels: {}".format(arr.shape[1]))
        print("# pulses per train: {}".format(arr.shape[2]))
        print("# trains: {}".format(arr.shape[3]))
    
    return arr




if __name__ == '__main__':
    pass