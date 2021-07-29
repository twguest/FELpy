# -*- coding: utf-8 -*-

from felpy.analysis.dataset import Dataset
from felpy.exp.shimadzu.preprocess import shimadzu_test_data



if __name__ == '__main__':
    
    data = shimadzu_test_data(512, 512, 5, 10)
    
    Dataset(data, data_order = (512,512), data_depth = (5,11))
    
    import numpy as np
    import h5py
    
    
 