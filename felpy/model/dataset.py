# -*- coding: utf-8 -*-
import os
import numpy as np
import h5py




class Dataset:
    """
    """
    
    def __init__(self):
        setattr(self, "test", 100)
    
    def load_hdf5(self):
        pass
    
    def store_hdf5(self, filename):
        ### in developmnt 
                    
        if filename.endswith(".h5") or filename.endswith(".hdf5"): 
            pass
        else:
            filename += ".h5"
            
            
        if os.path.exists(filename):
            raise Warning("File already exists.")
                   
        with h5py.File(filename, 'w') as h5file:
            
            for key, item in self.__dict__.items():
                    if self.__dict__.item is dict:
                        pass
                        
                    else:
                        h5file[filename + key] = item        
            

    
    @property
    def __nested__(self):
        return False
        
        
def Shimadzu(Dataset):
    

    
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    @property
    def __type__(self):
        pass

if __name__ == '__main__':
    
    ds = Dataset()
    a = dict()
    setattr(ds, "test", a)
    ds.store_hdf5("a")
    