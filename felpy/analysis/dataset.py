#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 17:53:52 2021

A generalised dataset class for storing and manipulating XFEL data, both 
simulated and experimental.

@author: twguest
"""
import datetime
import os
import numpy as np
import h5py

from felpy.utils.np_utils import get_mesh

def create_random_h5_file(folder, index):
    """create one random file"""
    name = os.path.join(folder, 'myfile_' + str(index) + ".h5")
    with h5py.File(name=name, mode='w') as f:
        d = f.create_dataset('data', (5, 10), 'i4')
        data = np.random.randint(low=0, high=100, size=(5*10))
        data = data.reshape(5, 10)
        d[:] = data
        f.close()
        
    return name

class Dataset:
    
    def __init__(self, **kwargs):
        """
        
        :param dataset: the data for analysis/storage. 
                        Here we expec
                
        :param data_dim: [tuple] dimensionality of data (note: this does notequal data.ndims, this is an explanation of
        the structure of the data. ie, an array of
        shape [nx,ny,a,b] would have an order [nx,ny])
                
        :param set_dim: [tuple] set_dim refers to the inverse of data_dimensions,
                    i.e. what is the order of the data_set. e.g
                    an array of shape [nx,ny,a,b] would have a data_depth
                    [a,b]
        """
        
        for key in kwargs:
            self.key = kwargs[key]
 
    
    @staticmethod
    def load_dataset(fdir, overwrite = True):
        """
        load a h5 dataset to python object struture w/ methods etc
        
        :param fdir: directory of h5 file
        """
        
 
        
        ds = loaded_dataset()
        
        if overwrite:
            ds.fdir = fdir
            ds.write_h5(fdir)

        with h5py.File(fdir, 'r') as f:
            
            for k in f.keys():
                
                exec(f'ds.{k} = f[k][()]')
       
        

        return ds 
    
    
    
    @staticmethod
    def create_virtual_dataset(fdir, files, key):
        """ 
        construct a virtual dataset, containing multiple h5 files
        
        :param fdir: location of virtual dataset
        :param files: list of files (strings) to be added to virtual dataset
        :param key: key from .h5 files to add to virtual dataset
        """
        sh = h5py.File(files[0], 'r')[key].shape  # get the first ones shape.
        layout = h5py.VirtualLayout(shape=(len(files),) + sh,
                                    dtype=np.float64)
        with h5py.File(fdir, 'w', libver='latest') as f:
            f.create_dataset("index", files)
            for i, filename in enumerate(files):
                vsource = h5py.VirtualSource(filename, key, shape=sh)
                layout[i,] = vsource
            
            f.create_virtual_dataset(key, layout, fillvalue=0)
            f.close()
            
    def add_virtual_dataset(self, files, key = None):
        pass
    
    def create_group(self, group_name):
        """
        add a subset dataset to this dataset 
        """
        self.f.create_group(group_name)
    
    def create_subgroup(self, group_name, sub_group):
        self.f[group_name].create_group(sub_group)
        
         
    def label_dataset(self, data_labels, set_labels):
        """
        add dataset labels for plotting purposes etc
        
        :param data_labels: data labels (list of strings)
        :param set_labels: set labels (list of strings)
        """
        
        assert len(data_labels) == len(self.data_dim)
        assert len(set_labels) == len(self.set_dim)
        
        self.data_labels = data_labels
        self.set_labels = set_labels
        
        
    def write_h5(self, fdir):
        """
        :param fdir: location + filename of the h5 file
        """
        
        self.fdir = fdir
        
        with h5py.File(fdir, 'w') as f:
              
            for item in self.__dict__:

                    f.create_dataset(item, data = self.__dict__[item])
        
        self.f = h5py.File(name=self.fdir, mode='w') 

    def close(self):
        self.f.close()
        
        
    def warnings(self):
        assert self.data_dim + self.set_dim == self.dataset.shape, "The dimensions of the data and its containing set should be representative of the data shape. i.e. data order + data depth = data shape.\n\nSee Dataset docstring for more info"

    def enable_analysis(self):
        """
        creates a group in the global h5 folder to store analysis files
        """
        assert fdir is not None, "please write h5 to file before enabling analysis"
        
        with h5py.File(self.fdir) as f:
            f.create_group("analysis")
    
    def set_mesh(self, px,py):
        self.mesh = get_mesh(self.dataset, px, py)
    
        
class loaded_dataset(Dataset):
    """ 
    A lazy Dataset object for loading h5 arrays
    """
    def __init__(self):
        self.last_loaded = str(datetime.datetime.now())


class Shimadzu(Dataset):
    """
    a data storage class specific to the shimadzu detector
    """

    def __init__(self, dataset, fdir = None):
        
        data_dim = dataset.shape[:2]
        set_dim = dataset.shape[2:]
        enable_analysis = False
    
        
        super().__init__(dataset, data_dim, set_dim, enable_analysis, fdir)
        
        self.label_dataset(["x", "y"],["Pulse", "Train"])
        
        
        
 
if __name__ == '__main__':
    
    fdir = "../data/tmp/dataset_tester.h5"

    
    from felpy.analysis.dataset import Dataset
    from felpy.exp.shimadzu.preprocess import shimadzu_test_data
 
# =============================================================================
#     
#     data = shimadzu_test_data(512, 512, 5, 10)
#     
#     ds = Dataset()
#     
#     import numpy as np
#     import h5py
#     
#     ### write dataset object to h5 file.
#     ds.write_h5(fdir)
#     
#     ### read dataset object from h5 file.
#     dr = Dataset.load_dataset(fdir)
#     
#     sh = Shimadzu(data)
#  
# =============================================================================
