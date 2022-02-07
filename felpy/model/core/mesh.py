# -*- coding: utf-8 -*-
import numpy as np

D2 = ["nx", "ny", "xMin", "xMax", "yMin", "yMax"]
D3 = ["nz", "zMin", "zMax"]
      
class Mesh:
    """
    generalised Mesh class
    
    is constructed via two-methods: from single-value field definitions or from
    numpy arrays defining exact coordinates.
    
    **kwargs
    :params nx: (optional) [int]
    :params ny: (optional) [int]
    :params nz: (optional) [int]

    :params xMin: (optional) [float64]
    :params xMax: (optional) [float64]

    :params yMin: (optional) [float64]
    :params yMax: (optional) [float64]

    :params zMin: (optional) [float64]
    :params zMax: (optional) [float64]

    :params x: (optional) [np.ndarray]
    :params y: (optional) [np.ndarray]
    :params z: (optional) [np.ndarray]
    """
  
    def __init__(self, **kwargs):

        self.ndims = 0
        self.__dict__.update(kwargs)
        
               
        
        if all(hasattr(self, attr) for attr in D2):
            
            self.ndims = 2
        
            if all(hasattr(self, attr) for attr in D3):
                self.ndims = 3
               
        if all(hasattr(self, attr) for attr in ["x", "y"]):
            
            self.ndims = 2
         
            if all(hasattr(self, attr) for attr in ["z"]):
                
                self.ndims = 3
        
            self.build_from_array()
        
            
        self.set_attributes()
        self.check_attributes()
    
    def __str__(self):
        print("")
        print("Mesh Type Object: ")
        print("")
        print("Spatial Parameters:")
        for key in self.D:
            print("{}: {}".format(key, self.__dict__[key]))
        return ""
 
    
    def check_attributes(self):
        missing = []
        
        for attr in self.D:
            
            if hasattr(self, attr) == False:
                missing.append(attr)
        

        if len(missing) > 0:
            raise Warning("You do not have the required attributes: {}".format(missing))
        
        
        del missing
        
    def set_attributes(self):
        """
        this sets the dimensions labels for attribute calling later
        
        note: that for now, the minumum number of dimensions is 2 -
        so D2 is set to be the default value for D
        """
        
        self.D = D2
        
        if self.ndims >= 2:
            self.D = D2
            
            if self.ndims == 3:
                self.D += D3
    
    def sampling(self):
        """
        this returns and records the field sampling in each of the spatial
        dimensions
        """
        if self.ndims >= 2:
            self.dx = (self.xMax-self.xMin)/self.nx
            self.dy = (self.yMax-self.yMin)/self.ny
            
            if self.ndims == 3:
                self.dz = (self.zMax-self.zMin)/self.nz
        
    def get_sampling(self):
        """ 
        return the two or three dimensioanl sampling of the field.
        """
        if self.ndims >= 2:
            return self.dx, self.dy
        elif self.ndims == 3:
            return self.dx, self.dy, self.dz
        
    def build_from_array(self):
        """ 
        construct a mesh object from linear-space variables
        """
        
        if self.ndims >= 2:
                
    
            self.nx = len(self.x)
            self.xMin = np.min(self.x)
            self.xMax = np.max(self.x)
            self.ny = len(self.y)
            self.yMin = np.min(self.y)
            self.yMax = np.max(self.y)
            
            if self.ndims == 3:
                self.nz = len(self.z)
                self.zMin = np.min(self.z)
                self.zMax = np.max(self.z)
    

    def get_attributes(self):
        """
        returns the necessary spatial attributes in dictionary form
        """
        return {key: self.__dict__[key] for key in self.D}
    
    def get_extent(self):
        """
        returns wpg style extent with respect to the transverse axis
        """
        return [self.xMin, self.xMax, self.yMin, self.yMax]
    
    def get_array(self, axes = 0):
        """
        return the mesh in array form
        """
        self.check_axes_options(axes)
        
        if axes in ['x', 0]:
            return np.linspace(self.xMin, self.xMax, self.nx)
        if axes in ['y', 1]:
            return np.linspace(self.yMin, self.yMax, self.ny)
        if axes in ['z', 2]:
            return np.linspace(self.zMin, self.zMax, self.nz)
    
    
    def check_axes_options(self, axes):
        """
        a sanity function to check the axis options.
        """
        options = ['x','y','z',0,1,2]
    
        if axes not in options:
            raise Warning("We're not sure what axes you require")
    
    def __update__(self, **kwargs):
        """
        update the attributes of the mesh via various methods
        """

        if "wfr" in kwargs:
    
            for o in [x for x in list(self.__dict__.keys()) if x in dir(kwargs['wfr'].params.Mesh)]:
                setattr(self, o, getattr(kwargs['wfr'].params.Mesh, o))
                
    
if __name__ == '__main__':
    m = Mesh(x = np.linspace(10,10,10), y = np.linspace(10,10,10))
    print(m)