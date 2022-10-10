# -*- coding: utf-8 -*-
from felpy.model.instrument import Instrument

def plot_image_surfaces():
    
    bl = Instrument()
 
    
    mirrors = ["HOM1", "HOM2", "NVE", "NHE"]
    for mirror in mirrors:
        bl.plot_mirror_profiles(mirror_name = mirror)
        
if __name__ == '__main__':
    plot_image_surfaces()