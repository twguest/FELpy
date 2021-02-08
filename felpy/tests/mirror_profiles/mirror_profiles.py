# -*- coding: utf-8 -*-
from felpy.model.beamlines.structure import BeamlineModel

def core():
    
    bl = BeamlineModel()
 
    
    mirrors = ["HOM1", "HOM2", "NVE", "NHE"]
    for mirror in mirrors:
        bl.plot_mirror_profiles(mirror_name = mirror, sdir = "")
        
if __name__ == '__main__':
    core()