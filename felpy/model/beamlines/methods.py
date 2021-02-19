# -*- coding: utf-8 -*-
""" 
some methods for calling beamlines and parameters etc.
"""

from felpy.model.beamlines.structure import BeamlineModel

def get_beamline_object(params = "", options = 'nano', ekev = 5.0):
    
    """ 
    return desired beamline

    note, params var currently has no use, at later date, it would be nice to
    be able to build the beamline from this file.
    
    this may or may not be worth the time, for now it is almost certainly not.
    """
 

    spb = BeamlineModel()
    params = spb.params
    
    
    spb.mirror_profiles(toggle = "on", aperture = True, overwrite = False)
    
    mirrors = ['HOM1', 'HOM2', 'NHE', 'NVE']
    
    for mirror in mirrors:
    
        if mirror in ['HOM1', 'HOM2']:    
            spb.adjust_mirror(mirror,
                              ekev,
                              params[mirror]['ang_max'])
        else:    
            spb.adjust_mirror(mirror,
                              ekev,
                              3.1e-03)
        
    spb.build_elements(focus = 'nano')
    spb.build_beamline(focus = 'nano')


    return spb.get_beamline()
    

def unit_test():
    get_beamline_object()

if __name__ == '__main__':
    unit_test()