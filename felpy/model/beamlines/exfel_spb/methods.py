#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

from felpy.model.beamlines.exfel_spb.exfel_spb import Instrument

 
def get_beamline_object(params = "", options = 'nano', ekev = 5.0,
                        apertures = True, surface = True, crop = None,
                        theta_HOM = 2.3e-03, theta_KB = 3.5e-03):
    
    """ 
    return desired beamline

    note, params var currently has no use, at later date, it would be nice to
    be able to build the beamline from this file.
    
    this may or may not be worth the time, for now it is almost certainly not.
    """
 

    spb = Instrument()
    params = spb.params
    
    
    
    mirrors = ['HOM1', 'HOM2', 'NHE', 'NVE']
    
    for mirror in mirrors:
    
        if mirror in ['HOM1', 'HOM2']:    
            spb.adjust_mirror(mirror,
                              ekev,
                              theta_HOM)#params[mirror]['ang_max'])
        else:    
            spb.adjust_mirror(mirror,
                              ekev,
                              theta_KB)
    ## bad python / quick fix
    if surface == True: 
        surface = "on"
    else: 
        surface = "off"
        
    spb.mirror_profiles(toggle = surface, aperture = apertures, overwrite = False)

        
    spb.build_elements(focus = 'nano')
    spb.build_beamline(focus = 'nano')
    
    if crop is not None:
        spb.crop_beamline(crop[0], crop[1])


    return spb.get_beamline()
    

def unit_test():
    get_beamline_object()

if __name__ == '__main__':
    unit_test()