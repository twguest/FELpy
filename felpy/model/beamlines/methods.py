#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "Apache"
__version__ = "1.0.0"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

from felpy.model.beamlines.structure import BeamlineModel

def get_beamline_object(params = "", options = 'nano', ekev = 5.0,
                        theta_hom = 2.2e-03, theta_kb = 3.5e-03):
    
    """ 
    return desired beamline

    note, params var currently has no use, at later date, it would be nice to
    be able to build the beamline from this file.
    
    this may or may not be worth the time, for now it is almost certainly not.
    """
    bl = BeamlineModel()
    
    bl.setupHOMs(ekev, theta_hom)
    bl.setupKBs(ekev, theta_kb)
    
    bl.mirrorProfiles(toggle = "on", aperture = True, overwrite = False)
    bl.buildElements(focus = options)
    bl.buildBeamline(focus = options)
    return bl.get_beamline()
    

def unit_test():
    get_beamline_object()

if __name__ == '__main__':
    unit_test()