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

import numpy as np ### interesting, first numpy import as applied in [A]

from felpy.model.beamlines.exfel_spb.exfel_spb import Instrument

 
def get_beamline_object(params = "", options = 'nano', ekev = 5.0,
                        apertures = True, surface = 'real', crop = None,
                        theta_HOM = 2.3e-03, theta_KB = 3.5e-03,
                        save_params = True):
    
    """ 
    return desired beamline

    note, params var currently has no use, at later date, it would be nice to
    be able to build the beamline from this file.
    
    this may or may not be worth the time, for now it is almost certainly not.
    """
 

    spb = Instrument(parameter_file = "spb-sfx_nkb_FAST")
    params = spb.params
    
    if apertures == False:
        theta_HOM = theta_KB = np.pi/2 ### is true. no aperture / mirror edges == 
    
    
    mirrors = ['HOM1', 'HOM2', 'NHE', 'NVE', 'MHE', 'MHP', 'MVE', 'MVP']

    spb.mirror_profiles(surface = surface, aperture = apertures)

    for mirror in mirrors:
    
        if mirror in ['HOM1', 'HOM2']:    
            spb.adjust_mirror(mirror,
                              ekev,
                              theta_HOM)
            
        else:    
            spb.adjust_mirror(mirror,
                              ekev,
                              theta_KB)

    
        
    if apertures == False: ## aperture claus
       
        for mirror in mirrors:
            spb.params[mirror]["design angle"] = np.pi/3 ### [A]s
            spb.params[mirror]["incidence angle"] = np.pi/3
            

    if apertures == False:
        
        for focus in spb.focus:
        
            spb.params['NVE_error']["design angle"] = np.pi/2 ### [A]s ## fix for now, should onkly accept elliptical mirrors 
            spb.params['NVE_error']["incidence angle"] = np.pi/2 ### should be for elliptical mirror surfaces

            spb.params['NHE_error']["design angle"] = np.pi/2 ### [A]s ## fix for now, should onkly accept elliptical mirrors 
            spb.params['NHE_error']["incidence angle"] = np.pi/2 ### should be for elliptical mirror surfaces
            

    spb.build_elements(focus = options)
    


    spb.build_beamline(focus = options)
    
    
    if save_params:
        spb.export_params()

        
    if crop is not None:
        if type(crop) == list:
            spb.crop_beamline(*crop)
        else:
            spb.crop_beamline(crop)


    return spb.get_beamline()
    

def setup_spb(params = "", options = 'nano', ekev = 5.0,
              apertures = True, surface = 'real', crop = None,
              theta_HOM = 2.3e-03, theta_KB = 3.5e-03,
              save_params = False):

    """ 
    return desired beamline

    note, params var currently has no use, at later date, it would be nice to
    be able to build the beamline from this file.

    this may or may not be worth the time, for now it is almost certainly not.
    """


    spb = Instrument()
    params = spb.params

    if apertures == False:
        theta_HOM = theta_KB = np.pi/2 ### is true. no aperture / mirror edges == 

    mirrors = spb.mirrors

    spb.mirror_profiles(surface = surface, aperture = apertures)

    for mirror in mirrors:

        if mirror in spb.focus:    
            spb.adjust_mirror(mirror,
            ekev,
            theta_KB)

        else:    
            spb.adjust_mirror(mirror,
            ekev,
            theta_HOM)



    if apertures == False: ## aperture claus

        for mirror in mirrors:
            spb.params[mirror]["design angle"] = np.pi/3 ### [A]s
            spb.params[mirror]["incidence angle"] = np.pi/3


    ### TO-DO: make a choice as to wether edits to the beamline will be via params or via beamline object 
    ### TO-DO: eventually, exfel_spb should be a sub-class of the Instrument class (formally)
    ### TO-DO: a front-end script to load and label mirrors that ships to a json file would be useful.

    if apertures == False:

        for focus in spb.focus:

            spb.params['NVE_error']["design angle"] = np.pi/2 ### [A]s ## fix for now, should onkly accept elliptical mirrors 
            spb.params['NVE_error']["incidence angle"] = np.pi/2 ### should be for elliptical mirror surfaces

            spb.params['NHE_error']["design angle"] = np.pi/2 ### [A]s ## fix for now, should onkly accept elliptical mirrors 
            spb.params['NHE_error']["incidence angle"] = np.pi/2 ### should be for elliptical mirror surfaces


    spb.build_elements(focus = options)
    spb.build_beamline(focus = options)


    if save_params:
        spb.export_params()


    if crop is not None:
        if type(crop) == list:
            spb.crop_beamline(*crop)
        else:
            spb.crop_beamline(crop)


    return spb

    

if __name__ == '__main__':
    from wpg.wpg_uti_wf import plot_intensity_map
    from felpy.model.src.coherent import construct_SA1_wavefront

    spb = setup_spb()
    bl = spb.bl
    wfr = construct_SA1_wavefront(512,512, 5, 0.25)
    bl.propagate_sequential(wfr)
    plot_intensity_map(wfr)