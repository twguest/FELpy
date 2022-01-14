#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

NOTE: the final goal of this class is to inherit from core.instrument.

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "1.0.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"


### TO-DO: Add a smart front-end that reads elements into the Instrument and saves names/properties etc
### TO-DO: each optical element should be an object with properties describing how to build it etc.
"""


import os
from os.path import exists
import json
import numpy as np
import sys
from matplotlib import pyplot as plt
import os
from felpy.model.core.beamline import Beamline
from wpg import srwlib
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture
from wpg.srwlib import SRWLOptT
from wpg.optical_elements import Mirror_elliptical as MirEl

from felpy.utils.os_utils import add_path, felpy_path
from felpy.model.materials.mirror_surface import genMirrorSurface, generate_mirror_surface
from felpy.model.materials.load_refl import get_refl, load_refl
from felpy.model.source.coherent import construct_SA1_wavefront
from wpg.optical_elements import calculateOPD
from felpy.model.beamlines.exfel_spb.params import get_params
from felpy.model.tools import propagation_parameters
from wpg.srwlib import srwl_opt_setup_surf_height_2d as MirPl

import os
import sys
 
import seaborn as sns

class Instrument:
    
    """
    A container for loading and modifying beamline data
    """
    
    
    def __init__(self, parameter_file = None, VERBOSE = True):
        
        self.VERBOSE = VERBOSE
        
        if VERBOSE:
            print("Initialising Single Particle Beamline")
        
        self.load_params(file = parameter_file)
        self.fpath = felpy_path() ### felpy path (for dev. purposes)
        
        self.mirrors = ["HOM1", "HOM2", "MHE", "MVE", "MHP", "MVP", "NVE", "NHE"]
        self.nano = ["HOM1", "HOM2", "NVE", "NHE"]
        self.focus = ["MHE", "MVE", "NVE", "NHE"]
        add_path()
        
    def load_params(self, file = None):
        """
        load beamline parameters from /data/spb/
        """
        
        if file is not None:
            with open(file, "r") as read_file:
                self.params = json.load(read_file)
        else:
            self.params = get_params()
            
    def export_params(self, sdir = None):
        """
        save the current set of beamline parameters in json format
        
        :param sdir: save directory
        """
        if sdir is None:
            with open(self.fpath + '/data/spb/parameters.json', 'w') as f:
                json.dump(self.params, f)
        else:
            with open(sdir + 'parameters.json', 'w') as f:
                json.dump(self.params, f)

        
        
    def adjust_mirror(self, mirror_name, ekev, new_ang, mirror_refl = None):
     
        if mirror_refl is None: 
            if ekev >= 7.5:
                material = "B4C"
            else: 
                material = "Ru"    
        
            refl = get_refl(load_refl(material), ekev, new_ang)
       
        new_ang = new_ang + np.tan(self.params[mirror_name]["xc"]/self.params[self.params[mirror_name]['next_drift']]['distance'])
       
        self.params[mirror_name]["design angle"] = new_ang
        self.params[mirror_name]["incidence angle"] = new_ang
        self.params[mirror_name]['reflectivity'] = (1-refl)#**2 ### note, this maps to an absorption parameter (hence 1-refl)



    def load_mirror_profiles(self, surface = "flat", aperture = True):
        
        
        for mirror in self.nano:
            
            if mirror in ["MHE", "MVP","MVE","MHP"]:
                fdir =  self.fpath + "/data/spb/mirror_surface/{}_mir_{}.dat".format(mirror,surface)
            else:
                fdir =  self.fpath + "/data/spb/mirror_surface/{}_mir_real.dat".format(mirror)
                
            if aperture:
                
                if os.path.exists(fdir):
                    if mirror in self.focus:
                        self.params[mirror+"_error"]['mirror profile'] = fdir 

                    else:
                        self.params[mirror]['mirror profile'] = fdir 

                else:
                    print("The mirror path could not be found")
                    print("Generating Random Mirror File")
                    generate_mirror_surface(512, 512,
                                           dx = self.params[mirror]['dx'],
                                           dy = self.params[mirror]['dy'],
                                           savedir = "../../data/spb/mirror_surface/",
                                           mode = surface,
                                           mirror_name = mirror)

                    self.params[mirror]['mirror profile'] = fdir 
                    
            elif aperture == False:
                
                ## ask liuba "RuntimeError: Failed to determine optical axis after reflection from mirror." for large el-mirror lengths

                if os.path.exists(fdir):
                    if mirror in self.focus:
                        self.params[mirror+"_error"]['mirror profile'] = fdir 
                    else:
                        self.params[mirror]['mirror profile'] = fdir 
                else:
                    generate_mirror_surface(512, 512,
                                           dx = 100,
                                           dy = 100,
                                           savedir = "../../data/spb/mirror_surface/",
                                           mode = surface,
                                           mirror_name = mirror)

                
 
                    self.params[mirror]['mirror profile'] = fdir 
 
    def build_elements(self, focus = "nano"):
        
        self.d1 =  Drift(self.params["HOM1"]['distance from source'])
        self.d1.name = self.params["d1"]['name']
        
        self.HOM1 = MirPl(np.loadtxt(self.params['HOM1']['mirror profile'].replace("../../","")),
                     _dim = self.params['HOM1']['orientation'],
                     _ang = self.params['HOM1']['incidence angle'], 
                     _amp_coef = 1,
                     _x = self.params['HOM1']['xc'], _y = self.params['HOM1']['yc'],
                     _refl = self.params['HOM1']['reflectivity']) 
        
        self.HOM1.name = self.params['HOM1']['name']

        self.d2 =  Drift(self.params["HOM2"]['distance from source']-self.params["HOM1"]['distance from source'])
        self.d2.name = self.params["d2"]['name']
        
        self.HOM2 = MirPl(np.loadtxt(self.params['HOM2']['mirror profile'].replace("../../","")),
                     _dim = self.params['HOM2']['orientation'],
                     _ang = self.params['HOM2']['incidence angle'], 
                     _amp_coef = 1,
                     _x = self.params['HOM2']['xc'], _y = self.params['HOM2']['yc'],
                     _refl = self.params['HOM2']['reflectivity']) 
        
        self.HOM2.name = self.params['HOM2']['name']
        
        if focus == "direct":
                    
            self.params["d3"]["distance"] = 656.424
            self.params["d4"]["distance"] = 1.20
            self.params["d5"]["distance"] = 1.00
            self.params["df"]["distance"] = 2.2
            
            self.d3 =  Drift(self.params["d3"]['distance'])
            self.d3.name = self.params["d3"]['name']
            
            self.d4 =  Drift(self.params["d4"]['distance'])
            self.d4.name = self.params["d4"]['name']
            
            self.d5 =  Drift(self.params["d5"]['distance'])
            self.d5.name = self.params["d5"]['name']
            
            self.df =  Drift(self.params["df"]['distance'])
            self.df.name =  'df'
            
            self.NKB_PSlit = Aperture(_shape=self.params["NKB_PSlit"]['shape'],
                                      _ap_or_ob=self.params["NKB_PSlit"]['type'],
                                      _Dx= self.params["NKB_PSlit"]['dx'],
                                      _Dy= self.params["NKB_PSlit"]['dy'],
                                      _x=self.params["NKB_PSlit"]['xc'],
                                      _y=self.params["NKB_PSlit"]['yc'])
            self.NKB_PSlit.name = self.params["NKB_PSlit"]['name']
            
        if focus == "micron":
            self.d3 =  Drift(self.params["d3"]['distance'])
            self.d3.name = self.params["d3"]['name']
                
            self.MKB_pslit = Aperture(_shape=self.params["MKB_pslit"]['shape'],
                                 _ap_or_ob=self.params["MKB_pslit"]['type'],
                                 _Dx= self.params["MKB_pslit"]['dx'],
                                 _Dy= self.params["MKB_pslit"]['dy'],
                                 _x=self.params["MKB_pslit"]['xc'],
                                 _y=self.params["MKB_pslit"]['yc'])
            self.MKB_pslit.name = self.params["MKB_pslit"]['name']
            
            self.d4 =  Drift(self.params["d4"]['distance'])
            self.d4.name = self.params["d4"]['name']
            

            self.MHP = MirPl(np.loadtxt(self.params['MHP']['mirror profile'].replace("../../","")),
                         _dim = self.params['MHP']['orientation'],
                         _ang = self.params['MHP']['incidence angle'], 
                         _refl = self.params['MHP']['transmission'],
                         _x = self.params['MHP']['xc'], _y = self.params['MHP']['yc']) 
            self.MHP.name = self.params['MHP']['name']
            
            self.d5 =  Drift(self.params["d5"]['distance'])
            self.d5.name = self.params["d5"]['name']
            
            self.MHE_error = MirPl(np.loadtxt(self.params['MHE_error']['mirror profile'].replace("../../","")),
             _dim = self.params['MHE_error']['orientation'],
             _ang = self.params['MHE_error']['incidence angle'], 
             _refl = self.params['MHE_error']['transmission'],
             _x = self.params['MHE_error']['xc'], _y = self.params['MHE_error']['yc']) 
            
            self.MHE_error.name = self.params['MHE_error']['name']
            

            self.MVE_error = MirPl(np.loadtxt(self.params['MVE_error']['mirror profile'].replace("../../","")),
             _dim = self.params['MVE_error']['orientation'],
             _ang = self.params['MVE_error']['incidence angle'], 
             _refl = self.params['MVE_error']['transmission'],
             _x = self.params['MVE_error']['xc'], _y = self.params['MVE_error']['yc']) 
            
            self.MVE_error.name = self.params['MVE_error']['name']
            
            self.MHE = MirEl(orient = self.params['MHE']["orientation"], p = self.params['MHE']["distance from source"], q = self.params['MHE']["distance to focus"],
                        thetaE = self.params['MHE']["design angle"], theta0 = self.params['MHE']["incidence angle"],
                        _x = self.params["MHE"]["xc"],
                        _y = self.params["MHE"]["yc"],
                        length = self.params['MHE']["length"],
                        roll = self.params['MHE']["roll"],
                        yaw = self.params['MHE']["yaw"],
                        _refl = self.params['MHE']["reflectivity"],
                        _ext_in = self.params['MHE']["_ext_in"], _ext_out = self.params['MHE']["_ext_out"]) 
            
            self.MHE.name = self.params['MHE']['name']
            
            self.d6 =  Drift(self.params["d7"]['distance'])
            self.d6.name = self.params["d7"]['name']
            
            self.MVE = MirEl(orient = self.params['MVE']["orientation"], p = self.params['MVE']["distance from source"], q = self.params['MVE']["distance to focus"],
                        thetaE = self.params['MVE']["design angle"], theta0 = self.params['MVE']["incidence angle"],
                        _x = self.params["MVE"]["xc"],
                        _y = self.params["MVE"]["yc"],
                        length = self.params['MVE']["length"],
                        roll = self.params['MVE']["roll"],
                        yaw = self.params['MVE']["yaw"],
                        _refl = self.params['MVE']["reflectivity"],
                        _ext_in = self.params['MVE']["_ext_in"], _ext_out = self.params['MVE']["_ext_out"]) 
                        
            self.MVE.name = self.params['MVE']['name']


            self.d8 =  Drift(self.params["d8"]['distance'])
            self.d8.name = self.params["d8"]['name']
            
            self.MVP = MirPl(np.loadtxt(self.params['MVP']['mirror profile'].replace("../../","")),
                     _dim = self.params['MVP']['orientation'],
                     _ang = self.params['MVP']['incidence angle'], 
                     _refl = self.params['MVP']['transmission'],
                     _x = self.params['MVP']['xc'], _y = self.params['MVP']['yc']) 
            self.MVP.name = self.params['MVP']['name']
            
            self.df =  Drift(self.params["df"]['distance'])
            self.df.name = self.params["df"]['name']
        
        elif focus == "nano":
            
            self.NKB_PSlit = Aperture(_shape=self.params["NKB_PSlit"]['shape'],
                                      _ap_or_ob=self.params["NKB_PSlit"]['type'],
                                      _Dx= self.params["NKB_PSlit"]['dx'],
                                      _Dy= self.params["NKB_PSlit"]['dy'],
                                      _x=self.params["NKB_PSlit"]['xc'],
                                      _y=self.params["NKB_PSlit"]['yc'])
            self.NKB_PSlit.name = self.params["NKB_PSlit"]['name']
            
            
            self.NHE = MirEl(orient = self.params['NHE']["orientation"], p = self.params['NHE']["distance from source"], q = self.params['NHE']["distance to focus"],
                    thetaE = self.params['NHE']["design angle"], theta0 = self.params['NHE']["incidence angle"],
                    _x = self.params["NHE"]["xc"],
                    _y = self.params["NHE"]["yc"],
                    length = self.params['NHE']["length"],
                    roll = self.params['NHE']["roll"],
                    yaw = self.params['NHE']["yaw"],
                    _refl = self.params['NHE']["reflectivity"],
                    _ext_in = self.params['NHE']["_ext_in"], _ext_out = self.params['NHE']["_ext_out"]) 
            
            self.NHE.name = self.params["NHE"]["name"]
            
            self.NVE = MirEl(orient = self.params['NVE']["orientation"], p = self.params['NVE']["distance from source"], q = self.params['NVE']["distance to focus"],
                    thetaE = self.params['NVE']["design angle"], theta0 = self.params['NVE']["incidence angle"],
                    _x = self.params["NVE"]["xc"],
                    _y = self.params["NVE"]["yc"],
                    length = self.params['NVE']["length"],
                    roll = self.params['NVE']["roll"],
                    yaw = self.params['NVE']["yaw"],
                    _refl = self.params['NVE']["reflectivity"],
                    _ext_in = self.params['NVE']["_ext_in"], _ext_out = self.params['NVE']["_ext_out"]) 
            
            self.NVE.name = self.params["NVE"]["name"]
            
            self.NVE_error = MirPl(np.loadtxt(self.params['NVE_error']['mirror profile'].replace("../../","")),
                            _dim = self.params['NVE_error']['orientation'],
                            _ang = self.params['NVE_error']['incidence angle'], ### + self.params['NVE']['incidence angle'], 
                            _refl = self.params['NVE_error']['transmission'],
                            _x = self.params['NVE_error']['xc'], _y = self.params['NVE_error']['yc']) 
            
            self.NVE_error.name = self.params['NVE_error']['name']
            
            self.NHE_error = MirPl(np.loadtxt(self.params['NHE_error']['mirror profile'].replace("../../","")),
                _dim = self.params['NHE_error']['orientation'],
                _ang = self.params['NHE_error']['incidence angle'], ###+self.params['NHE']['incidence angle'], 
                _refl = self.params['NHE_error']['transmission'],
                _x = self.params['NHE_error']['xc'], _y = self.params['NHE_error']['yc']) 

            
                  
            self.NHE_error.name = self.params['NHE_error']['name']
            
            self.params["d3"]["distance"] = 656.424
            self.params["d4"]["distance"] = 1.20
            self.params["d5"]["distance"] = 1.00
            self.params["df"]["distance"] = 2.2
            
            self.d3 =  Drift(self.params["d3"]['distance'])
            self.d3.name = self.params["d3"]['name']
            
            self.d4 =  Drift(self.params["d4"]['distance'])
            self.d4.name = self.params["d4"]['name']
            
            self.d5 =  Drift(self.params["d5"]['distance'])
            self.d5.name = self.params["d5"]['name']
            
            self.df =  Drift(self.params["df"]['distance'])
            self.df.name =  'df'
            
    def build_beamline(self, focus = "nano", screens = "false"):
        """
        Construct the beamline object
        
        :param focus: what beamline configuration (micron or nano)
        """
        
        self.bl = Beamline()
        
        if focus == 'direct':
            
            self.bl.append(self.d1, propagation_parameters(1,1,1,1, mode = "quadratic"))
    
            self.bl.append(self.HOM1, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d2, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.HOM2,  propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d3, propagation_parameters(2,1,2,1, mode = 'fraunhofer'))
            
            self.bl.append(self.NKB_PSlit, propagation_parameters(1/10, 1, 1/10,  1, mode = 'fresnel'))
            self.bl.append(self.d4, propagation_parameters(1, 1, 1, 1, mode = 'fraunhofer'))
            self.bl.append(self.d5, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.df, propagation_parameters(1,1,1,1, mode = 'quadratic'))

            
       
        if focus == "micron":
            
            
            self.bl.append(self.d1, propagation_parameters(1,1,1,1, mode = "fraunhofer"))
    
            self.bl.append(self.HOM1, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d2, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.HOM2,  propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d3, propagation_parameters(1,1,1,1, mode = 'fraunhofer'))
            
            self.bl.append(self.d4, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.MHP, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d5, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.MHE, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.MHE_error, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d6, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.MVE, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.MVE_error, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d8, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
        
            self.bl.append(self.MVP, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.df, propagation_parameters(5,1,5,1, mode = 'converge'))
       
        elif focus == "nano":
            
            self.bl.append(self.d1, propagation_parameters(*self.params['d1']['pp']))
    
            self.bl.append(self.HOM1, propagation_parameters(*self.params['HOM1']['pp']))
            self.bl.append(self.d2, propagation_parameters(*self.params['d2']['pp']))
            self.bl.append(self.HOM2,  propagation_parameters(*self.params['HOM2']['pp']))
            self.bl.append(self.d3, propagation_parameters(*self.params['d3']['pp']))
            
            self.bl.append(self.NKB_PSlit, propagation_parameters(*self.params['NKB_PSlit']['pp']))
            self.bl.append(self.d4, propagation_parameters(*self.params['d4']['pp']))
            self.bl.append(self.NHE_error, propagation_parameters(*self.params['NHE_error']['pp']))
            self.bl.append(self.NHE, propagation_parameters(*self.params['NHE']['pp']))
            
            
            self.bl.append(self.d5, propagation_parameters(*self.params['d5']['pp']))
            self.bl.append(self.NVE_error, propagation_parameters(*self.params['NVE_error']['pp']))
            self.bl.append(self.NVE, propagation_parameters(*self.params['NVE']['pp']))
            
            #self.bl.append(self.df, propagation_parameters(1,1,1,1, mode = 'converge'))
            #self.bl.append(self.df, propagation_parameters(15,1,15,1, mode = 'converge')) ### perfect mirrors
            
        self.bl.params = self.params
        
    def get_beamline(self):
        return self.bl
    
    
    def crop_beamline(self, element1 = None, element2 = None):
        """
        crop the beamline to some position as defined by the name of the optical element
        
        :param bl: beamline object
        :param position: position to be cropped to
        """
        
        if element1 is not None:
            names = [el.name for el in self.bl.propagation_options[0]['optical_elements']]
            idx1 = names.index(element1)
        else: 
            idx1 = 0
            
        if element2 is not None:
            names = [el.name for el in self.bl.propagation_options[0]['optical_elements']]
            idx2 = names.index(element2)
            
             
        elif element2 is not None:
            self.bl.propagation_options[0]['optical_elements'] = self.bl.propagation_options[0]['optical_elements'][idx1:idx2+1]
            self.bl.propagation_options[0]['propagation_parameters'] = self.bl.propagation_options[0]['propagation_parameters'][idx1:idx2+1]
        else:
            self.bl.propagation_options[0]['optical_elements'] = self.bl.propagation_options[0]['optical_elements'][:idx1+1]
            self.bl.propagation_options[0]['propagation_parameters'] = self.bl.propagation_options[0]['propagation_parameters'][:idx1+1]

    def add_screen(self,position, distance, screenName = None):
        """
        add a screening plane at drift beyond some element posiFalsetion
        
        :param position: last optical element before screen
        :param distance: position from last optical element to screen
        :param screenName: name of the screen element (ie., MKB-scr etc) [str]
        """
        self.crop_beamline(element1 = position)
        
        drift2screen = Drift(distance)
        if screenName is not None:
            drift2screen.name = "screen"
        else:
            drift2screen.name = screenName
        self.bl.append(Drift(distance), propagation_parameters(1, 1, 1, 1, m = 'quadratic'))
    
    def mirror_profiles(self, surface, aperture, overwrite = True):
        """
        toggle for mirror surfaces
        """ 
        self.load_mirror_profiles(aperture = aperture, surface = surface)
            
        
    def get_mirror_profile(self, mirror_name, sdir = None):
        
        
        surface = np.loadtxt(self.params[mirror_name]['mirror profile'])
                
        x = surface[1:, 0]*1e3
        y = surface[0, 1:]*1e3
        
        surface = surface[1:,1:]*1e9
        
        return surface, x, y
    
    def list_elements(self):
        return [el.name for el in self.bl.propagation_options[0]['optical_elements']]
        
    def get_index(self, element_name):
        """
        get the index of an element in a beamline by name
        """
        ### get index
        names = self.list_elements()
    
        if self.VERBOSE:
            print("List of Elements: {}".format(names))
        try:
            index = names.index(element_name)
        except(ValueError):
            print("Beamline does not contain optical element: {}".format(element_name))

        return index
        
    def remove_element(self, element_name):
        """ 
        remove an element from the beamline by name
        """
        
        index = self.get_index(element_name)

        del self.bl.propagation_options[0]['optical_elements'][index]
        del self.bl.propagation_options[0]['propagation_parameters'][index]

    def edit_propagation_parameters(self, element_name, new_parameters):
        """
        edit the propagation parameters of an element by name
        """
        self.bl.propagation_options[0]['propagation_parameters'][self.get_index(element_name)] = new_parameters

        
    def rebuild(self, focus = 'nano'):
        self.build_elements(focus)
        self.build_beamline(focus)