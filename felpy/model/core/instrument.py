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


import os
from os.path import exists
import json
import numpy as np
import sys
from matplotlib import pyplot as plt

from felpy.model.core.beamline import Beamline
from wpg import srwlib
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture
from wpg.srwlib import SRWLOptT
from wpg.optical_elements import Mirror_elliptical as MirEl

from felpy.utils.os_utils import add_path, felpy_path
from felpy.model.materials.mirror_surface import genMirrorSurface
from felpy.model.materials.load_refl import get_refl, load_refl
from felpy.model.src.coherent import construct_SA1_wavefront
from wpg.optical_elements import calculateOPD
from felpy.model.beamlines.exfel_spb.params import get_params
from felpy.model.tools import propagation_parameters
from wpg.srwlib import srwl_opt_setup_surf_height_2d as MirPl

from felpy.utils.vis_utils import colorbar_plot
from felpy.utils.np_utils import get_mesh

import seaborn as sns

class Instrument:
    
    """
    A container for loading and modifying beamline data
    """
    
    
    def __init__(self, VERBOSE = True):
        
        self.VERBOSE = VERBOSE
        
        if VERBOSE:
            print("Initialising Single Particle Beamline")
        
        self.load_params()
        self.fpath = felpy_path() ### felpy path (for dev. purposes)
        
        add_path()
        
    def load_params(self, fromFile = False):
        """
        load beamline parameters from /data/spb/
        """
        
        if fromFile:
            with open("../../data/params/exfel_spb.json", "r") as read_file:
                self.params = json.load(read_file)
        else:
            self.params = get_params()
            
    def export_params(self, sdir = None):
        """
        save the current set of beamline parameters in json format
        
        :param sdir: save directory
        """
        if sdir is None:
            with open('../../data/spb/parameters.json', 'w') as f:
                json.dump(self.params, f)
        else:
            with open(sdir + 'parameters.json', 'w') as f:
                json.dump(self.params, f)

        
        
    def adjust_mirror(self, mirror_name, ekev, new_ang, mirror_refl = None):
     
        if mirror_refl == None: 
            if ekev >= 7.5:
                material = "B4C"
            else: 
                material = "Ru"    

            refl = get_refl(load_refl(material), ekev, new_ang)
        
        else:
            refl = mirror_refl 
        
        new_ang = new_ang + np.tan(self.params[mirror_name]["xc"]/self.params[self.params[mirror_name]['next_drift']]['distance'])
        self.params[mirror_name]["incidence angle"] = new_ang
        self.params[mirror_name]['reflectivity'] = refl


    def define_mirror_profiles(self, overwrite = False, surface = 'flat', plot = False, aperture = True):
        """
        Define the plane mirror profiles by loading from /data/. 
        If mirror profiles not defined (unlikely), generate profiles via genMirrorSurface
        
        :param overwrite: bool to overwrite current file.
        :param surface: surface type (flat, random or real)
        """
    
        
        if surface == 'real':
            
            mm = 'random' ### fix for no 'real' surface MHE, MVE, MVP surfaces etc.
            
            if overwrite == True:
                
                if aperture == True:
                    genMirrorSurface(500, 500, [self.params["MHP"]['dx'],self.params["MHP"]['dy']], "../../data/spb/mirror_surface/mhp_", mode = mm, plot = plot, mirrorName = "MHP") 
                    genMirrorSurface(500, 500, [self.params["MVP"]['dx'],self.params["MVP"]['dy']], "../../data/spb/mirror_surface/mvp_", mode = mm, plot = plot, mirrorName = "MVP")  
                    genMirrorSurface(500, 500, [self.params["MHE"]['dx'],self.params["MHE"]['dy']], "../../data/spb/mirror_surface/mhe_", mode = mm, plot = plot, mirrorName = "MHE")
                    genMirrorSurface(500, 500, [self.params["MVE"]['dx'],self.params["MVE"]['dy']], "../../data/spb/mirror_surface/mve_", mode = mm, plot = plot, mirrorName = "MVE")  
                    genMirrorSurface(500, 500, [self.params["NHE"]['dx'],self.params["NHE"]['dy']], "../../data/spb/mirror_surface/nhe_", mode = surface, plot = plot, mirrorName = "NHE")
                    genMirrorSurface(500, 500, [self.params["NVE"]['dx'],self.params["NVE"]['dy']], "../../data/spb/mirror_surface/nve_", mode = surface, plot = plot, mirrorName = "NVE")  
                    
                elif aperture == False:
                    genMirrorSurface(500, 500, [100,100], "data/spb/mirror_surface/hom1_", mode = mm, plot = plot, mirrorName = "HOM1") 
                    genMirrorSurface(500, 500, [100,100], "data/spb/mirror_surface/hom2_", mode = mm, plot = plot, mirrorName = "HOM2")  
                    genMirrorSurface(500, 500, [100,100], "data/spb/mirror_surface/mhp_", mode = mm, plot = plot, mirrorName = "MHP") 
                    genMirrorSurface(500, 500, [100,100], "data/spb/mirror_surface/mvp_", mode = mm, plot = plot, mirrorName = "MVP")  
                    genMirrorSurface(500, 500, [100,100], "data/spb/mirror_surface/mhe_", mode = mm, plot = plot, mirrorName = "MHE")
                    genMirrorSurface(500, 500, [100,100], "../../data/spb/mirror_surface/mve_", mode = mm, plot = plot, mirrorName = "MVE")  
                    genMirrorSurface(500, 500, [100,100], "../../data/spb/mirror_surface/nhe_", mode = mm, plot = plot, mirrorName = "NHE")
                    genMirrorSurface(500, 500, [100,100], "../../data/spb/mirror_surface/nve_", mode = mm, plot = plot, mirrorName = "NVE")  
                    
        
        elif surface == 'flat':
            
            mm = 'flat'
        
            if overwrite == True:
                
                if aperture == True:
                    genMirrorSurface(500, 500, [self.params["MHP"]['dx'],self.params["MHP"]['dy']], "../../data/spb/mirror_surface/mhp_", mode = surface, plot = plot, mirrorName = "MHP") 
                    genMirrorSurface(500, 500, [self.params["MVP"]['dx'],self.params["MVP"]['dy']], "../../data/spb/mirror_surface/mvp_", mode = surface, plot = plot, mirrorName = "MVP")  
                    genMirrorSurface(500, 500, [self.params["MHE"]['dx'],self.params["MHE"]['dy']], "../../data/spb/mirror_surface/mhe_", mode = surface, plot = plot, mirrorName = "MHE")
                    genMirrorSurface(500, 500, [self.params["MVE"]['dx'],self.params["MVE"]['dy']], "../../data/spb/mirror_surface/mve_", mode = surface, plot = plot, mirrorName = "MVE")  
                    genMirrorSurface(500, 500, [self.params["NHE"]['dx'],self.params["NHE"]['dy']], "../../data/spb/mirror_surface/nhe_", mode = surface, plot = plot, mirrorName = "NHE")
                    genMirrorSurface(500, 500, [self.params["NVE"]['dx'],self.params["NVE"]['dy']], "../../data/spb/mirror_surface/nve_", mode = surface, plot = plot, mirrorName = "NVE")  
                
                if aperture == False:
                    genMirrorSurface(500, 500, [100,100], "../../data/spb/mirror_surface/hom1_", mode = surface, plot = plot, mirrorName = "HOM1") 
                    genMirrorSurface(500, 500, [100,100], "../../data/spb/mirror_surface/hom2_", mode = surface, plot = plot, mirrorName = "HOM2")   
                    genMirrorSurface(500, 500, [100,100], "../../data/spb/mirror_surface/mhp_", mode = surface, plot = plot, mirrorName = "MHP") 
                    genMirrorSurface(500, 500, [100,100], "../../data/spb/mirror_surface/mvp_", mode = surface, plot = plot, mirrorName = "MVP")  
                    genMirrorSurface(500, 500, [1000,1000], "../../data/spb/mirror_surface/mhe_", mode = surface, plot = plot, mirrorName = "MHE")
                    genMirrorSurface(500, 500, [1000,1000], "../../data/spb/mirror_surface/mve_", mode = surface, plot = plot, mirrorName = "MVE")  
                    genMirrorSurface(500, 500, [100,100], "../../data/spb/mirror_surface/nhe_", mode = surface, plot = plot, mirrorName = "NHE")
                    genMirrorSurface(500, 500, [100,100], "../../data/spb/mirror_surface/nve_", mode = surface, plot = plot, mirrorName = "NVE")  
                
        
        if aperture == True:
            self.params['HOM1']['mirror profile'] = "../../data/spb/mirror_surface/hom1_mir_{}.dat".format(surface)
            self.params['HOM2']['mirror profile'] = "../../data/spb/mirror_surface/hom2_mir_{}.dat".format(surface)
        else:
            self.params['HOM1']['mirror profile'] = "../../data/spb/mirror_surface/hom1_mir_{}.dat".format(mm)
            self.params['HOM2']['mirror profile'] = "../../data/spb/mirror_surface/hom2_mir_{}.dat".format(mm)
            self.params['MHE']["length"] = 10
            self.params['MVE']["length"] = 10
            
        self.params['MHP']['mirror profile'] = "../../data/spb/mirror_surface/mhp_mir_{}.dat".format(mm)
        self.params['MVP']['mirror profile'] = "../../data/spb/mirror_surface/mvp_mir_{}.dat".format(mm)
        self.params['MHE_error']['mirror profile'] = "../../data/spb/mirror_surface/mhe_mir_{}.dat".format(mm)
        self.params['MVE_error']['mirror profile'] = "../../data/spb/mirror_surface/mve_mir_{}.dat".format(mm)
        self.params['NHE_error']['mirror profile'] = "../../data/spb/mirror_surface/nhe_mir_{}.dat".format(mm)
        self.params['NVE_error']['mirror profile'] = "../../data/spb/mirror_surface/nve_mir_{}.dat".format(mm)
        
        
 
    def build_elements(self, focus = "nano"):
        

        self.d1 =  Drift(self.params["HOM1"]['distance from source'])
        self.d1.name = self.params["d1"]['name']
        
        self.HOM1 = MirPl(np.loadtxt(self.fpath + self.params['HOM1']['mirror profile'].replace("../../","")),
                     _dim = self.params['HOM1']['orientation'],
                     _ang = self.params['HOM1']['incidence angle'], 
                     _amp_coef = 1,
                     _x = self.params['HOM1']['xc'], _y = self.params['HOM1']['yc'],
                     _refl = self.params['HOM1']['transmission']) 
        
        self.HOM1.name = self.params['HOM1']['name']

        self.d2 =  Drift(self.params["HOM2"]['distance from source']-self.params["HOM1"]['distance from source'])
        self.d2.name = self.params["d2"]['name']
        
        self.HOM2 = MirPl(np.loadtxt(self.fpath + self.params['HOM2']['mirror profile'].replace("../../","")),
                     _dim = self.params['HOM2']['orientation'],
                     _ang = self.params['HOM2']['incidence angle'], 
                     _amp_coef = 1,
                     _x = self.params['HOM2']['xc'], _y = self.params['HOM2']['yc'],
                     _refl = self.params['HOM2']['transmission']) 
        
        self.HOM2.name = self.params['HOM2']['name']
        

        
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
            

            self.MHP = MirPl(np.loadtxt(self.fpath + self.params['MHP']['mirror profile']),
                         _dim = self.params['MHP']['orientation'],
                         _ang = self.params['MHP']['incidence angle'], 
                         _refl = self.params['MHP']['transmission'],
                         _x = self.params['MHP']['xc'], _y = self.params['MHP']['yc']) 
            self.MHP.name = self.params['MHP']['name']
            
            self.d5 =  Drift(self.params["d5"]['distance'])
            self.d5.name = self.params["d5"]['name']
            
            self.MHE_error = MirPl(np.loadtxt(self.fpath + self.params['MHE_error']['mirror profile']),
             _dim = self.params['MHE_error']['orientation'],
             _ang = self.params['MHE_error']['incidence angle'], 
             _refl = self.params['MHE_error']['transmission'],
             _x = self.params['MHE_error']['xc'], _y = self.params['MHE_error']['yc']) 
            
            self.MHE_error.name = self.params['MHE_error']['name']
            

            self.MVE_error = MirPl(np.loadtxt(self.fpath + self.params['MVE_error']['mirror profile']),
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
            
            self.MVP = MirPl(np.loadtxt(self.fpath + self.params['MVP']['mirror profile']),
                     _dim = self.params['MVP']['orientation'],
                     _ang = self.params['MVP']['incidence angle'], 
                     _refl = self.params['MVP']['transmission'],
                     _x = self.params['MVP']['xc'], _y = self.params['MVP']['yc']) 
            self.MVP.name = self.params['MVP']['name']
            
            self.df =  Drift(self.params["df"]['distance'])
            self.df.name = self.params["df"]['name']
        
        elif focus == "nano":
            
            self.NKB_pslit = Aperture(_shape=self.params["NKB_pslit"]['shape'],
                                      _ap_or_ob=self.params["NKB_pslit"]['type'],
                                      _Dx= self.params["NKB_pslit"]['dx'],
                                      _Dy= self.params["NKB_pslit"]['dy'],
                                      _x=self.params["NKB_pslit"]['xc'],
                                      _y=self.params["NKB_pslit"]['yc'])
            self.NKB_pslit.name = self.params["NKB_pslit"]['name']
            
            
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
            
            self.NVE_error = MirPl(np.loadtxt(self.fpath + self.params['NVE_error']['mirror profile'].replace("../../","")),
                            _dim = self.params['NVE_error']['orientation'],
                            _ang = self.params['NVE_error']['incidence angle']+self.params['NVE']['incidence angle'], 
                            _refl = self.params['NVE_error']['transmission'],
                            _x = self.params['NVE_error']['xc'], _y = self.params['NVE_error']['yc']) 
            
            self.NVE_error.name = self.params['NVE_error']['name']
            
            self.NHE_error = MirPl(np.loadtxt(self.fpath + self.params['NHE_error']['mirror profile'].replace("../../","")),
                _dim = self.params['NHE_error']['orientation'],
                _ang = self.params['NHE_error']['incidence angle']+self.params['NHE']['incidence angle'], 
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
            self.df.name = self.params["df"]['name']
            
    def build_beamline(self, focus = "nano", screens = "false"):
        """
        Construct the beamline object
        
        :param focus: what beamline configuration (micron or nano)
        """
        
        self.bl = Beamline()
        
        if focus == "micron":
            
            
            self.bl.append(self.d1, propagation_parameters(1,1,1,1, mode = "fraunhofer"))
    
            self.bl.append(self.HOM1, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d2, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.HOM2,  propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d3, propagation_parameters(1,1,1,1, mode = 'fraunhofer'))
            
            self.bl.append(self.MKB_pslit, propagation_parameters(1/5, 1, 1/5, 1, mode = 'fresnel'))
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
            
            self.bl.append(self.d1, propagation_parameters(1,1,1,1, mode = "fraunhofer"))
    
            self.bl.append(self.HOM1, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d2, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.HOM2,  propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.d3, propagation_parameters(1,1,1,1, mode = 'fraunhofer'))
            
            self.bl.append(self.NKB_pslit, propagation_parameters(1/10, 1, 1/10,  1, mode = 'fresnel'))
            self.bl.append(self.d4, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.NHE_error, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.NHE, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            
            
            self.bl.append(self.d5, propagation_parameters(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.NVE_error, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            self.bl.append(self.NVE, propagation_parameters(1, 1, 1, 1, mode = 'fresnel'))
            
            
            self.bl.append(self.df, propagation_parameters(1/20,1,1/20,1, mode = 'converge'))

            
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
    
        if element2 is not None:
            names = [el.name for el in self.bl.propagation_options[0]['optical_elements']]
            idx2 = names.index(element2)
            
            
        if element1 not in names:
            pass 
        elif element2 is not None:
            self.bl.propagation_options[0]['optical_elements'] = self.bl.propagation_options[0]['optical_elements'][idx1:idx2+1]
            self.bl.propagation_options[0]['propagation_parameters'] = self.bl.propagation_options[0]['propagation_parameters'][idx1:idx2+1]
        else:
            self.bl.propagation_options[0]['optical_elements'] = self.bl.propagation_options[0]['optical_elements'][:idx1+1]
            self.bl.propagation_options[0]['propagation_parameters'] = self.bl.propagation_options[0]['propagation_parameters'][:idx1+1]

    def add_screen(self,position, distance, screenName = None):
        """
        add a screening plane at drift beyond some element position
        
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
    
    def mirror_profiles(self, toggle = "on", aperture = True, overwrite = False):
        """
        toggle for mirror surfaces
        """
        if toggle == "on":
            self.define_mirror_profiles(overwrite = overwrite, aperture = aperture, surface = 'real')
        if toggle == "off":
            self.define_mirror_profiles(overwrite = overwrite, aperture = aperture, surface = 'flat')
            
    def get_mirror_profiles(self, mirror_name, context = 'talk', sdir = None):
        
        sns.set()
        
        surface = np.loadtxt(self.params[mirror_name]['mirror profile'])
        
        x = surface[1:, 0]
        y = surface[0, 1:]
        
        surface = surface[1:,1:]
        
        mesh = get_mesh(surface, abs(x[1]-x[0]), abs(y[1]-y[0]))

        return surface, mesh
 
 
 