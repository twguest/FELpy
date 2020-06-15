#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 14:43:15 2020

@author: twguest
"""

###############################################################################
import sys
sys.path.append("/opt/WPG/") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/WPG") # DESY MAXWELL PATH

sys.path.append("/opt/spb_model") # LOCAL PATH
sys.path.append("/gpfs/exfel/data/user/guestt/spb_model") # DESY MAXWELL PATH
###############################################################################
###############################################################################
import os
from os.path import exists
import json
import numpy as np

from matplotlib import pyplot as plt

from wpg.beamline import Beamline
from wpg import srwlib
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture
from wpg.srwlib import SRWLOptT
from wpg.optical_elements import Mirror_elliptical as MirEl

from model.materials.mirrorSurface import genMirrorSurface
from model.materials.load_refl import get_refl, load_refl
from model.src.coherent import coherentSource
from wpg.optical_elements import calculateOPD

from wpg.srwlib import srwl_opt_setup_surf_height_2d as MirPl

def propParams(sx, zx, sy, zy, mode = "normal"):
    """
    wrapper for propagation parameters
    
    :param zx: horizontal scaling factor 
    :param rx: horizontal zoom factor
    :param zy: vertical scaling factor
    :param ry: vertical zoom factor
    :param mode: normal, semi-analytical, converge or diverge
    
    :return propagation parameters:
    """
    
    if mode == "normal":
        m = 0
    elif mode == "quadratic":
        m = 1
    elif mode == "farfield":
        m = 2
    elif mode == "diverge":
        m = 3
    elif mode == "converge":
        m = 4
    
    return [0,0,1,m,0,sx,zx/sx,sy,zy/sy,0,0,0]



class BeamlineModel:
    
    """
    A container for loading and modifying beamline data
    """
    
    
    def __init__(self):
        print("Initialising Single Particle Beamline")
        self.load_params()

    def load_params(self):
        """
        load beamline parameters from /data/input/
        """
        with open("../../data/input/parameters.json", "r") as read_file:
            self.params = json.load(read_file)
    
    def export_params(self, outdir = None):
        """
        save the current set of beamline parameters in json format
        
        :param outdir: save directory
        """
        if outdir is None:
            with open('../../data/input/parameters.json', 'w') as f:
                json.dump(self.params, f)
        else:
            with open(outdir + 'parameters.json', 'w') as f:
                json.dump(self.params, f)
    
    def setupHOMs(self, ekev, ang = 2.2e-03):
    
        refl_data = load_refl()
        refl, ang = get_refl(refl_data, ekev, ang, limits = [1.1e-03, 3.6e-03])
        
        self.adjustHOMs(refl, ang)           
        
    def setupKBs(self, ekev, ang = 3.5e-03, misorientation = 0, focus = "nano"):
        
        if ekev >= 7.5:
            material = "B4C"
        else: 
            material = "Ru"
            
        refl_data = load_refl(material, indir = "../../data/kb_refl")
        refl, ang = get_refl(refl_data, ekev, ang, limits = [0, 5.5e-03])
        
        if focus == 'nano':
            self.adjustNHE(refl, ang)
            
            ang_e = ang + np.tan(self.params["NVE"]["xc"]/self.params["df"]["distance"])
            refl_e, ang_e = get_refl(refl_data, ekev, ang_e, limits = [0, 5.5e-03])
            
            self.adjustNVE(refl_e, ang_e)

    def adjustNHE(self, refl = None, ang = None, misorientation = 0):
        
        if refl is not None:
            self.params["NHE"]['reflectivity'] = refl
        if ang is not None:
            self.params["NHE"]["design angle"] = ang
        if ang or misorientation is not None:
            self.params["NHE"]["incidence angle"] = self.params["NHE"]["design angle"] + misorientation
            
    def adjustNVE(self, refl, ang, misorientation = 0):
        
        if refl is not None:
            self.params["NVE"]['reflectivity'] = refl
        if ang is not None:
            self.params["NVE"]["design angle"] = ang
        if ang or misorientation is not None:
            self.params["NVE"]["incidence angle"] = self.params["NVE"]["design angle"] + misorientation
            
    def adjustKBs(self, refl = None, ang = None, misalignment = 0):
        """
        Wrapper to adjust the components of the kbs via editing the current 
        parameters file
        """
        
        if refl is not None:
            self.params["MHE"]['reflectivity'] = refl
            self.params["MVE"]['reflectivity'] = refl
            self.params["NHE"]['reflectivity'] = refl
            
        if ang is not None:
            self.params["MHE"]["design angle"] = ang
            self.params["MVE"]["design angle"] = ang
            self.params["NHE"]["design angle"] = ang
            
            self.params["MHE"]["incidence angle"] = ang + misalignment
            self.params["MVE"]["incidence angle"] = ang + misalignment
            self.params["NHE"]["incidence angle"] = ang + misalignment
            
    
    def adjustHOMs(self, refl = None, ang = None):
        """
        Wrapper to adjust the components of the HOMs via editing the current 
        parameters file
        """
        
        if refl is not None:
            self.params["HOM1"]['transmission'] = refl
            self.params["HOM2"]['transmission'] = refl
        if ang is not None:
            self.params["HOM1"]["incidence angle"] = ang
            self.params["HOM2"]["incidence angle"] = ang
    
    def adjustMHE(self, refl = None, ang = None):
        """
        Wrapper to adjust the components of the MHE via editing the current 
        parameters file
        """
        
        if refl is not None:
            self.params["MHE"]['transmission'] = refl
            self.params["MHP"]['transmission'] = refl
        if ang is not None:
            self.params["MHE"]["incidence angle"] = ang
            self.params["MHP"]["incidence angle"] = ang

    def adjustMVE(self, refl = None, ang = None):
        """
        Wrapper to adjust the components of the MVE via editing the current 
        parameters file
        """
        if refl is not None:
            self.params["MVE"]['transmission'] = refl
            self.params["MVP"]['transmission'] = refl
        if ang is not None:
            self.params["MVE"]["incidence angle"] = ang
            self.params["MVP"]["incidence angle"] = ang
    

    def defineMirrorProfiles(self, overwrite = False, surface = 'flat', plot = False, aperture = True):
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
                    genMirrorSurface(500, 500, [self.params["MHP"]['dx'],self.params["MHP"]['dy']], "../../data/input/mhp_", mode = mm, plot = plot, mirrorName = "MHP") 
                    genMirrorSurface(500, 500, [self.params["MVP"]['dx'],self.params["MVP"]['dy']], "../../data/input/mvp_", mode = mm, plot = plot, mirrorName = "MVP")  
                    genMirrorSurface(500, 500, [self.params["MHE"]['dx'],self.params["MHE"]['dy']], "../../data/input/mhe_", mode = mm, plot = plot, mirrorName = "MHE")
                    genMirrorSurface(500, 500, [self.params["MVE"]['dx'],self.params["MVE"]['dy']], "../../data/input/mve_", mode = mm, plot = plot, mirrorName = "MVE")  
                    genMirrorSurface(500, 500, [self.params["NHE"]['dx'],self.params["NHE"]['dy']], "../../data/input/nhe_", mode = mm, plot = plot, mirrorName = "NHE")
                    genMirrorSurface(500, 500, [self.params["NVE"]['dx'],self.params["NVE"]['dy']], "../../data/input/nve_", mode = mm, plot = plot, mirrorName = "NVE")  
                    
                elif aperture == False:
                    genMirrorSurface(500, 500, [100,100], "../../data/input/hom1_", mode = mm, plot = plot, mirrorName = "HOM1") 
                    genMirrorSurface(500, 500, [100,100], "../../data/input/hom2_", mode = mm, plot = plot, mirrorName = "HOM2")  
                    genMirrorSurface(500, 500, [100,100], "../../data/input/mhp_", mode = mm, plot = plot, mirrorName = "MHP") 
                    genMirrorSurface(500, 500, [100,100], "../../data/input/mvp_", mode = mm, plot = plot, mirrorName = "MVP")  
                    genMirrorSurface(500, 500, [100,100], "../../data/input/mhe_", mode = mm, plot = plot, mirrorName = "MHE")
                    genMirrorSurface(500, 500, [100,100], "../../data/input/mve_", mode = mm, plot = plot, mirrorName = "MVE")  
                    genMirrorSurface(500, 500, [100,100], "../../data/input/nhe_", mode = mm, plot = plot, mirrorName = "NHE")
                    genMirrorSurface(500, 500, [100,100], "../../data/input/nve_", mode = mm, plot = plot, mirrorName = "NVE")  
                    
        
        elif surface == 'flat':
            
            mm = 'flat'
        
            if overwrite == True:
                
                if aperture == True:
                    genMirrorSurface(500, 500, [self.params["MHP"]['dx'],self.params["MHP"]['dy']], "../../data/input/mhp_", mode = surface, plot = plot, mirrorName = "MHP") 
                    genMirrorSurface(500, 500, [self.params["MVP"]['dx'],self.params["MVP"]['dy']], "../../data/input/mvp_", mode = surface, plot = plot, mirrorName = "MVP")  
                    genMirrorSurface(500, 500, [self.params["MHE"]['dx'],self.params["MHE"]['dy']], "../../data/input/mhe_", mode = surface, plot = plot, mirrorName = "MHE")
                    genMirrorSurface(500, 500, [self.params["MVE"]['dx'],self.params["MVE"]['dy']], "../../data/input/mve_", mode = surface, plot = plot, mirrorName = "MVE")  
                    genMirrorSurface(500, 500, [self.params["NHE"]['dx'],self.params["NHE"]['dy']], "../../data/input/nhe_", mode = surface, plot = plot, mirrorName = "NHE")
                    genMirrorSurface(500, 500, [self.params["NVE"]['dx'],self.params["NVE"]['dy']], "../../data/input/nve_", mode = surface, plot = plot, mirrorName = "NVE")  
                
                if aperture == False:
                    genMirrorSurface(500, 500, [100,100], "../../data/input/hom1_", mode = surface, plot = plot, mirrorName = "HOM1") 
                    genMirrorSurface(500, 500, [100,100], "../../data/input/hom2_", mode = surface, plot = plot, mirrorName = "HOM2")   
                    genMirrorSurface(500, 500, [100,100], "../../data/input/mhp_", mode = surface, plot = plot, mirrorName = "MHP") 
                    genMirrorSurface(500, 500, [100,100], "../../data/input/mvp_", mode = surface, plot = plot, mirrorName = "MVP")  
                    genMirrorSurface(500, 500, [1000,1000], "../../data/input/mhe_", mode = surface, plot = plot, mirrorName = "MHE")
                    genMirrorSurface(500, 500, [1000,1000], "../../data/input/mve_", mode = surface, plot = plot, mirrorName = "MVE")  
                    genMirrorSurface(500, 500, [100,100], "../../data/input/nhe_", mode = surface, plot = plot, mirrorName = "NHE")
                    genMirrorSurface(500, 500, [100,100], "../../data/input/nve_", mode = surface, plot = plot, mirrorName = "NVE")  
                
        
        if aperture == True:
            self.params['HOM1']['mirror profile'] = "../../data/input/hom1_mir_{}.dat".format(surface)
            self.params['HOM2']['mirror profile'] = "../../data/input/hom2_mir_{}.dat".format(surface)
        else:
            self.params['HOM1']['mirror profile'] = "../../data/input/hom1_mir_{}.dat".format(mm)
            self.params['HOM2']['mirror profile'] = "../../data/input/hom2_mir_{}.dat".format(mm)
            self.params['MHE']["length"] = 10
            self.params['MVE']["length"] = 10
            
        self.params['MHP']['mirror profile'] = "../../data/input/mhp_mir_{}.dat".format(mm)
        self.params['MVP']['mirror profile'] = "../../data/input/mvp_mir_{}.dat".format(mm)
        self.params['MHE_error']['mirror profile'] = "../../data/input/mhe_mir_{}.dat".format(mm)
        self.params['MVE_error']['mirror profile'] = "../../data/input/mve_mir_{}.dat".format(mm)
        self.params['NHE_error']['mirror profile'] = "../../data/input/nhe_mir_{}.dat".format(mm)
        self.params['NVE_error']['mirror profile'] = "../../data/input/nve_mir_{}.dat".format(mm)
        
        
 
    def buildElements(self, focus = "micron"):
       
        self.d1 =  Drift(self.params["HOM1"]['distance from source'])
        self.d1.name = self.params["d1"]['name']
        
        self.HOM1 = MirPl(np.loadtxt(self.params['HOM1']['mirror profile']),
                     _dim = self.params['HOM1']['orientation'],
                     _ang = self.params['HOM1']['incidence angle'], 
                     _amp_coef = 1,
                     _x = self.params['HOM1']['xc'], _y = self.params['HOM1']['yc'],
                     _refl = self.params['HOM1']['transmission']) 
        
        self.HOM1.name = self.params['HOM1']['name']

        self.d2 =  Drift(self.params["HOM2"]['distance from source']-self.params["HOM1"]['distance from source'])
        self.d2.name = self.params["d2"]['name']
        
        self.HOM2 = MirPl(np.loadtxt(self.params['HOM2']['mirror profile']),
                     _dim = self.params['HOM2']['orientation'],
                     _ang = self.params['HOM2']['incidence angle'], 
                     _amp_coef = 1,
                     _x = self.params['HOM2']['xc'], _y = self.params['HOM2']['yc'],
                     _refl = self.params['HOM2']['transmission']) 
        
        self.HOM2.name = self.params['HOM2']['name']
        

        
        if focus == "micron":
            self.d3 =  Drift(self.params["MKB_pslit"]['distance from source']-self.params["HOM2"]['distance from source'])
            self.d3.name = self.params["d3"]['name']
                
            self.MKB_pslit = Aperture(_shape=self.params["MKB_pslit"]['shape'],
                                 _ap_or_ob=self.params["MKB_pslit"]['type'],
                                 _Dx= self.params["MKB_pslit"]['dx'],
                                 _Dy= self.params["MKB_pslit"]['dy'],
                                 _x=self.params["MKB_pslit"]['xc'],
                                 _y=self.params["MKB_pslit"]['yc'])
            self.MKB_pslit.name = self.params["MKB_pslit"]['name']
            
            self.d4 =  Drift(self.params["MHP"]['distance from source']-self.params["MKB_pslit"]['distance from source'])
            self.d4.name = self.params["d4"]['name']
            

            self.MHP = MirPl(np.loadtxt(self.params['MHP']['mirror profile']),
                         _dim = self.params['MHP']['orientation'],
                         _ang = self.params['MHP']['incidence angle'], 
                         _refl = self.params['MHP']['transmission'],
                         _x = self.params['MHP']['xc'], _y = self.params['MHP']['yc']) 
            self.MHP.name = self.params['MHP']['name']
            
            self.d5 =  Drift(self.params["MHE"]['distance from source']-self.params["MHP"]['distance from source'])
            self.d5.name = self.params["d5"]['name']
            
            self.MHE_error = MirPl(np.loadtxt(self.params['MHE_error']['mirror profile']),
             _dim = self.params['MHE_error']['orientation'],
             _ang = self.params['MHE_error']['incidence angle'], 
             _refl = self.params['MHE_error']['transmission'],
             _x = self.params['MHE_error']['xc'], _y = self.params['MHE_error']['yc']) 
            
            self.MHE_error.name = self.params['MHE_error']['name']
            

            self.MVE_error = MirPl(np.loadtxt(self.params['MVE_error']['mirror profile']),
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
            
            self.d6 =  Drift(self.params["MHE"]['distance from source']-self.params["MVE"]['distance from source'])
            self.d6.name = self.params["d6"]['name']
            
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


            self.d8 =  Drift(self.params["MVP"]['distance from source']-self.params["MVE"]['distance from source'])
            self.d8.name = self.params["d8"]['name']
            
            self.MVP = MirPl(np.loadtxt(self.params['MVP']['mirror profile']),
                     _dim = self.params['MVP']['orientation'],
                     _ang = self.params['MVP']['incidence angle'], 
                     _refl = self.params['MVP']['transmission'],
                     _x = self.params['MVP']['xc'], _y = self.params['MVP']['yc']) 
            self.MVP.name = self.params['MVP']['name']
            
            self.df =  Drift(self.params["df"]['distance from source']-self.params["MVP"]['distance from source'])
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
            
            self.NVE_error = MirPl(np.loadtxt(self.params['NVE_error']['mirror profile']),
                            _dim = self.params['NVE_error']['orientation'],
                            _ang = self.params['NVE_error']['incidence angle'], 
                            _refl = self.params['NVE_error']['transmission'],
                            _x = self.params['NVE_error']['xc'], _y = self.params['NVE_error']['yc']) 
            
            self.NVE_error.name = self.params['NVE_error']['name']
            
            self.NHE_error = MirPl(np.loadtxt(self.params['NHE_error']['mirror profile']),
                _dim = self.params['NHE_error']['orientation'],
                _ang = self.params['NHE_error']['incidence angle'], 
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
            
    def buildBeamline(self, focus = "micron", screens = "false"):
        """
        Construct the beamline object
        
        :param focus: what beamline configuration (micron or nano)
        """
        
        self.bl = Beamline()
        
        if focus == "micron":
            
            
            self.bl.append(self.d1, propParams(1,1,1,1, mode = "farfield"))
    
            self.bl.append(self.HOM1, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d2, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.HOM2,  propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d3, propParams(1,1,1,1, mode = 'farfield'))
            
            self.bl.append(self.MKB_pslit, propParams(1/5, 1, 1/5, 1, mode = 'normal'))
            self.bl.append(self.d4, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.MHP, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d5, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.MHE, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.MHE_error, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d6, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.MVE, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.MVE_error, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d8, propParams(1, 1, 1, 1, mode = 'quadratic'))
        
            self.bl.append(self.MVP, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.df, propParams(20, 10, 20, 10, mode = 'converge'))
       
        elif focus == "nano":
            
            self.bl.append(self.d1, propParams(1,1,1,1, mode = "farfield"))
    
            self.bl.append(self.HOM1, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d2, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.HOM2,  propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d3, propParams(1,1,1,1, mode = 'farfield'))
            
            self.bl.append(self.NKB_pslit, propParams(1/5, 1, 1/5, 1, mode = 'normal'))
            self.bl.append(self.d4, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.NHE_error, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.NHE, propParams(1, 1, 1, 1, mode = 'normal'))
            
            
            self.bl.append(self.d5, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.NVE_error, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.NVE, propParams(1, 1, 1, 1, mode = 'normal'))
            
            
            self.bl.append(self.df, propParams(20, 10, 20, 10, mode = 'converge'))

            
        self.bl.params = self.params
        
    def get_beamline(self):
        return self.bl
    
    
    def cropBeamline(self, element1 = None, element2 = None):
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

    def addScreen(self,position, distance, screenName = None):
        """
        add a screening plane at drift beyond some element position
        
        :param position: last optical element before screen
        :param distance: position from last optical element to screen
        :param screenName: name of the screen element (ie., MKB-scr etc) [str]
        """
        self.cropBeamline(element1 = position)
        
        drift2screen = SRWLOptD(distance)
        if screenName is not None:
            drift2screen.name = "screen"
        else:
            drift2screen.name = screenName
        self.bl.append(Drift(distance), propParams(1, 1, 1, 1, m = 'quadratic'))
    
    def mirrorProfiles(self, toggle = "on", aperture = True, overwrite = False):
        """
        toggle for mirror surfaces
        """
        if toggle == "on":
            self.defineMirrorProfiles(overwrite = overwrite, aperture = aperture, surface = 'real')
        if toggle == "off":
            self.defineMirrorProfiles(overwrite = overwrite, aperture = aperture, surface = 'flat')
            
    def plotMirrorProfile(self, mirror, outdir = None):
        
        surface = np.loadtxt(self.params[mirror]['mirror profile'])
        
        x = surface[1:, 0]
        y = surface[0, 1:]
        
        surface = surface[1:,1:]
        
        extent = [x.min()*1e03, x.max()*1e03, y.min()*1e03, y.max()*1e03]
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        

        ax.set_title(mirror + " Surface")

        
        img = ax.imshow(surface*1e9,
                        extent = extent,
                        aspect = 'auto')
        
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        
        cb = plt.colorbar(img, ax = ax)
        cb.ax.get_yaxis().labelpad = 15
        cb.ax.set_ylabel("Height Error (nm)", rotation = 270)
        
        if outdir is not None:
            fig.savefig(outdir + mirror + "_surface.png")
        else:
            plt.show()
    
