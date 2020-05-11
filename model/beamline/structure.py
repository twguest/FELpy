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

from wpg.beamline import Beamline
from wpg import srwlib
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptA as Aperture

from wpg.optical_elements import Mirror_elliptical as MirEl

from model.materials.load_refl import get_refl, load_refl
from model.src.coherent import coherentSource


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


def genMirrorSurface(nx, ny, mirDim, outdir, mode = 'Flat'):
    """
    Generate a plane mirror surface
    
    :param nx: number of horizontal pixels [int]
    :param ny: number of vertical pixels [int]
    :param mirDim: list of mirror dimensions [dx,dy] [m]
    :param outdir: save directory
    :param mode: type of mirror surface to be generated
    """
    if mode == 'Flat':
        
        mirLen = mirDim[0]
        mirWid = mirDim[1]
        
        surface = np.zeros((nx,ny))
        
        surface[1:, 0] = np.linspace(-mirLen/2, mirLen/2, nx-1) 
        surface[0, 1:] = np.linspace(-mirWid/2, mirWid/2, ny-1) 
        
    np.savetxt(outdir+"mir_"+ mode +".dat", surface, delimiter='\t')


class BeamlineModel:
    
    """
    A container for loading and modifying beamline data
    """
    
    
    def __init__(self, overwrite_mirrors = False):
        print("Initialising Single Particle Beamline")
        self.load_params()
        self.defineMirrorProfiles(overwrite = overwrite_mirrors)
        
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
            
    def defineMirrorProfiles(self, overwrite = False):
        """
        Define the plane mirror profiles by loading from /data/. 
        If mirror profiles not defined (unlikely), generate profiles via genMirrorSurface
        
        :param overwrite: bool to overwrite current file.
        """
        
 
        if exists("../../data/hom1_mir_Flat.dat") and overwrite == False:
            self.hom1_profile = "../../data/hom1_mir_Flat.dat"
        elif overwrite == True:
            genMirrorSurface(500, 500, [0.010, 0.800], "../../data/hom1_", mode = 'Flat')
            self.hom1_profile = "../../data/hom1_mir_Flat.dat"
        else:
            genMirrorSurface(500, 500, [0.010, 0.800], "../../data/hom1_", mode = 'Flat')
            self.hom1_profile = "../../data/hom1_mir_Flat.dat"
            
        ### HOM2    
        if exists("../../data/hom2_mir_Flat.dat") and overwrite == False:
            self.hom2_profile = "../../data/hom2_mir_Flat.dat"
        elif overwrite == True:
            genMirrorSurface(500, 500, [0.010, 0.800], "../../data/hom1_", mode = 'Flat')
            self.hom1_profile = "../../data/hom1_mir_Flat.dat"
        else:
            genMirrorSurface(500, 500, [0.010, 0.800], "../../data/hom2_", mode = 'Flat')
            self.hom2_profile = "../../data/hom2_mir_Flat.dat"
        
        ### MHP
        if exists("../../data/mhp_mir_Flat.dat") and overwrite == False:
            self.mhp_profile = "../../data/mhp_mir_Flat.dat"
        else:
            genMirrorSurface(500, 500, [0.025, 0.950 ], "../../data/mhp_", mode = 'Flat')
            self.mhp_profile = "../../data/mhp_mir_Flat.dat"
        
        ### MHE                
        if exists("../../data/mvp_mir_Flat.dat") and overwrite == False:
            self.mvp_profile = "../../data/mvp_mir_Flat.dat"
        else:
           genMirrorSurface(500, 500, [0.025, 0.950 ], "../../data/mvp_", mode = 'Flat')
           self.mvp_profile = "../../data/mvp_mir_Flat.dat"
 

    def buildElements(self, focus = "micron", screens = False):
        
        if focus == "micron" and screens == False:
            self.d1 =  Drift(self.params["d1"]['distance'])
            self.d1.name = self.params["d1"]['name']
            
            self.HOM1 = MirPl(np.loadtxt(self.params['HOM1']['mirror profile']),
                         _dim = self.params['HOM1']['orientation'],
                         _ang = self.params['HOM1']['incidence angle'], 
                         _amp_coef = 1,
                         _x = self.params['HOM1']['xc'], _y = self.params['HOM1']['yc'],
                         _refl = self.params['HOM1']['transmission']) 
            
            self.HOM1.name = self.params['HOM1']['name']
            
            self.d2 =  Drift(self.params["d2"]['distance'])
            self.d2.name = self.params["d2"]['name']
            
            self.HOM2 = MirPl(np.loadtxt(self.params['HOM2']['mirror profile']),
                         _dim = self.params['HOM2']['orientation'],
                         _ang = self.params['HOM2']['incidence angle'], 
                         _amp_coef = 1,
                         _x = self.params['HOM2']['xc'], _y = self.params['HOM2']['yc'],
                         _refl = self.params['HOM2']['transmission']) 
            
            self.HOM2.name = self.params['HOM2']['name']
            
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
            
            self.MHE_ap = Aperture(_shape=self.params["MHE_ap"]['shape'],
                                 _ap_or_ob=self.params["MHE_ap"]['type'],
                                 _Dx= self.params["MHE_ap"]['dx'],
                                 _Dy= self.params["MHE_ap"]['dy'],
                                 _x=self.params["MHE_ap"]['xc'],
                                 _y=self.params["MHE_ap"]['yc'])
            self.MHE_ap.name = self.params["MHE_ap"]['name']
            
            self.MHP = MirPl(np.loadtxt(self.params['MHP']['mirror profile']),
                         _dim = self.params['MHP']['orientation'],
                         _ang = self.params['MHP']['incidence angle'], 
                         _refl = self.params['MHP']['transmission'],
                         _x = self.params['MHP']['xc'], _y = self.params['MHP']['yc']) 
            self.MHP.name = self.params['MHP']['name']
            
            self.d5 =  Drift(self.params["d5"]['distance'])
            self.d5.name = self.params["d5"]['name']
            
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
            
            self.d6 =  Drift(self.params["d6"]['distance']+self.params["d7"]['distance'])
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

            self.MVE_ap = Aperture(_shape=self.params["MVE_ap"]['shape'],
                     _ap_or_ob=self.params["MVE_ap"]['type'],
                     _Dx= self.params["MVE_ap"]['dx'],
                     _Dy= self.params["MVE_ap"]['dy'],
                     _x=self.params["MVE_ap"]['xc'],
                     _y=self.params["MVE_ap"]['yc'])
            self.MVE_ap.name = self.params["MVE_ap"]['name']
    
            self.d8 =  Drift(self.params["d8"]['distance'])
            self.d8.name = self.params["d8"]['name']
            
            self.MVP = MirPl(np.loadtxt(self.params['MVP']['mirror profile']),
                     _dim = self.params['MVP']['orientation'],
                     _ang = self.params['MVP']['incidence angle'], 
                     _refl = self.params['MVP']['transmission'],
                     _x = self.params['MVP']['xc'], _y = self.params['MVP']['yc']) 
            self.MVP.name = self.params['MVP']['name']
            
            self.df =  Drift(self.params["df"]['distance'])
            self.df.name = self.params["df"]['name']
            
    def buildBeamline(self, focus = "micron", screens = "false"):
        """
        Construct the beamline object
        
        :param focus: what beamline configuration (micron or nano)
        :param screens: include sensing planes [bool]
        """
        if focus == "micron":
            
            self.bl = Beamline()
            self.bl.append(self.d1, propParams(1,1,1,1, mode = "farfield"))
    
            self.bl.append(self.HOM1, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d2, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.HOM2,  propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d3, propParams(1,1,1,1, mode = 'farfield'))
            
            self.bl.append(self.MKB_pslit, propParams(1/5, 1, 1/5, 1, mode = 'normal'))
            self.bl.append(self.d4, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.MHP, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d5, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.MHE_ap, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.MHE, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d6, propParams(1, 1, 1, 1, mode = 'quadratic'))
            self.bl.append(self.MVE_ap, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.MVE, propParams(1, 1, 1, 1, mode = 'normal'))
            self.bl.append(self.d8, propParams(1, 1, 1, 1, mode = 'quadratic'))
        
            self.bl.append(self.MVP, propParams(1, 1, 1, 1, mode = 'normal'))
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

    def addScreen(position, distance, screenName = None):
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
        self.bl.append(SRWLOptD(distance), propParams(1, 1, 1, 1, m = 'quadratic'))

if __name__ == '__main__':

    wfr = coherentSource(1048, 1048, 16, 1)

    spb = BeamlineModel()
    spb.buildElements(focus = "micron", screens = False)
    
    bl = spb.buildBeamline(focus = "micron", screens = False)

     