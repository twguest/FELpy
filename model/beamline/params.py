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
import datetime

import numpy as np

def export(params):
    with open('../../data/input/parameters.json', 'w') as f:
        json.dump(params, f)
        
if __name__ == "__main__":

    date = datetime.datetime.now()
    mn,h,d,m,y = date.minute, date.hour, date.day, date.month, date.year

    params = {}
    
    params["d1"] = {'name': "drift1",
                    "description": "Drift from Source to HOM1",
                    'distance': 246.5}
    
    params["HOM1"] = {"name": "HOM1",
                      "description": "First Horizontal Offset Mirror",
                      "mirror profile": "../../data/hom1_mir_Flat.dat",
                      "orientation": 'x',
                      "incidence angle": 2.1e-03,
                      "xc": 0,
                      "yc": 0,
                      "transmission": 1}
    
    params["d2"] = {'name': "drift2",
                    "description": "Drift from HOM1 to HOM2",
                    'distance': 11.36}
    
    params["HOM2"] = {"name": "HOM2",
                      "description": "Second Horizontal Offset Mirror",
                      "mirror profile": "../../data/hom2_mir_Flat.dat",
                      "orientation": 'x',
                      "incidence angle": 2.4e-03,
                      "transmission": 1,
                      "xc": 0,
                      "yc": 0}
    
    params["d3"] = {'name': "drift3",
                    "description": "Drift from HOM2 to Effective Tunnel Entrance",
                    'distance': 634.669}

    params["MKB_pslit"] = {"name": "MKB-Pslit",
                  "description": "Power Slit Aperture Prior to MKB",
                  "shape": 'r',
                  "type": 'a',
                  "dx": 0.0038,
                  "dy": 0.0038,
                  "xc": 0,
                  "yc": 0}

    params["d4"] = {'name': "drift4",
                "description": "Drift MKB-PSLIT to MHP",
                'distance': 1.200}
    
    
    
    params["MHP"] = {"name": "MHP",
                     "mirror profile": "../../data/mhp_mir_Flat.dat",
                     "orientation": 'x',
                     "incidence angle": 1.1e-03,
                     "transmission": 1,
                     "xc": 0,
                     "yc": 0}
    
    params["d5"] = {'name': "drift5",
                    "description": "Drift from MHP to MHE",
                    'distance': 1.050}
    
    params["MHE_ap"] = {"name": "MHE_ap",
              "description": "Aperture of the Horizontal Elliptical Mirror",
              "shape": 'r',
              "type": 'a',
              "dx": 0.950,
              "dy": 0.0250,
              "xc": 0,
              "yc": 0}

    params["MHE"] = {"name": "MHE",
                     "orientation": 'x',
                     "description": "Micron Focus Horizontal Elliptical Mirror",
                     "distance from source":  894.779,
                     "distance to focus": 23.905,
                     "design angle": 1e-03,
                     "incidence angle": 1e-03,
                     "xc": 0,
                     "yc": 0,
                     "length": 1,
                     "roll": 0,
                     "yaw": 0,
                     "reflectivity": 1,
                     "_ext_in": 0.5,
                     "_ext_out": 0.5}
    
    params["d6"] = {'name': "MKB-SCR",
                    "description": "Drift from MHE tp MKB-SCR",
                    'distance': 0.360}
    
    params["d7"] = {'name': "drift7",
                    "description": "Drift from MKB-SCR to MVE",
                    'distance': 1.320}

    params["MVE_ap"] = {"name": "MVE_ap",
          "description": "Aperture of the Horizontal Elliptical Mirror",
          "shape": 'r',
          "type": 'a',
          "dy": 0.950,
          "dx": 0.0250,
          "xc": 0,
          "yc": 0}

    
    params["MVE"] = {"name": "MVE",
                     "orientation": 'y',
                     "description": "Micron Focus Vertica; Elliptical Mirror",
                     "distance from source":  896.459,
                     "distance to focus": 22.225,
                     "design angle": 1e-03,
                     "incidence angle": 1e-03,
                     "xc": 0,
                     "yc": 0,
                     "length": 1,
                     "roll": 0,
                     "yaw": 0,
                     "reflectivity": 1,
                     "_ext_in": 0.5,
                     "_ext_out": 0.5}
    
    params["d8"] = {'name': "drift8",
                    "description": "Drift from MVE to MVP",
                    'distance': 1.050}
    
    params["MVP"] = {"name": "MVP",
                     "description": "Vertical Plane Mirror of MKB",
                     "mirror profile": "../../data/mvp_mir_Flat.dat",
                     "orientation": 'y',
                     "incidence angle": 1.1e-03,
                     "transmission": 1,
                     "xc": 0,
                     "yc": 0}

    params["df"] = {'name': "focus",
                    "description": "Drift to Focus",
                    'distance': 21.175}
    

    
    export(params)
    
        
    
    