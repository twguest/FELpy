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
import datetime

import numpy as np

def export(params, filename):
    with open('../../../data/params/{}.json'.format(filename), 'w') as f:
        json.dump(params, f)

def get_params():
    
    date = datetime.datetime.now()
    mn,h,d,m,y = date.minute, date.hour, date.day, date.month, date.year

    params = {}
    
    params["d1"] = {'name': "d1",
                    "description": "Drift from Source to HOM1",
                    'distance': 246.5,
                    'pp': [1,2.5,1,2.5, "quadratic"]}
    
    params["HOM1"] = {"name": "HOM1",
                      "distance from source": 246.5,
                      "description": "First Horizontal Offset Mirror",
                      "mirror profile": "../../data/spb/mirror_surface/HOM1_mir_real.dat",
                      "orientation": 'x',
                      "incidence angle": 2.1e-03,
                      "xc": 0,
                      "yc": 0,
                      "reflectivity": 1,
                      'next_drift' : 'd2',
                      'ang_min': 1.1e-03,
                      'ang_max': 3.6e-03,
                      'dx': 0.525,
                      'dy': 0.010,
                      'pp': [1, 1, 1, 1,'fresnel']}
    
    params["d2"] = {'name': "d2",
                    "description": "Drift from HOM1 to HOM2",
                    'distance': 11.36,
                    'pp': [1, 1, 1, 1,'quadratic']}
    
    params["HOM2"] = {"name": "HOM2",
                      "description": "Second Horizontal Offset Mirror",
                      "distance from source": 257.86,
                      "mirror profile": "../../data/spb/mirror_surface/HOM2_mir_real.dat",
                      "orientation": 'x',
                      "incidence angle": 2.4e-03,
                      "reflectivity": 1,
                      "xc": 0,
                      "yc": 0,
                      'next_drift' : 'd3',
                      'ang_min': 1.1e-03,
                      'ang_max': 3.6e-03,
                      'dx': 0.525,
                      'dy': 0.010,
                      'pp': [1, 1, 1, 1,'fresnel']}
    
    params["d3"] = {'name': "d3",
                    "description": "Drift from HOM2 to Effective Tunnel Entrance",
                    'distance': 634.669,
                      'pp': [2, 1, 2, 1,'fraunhofer']}

    params["MKB_pslit"] = {"name": "MKB-Pslit",
                  "description": "Power Slit Aperture Prior to MKB",
                  "distance from source": 892.529,
                  "shape": 'r',
                  "type": 'a',
                  "dx": 0.0038,
                  "dy": 0.0038,
                  "xc": 0,
                  "yc": 0,
                  'next_drift': 'd4',
                      'pp': [1, 1, 1, 1,'fresnel']}

    params["d4"] = {'name': "d4",
                "description": "Drift MKB-PSLIT to MHP",
                'distance': 1.200,
                      'pp': [1, 1, 1, 1,'fraunhofer']}
    
    
    
    params["MHP"] = {"name": "MHP",
                     "mirror profile": "../../data/spb/mirror_surface/MHP_mir_flat.dat",
                     "distance from source": 893.729,
                     "orientation": 'x',
                     "incidence angle": 1.1e-03,
                     "transmission": 1,
                     "dx": 0.0250,
                     "dy": 0.950,
                     "xc": 0,
                     "yc": 0,
                     'next_drift': 'd5',
                      'ang_min': -0.5e-03,
                      'ang_max': 5.5e-03,
                      'pp': [1, 1, 1, 1,'fresnel']}
    
    params["d5"] = {'name': "d5",
                    "description": "Drift from MHP to MHE",
                    'distance': 1.050,
                      'pp': [1, 1, 1, 1,'quadratic']}
    




    params["MHE"] = {"name": "MHE",
                     "orientation": 'x',
                     "distance from source": 894.779,
                     "description": "Micron Focus Horizontal Elliptical Mirror",
                     "distance from source":  894.779,
                     "distance to focus": 23.905,
                     "design angle": 1e-03,
                     "incidence angle": 1e-03,
                     "dx": 0.950,
                     "dy": 0.0250,
                     "xc": 0,
                     "yc": 0,
                     "length": 1,
                     "roll": 0,
                     "yaw": 0,
                     "reflectivity": 1,
                     "_ext_in": 0.5,
                     "_ext_out": 0.5,
                     'next_drift' : 'd6',
                     'ang_min': -0.5e-03,
                     'ang_max': 5.5e-03,
                      'pp': [1, 1, 1, 1,'fresnel']}
    
    params["MHE_error"] = {"name": "MHE_error",
                           "orientation": 'x',
                           "description": "Micron Focus Horizontal Elliptical Mirror Surface Height Error",
                           "mirror profile": "../../data/spb/mirror_surface/MHE_mir_flat.dat",
                           "incidence angle": np.pi/2,
                           "transmission": 1,
                           "xc": 0,
                           "yc": 0,
                      'pp': [1, 1, 1, 1,'fresnel']}

    params["d6"] = {'name': "d6",
                    "description": "Drift from MKB-SCR to MVE",
                    'distance': 1.680,
                      'pp': [1, 1, 1, 1,'fresnel']}


    
    params["MVE"] = {"name": "MVE",
                     "orientation": 'y',
                     "distance from source": 896.459,
                     "description": "Micron Focus Vertica; Elliptical Mirror",
                     "distance from source":  896.459,
                     "distance to focus": 22.225,
                     "design angle": 1e-03,
                     "incidence angle": 1e-03,
                     "dx": 0.0250,
                     "dy": 0.950,
                     "xc": 0,
                     "yc": 0,
                     "length": 1,
                     "roll": 0,
                     "yaw": 0,
                     "reflectivity": 1,
                     "_ext_in": 0.5,
                     "_ext_out": 0.5,
                     'next_drift' : 'd7',
                     'ang_min': -5.0e-03,
                     'ang_max': 5.0e-03,
                      'pp': [1, 1, 1, 1,'fresnel']}
    
        
    params["MVE_error"] = {"name": "MVE_error",
                       "description": "Micron Focus Horizontal Elliptical Mirror Surface Height Error",
                       "mirror profile": "../../data/spb/mirror_surface/MVE_mir_flat.dat",
                       "orientation": 'y',
                       "incidence angle": np.pi/2,
                       "transmission": 1,
                       "xc": 0,
                       "yc": 0,
                      'pp': [1, 1, 1, 1,'fresnel']}

    
    params["d7"] = {'name': "d7",
                    "description": "Drift from MVE to MVP",
                    'distance': 1.050,
                      'pp': [1, 1, 1, 1,'fresnel']}
    
    params["MVP"] = {"name": "MVP",
                     "description": "Vertical Plane Mirror of MKB",
                     "distance from source": 897.509,
                     "mirror profile": "../../data/spb/mirror_surface/MVP_mir_flat.dat",
                     "orientation": 'y',
                     "incidence angle": 1.1e-03,
                     "transmission": 1,
                     "dx": 0.0250,
                     "dy": 0.950,
                     "xc": 0,
                     "yc": 0,
                     'next_drift': 'df',
                     'ang_min': -2.0e-03,
                     'ang_max': -2.0e-03,
                      'pp': [1, 1, 1, 1,'fresnel']}

    params["df"] = {'name': "focus",
                    "description": "Drift to Focus",
                    "distance from source": 918.684,
                    'distance': 21.175,
                      'pp': [1, 1, 1, 1,'fresnel']}
    
    params["NHE"] = {"name": "NHE",
                     "orientation": 'x',
                     "distance from source": 915.484,
                     "description": "Nano Focus Horizontal; Elliptical Mirror",
                     "mirror profile": "../../data/spb/mirror_surface/NHE_mir_real.dat",
                     "distance to focus": 3.2,
                     "design angle": 1e-03,
                     "incidence angle": 1e-03,
                     "dx": 0.950,
                     "dy": 0.0250,
                     "xc": 0,
                     "yc": 0,
                     "length": 1,
                     "roll": 0,
                     "yaw": 0,
                     "reflectivity": 1,
                     "_ext_in": 0.5,
                     "_ext_out": 0.5,
                     'next_drift': 'd5',
                     'ang_min': 0.5e-03,
                     'ang_max': 5.5e-03,
                      'pp': [1, 1, 1, 1,'fresnel']}
    
    

    params["NHE_error"] = {"name": "NHE_error",
                       "orientation": 'x',
                       "description": "Nano Focus Horizontal Elliptical Mirror Surface Height Error",
                       "mirror profile": "../../data/spb/mirror_surface/NHE_mir_flat.dat",
                       "orientation": 'x',
                       "incidence angle": np.pi/2,
                       "transmission": 1,
                       "xc": 0,
                       "yc": 0,
                      'pp': [1, 1, 1, 1,'fresnel']}

    params["NVE"] = {"name": "NVE",
                     "orientation": 'y',
                     "distance from source": 916.484,
                     "description": "Nano Focus Vertical; Elliptical Mirror",
                     "mirror profile": "../../data/spb/mirror_surface/NVE_mir_real.dat",
                     "distance from source":  916.484,
                     "distance to focus": 2.2,
                     "design angle": 1e-03,
                     "incidence angle": 1e-03,
                     "dx": 0.0250,
                     "dy": 0.950,
                     "xc": 0,
                     "yc": 0,
                     "length": 1,
                     "roll": 0,
                     "yaw": 0,
                     "reflectivity": 1,
                     "_ext_in": 0.5,
                     "_ext_out": 0.5,
                     'next_drift': 'df',
                     'ang_min': 0.5e-03,
                     'ang_max': 5.0e-03,
                      'pp': [1, 1, 1, 1,'fresnel']}
    
    
    params["NVE_error"] = {"name": "NVE_error",
                       "orientation": 'x',
                       "description": "Nano Focus Horizontal Elliptical Mirror Surface Height Error",
                       "mirror profile": "../../data/spb/mirror_surface/NVE_mir_flat.dat",
                       "orientation": 'x',
                       "incidence angle": np.pi/2,
                       "transmission": 1,
                       "xc": 0,
                       "yc": 0,
                      'pp': [1, 1, 1, 1,'fresnel']}

    
    params["NKB_PSlit"] = {"name": "NKB_PSlit",
                           "description": "Power Slit Aperture Prior to NKB",
                           "shape": 'r',
                           "type": 'a',
                           "distance from source":  914.284,
                           "dx": 0.0038,
                           "dy": 0.0038,
                           "xc": 0,
                           "yc": 0,
                           'next_drift': 'd4',
                      'pp': [1/15, 2, 1/15, 2,'fresnel']}

    return params
                 
if __name__ == "__main__":
    params = get_params()
    export(params, filename = "spb-sfx_nkb_FAST")
    
        
    
    