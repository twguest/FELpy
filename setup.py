#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 15:08:38 2020

@author: twguest
"""

from os import system

from distutils.core import setup

from distutils.command.install import install as DistutilsInstall


class MyInstall(DistutilsInstall):
    def run(self):
        
        DistutilsInstall.run(self)
        system("echo making WPG")
        system("rm develop")
        system("git submodule update --init")
        system("cd FELpy/WPG; make")
        system("pip install -e FELpy/WPG")
        system("cd ../; pip install -e FELpy")
        
setup(name='FELpy',
      version='0.1.1',
      ext_modules=[],
      cmdclass={'install': MyInstall},
      install_requires= ['h5py>=2.10.0',
                         'imageio>=2.8.0',
                         'ipython>=7.13.0',
                         'joblib>=0.14.1',
                         'jupyter_client>=6.1.2',
                         'jupyter_core>=4.6.3',
                         'matplotlib>=3.1.1',
                         'numpy>=1.18.1',
                         'pandas>=1.0.3',
                         'pillow>=7.0.0',
                         'pip>=20.0.2',
                         'scikit-learn>=0.22.1',
                         'scipy>=1.4.1',
                         'spyder>=4.1.0',
                         'tqdm>=4.46.0',
                         'wheel>=0.34.2'
                         ],
      
      packages=[ 'felpy',
                 'felpy.data',
                 'felpy.exp',
                 'felpy.utils',
                 'felpy.model',
                 'felpy.data.spb-sfx',
                 'felpy.data.samples',
                 'felpy.data.spb-sfx.hom_refl',
                 'felpy.data.spb-sfx.kb_refl',
                 'felpy.data.spb-sfx.mirror_surface',
                 'felpy.analysis',
                 'felpy.model.src',
                 'felpy.model.materials',
                 'felpy.model.beamline'
                ]
      )