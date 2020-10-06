#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  6 15:08:38 2020

@author: twguest
"""


from setuptools import Extension
from distutils.core import setup, Extension

from distutils.command.install import install as DistutilsInstall
from os import system

class MyInstall(DistutilsInstall):
    def run(self):
        
        DistutilsInstall.run(self)
        system("rm develop")
 
        system("wget --no-check-certificate -nc https://github.com/twguest/WPG/tarball/develop")
        system("mkdir sources/WPG")
        system("tar -xvf  develop -C sources/WPG --strip-components=1")
        system("cd sources/WPG")
        system("pip install -e sources/WPG")
        
        
setup(name='foo',
      version='1.0',
      ext_modules=[],
      cmdclass={'install': MyInstall},
      install_requires= [''],
      )