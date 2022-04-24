#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.2.1"
__maintainer__ = "Trey Guest"
__email__ = "trey.guest@xfel.eu"
__status__ = "Developement"
"""

from felpy.model.wavefront import Wavefront


def wavefront_test():
    
    assert Wavefront()
    
def wavefront_complex_test():
    pass

if __name__ == '__main__':
    
    wavefront_test()
    wavefront_complex_test()