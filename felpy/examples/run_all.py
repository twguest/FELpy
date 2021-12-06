#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
FELPY

__author__ = "Trey Guest"
__credits__ = ["Trey Guest"]
__license__ = "EuXFEL"
__version__ = "0.1.1"
__maintainer__ = "Trey Guest"
__email__ = "twguest@students.latrobe.edu.au"
__status__ = "Developement"
"""

from felpy.tests.coherent_source.coherent_source import core as coherent_source_test
from felpy.tests.mirror_angle_optimisation.mirror_angle_optimisation import core as mirror_angle_optimisation
from felpy.tests.mirror_profiles.mirror_profiles import core as test_mirror_profiles

if __name__ == '__main__':
    coherent_source_test()
    mirror_angle_optimisation()
    test_mirror_profiles()