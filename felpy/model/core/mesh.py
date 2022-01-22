# -*- coding: utf-8 -*-
import numpy as np

class Mesh:
    """
    generalised Mesh class
    
    WIP
    """
    def __init__(self, nx, ny, xMin, xMax, yMin, yMax):
    
        self.nx = nx
        self.ny = ny
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        
        self.xRange = xMax - xMin
        self.yRange = yMax - yMin 