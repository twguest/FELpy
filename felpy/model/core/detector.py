# -*- coding: utf-8 -*-

"""
hybrid wpg-felpy detector module
"""

from felpy.model.source.coherent import construct_SA1_pulse
from wpg.srwlib import SRWLDet
from felpy.model.core.beamline import Beamline
from felpy.model.tools import propagation_parameters
from wpg.optical_elements import Drift

class Detector:
    
    def __init__(self, dx, dy, nx, ny):
        """
        :param dx: detector pixel size
        :param dy: detector pixel size
        :param nx: number of pixels
        :param ny: number of pixels
        """
        
        self.dx = dx
        self.dy = dy
        self.nx = nx
        self.ny = ny
    
    def detect(self, wfr):
        
        xMin, xMax, yMax, yMin = wfr.get_limits()
        
        idx, idy = xMax-xMin, yMax-yMin
        inx, iny = wfr.params.Mesh.nx, wfr.params.Mesh.ny
        print(idx,idy)
        print((self.dy*self.ny)/idy)
        bl = Beamline()
        bl.append(Drift(0),
                  propagation_parameters((self.dx*self.nx)/idx, self.nx/inx, (self.dy*self.ny)/idy, self.ny/iny))
    
        bl.propagate(wfr)
        

if __name__ == '__main__':
    
    wfr = construct_SA1_pulse(512, 512, 5, 5.0, 0.25)
    
    
    det = Detector(1e-06, 1e-06, 1024, 1024)
    det.detect(wfr)