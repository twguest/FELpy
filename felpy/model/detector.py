# -*- coding: utf-8 -*-

""" 
Small detector module for scaling

"""

from felpy.model.mesh import Mesh
from felpy.model.beamline import Beamline
from wpg.optical_elements import Drift
from felpy.model.tools import propagation_parameters

class Detector:
    
    def __init__(self, **kwargs):
        
        self.__dict__.update(kwargs)

        ### this means that we first defer to nx/ny/nz definitions of the field
        if all(hasattr(self, attr) for attr in ["nx", "ny", "dx", "dy"]):
            pass
        else:
            if hasattr(self, "mesh"):
                self.__dict__.update(self.mesh.get_attributes())
        
    def detect(self, wfr):
        """
        return the detected intensity distribution
        """
 
        
        bl = Beamline()
        bl.append(Drift(0), propagation_parameters(wfr.dx/self.mesh.dx, self.mesh.nx/wfr.nx,
                                                  wfr.dy/self.mesh.dy, self.mesh.ny/wfr.ny,
                                                  mode = 'fresnel')
                  )
        bl.propagate(wfr)
        
        return wfr.get_intensity().sum(-1)
    

if __name__ == '__main__':
    
    from felpy.model.src.coherent import construct_SA1_pulse
    
    wfr = construct_SA1_pulse(512,512,3,10,0.25)
    [wxMin, wxMax, wyMax, wyMin] = wfr.get_limits()
    wnx, wny = wfr.nx, wfr.ny
    wdx, wdy = wfr.dx, wfr.dy
    

    
    detector_mesh = Mesh(nx = 1024, ny = 1024, dx = 1e-06, dy = 1e-06)
    
    d = Detector(mesh = detector_mesh)
    i = d.detect(wfr)
    
    print("Desired Properties/Actual Values")
    print("nx: {}/{}".format(detector_mesh.nx, wfr.params.Mesh.nx))
    print("ny: {}/{}".format(detector_mesh.ny, wfr.ny))
    print("dx: {}/{}".format(detector_mesh.dx, wfr.dx))
    print("dy: {}/{}".format(detector_mesh.dy, wfr.dy))