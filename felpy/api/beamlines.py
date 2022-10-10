from felpy.model.beamline import Beamline
from wpg.optical_elements import Drift
from felpy.model.tools import propagation_parameters

def get_drift_beamline(z = 1, propagation_parameters = propagation_parameters(1,1,1,1, mode = 'quadratic'), drift_name = 'free-space'):
    """
    returns a beamline w/ drift of distance z

    :param z: beamline drift distance
    :param propagation_parameters: wpg propagation parameters

    :returns bl: free space drift beamline
    """

    bl = Beamline()
    
    D = Drift(z)
    D.name = drift_name

    bl.append(D, propagation_parameters)
    
    return bl