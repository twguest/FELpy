from felpy.model.beamline import Beamline
from wpg.optical_elements import Drift
from felpy.model.tools import propagation_parameters
from felpy.model.beamlines.exfel_spb.methods import setup_spb
 


def spb_to_m2(ekev, theta = 5e-03, z = 0):
    """
    api function to return a beamline cropped to the M2 mirror of the SPB/SFX instrument of the European XFEL

    :param ekev: photon beam wavelength - defines mirror transmission
    :param theta: angle of the horizontal mirrors
    :param z: post mirror drift (m)
    """

    spb = setup_spb(parameter_file = "/gpfs/exfel/data/user/guestt/FELpy/felpy/data/params/spb-sfx_nkb_GM_4.98.json",
                            theta_KB = 5e-03, theta_HOM = theta, crop = ["HOM2"],
                            apertures = True, surface = True, ekev =  ekev)


    spb.edit_propagation_parameters("d1", propagation_parameters(1/4, 1, 1/4, 1 ,mode = 'fraunhofer'))
                        
    spb.edit_propagation_parameters("HOM1", propagation_parameters(1/5, 1, 1/5, 1 ,mode = 'fresnel'))
    spb.edit_propagation_parameters("d2", propagation_parameters(1/1, 1, 1/1, 1 ,mode = 'quadratic'))

    spb.edit_propagation_parameters("HOM2", propagation_parameters(1, 1, 1, 1 ,mode = 'fresnel'))


    D = Drift(z)

    D.name = 'post-M2-drift'
    
    spb.bl.append(D, propagation_parameters(1,1,1,1,mode = 'quadratic'))
    
    return spb.bl

def beamline_to_IMPII45(propagation_parameters = propagation_parameters(1/5,1,1/5,1,mode = 'fraunhofer')):
    """ 
    return a wpg beamline that propagates the the IMPII45 imaging screen in the SASE1 undulator
    """
    return get_drift_beamline(z = 230.3, propagation_parameters = propagation_parameters)

    
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