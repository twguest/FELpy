from felpy.model.beamline import Beamline
from wpg.optical_elements import Drift
from felpy.model.tools import propagation_parameters
from felpy.model.beamlines.exfel_spb.methods import setup_spb
from felpy.model.beamlines.exfel_spb.exfel_spb import Instrument
from felpy.model.materials.load_refl import load_refl, get_refl

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



def direct_beamline(ekev = 5.0, z = 15):
    
    print("Setting up SPB/SFX Direct Beamline Configuration")
    print("Observation Point: {} m downstream of NanoKB Powerslits".format(z))
    
    spb = setup_spb(parameter_file = "/gpfs/exfel/data/user/guestt/FELpy/felpy/data/params/spb-sfx_nkb_GM_4.98.json",
                              theta_KB = 3.5e-03, theta_HOM = 2.5e-03, crop = ["NKB_PSlit"],
                             apertures = True, surface = True, ekev =  ekev)
    
    print("Instrument Setup")
    
    spb.edit_propagation_parameters("d1", propagation_parameters(1/4, 1, 1/4, 1 ,mode = 'fraunhofer'))
    spb.edit_propagation_parameters("HOM1", propagation_parameters(1/2, 1, 1/2, 1 ,mode = 'fresnel'))
    spb.edit_propagation_parameters("d2", propagation_parameters(1/1, 1, 1/5, 1 ,mode = 'quadratic'))
    spb.edit_propagation_parameters("HOM2", propagation_parameters(1, 1, 1, 1 ,mode = 'fresnel'))
    spb.edit_propagation_parameters("d3", propagation_parameters(1.25, 1, 2, 1 ,mode = 'quadratic'))
    spb.edit_propagation_parameters("NKB_PSlit", propagation_parameters(1, 1, 1, 1 ,mode = 'fresnel'))
    
    D = Drift(z)
    spb.bl.append(D, propagation_parameters(1,1,1,1,mode = 'quadratic'))
    
    return spb.bl




def setup_spb(parameter_file = "../../../data/params/spb-sfx_nkb_FAST.json", options = 'nano', ekev = 5.0,
              apertures = True, surface = 'real', crop = None,
              theta_HOM = 2.3e-03, theta_KB = 3.5e-03,
              save_params = False):

    """
    return desired beamline
    """


    spb = Instrument(parameter_file = parameter_file)



    mirrors = spb.mirrors

    spb.mirror_profiles(surface = surface, aperture = apertures)

    for mirror in mirrors:

        if mirror in spb.focus:
            spb.adjust_mirror(mirror,
            ekev,
            theta_KB)

        else:
            spb.adjust_mirror(mirror,
            ekev,
            theta_HOM)


    if ekev <= 7.5:
        material = "B4C"
    else:
        material = "Ru"

    if apertures == False: ## aperture claus

        for mirror in mirrors:

            spb.params[mirror]["dx"] = 5
            spb.params[mirror]["dy"] = 5
            spb.params[mirror]["mirror profile"] = generate_infinite_mirror()
    
        for focus in spb.focus:
            
            spb.params[focus]["design angle"] = np.pi/3 ### [A]s ## fix for now, should onkly accept elliptical mirrors
            spb.params[focus]["incidence angle"] = np.pi/3 ### should be for elliptical mirror surfaces


            
            
    for mirror in mirrors:
        if mirror in spb.focus:
            spb.params[mirror]['reflectivity'] = get_refl(load_refl(material), ekev, theta_KB)
        else:
            spb.params[mirror]['reflectivity'] = get_refl(load_refl(material), ekev, theta_HOM)

    ### TO-DO: make a choice as to wether edits to the beamline will be via params or via beamline object
    ### TO-DO: eventually, exfel_spb should be a sub-class of the Instrument class (formally)
    ### TO-DO: a front-end script to load and label mirrors that ships to a json file would be useful.



    spb.build_elements(focus = options)
    spb.build_beamline(focus = options)
    
    if apertures == False and surface == 'flat':
        spb.remove_element("NVE_error")
        spb.remove_element("NHE_error")
        spb.remove_element("NKB_PSlit")
        

    if save_params:
        spb.export_params()


    if crop is not None:
        if type(crop) == list:
            spb.crop_beamline(*crop)
        else:
            spb.crop_beamline(crop)


    return spb
