from felpy.model.beamlines.exfel_spb.methods import get_beamline_object
from felpy.model.core.beamline import Beamline
from wpg.optical_elements import Drift, Aperture
from felpy.model.tools import propagation_parameters, scale
from felpy.model.core.wavefront import Wavefront
from labwork.about import dCache
from wpg.wpg_uti_wf import plot_intensity_map as plt
    
import sys 


def execute(wfr_directory, sdir, focus, analysis = False, crop = [], append = [], VERBOSE = True):
    
    wfr = Wavefront()
    wfr.load_dir(wdir)
    sdir = sdir + "/{}/".format(focus)
    
    bl = get_beamline_object(ekev = 4.96, options = focus, crop = crop)
    
    for item in append:
        bl.append(item[0],item[1])
    
    wfr.log(bl)
    
    bl.propagate(wfr)
    
    if analysis: 
        wfr.analysis()
        
    wfr.store_hdf5(sdir)
    

    
if __name__ == '__main__':
    
    print("workling")
    sys.argv[1] = wdir
    sys.argv[2] = sdir 
    sys.argv[3] = focus
    sys.argv[4] = analysis
    sys.argv[5] = crop
    sys.argv[6] = append
    execute(wdir, sdir, focus, analysis, crop, append)
    