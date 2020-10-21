# -*- coding: utf-8 -*-

""" 
a simulation of the nanoprobe beamline at the AS
"""

import numpy as np 

from wpg.wavefront import Wavefront
from wpg.beamline import Beamline
from wpg.optical_elements import Drift, Aperture
from wpg.srwlib import SRWLOptMirEl, SRWLOptMirTor

from felpy.model.tools import constructPulse, propParams

from wpg.wpg_uti_wf import plot_intensity_map as plotIntensity

# np.loadtxt("../../data/nanoprobe/B4CRu_mlayer.dat")

wfr = constructPulse(d2waist = 25)

plotIntensity(wfr)


drift_1 = Drift(14)

ap_1 = Aperture(shape = 'r',
                   ap_or_ob = 'ap',
                   Dx = 0.5e-03,
                   Dy = 0.3e-03)


drift_2 = Drift(0.5)
drift_3 = Drift(1.0)

VFM =  SRWLOptMirEl(_p = 15.5,
                    _q = 15.501,
                    _ang_graz = 1)

drift_4 = Drift(1.0)

HFM =  SRWLOptMirEl(_p = 16.5,
                    _q = 14.5,
                    _ang_graz = 1)

drift_5 = Drift(14.5)
drift_6 = Drift(61.0)

VKB =  SRWLOptMirEl(_p = 16.5,
                    _q = 14.5,
                    _ang_graz = 1)

HKB =  SRWLOptMirEl(_p = 16.5,
                    _q = 14.5,
                    _ang_graz = 1)






bl = Beamline()
bl.append(drift_1, propParams(1,1,1,1, mode = 'quadratic'))
bl.append(ap_1, propParams(1,1,1,1, mode = 'fresnel'))
bl.append(drift_2, propParams(1,1,1,1, mode = 'quadratic'))
bl.append(drift_3, propParams(1,1,1,1, mode = 'quadratic'))


bl.append(VFM, propParams(1,1,1,1, mode = 'fresnel'))


bl.append(drift_4, propParams(1,1,1,1, mode = 'quadratic'))
bl.append(HFM, propParams(1,1,1,1, mode = 'fresnel'))

bl.append(drift_5, propParams(1,1,1,1, mode = 'converge'))

bl.append(drift_6, propParams(1,1,1,1, mode = 'quadratic'))

bl.propagate(wfr)

plotIntensity(wfr)
