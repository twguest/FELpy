# -*- coding: utf-8 -*-
import os
from wpg.beamline import Beamline as WPG_Beamline
from wpg import srwlib
from wpg.srw import srwlpy as srwl
from wpg.wpg_uti_wf import plot_intensity_map

import pandas as pd

class Beamline(WPG_Beamline):
    
    def __init__(self, _srwl_wf = None):
        super().__init__(_srwl_wf)
      
    @property
    def index(self):
        index = {}
        
        for itr, a in enumerate(self.propagation_options[0]['optical_elements']):
            
            try:
                index[a.name] = itr
            except(AttributeError):
                index['Element {}'.format(itr)] = itr
            
        return index
    
    def print_element(self, oe_name):
        print(self.propagation_options[0]['optical_elements'][self.index[oe_name]].__dict__)
    
    def print_propagation_parameters(self, oe_name):
        print(self.propagation_options[0]['propagation_parameters'][self.index[oe_name]])
        
    def replace_element(self,oe_name, el):
        """ 
        replace an optical element labelled oe_name w/ a pre-defined element oe.
        """
        self.propagation_options[0]['propagation_parameters'][self.index[oe_name]] = el
    
    def edit_element_property(self, oe_name, prop, value):
        self.propagation_options[0]['optical_elements'][self.index[oe_name]].__dict__[prop] = value
    
    def propagate_sequential(self, wfr, return_intensity = False, return_mesh = False, savedir = "", checkpoints = []):
        """
        Propagate sequentially through each optical element in beamline.

        :param wfr: Input wavefront (will be re-writed after propagation)
        :param outdir: save directory
        """
               
        
        if return_intensity:
            intensity = []
        if return_mesh:
            mesh = []
            
        for itr in range(len(self.propagation_options[0]['optical_elements'])):
            oe = self.propagation_options[0]['optical_elements'][itr]
            pp = self.propagation_options[0]['propagation_parameters'][itr]
            
            bl = srwlib.SRWLOptC([oe],[pp])
            
            try:
                print(oe.name)
            except:
                print(oe)
             
            
            srwl.PropagElecField(wfr._srwl_wf, bl)
            
            if return_intensity:
                intensity.append(wfr.get_intensity().sum(-1))
            if return_mesh:
                mesh.append(wfr.get_mesh())
                    
            plot_intensity_map(wfr)
            
        if return_intensity:
            
            if return_mesh:
                return intensity, mesh
            else:
                return intensity
            
            
if __name__ == '__main__':

    from felpy.model.src.coherent import construct_SA1_pulse
    from felpy.model.tools import propagation_parameters
    from wpg.optical_elements import Drift
    wfr = construct_SA1_pulse(512,512,2,5.0,1)

    bl = Beamline()
    D = Drift(10)
    D.name = 'd'
    bl.append(D, propagation_parameters(2,1,2,1,'quadratic'))
    bl.append(Drift(1), propagation_parameters(1,1,1,1,'quadratic'))
    #bl.propagate_sequential(wfr)
