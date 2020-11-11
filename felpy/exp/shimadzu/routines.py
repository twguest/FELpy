# -*- coding: utf-8 -*-
"""
some routines that run successively for shimadzu whitefield data analysis
"""
import os 
import numpy as np
from matplotlib import pyplot as plt
from felpy.utils.np_utils import get_mesh
from felpy.exp.shimadzu.preprocess import preprocess_shimadzu_data, shimadzu_test_data
from felpy.utils.vis_utils import extract_animation

from felpy.exp.shimadzu.methods import correlation_method_a, correlation_method_b, correlation_method_c

def burger_w_the_lot():
    pass


def extract_whitefield_animation(ii, mesh, sdir, run, modes = ['all']):
       
    for m in modes:
        if os.path.exists(sdir + run + "_{}.gif".format(m)):
            pass
        else:
            print("extracting all images for .gif")
            print("saving to directory: {}{}".format(sdir, run + "_{}.gif".format(m)))
            
            extract_animation(ii, mesh, fname = run + "_{}".format(m), sdir = sdir,
                              mode = m)
            
# =============================================================================
# def correlation_analysis(ii, mesh, method, sdir, desc, mpi = False):
#     
#     
#     corr = method(ii, mpi = False)
#         
#     for train in range(arr.shape[-1]):
#            
#         tdir = sdir + "/tmp/"
#         mkdir_p(tdir)
#         
#         for pulse in range(tcorr.shape[-2]):
# 
#             correlation_plot(tcorr[:,:,pulse,train], mesh,
#                              sdir = tdir + "train_{:04d}_pulse_{:04d}.png".format(train, pulse),
#                              label = "train_{}".format(train))
#             
#     
#         animate(indir = tdir, outdir = sdir,
#                 fname = fname + "_sequential_correlation",
#                 delay = 0.06,
#                 rmdir = True)
#         
#     np.save(sdir + "/sequential_correlation", tcorr)
# =============================================================================

if __name__ == '__main__':
    
    px = 1e-06
    
    ii = shimadzu_test_data(250,250,2,5)
    mesh = get_mesh(ii, px,px)
    extract_whitefield_animation(ii, mesh, sdir = '/opt/labwork/scratch/',
                                 run = 'r0069')
    
    