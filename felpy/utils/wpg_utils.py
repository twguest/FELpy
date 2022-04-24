# -*- coding: utf-8 -*-
import numpy as np


def complex_to_wpg(arr): ### converter
    new_arr = np.zeros([arr.shape[0], arr.shape[1], arr.shape[2], 2])
    new_arr[:,:,:,0] = arr.real
    new_arr[:,:,:,1] = arr.imag
    return new_arr
