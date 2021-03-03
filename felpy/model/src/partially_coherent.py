#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 14:21:55 2021

@author: twguest
"""

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from felpy.utils.np_utils import gaussian_2d

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

nSamples = 200
time = np.linspace(-100,100,nSamples)
pulse_time = 40e-15
sampling_interval = (2*np.pi)//pulse_time


random_phases = np.ones(nSamples, dtype = 'complex64')*np.random.rand(nSamples)*(2*np.pi)-np.pi
plt.plot(random_phases.imag)
plt.show()

envelope = gaussian(time, 0, 35)
plt.plot(envelope)
spectral_env = envelope*np.exp(1j*random_phases)
temporal_env = np.fft.fft(spectral_env)*envelope

plt.plot(time, abs(temporal_env)**2)
plt.show()

profile = gaussian_2d(100, 100)+1j*gaussian_2d(100, 100)
plt.imshow(abs(profile**2))
plt.show()

wfr = np.zeros([100,100,nSamples], dtype = 'complex64')

for itr in range(nSamples):
    wfr[:,:,itr] = profile
plt.imshow(wfr.mean(-1).imag)
plt.show()

wfr = wfr*temporal_env

plt.plot(abs(wfr.mean(0).mean(0))**2)
plt.show()

plt.imshow(wfr.mean(-1).imag)
plt.imshow(abs(wfr.mean(-1)**2))
plt.show()
