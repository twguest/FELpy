#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 15:55:39 2020

@author: twguest
"""
import sys
sys.path.append("/opt/spb_model/")

import os

import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

import pandas as pd

def bsi():
    """
    work with integrated beam size
    """
    
    bs = np.load("../../data/NanoKB/bsi.npy")[:1,:2,:].reshape([70,2])*1e6
    
    ba = np.pi * bs[:,0] * bs[:,1]
    print("average beam area: {} um^2".format(np.mean(ba)))
    print("std. dev in beam area: {} um^2".format(np.std(ba)))
    return bs

def bss():
    """
    work with pulsed beam size
    """
    
    bs = np.load("../../data/NanoKB/bss.npy")

    for pulse in range(6,7):    
        plt.scatter(bs[:, 1,pulse], bs[:,3 ,pulse])
    return bs

def cni():
    """
    work with pulse centroid integrated
    """
    bs = np.load("../../data/NanoKB/cni.npy")

    print(np.std(bs[:,0,:]))
    print(np.std(bs[:,1,:]))
    
    return bs

def cns():
    """
    work with pulse centroid integrated
    """
    bs = np.load("../../data/NanoKB/cnp.npy")

    print()
    print(np.std(bs[:,1,], axis = 1))
    
    plt.plot(np.std(bs[:,0,], axis = 1))
    plt.show()
    sns.set()
    plt.scatter(bs[:,0,:], bs[:,1,:])
    return bs

def coh():
    """
    work with pulse coherence integrated
    """
    bs = np.load("../../data/NanoKB/coh.npy")
    bs[:,0] *= (1e15)**2
    
    print("Mean Coherence Time (fs)", np.mean(bs[:,0]))
    print("Std. Dev in Coherence Time (fs)", np.std(bs[:,0]))
    return bs


def esi():
    """
    work with pulse coherence integrated
    """
    bs = np.load("../../data/NanoKB/esi.npy")
    
    return bs



def esp():
    """
    work with pulse coherence integrated
    """
    bs = np.load("../../data/NanoKB/esp.npy")
    
    return bs    



if __name__ == '__main__':
    
    ### load data
    beamSize = np.load("../../data/NanoKB/bsi.npy")[:1,:2,:].reshape([70,2])
    beamArea = np.pi * beamSize[:,0] * beamSize[:,1]
    beamEnergy = np.load("../../data/NanoKB/esi.npy")
    coherence = np.load("../../data/NanoKB/coh.npy")
    coherence[:,0] *= (1e15)    
    beamCentroid = np.load("../../data/NanoKB/cni.npy")
    
    dfI = pd.DataFrame()
    dfI['dx'] = beamSize[:,0]
    dfI['dy'] = beamSize[:,1]
    dfI['area'] = beamArea
    dfI['energy'] = beamEnergy[0,0,:]
    dfI['nPhotons'] = beamEnergy[0,1,:]
    dfI['flux'] = beamEnergy[0,2,:]
    dfI['cohtime'] = coherence[:,0]
    dfI['cohlen'] = coherence[:,1]
    dfI['tdoc'] = coherence[:,2]
    dfI['cx'] = beamCentroid[0,0,:]
    dfI['cy'] = beamCentroid[0,1,:]


    pulseSize = np.load("../../data/NanoKB/bss.npy")
    pulseCentroid = np.load("../../data/NanoKB/cnp.npy") 
    pulseEnergy = np.load("../../data/NanoKB/esp.npy")
    
    dfpulseArea = pd.DataFrame()
    dfpulseFlux = pd.DataFrame()
    dfpulseCentroidX = pd.DataFrame()
    dfpulseCentroidY = pd.DataFrame()
    
    dfpulseFlux['sliceNumber'] = pulseSize[:,0,0]
    dfpulseFlux['time'] = pulseCentroid[:,2,0]
    
    dfpulseArea['sliceNumber'] = pulseSize[:,0,0]
    dfpulseArea['time'] = pulseCentroid[:,2,0]

    dfpulseCentroidX['sliceNumber'] = pulseSize[:,0,0]
    dfpulseCentroidX['time'] = pulseCentroid[:,2,0]
    
    dfpulseCentroidY['sliceNumber'] = pulseSize[:,0,0]
    dfpulseCentroidY['time'] = pulseCentroid[:,2,0]
        
    
    for i in range(pulseSize.shape[-1]):
        dfpulseFlux[str(i)] = pulseEnergy[:,2,i]
        dfpulseArea[str(i)] = np.pi * pulseSize[:,2,i] * pulseSize[:,3,i]
        dfpulseCentroidX[str(i)] = pulseCentroid[:,0,i]
        dfpulseCentroidY[str(i)] = pulseCentroid[:,1,i]
        
    dfpulseFlux.to_csv("/opt/spb_model/data/NanoKB/pulseFlux.ots")
    dfpulseArea.to_csv("/opt/spb_model/data/NanoKB/pulseArea.ots")
    dfpulseCentroidX.to_csv("/opt/spb_model/data/NanoKB/pulseCentroidX.ots")
    dfpulseCentroidY.to_csv("/opt/spb_model/data/NanoKB/pulseCentroidy.ots")
    
    dfI.to_csv("/opt/spb_model/data/NanoKB/integratedData.ots")

    
    
    plt.plot(dfpulseArea['time'], np.mean(pulseEnergy[:,2], axis = 1))
    plt.plot(dfpulseArea['time'], np.std(pulseEnergy[:,2], axis = 1))
    plt.show()
    
    plt.plot(dfpulseArea['time'], np.mean(pulseCentroid[:,0,:], axis = 1))
    plt.plot(dfpulseArea['time'], np.mean(pulseCentroid[:,1,:], axis = 1))