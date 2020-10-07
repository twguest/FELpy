#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script houses the analysis pipepline for a set of pulses in a directory.
It primarily concerns data extraction and cleaning; the structure of the main()
function of the script is as follows:
    
    * - yields plot
    ** - yields data of merit
    
    - Data handling methods:
        -Saving methods (ie building memmaps and assigning a location to each
                         pulse) (DONE)

    - Pulse handling methods:
        - Extraction of intensity * (DONE)
        - Construction of multi-pulse profiles etc. ** (DONE)
        
    ## PULSE STATISTICS
    
    - Basic statistics [individual pulse]:
        - Pulse Fluence * (DONE)
        - System Efficiency * 
        - Pulse Centroid Intensity and Location **Centroid Density Plot (3D - X,Y,INTENSITY) (DONE)
            see(https://python-graph-gallery.com/371-surface-plot/)
        - Pulse Energy and Number of Photons. *(DONE)
    
    - Coherence analysis [individual pulse]
        - Coherence Length Measurement * (DONE)
        - Coherence Time Measurement * (DONE)
        - Transverse Degree of Coherence * (DONE)
    
    - Beam size analysis [Individual Pulse]:
        - Enclosed energy of a single pulse *Comparison w/ Legacy Values (DONE)

    ## INTRA-PULSE STATISTICS:
            
    - Basic statistics [Intra-Pulse]:
        - Slice Centroid Intensity and Location vs. Time **(4D Scatter - X,Y, Slice Intensity, Time) (DONE)

    - Beam size analysis [Intra-Pulse]:
        - Intra-Pulse Enclosed Energy **PLOT BEAM SIZE v. Pulse Time (DONE)

    ## PULSE TRAIN STATISTICS:
    
    - Coherence analysis [multi pulse]
        - Coherence Length Measurement ** (DONE)
        - Coherence Time Measurement ** (DONE)
        - Transverse Degree of Coherence** (DONE)
    
    - Beam size analysis [Multi-Pulse]:
        - Enclosed energy of a integrated pulse** (DONE)
        
    USAGE:
        >>> python -c "from pipeline import launch; launch(); launch(multi = True)"
        
    
@author: twguest
@version: 0.0.0
"""


import sys
import os

sys.path.append("../../")

import numpy as np

from multiprocessing.pool import ThreadPool as Pool

from model.tools import mkdir_p, memoryMap, readMap

from model.tools import constructPulse ## for testing

from model.analysis.coherence import run as coherence
from model.analysis.enclosedEnergy import run as beamSize
from model.analysis.energyStatistics import getPulseEnergy

from wpg.wavefront import Wavefront

from wpg.wpg_uti_wf import getCentroid

from utils.job_utils import JobScheduler

def setup(VERBOSE = False):
    """
    These are the relevant parameters for loading and saving from the analysis
    pipeline. a dictionary is constructed containing all the relevant
    parameters and directories. Due to the large size of some of the data
    generated, it is necessary to send each of the pulses to a worker. 
    Consequently, we must save data to a location on disk, rather than in mem
    using np.memmap. 
    
    The following definitions are not limiting. However it is recommended that
    we only change the memmap directory name, and not the parameter key.
    
    In future, adding any analysis requires new memmaps to be generated for 
    the new data, and as a consequence a new dictionary entry will need to be 
    added. The workflow for this should be as such:
        
        - create key entry for memMap 
        - create memMap
        - read memMap into required function
        - add memMap into flush()
        
    ### TODO: in the future a helper function might be able to ease this process by 
    identifying all memmaps in the params file or other.
    
    The utility of each dictionary entry is included as a comment next to its
    definition
    
    :param VERBOSE: [bool.] print dictionary definition to console
    
    :return params: parameter dictionary
    """
    
    
    params = {}
    
    params['nProc'] = 8 ## number of cpus for MPI 
    
    params['indir'] = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/out/" ## input directory
    params['global'] = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/" ## global data directory
    
    ### train 
    params ['traindir'] = "/gpfs/exfel/data/group/spb-sfx/user/guestt/h5/NanoKB-Pulse/data/pulseTrain/" ## directory to save pulse train
    params['trainwfr'] = "pulseTrain.h5" ## pulse train name 
    
    params['train'] = 'train' ## pulse train e-field memmap
    params['tsi'] = "tsi" ## pulse train integrated beam size memmap
    params['tss'] = "tss" ## pulse train slice-to-slice beam size memmap
    params['tcoh'] = 'tcoh' ## pulse train coherence memmap 


    params['tesi'] = "tesi" ## pulse train integrated energy statistics memmap
    params['tesp'] = "tesp" ## pulse train slice-to-slice energy stats memmap
    params['tcni'] = 'tcni' ## pulse train integrated centroid memmap
    params['tcnp'] = 'tcnp' ## pulse train slice-to-slice centroid memmap
    
    
    ### puse-to-pulse
    params['esi'] = "esi" ## energy statistics integrated memmap 
    params['esp'] = "esp" ## energy statistics pulsed memmap 
    params['cni'] = 'cni' ## centroid integrated memmap 
    params['cnp'] = "cnp" ## centroid pulse memmap 
    
    params['bsi'] = "bsi" ## integrated beam size memmap 
    params['bss'] = "bss" ## pulsed beam size memmap 
    params['coh'] = 'coh' ## single-coherence memmap

    

    mkdir_p(params['global'])
    mkdir_p(params['traindir'])
    
    ### hardcoded pulse dimensions
    nx = 1024
    ny = 1024
    nSlice = 2250 ### number of slices per pulse
    
    
    
    ### construct memmaps
    train = memoryMap(params['global'], params['train'],
          shape = (nx,ny,nSlice,2))  

    tsi = memoryMap(params['global'], params['tsi'],
              shape = 2)
    
    tss = memoryMap(params['global'], params['tss'], 
          shape=(nSlice, 4, 1))
    
    
    tcoh = memoryMap(params['global'], params['tcoh'], 
          shape=(1, 3))


    tesi = memoryMap(params['global'], params['tesi'],
              shape = (1,3,1)) 
    
 
    tesp = memoryMap(params['global'], params['tesp'], 
              shape=(nSlice,3,1))
    
    tcni = memoryMap(params['global'], params['tcni'], 
              shape=(1,2,1))
    
    
    tcnp = memoryMap(params['global'], params['tcnp'],
          shape = (nSlice,4,1)) 

    
    esi = memoryMap(params['global'], params['esi'],
              shape = (1,3,len(os.listdir(params['indir'])))) 
    
 
    esp = memoryMap(params['global'], params['esp'], 
              shape=(nSlice,3,len(os.listdir(params['indir']))))
    
    cni = memoryMap(params['global'], params['cni'], 
              shape=(1,2,len(os.listdir(params['indir']))))
    
    
    cnp = memoryMap(params['global'], params['cnp'],
          shape = (nSlice,4,len(os.listdir(params['indir']))))  


    bsi = memoryMap(params['global'], params['bsi'],
              shape = (len(os.listdir(params['indir'])),2)) 
    
 
    bss = memoryMap(params['global'], params['bss'], 
              shape=(nSlice, 4, len(os.listdir(params['indir']))))
    
    coh = memoryMap(params['global'], params['coh'], 
              shape=(len(os.listdir(params['indir'])), 3))
    
    
    ### define memorymap shapes
    params["train_shape"] = train.shape
    
    params["tsi_shape"] = tsi.shape
    params["tss_shape"] = tss.shape
    params["tcoh_shape"] = tcoh.shape
    

    params["tesi_shape"] = tesi.shape
    params["tesp_shape"] = tesp.shape
    params["tcni_shape"] = tcni.shape
    params["tcnp_shape"] = tcnp.shape
    

    params["esi_shape"] = esi.shape
    params["esp_shape"] = esp.shape
    params["cni_shape"] = cni.shape
    params["cnp_shape"] = cnp.shape
    

    params["bsi_shape"] = bsi.shape
    params["bss_shape"] = bss.shape
    params["coh_shape"] = coh.shape
    

    if VERBOSE:
        
        print("Analysis Pipeline Parameters")
        
        for key in params:
            print(key, params[key])
    
        print("\n")
    return params

def generateTestPulses(savedir, nx = 1024, ny = 1024, N = 5):
    """
    generate a set of test pulses
    
    :param savedir: directory for test pulses to be saved
    :param N: number of pulses to generate
    """
    print("Constructing {} Test Pulses".format(N))
    
    for n in range(N):
        
        wfr = constructPulse(nx, ny, nz = 6, tau = 1e-12)
        
        wfr.data.arrEhor*= np.random.uniform(0.75, 1, size = wfr.data.arrEhor.shape)
        
        wfr.store_hdf5(savedir + "testWavefront_{}.h5".format(n))
        print("Storing Wavefront @" + savedir + "testWavefront_{}.h5".format(n))    
    
def BeamSize(wfr, mode, memMap, ID, VERBOSE = False):
        
    if mode == 'integrated':
        memMap[ID, :] = beamSize(wfr, mode = mode, VERBOSE = VERBOSE)
    if mode == 'pulse':
        memMap[:,:,ID] = beamSize(wfr, mode = mode, VERBOSE = VERBOSE, nSlc = 250) 
       
def Coherence(wfr, memMap, ID, VERBOSE):
    
    memMap[ID,:] = coherence(wfr)

def superimpose(fname, params):    
    """ 
    function to add wavefronts
    """
    
    
    ### TODO: FIX HARDCODING OF FILENAMES ETC.
    memMap = readMap(params['global'] + params['train'],
                     shape = params['train_shape'])
    wfr = Wavefront()
    wfr.load_hdf5(params['indir'] + fname)
    
    memMap += wfr.data.arrEhor 

def getTrain(indir, params, mpi = True, nProc = 1):
    """ 
    wrapper function linear addition of a set of wavefronts parsed through from some input
    directory
    """
    ### TODO: FIX HARDCODING OF FILENAMES ETC.
    
    from functools import partial
    
    if mpi:
        
        p = partial(superimpose, params = params)
        pool = Pool(processes = nProc)
        pool.map(p, iterable =  os.listdir(indir))
    
    wfr = Wavefront()
    wfr.load_hdf5(params['indir'] + os.listdir(params['indir'])[
        np.random.randint(0,len(os.listdir(params['indir'])))])

    wfr.data.arrEhor = readMap(params['global'] + params['train'],
                     shape = params['train_shape'])
    
    wfr.store_hdf5(params['traindir'] + params['trainwfr'])

def trainAnalysis():
    
    params = setup(VERBOSE = False)
    
    ID = 0

    ## get multipulse profiles
    getTrain(params['indir'], params = params, mpi = True, nProc = params['nProc'])
    print("Pulse Train Succesfully Added")
    ## load relevant train memmap
    tsi = readMap(params['global'] + params['tsi'], shape = params['tsi_shape'])
    tss = readMap(params['global'] + params['tss'], shape = params['tss_shape'])
    tcoh = readMap(params['global'] + params['tcoh'], shape = params['tcoh_shape'])

    tesi = readMap(params['global'] + params['tesi'], shape = params['tesi_shape'])
    tesp = readMap(params['global'] + params['tesp'], shape = params['tesp_shape'])
    tcni = readMap(params['global'] + params['tcni'], shape = params['tcni_shape'])
    tcnp = readMap(params['global'] + params['tcnp'], shape = params['tcnp_shape'])
    
    ## load multipulse for analysis prior to launch
    tfr = Wavefront()
    tfr.load_hdf5(params['traindir'] + params['trainwfr'])
    
                
    pulseEnergy(tfr, tesi, ID, mode = 'integrated')
    pulseEnergy(tfr, tesp, ID, mode = 'pulse')
    print("Train Energy Calculated")
    centroid(tfr, tcni, ID, mode = 'integrated')
    centroid(tfr, tcnp, ID, mode = 'pulse')
    print("Train Centrod Calculated")
    
    BeamSize(tfr, mode = 'integrated', memMap = tsi, ID = 0)
    BeamSize(tfr, mode = 'pulse', memMap = tss, ID = 0)
    print("Train Beam Size Calculated")
    Coherence(tfr, tcoh, ID = 0, VERBOSE = True)
    print("Train Coherence Calculated")

    del tsi, tss, tcoh, tesi, tesp, tcni, tcnp
    
def pulseEnergy(wfr, memMap, ID, mode = "integrated"):
    
    memMap[:,:,ID] = getPulseEnergy(wfr, mode = mode)
    
def centroid(wfr, memMap, ID, mode = "integrated"):

    memMap[:,:,ID] = getCentroid(wfr, mode = mode, idx = False)

def flush(params):
    """
    save memMaps to final npy file
    """    
    
  
    bsi = readMap(params['global'] + params['bsi'], shape = params['bss_shape'])
    np.save(params['global'] + params['bsi'], bsi)
    
    bss = readMap(params['global'] + params['bss'], shape = params['bss_shape'])
    np.save(params['global'] + params['bss'], bss)
    
    coh = readMap(params['global'] + params['coh'], shape = params['coh_shape'])
    np.save(params['global'] + params['coh'], coh)
    
    esi = readMap(params['global'] + params['esi'], shape = params['esi_shape'])
    np.save(params['global'] + params['esi'], esi)
    
    esp = readMap(params['global'] + params['esp'], shape = params['esp_shape'])
    np.save(params['global'] + params['esp'], esp)
    
    cni = readMap(params['global'] + params['cni'], shape = params['cni_shape'])
    np.save(params['global'] + params['cni'],cni)
    
    cnp = readMap(params['global'] + params['cnp'], shape = params['cnp_shape'])
    np.save(params['global'] + params['cnp'], cnp)
    
    tsi = readMap(params['global'] + params['tsi'], shape = params['tsi_shape'])
    np.save(params['global'] + params['tsi'], tsi)
    
    tss = readMap(params['global'] + params['tss'], shape = params['tss_shape'])
    np.save(params['global'] + params['tss'], tss)
    
    tcoh = readMap(params['global'] + params['tcoh'], shape = params['tcoh_shape'])
    np.save(params['global'] + params['tcoh'], tcoh)
    
    tesi = readMap(params['global'] + params['tesi'], shape = params['tesi_shape'])
    np.save(params['global'] + params['tesi'], tesi)
    
    tesp = readMap(params['global'] + params['tesp'], shape = params['tesp_shape'])
    np.save(params['global'] + params['tesi'], tesp)
    
    tcni = readMap(params['global'] + params['tcni'], shape = params['tcni_shape'])
    np.save(params['global'] + params['tcni'], tcni)
    
    tcnp = readMap(params['global'] + params['tcnp'], shape = params['tcnp_shape'])
    np.save(params['global'] + params['tcnp'], tcnp)
    
    print("Items Saved, Analysis Complete")
 
def launch(multi = False):
    """ 
    job scheduler function
    """
    
    params = setup()
    

    
    if multi:
        
        js = JobScheduler(pycmd = " -c 'import pipeline; pipeline.trainAnalysis()' ",
                          jobName = "pulseTrainAnalysis", logDir = "../../logs/",
                          jobType = 'single', nodes = 4)
        
        js.run(test = False)
        
    else:
        js = JobScheduler(pycmd = os.getcwd() + "/pipeline.py", 
                          jobName = "FEL_PulseAnalysis",
                          logDir = "../../logs/",
                          jobType = 'array',
                          nodes = 1,
                          jobArray = range(len(os.listdir(params['indir']))))
        
        js.run(test = False)
   

if __name__ == '__main__':
    
    ID = int(sys.argv[1])
    print("Job ID ", ID)
    
    params = setup(VERBOSE = True)


    
    bsi = readMap(params['global'] + params['bsi'], shape = params['bsi_shape'])
    bss = readMap(params['global'] + params['bss'], shape = params['bss_shape'])
    coh = readMap(params['global'] + params['coh'], shape = params['coh_shape'])

    esi = readMap(params['global'] + params['esi'], shape = params['esi_shape'])
    esp = readMap(params['global'] + params['esp'], shape = params['esp_shape'])
    cni = readMap(params['global'] + params['cni'], shape = params['cni_shape'])
    cnp = readMap(params['global'] + params['cnp'], shape = params['cnp_shape'])
    
    fname = os.listdir(params['indir'])[ID]
    print("loading wavefront: {}".format(fname))
    wfr = Wavefront()
    wfr.load_hdf5(params['indir']+fname)
    
    
    pulseEnergy(wfr, esi, ID, mode = 'integrated')
    pulseEnergy(wfr, esp, ID, mode = 'pulse')
    print("Pulse Energy Calculated")
    centroid(wfr, cni, ID, mode = 'integrated')
    centroid(wfr, cnp, ID, mode = 'pulse')
    print("Pulse Centroid Calculated")
    ### integrated beam profile
    print("Calculating Integrated Beam Profiles")
    BeamSize(wfr, mode = "integrated", memMap = bsi, ID = ID)
    print("Integrated Beam Profiles Calculated")
    ### pulsed beam profile
    print("calculating pulsed beam profiles")
    BeamSize(wfr, mode = "pulse", memMap = bss, ID = ID, VERBOSE = True)
    print("Beam Size Calculated")
    ### single-pulse coherence analysis
    Coherence(wfr, coh, ID, VERBOSE = True)
    print("Beam Coherence Calculated")
    
    del bsi, bss, coh, esi, esp, cni, cnp ## save to mem?!?
    
    print("Saving Results...")
    flush()
    print("Results Saved")
    
    