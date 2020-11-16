# -*- coding: utf-8 -*-
import numpy as np

from scipy.constants import c,h

def ekev2wav(ekev):
    
    return (h*c)/(ekev/1000)
