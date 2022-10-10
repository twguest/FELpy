import numpy as np


def gaussian(x, sigma, x0 = 0, a = 1):
    """
    one-dimensional gaussian function for fitting
    """

    return a * np.exp(-(x-x0)**2/(2*sigma**2))
