import numpy as np

def complex_gaussian_2d(nx, ny, sigma = 1, mu = 0):
    """
    generate a 2d gaussian-like beam in array form

    :param nx: number of columns
    :param ny: number of rows
    """
    x, y = np.meshgrid(np.linspace(-1,1,nx), np.linspace(-1,1,ny))
    d = np.sqrt(x*x+y*y)
    g = np.exp(-(1j* (d-mu)**2 / ( 2.0 * sigma**2 ) ) )
    return g
