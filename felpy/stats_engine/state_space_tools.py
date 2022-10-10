import numpy as np
from pandas import DataFrame
from felpy.utils.vis_utils import Grids
from matplotlib.pyplot import cm
from skimage.measure import EllipseModel
from matplotlib.patches import Ellipse
from numpy import sin, cos
import pandas as pd

def solution_for_y(x, h,w, cx, cy, theta):
    """
    :param a: width
    :param b: height
    :param cx: horizontal center
    :param cy: vertical center
    :param theta: rotation angle
    """
    x -= cx
    a = ((sin(theta)**2)/w**2 + (cos(theta)**2)/h**2)
    b = 2*sin(theta)*cos(theta)*((1/(w**2))-(1/(h**2)))
    c = sin(theta)**2/w**2 + cos(theta)**2/h**2
    return cy+((-b*x)-np.sqrt(b**2*x**2-4*a*(c*x**2-1)))/(2*a), cy+((-b*x)+np.sqrt(b**2*x**2-4*a*(c*x**2-1)))/(2*a)

def solution_from_df(df, pulse, val, method = None, dim = 'x'):
    """ 
    a wrapper function to return the (n+1)th pulse value from the value and pulse number of the nth pulse
    solution from the values of an ellipse defined by the dataframe df, previous value, and pulse number of that
    value.
    
    if method is none, assumes distribution to be normal
    
    :param df: dataframe mapping the 95% confidence ellipse parameter for each pulse
    :param pulse: the nth pulse [int]
    :param val: the value of the nth pulse
    
    :returns sol: value of the n+1th pulse
    """
    
    lower, upper = solution_for_y(val,
                                  h = df['Height'][pulse]/2,
                                  w = df['Width'][pulse]/2,
                                  cx = df['Centre'][pulse][0],
                                  cy = df['Centre'][pulse][1],
                                  theta = df['Angle'][pulse])
    
    
    while np.isnan(lower) and np.isnan(upper):
        val = df['Centre'][pulse][0]
     
        lower, upper = solution_for_y(val,
                                  h = df['Height'][pulse]/2,
                                  w = df['Width'][pulse]/2,
                                  cx = df['Centre'][pulse][0],
                                  cy = df['Centre'][pulse][1],
                                  theta = df['Angle'][pulse])

    
    std = abs(upper-lower)/(2*df['Sigma'][pulse])
     
    if method is None:
        sol = np.random.normal(loc = (upper+lower)/2, scale = std)

        #sol = np.random.normal(loc = (upper+lower)/2, scale = std)
        

    return sol

def statengine_covariance_ellipse(arr, M = 1, nstd = 2, constraint = None, initial_condition = 0):
    """
    dataframe containing covariance arrays in format of .extract_elliptical_stats output

    :param arr: variable data
    :param M: number of pulse trains to generate of length arr.shape[0]-1
    :param initial_condition: seed value for the n = 0, list of length M or float. 
    
    :returns df_input: definition of covariance ellipses over nxm input array
    :returns df_output: set of M statistical models of arr
    """
    
    df_input = extract_elliptical_stats(arr, nstd = 2)

    outputs = {}

    if constraint is None:

        for m in range(M):
            output = []
            
            if type(initial_condition) == float or int:
                output.append(initial_condition)
            elif type(initial_condition) == list:
                output.append(initial_condition[m])

            for n in range(arr.shape[0]-1):
                output.append(solution_from_df(df_input, pulse = n, val = output[n], method = None, dim = 'x'))

            outputs['{}'.format(m)] = output

    df_output = pd.DataFrame.from_dict(outputs)

    return df_input, df_output




def get_cov_ellipse_params(cov, centre, nstd, **kwargs):
    """
    Return a matplotlib Ellipse patch representing the covariance matrix
    cov centred at centre and scaled by the factor nstd.
    
    via: https://scipython.com/book/chapter-7-matplotlib/examples/bmi-data-with-confidence-ellipses/
    

    """

    # Find and sort eigenvalues and eigenvectors into descending order
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # The anti-clockwise angle to rotate our ellipse by 
    vx, vy = eigvecs[:,0][0], eigvecs[:,0][1]
    theta = np.arctan2(vy, vx)

    # Width and height of ellipse to draw
    width, height = 2 * nstd * np.sqrt(eigvals)
    params = [centre, width, height, theta]

    
    return params

def sequential_pairs(arr, n = 0):
    """
    :param n:
    :return N: length array as a set of pairs for n,n+1
    """
        
    M = arr.shape[1]
    
    return np.array([(arr[n,m], arr[n+1,m]) for m in range(M-1)])



def extract_elliptical_stats(arr, nstd = 2):
    """
    wrapper for convidence_ellipse to extract data to pandas dataframe
    """

    pre_ = []
    cov = confidence_ellipse(arr, nstd = 2)


    for n in range(len(cov)):
        pre_.append([n,*cov[n], nstd])


    df = DataFrame(pre_, columns = ['Pulse','Centre', 'Width', 'Height', 'Angle', 'Sigma'])
    df.set_index('Pulse')

    return df


def confidence_ellipse(arr, nstd = 2):
    """
    function to return a list of covariance (95% condidence interval) ellipses from an NxM dimensional dataset
    
    :param arr: variable arr of size NxM describing the evolution of N -> N+1 
    :returne cov: M covariance parameters lists
    """
    
    N = arr.shape[0]-1
    M = arr.shape[1]

    cov = []

    for n in range(N):

        coordinates = sequential_pairs(arr, n = n)
        N0 = coordinates[:, 0]
        N1 = coordinates[:, 1]

        cov_matrix = np.cov(N0,N1)
        cov.append(get_cov_ellipse_params(cov_matrix, (N0.mean(), N1.mean()), nstd = nstd))

    return cov 


def plot_covariance_ellipse(arr, cov = None, grid = None, N = 0, axis = 0, link_dots = False, add_colorbar = True):
    """
    plot the covariance ellipse of a given array

    :param arr:
    :param cov:
    :param grid:
    :param N:
    """
    
    color = cm.rainbow(np.linspace(0, 1, arr.shape[1]-1))

    if cov == None:
        cov = confidence_ellipse(arr)
    
    if grid is None:

        grid = Grids(scale = 2, global_aspect = 1)
        grid.create_grid(n = 1, m = 1, share_x = False, share_y = False)
        grid.pad(4)
    
    if add_colorbar:
        grid.add_global_colorbar(vmin = 0, vmax = arr.shape[-1], cmap = 'rainbow', clabel = "Train No.", fontsize = 22)
    
    
    
    if type(grid.axes) == list:
        ax = grid.axes[0]
    elif type(grid.axes) == dict:
        ax = grid.axes[list(grid.axes.keys())[axis]]
    else:
        ax = grid.axes
    
    grid.set_fontsize(22)
    
    
    for m in range(arr.shape[-1]-1):
        
        c = color[m]
            
        coordinates = sequential_pairs(arr, n = N)
        N0 = coordinates[:, 0]
        N1 = coordinates[:, 1]

        ax.scatter(N0[m], N1[m], color = c)
    
    if link_dots:
        ax.plot(coordinates[:,0],coordinates[:,1], linewidth = 0.25)


    ex = Ellipse(xy= cov[N][0], width= cov[N][1], height= cov[N][2],angle=np.degrees(cov[N][3]), fill = None)
    ax.add_artist(ex)

    return grid