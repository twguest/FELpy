from scipy import optimize
from felpy.utils.maths.fit_funcs import gaussian


def fit_gaussian(x, data, VERBOSE = True, **kwargs):
    """
    optimisation for fit of autocovariance function

    see scipy.optimize.curve_fit for keyword arguments

    :param x:
    :param data: data to be fit


    """

    params, params_covariance = optimize.curve_fit(gaussian, x, data, **kwargs)
    
    if VERBOSE:
        print("sigma: {}".format(params[0]))
        print("x0: {}".format(params[1]))
        print("a: {}".format(params[2]))
        

    
    return params, params_covariance



if __name__ == '__main__':

    ### test gaussian fit

    import numpy as np
    from matplotlib import pyplot as plt
    x = np.linspace(-100,100,100)
    
    f = gaussian(x, sigma = 15, x0 = 10)*np.arange(100)
    
    plt.plot(x,f)
    params, params_covariance = fit_gaussian(x, f)
    plt.scatter(x, gaussian(x, *params))
    
    