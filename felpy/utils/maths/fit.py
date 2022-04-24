from scipy import optimize
from felpy.utils.maths.fit_funcs import gaussian
from felpy.utils.maths.constants import sigma_to_fwhm


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
    ### usage
    
    import numpy as np
    from matplotlib import pyplot as plt
    
    
# =============================================================================
#     ### test gaussian fit
#     x = np.linspace(-100,100,100)
#     
#     f = gaussian(x, sigma = 15, x0 = 10)*np.arange(100)
#     
#     plt.plot(x,f)
#     params, params_covariance = fit_gaussian(x, f)
#     plt.scatter(x, gaussian(x, *params))
# =============================================================================
    
    

    ### test gaussian fit w/ wavefront
    
    
    from felpy.model.src.coherent import construct_SA1_pulse
    
    # setup
    wfr = construct_SA1_pulse(512, 512, 2, 5.0, 0.25)
    
    
    x = wfr.x_axis
    y = wfr.y_axis

    ix = wfr.x_profile
    iy = wfr.y_profile
    
    imax = wfr.peak_intensity
    
    ### fit x
    initial_guess = [wfr.get_fwhm()[0], wfr.com[0], imax]
    p, co = fit_gaussian(x, ix, p0 = initial_guess)
    
    fwhm_x = p[0]*sigma_to_fwhm
        
    ### fit y
    initial_guess = [wfr.get_fwhm()[1], wfr.com[1], imax]
    p, co = fit_gaussian(y, iy, p0 = initial_guess)
    
    fwhm_y = p[0]*sigma_to_fwhm
       
    