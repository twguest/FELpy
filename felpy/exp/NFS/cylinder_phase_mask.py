# -*- coding: utf-8 -*-

from felpy.model.materials.phase_mask import phase_mask
import numpy as np
from felpy.utils.vis_utils import basic_plot
from matplotlib import pyplot as plt
from felpy.utils.np_utils import get_mesh
from felpy.utils.np_utils import gaussian_2d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib import cm
import seaborn as sns


def hollow_cylinder_thickness(x,a,b):
    
    #assert(a>b)
    
    return np.sqrt(a**2-x**2)-np.sqrt(b**2-x**2)


def phase_gradient(x, a, b, k, delta):
    return 2*x*delta*((1/np.sqrt(a**2-x**2)) )#-(1/np.sqrt(b**2-x**2)))



def phase(x, a, b, k, delta):
    return k*delta*(2*np.sqrt(a**2-x**2)-2*np.sqrt(b**2-x**2))

if __name__ == '__main__':  
    
    cmap = [cm.viridis(a) for a in np.linspace(0,1,5)]
    exp_setup = {}
    exp_setup['z1'] = np.linspace(0.5, 2.5, 5)
    exp_setup['z2'] = np.linspace(2, 0.5, 5)
    
    fig, ax1 = plt.subplots()
    
    for itr in range(len(exp_setup['z1'])):
        
        a = 1e-03
        b = 0.5e-03
        
        M =(exp_setup['z1'][itr]+exp_setup['z2'][itr])/exp_setup['z1'][itr]        
        from felpy.utils.vis_utils import simple_line_plot  
        
        dx = 500
        
        X = np.linspace(-a*2 , a*2 , dx)
        
        
        ans = hollow_cylinder_thickness(X,a,b)
        #plt.plot(X*1e6,ans)
        #plt.show()
        
    
        ax1 = simple_line_plot(X*1e3, phase_gradient(X, a, b, np.pi*2*(1/5e-11), 1.5e-06)*1e6,
                         xlabel = "x(mm)", ylabel = "Phase Gradient ($\mu$rad)", parse_axes = ax1, return_axes = True,
                         label = "{:.0f} ({:.0f},{:.0f})".format(M, exp_setup['z1'][itr], exp_setup['z2'][itr] ),
                         color = cmap[itr])
        plt.legend()
        
        #plt.show()
    #ax1.set_xlim([-7.5, 7.5])
    #ax1.set_ylim([-20,20])
    plt.show()
        
    wrapped_phases = (np.pi * phase(X, a, b, np.pi*2*(1/5e-10), 2e-04)*1e6) % (2 * np.pi) - np.pi
    #plt.plot(X*1e6, wrapped_phases)
    ii = np.exp(X/50**2)
    ans *= dx
    ans = np.ones([dx, dx])
    
# =============================================================================
#     #plt.plot(X*1e6, wrapped_phases)
#     #plt.show()
#     
# =============================================================================
    ans *= wrapped_phases
    basic_plot(ans,
               mesh = get_mesh(ans, a/dx, a/dx),
               cmap = 'hsv')
        
# =============================================================================
#     # Cyclic colormap
#     my_cmap = ListedColormap(sns.color_palette("husl", 256))
#     # Alpha-blending grayscale
#     alpha_cmap = cm.gray(np.linspace(0, 1, 256))
#     alpha_cmap[:, -1] = 1 - np.kaiser(256, 15)
#     alpha_cmap = ListedColormap(alpha_cmap)
#     
#     def f(x, y): # arbitrary complex function with 1 pole and 3 zeros
#         z = x + 1j*y
#         return z**2 - 8/z
#     
#     x = y = np.linspace(-3, 3, 1000)
#     X, Y = np.meshgrid(x, y)
#     C = f(X, Y)
#     
#     gauss = gaussian_2d(dx, dx, dx/5)
#      
#     fig, ax = plt.subplots()
#     im = ax.imshow(np.angle(ans*1j*gauss), cmap=my_cmap, vmin=-1, vmax=1, extent=(-3, 3, -3, 3), origin='lower')
#     ax.imshow(gauss.real, cmap=alpha_cmap, norm=LogNorm(0.5, 1), extent=(-3, 3, -3, 3), origin='lower')
#     cbar = plt.colorbar(im)
#     cbar.set_label('phase ($\pi$)')
#     ax.set_xlabel('real')
# =============================================================================
# =============================================================================
#     ax.set_ylabel('imaginary')
#     fig.tight_layout()
#     plt.show()
# =============================================================================
        
        