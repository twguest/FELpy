# -*- coding: utf-8 -*-

 import numpy as np
spb = Instrument()

surfaces = []
x = []
y = []
 
for mirror in ["HOM1", "HOM2", "NHE", "NVE"]:
    surface, x1, y1 = spb.get_mirror_profile(mirror)
    surfaces.append(surface)
    x.append(x1)
    y.append(y1)
    
    
from felpy.utils.vis_utils import simple_line_plot


### plot HOM1

ax1 = simple_line_plot(y = surfaces[0][0,1:], x = x[0][0:], label = "HOM1", 
                       xlabel = "x (mm)", ylabel = "Height Error (nm)",
                       color = 'red', return_axes = True)  

ax1 = simple_line_plot(y = surfaces[1][0,0:], x = x[1][1:], label = "HOM2", 
                       xlabel = "x (mm)", ylabel = "Height Error (nm)",
                       color = 'blue', return_axes = True, parse_axes = ax1)  

ax1.legend()


ax2 = simple_line_plot(y = surfaces[2][0,0:]*1e-9, x = y[2][0:], label = "NHE", 
                       xlabel = "x (mm)", ylabel = "Height Error (nm)",
                       color = 'red', return_axes = True)  

ax2 = simple_line_plot(y = surfaces[3][0,0:]*1e-9-np.mean(surfaces[3][0,0:]*1e-9), x = x[3][0:], label = "NVE", 
                       xlabel = "x/y (mm)", ylabel = "Height Error (nm)",
                       color = 'blue', return_axes = True, parse_axes = ax2)  

ax2.legend()

