"""
================================
Example ams pyart course 10
================================

This is the 10th example in the AMS Radar conference 2015 pyart course.
Map many radars to a Cartesian grid

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import pyart
from matplotlib import pyplot as plt
import numpy as np
from time import time

radarpath='/data/pyart_examples/'
radarfiles=[radarpath+'KVNX20120501_040957_V06',
            radarpath+'KTLX20120501_040805_V03',
            radarpath+'KINX20120501_041024_V03']


#read files
radars = [pyart.io.read(radarfile) for radarfile in radarfiles]

print('file reading finished')

# grid data
t1 = time()
grids = pyart.map.grid_from_radars(
         radars, grid_shape=(52, 251, 251),
        grid_limits=((0, 17000.0),(-200000, 200000), (-200000, 200000)),
        fields=['reflectivity'], gridding_algo="map_gates_to_grid",
        weighting_function='BARNES')
print(time() - t1)

#display data
display = pyart.graph.GridMapDisplay(grids)
fig = plt.figure(figsize=[15, 7])

# panel sizes
map_panel_axes = [0.05, 0.05, .4, .80]
x_cut_panel_axes = [0.55, 0.10, .4, .25]
y_cut_panel_axes = [0.55, 0.50, .4, .25]

# parameters
level = 5
vmin = -8
vmax = 64
lat = 36.4
lon = -97.0

# panel 1, basemap, radar reflectivity and NARR overlay
ax1 = fig.add_axes(map_panel_axes)
display.plot_basemap(lon_lines = np.arange(-104, -93, 2) )
display.plot_grid('reflectivity', level=level, vmin=vmin, vmax=vmax,
                 cmap = pyart.graph.cm.NWSRef)
display.plot_crosshairs(lon=lon, lat=lat)

# panel 2, longitude slice.
ax2 = fig.add_axes(x_cut_panel_axes)
display.plot_longitude_slice('reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax,
                            cmap = pyart.graph.cm.NWSRef)
ax2.set_ylim([0,17])
ax2.set_xlim([-200,200])
ax2.set_xlabel('Distance from SGP CF (km)')

# panel 3, latitude slice
ax3 = fig.add_axes(y_cut_panel_axes)
ax3.set_ylim([0,17])
ax3.set_xlim([-200,200])
display.plot_latitude_slice('reflectivity', lon=lon, lat=lat, vmin=vmin, vmax=vmax,
                           cmap = pyart.graph.cm.NWSRef)
						   
plt.show()
