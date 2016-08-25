"""
================================
Example ams pyart course 7
================================

This is the 7th example in the AMS Radar conference 2015 pyart course.
Simple processing and adding a field

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import pyart
from matplotlib import pyplot as plt
import numpy as np
print(pyart.__version__)

radarpath='/data/pyart_examples/'
radarfile='3402338_KAMX_20140417_1056'

radar = pyart.io.read(radarpath+radarfile)

display = pyart.graph.RadarMapDisplay(radar)
f = plt.figure(figsize = [17,4])
plt.subplot(1, 3, 1) 
display.plot_ppi_map('differential_reflectivity', max_lat = 26.5, min_lat =25.4, min_lon = -81., max_lon = -79.5,
                     vmin = -7, vmax = 7, lat_lines = np.arange(20,28,.2), lon_lines = np.arange(-82, -79, .5),
                     resolution = 'l')
plt.subplot(1, 3, 2) 
display.plot_ppi_map('reflectivity', max_lat = 26.5, min_lat =25.4, min_lon = -81., max_lon = -79.5,
                     vmin = -8, vmax = 64, lat_lines = np.arange(20,28,.2), lon_lines = np.arange(-82, -79, .5),
                     resolution = 'l')
plt.subplot(1, 3, 3) 
display.plot_ppi_map('velocity', sweep = 1, max_lat = 26.5, min_lat =25.4, min_lon = -81., max_lon = -79.5,
                     vmin = -15, vmax = 15, lat_lines = np.arange(20,28,.2), lon_lines = np.arange(-82, -79, .5),
                     resolution = 'l')
plt.show()


# smooth zdr and add the smoothed field 
zdr = radar.fields['differential_reflectivity']['data']
smooth_zdr = np.zeros_like(zdr)
for i in range(smooth_zdr.shape[0]):
    smooth_zdr[i,:] = \
            pyart.correct.phase_proc.smooth_and_trim(zdr[i,:], 8)

radar.add_field_like('differential_reflectivity', 
                     'differential_reflectivity_smooth', 
                     smooth_zdr, replace_existing = True)

display = pyart.graph.RadarMapDisplay(radar)
f = plt.figure(figsize = [17,4])
plt.subplot(1, 3, 1) 
display.plot_ppi_map('differential_reflectivity_smooth', max_lat = 26.5, min_lat =25.4, min_lon = -81., max_lon = -79.5,
                     vmin = -7, vmax = 7, lat_lines = np.arange(20,28,.2), lon_lines = np.arange(-82, -79, .5),
                     resolution = 'i')
plt.subplot(1, 3, 2) 
display.plot_ppi_map('reflectivity', max_lat = 26.5, min_lat =25.4, min_lon = -81., max_lon = -79.5,
                     vmin = -8, vmax = 64, lat_lines = np.arange(20,28,.2), lon_lines = np.arange(-82, -79, .5),
                     resolution = 'i')
plt.subplot(1, 3, 3) 
display.plot_ppi_map('velocity', sweep = 1, max_lat = 26.5, min_lat =25.4, min_lon = -81., max_lon = -79.5,
                     vmin = -15, vmax = 15, lat_lines = np.arange(20,28,.2), lon_lines = np.arange(-82, -79, .5),
                     resolution = 'i')
					 
plt.show()


# detailed view of the smoothed field
display = pyart.graph.RadarMapDisplay(radar)
f = plt.figure(figsize = [17,4])
plt.subplot(1, 3, 1) 
display.plot_ppi_map('differential_reflectivity_smooth', max_lat = 26.4, min_lat =26, min_lon = -80.75, max_lon = -80.25,
                     vmin = -7, vmax = 7, lat_lines = np.arange(20,28,.1), lon_lines = np.arange(-82, -79, .2),
                     resolution = 'l')
plt.subplot(1, 3, 2) 
display.plot_ppi_map('reflectivity', max_lat = 26.4, min_lat =26, min_lon = -80.75, max_lon = -80.25,
                     vmin = -8, vmax = 64, lat_lines = np.arange(20,28,.1), lon_lines = np.arange(-82, -79, .2),
                     resolution = 'l')
plt.subplot(1, 3, 3) 
display.plot_ppi_map('velocity', sweep = 1, max_lat = 26.4, min_lat =26, min_lon = -80.75, max_lon = -80.25,
                     vmin = -15, vmax = 15, lat_lines = np.arange(20,28,.1), lon_lines = np.arange(-82, -79, .2),
                     resolution = 'l')
					 
plt.show()


# plot phidp instead of radial velocity
display = pyart.graph.RadarMapDisplay(radar)
f = plt.figure(figsize = [17,4])
plt.subplot(1, 3, 1) 
display.plot_ppi_map('differential_reflectivity_smooth', max_lat = 26.4, min_lat =26, min_lon = -80.75, max_lon = -80.25,
                     vmin = -7, vmax = 7, lat_lines = np.arange(20,28,.1), lon_lines = np.arange(-82, -79, .2),
                     resolution = 'l')
plt.subplot(1, 3, 2) 
display.plot_ppi_map('reflectivity', max_lat = 26.4, min_lat =26, min_lon = -80.75, max_lon = -80.25,
                     vmin = -8, vmax = 64, lat_lines = np.arange(20,28,.1), lon_lines = np.arange(-82, -79, .2),
                     resolution = 'l')
plt.subplot(1, 3, 3) 
display.plot_ppi_map('differential_phase', sweep = 0, max_lat = 26.4, min_lat =26, min_lon = -80.75, max_lon = -80.25,
                     vmin = 0, vmax = 360, lat_lines = np.arange(20,28,.1), lon_lines = np.arange(-82, -79, .2),
                     resolution = 'l')

plt.show()

# write modified radar instance in a file
pyart.io.write_cfradial(radarpath+'example_ams_7.nc', radar)

