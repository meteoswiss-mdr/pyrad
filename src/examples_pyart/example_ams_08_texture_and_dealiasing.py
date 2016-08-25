"""
================================
Example ams pyart course 8
================================

This is the 8th example in the AMS Radar conference 2015 pyart course.
Texture and dealiasing

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import pyart
from matplotlib import pyplot as plt
import numpy as np
from scipy import ndimage, signal
import time

radarpath='/data/pyart_examples/'
radarfile='csapr_test_case.nc'

radar = pyart.io.read(radarpath+radarfile)

# plot reflectivity
display = pyart.graph.RadarMapDisplay(radar)
fig = plt.figure(figsize = [10,8])
display.plot_ppi_map('reflectivity', sweep = 2, resolution = 'i',
                    vmin = -10, vmax = 64, mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef)
plt.show()

# plot velocity
nyq = radar.instrument_parameters['nyquist_velocity']['data'][0]

display = pyart.graph.RadarMapDisplay(radar)
fig = plt.figure(figsize = [10,8])
display.plot_ppi_map('velocity', sweep = 2, resolution = 'i',
                    vmin = -nyq, vmax = nyq, mask_outside = False,
                    cmap = pyart.graph.cm.NWSVel)
plt.show()

# compute and plot velocity texture
start_time = time.time()
data = ndimage.filters.generic_filter(radar.fields['velocity']['data'],
                                            pyart.util.interval_std, size = (4,4),
                                           extra_arguments = (-nyq, nyq))
total_time = time.time() - start_time
print(total_time)

filtered_data = ndimage.filters.median_filter(data, size = (4,4))
texture_field = pyart.config.get_metadata('velocity')
texture_field['data'] = filtered_data
radar.add_field('velocity_texture', texture_field, replace_existing = True)

display = pyart.graph.RadarMapDisplay(radar)
fig = plt.figure(figsize = [10,8])
display.plot_ppi_map('velocity_texture', sweep = 2, resolution = 'i',
                    mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    vmin = 0, vmax = 14)
plt.show()


# compute and plot histogram of velocity texture
n, bins = np.histogram(filtered_data.flatten(), bins = 150)

peaks = signal.find_peaks_cwt(n, np.array([10]))
centers = bins[0:-1] + (bins[1] - bins[0])
search_data = n[peaks[0]:peaks[1]]
search_centers = centers[peaks[0]:peaks[1]]
locs = search_data.argsort()
location_of_minima = locs[0]

fig = plt.figure(figsize = [10,6])
plt.plot(centers, n)
zmax = n.max()
plt.xlabel('Radial velocity texture')
plt.ylabel('Number of Gates')

plt.plot([centers[peaks[0]], centers[peaks[0]]], [0, zmax])
plt.plot([centers[peaks[1]], centers[peaks[1]]], [0, zmax])
plt.plot([search_centers[location_of_minima], search_centers[location_of_minima]], [0, zmax])
plt.show()


noise_threshold = search_centers[locs[0]]
print(noise_threshold)

# filter noise from reflectivity
likely_noise = filtered_data > noise_threshold
likely_signal = np.logical_not(likely_noise)

z_masked = np.ma.masked_where(likely_noise, radar.fields['reflectivity']['data'])
radar.add_field_like('reflectivity', 
                     'reflectivity_masked', 
                     z_masked, replace_existing = True)

					 
# plot filtered and unfiltered reflectivity
display = pyart.graph.RadarMapDisplay(radar)
fig = plt.figure(figsize = [15,7])
plt.subplot(1,2,1)
display.plot_ppi_map('reflectivity_masked', sweep = 2, resolution = 'i',
                    mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    vmin = -10, vmax = 64)
plt.subplot(1,2,2)
display.plot_ppi_map('reflectivity', sweep = 2, resolution = 'i',
                    mask_outside = False,
                    cmap = pyart.graph.cm.NWSRef,
                    vmin = -10, vmax = 64)
plt.show()

# filter noise from velocity
gatefilter = pyart.correct.GateFilter(radar)
gatefilter.exclude_masked('reflectivity_masked')

corr_vel = pyart.correct.dealias_region_based(
    radar, vel_field='velocity', keep_original=False, 
    gatefilter = gatefilter, nyquist_vel=nyq, centered = True)
radar.add_field('corrected_velocity', corr_vel, replace_existing = True)

# plot filtered and unfiltered velocity
display = pyart.graph.RadarMapDisplay(radar)
fig = plt.figure(figsize = [15,7])
plt.subplot(1,2,1)
display.plot_ppi_map('velocity', sweep = 2, resolution = 'i',
                    mask_outside = False,
                    cmap = pyart.graph.cm.NWSVel,
                    vmin = -nyq, vmax = nyq)
plt.subplot(1,2,2)
display.plot_ppi_map('corrected_velocity', sweep = 2, resolution = 'i',
                    mask_outside = False,
                    cmap = pyart.graph.cm.NWSVel,
                    vmin = -1.5*nyq, vmax = 1.5*nyq)
plt.show()







