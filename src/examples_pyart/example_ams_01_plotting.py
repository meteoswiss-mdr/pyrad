"""
================================
Example ams pyart course 1
================================

This is the first example in the AMS Radar conference 2015 pyart course. Simple figure plotting 

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import matplotlib.pyplot as plt
import pyart

radarpath='/data/pyart_examples/'
radarfile='sgpwsacrcwrhiC1.a1.20120820.204016.nc'

plt.rcParams['figure.figsize'] = [12.0, 9.0]

radar = pyart.io.read(radarpath+radarfile)
radar.info('compact')
# try out 'standard' or 'full' also

display = pyart.graph.RadarDisplay(radar)

fig = plt.figure()
display.plot('reflectivity', sweep=1, vmin=-32, vmax=20)
plt.show()

fig = plt.figure()
display.plot('reflectivity', sweep=1, vmin=-32, vmax=20, cmap='pyart_NWSRef')
plt.show()

fig = plt.figure()
display.plot('reflectivity', sweep=1, vmin=-32, vmax=20, cmap='pyart_NWSRef')
display.set_limits(ylim=(7, 12), xlim=(-18, -10))
plt.show()

fig = plt.figure()
display.plot('snr', sweep=1, vmin=-20, vmax=30)
plt.show()

fig = plt.figure()
display.plot('mean_doppler_velocity', sweep=1, vmin=-4, vmax=4, cmap='pyart_NWSVel')
plt.show()

fig = plt.figure()
display.plot('spectral_width', sweep=1, vmin=0, vmax=0.4)
plt.show()

fig = plt.figure()
display.plot('linear_depolarization_ratio', sweep=1, vmin=-30, vmax=0, cmap='pyart_Carbone17_r')
plt.show()







