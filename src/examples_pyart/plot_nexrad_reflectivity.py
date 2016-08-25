"""
====================================
Create a plot of NEXRAD reflectivity
====================================

An example which creates a plot containing the first collected scan from a
NEXRAD file.

"""
print (__doc__)

# Author: Jonathan J. Helmus (jhelmus@anl.gov)
# License: BSD 3 clause
# 2016.10.05 tested and works

import matplotlib.pyplot as plt
import pyart

# open the file, create the displays and figure
filepath='/data/pyart_examples/'
filename='KDDC_20150513_1940'

# original file name
#filename = 'Level2_KATX_20130717_1950.ar2v'

radar = pyart.io.read_nexrad_archive(filepath+filename)
display = pyart.graph.RadarDisplay(radar)
fig = plt.figure(figsize=(6, 5))

# plot super resolution reflectivity
ax = fig.add_subplot(111)
display.plot('reflectivity', 0, title='NEXRAD Reflectivity',
             vmin=-32, vmax=64, colorbar_label='', ax=ax)
display.plot_range_ring(radar.range['data'][-1]/1000., ax=ax)
display.set_limits(xlim=(-500, 500), ylim=(-500, 500), ax=ax)
plt.show()
