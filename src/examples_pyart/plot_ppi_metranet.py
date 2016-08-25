"""
======================================
Create a PPI plot from a METRANET file
======================================

An example which creates a PPI plot of a METRANET file using a RadarDisplay object.

"""
print(__doc__)

# Author: fvj
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

filepath='/data/rad4alp/rawdata/15279/PHA15279/'
filename = 'PHA1527913307U.003'

# create the plot using RadarDisplay
radar = pyart.aux_io.read_metranet(filepath+filename)
display = pyart.graph.RadarDisplay(radar)
fig = plt.figure(figsize=[5, 5])
ax = fig.add_subplot(111)
display.plot('reflectivity', 0, vmin=-32, vmax=64.)
display.set_limits(ylim=[-120, 120], xlim=[-120, 120])
display.plot_range_rings([10, 20, 30, 40])
display.plot_cross_hair(5.)
plt.show()
