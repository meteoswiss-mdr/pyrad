"""
=================================
Create a RHI plot from a MDV file
=================================

An example which creates a RHI plot of a MDV file using a RadarDisplay object.

"""
print(__doc__)

# Author: Jonathan J. Helmus (jhelmus@anl.gov)
# License: BSD 3 clause
# 2016.05.10 fvj tested. Works properly

import matplotlib.pyplot as plt
import pyart

filepath = '/data/pyart_examples/'
filename = '110635.mdv'

# create the plot using RadarDisplay
radar = pyart.io.read_mdv(filepath+filename)
display = pyart.graph.RadarDisplay(radar)
fig = plt.figure(figsize=[5, 5])
ax = fig.add_subplot(111)
display.plot('reflectivity', 0, vmin=-16, vmax=64.0)
plt.show()
