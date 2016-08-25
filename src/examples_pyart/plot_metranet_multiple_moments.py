"""
======================================================
Create a plot of multiple moments from a METRANET file
======================================================

An example which creates a plot containing multiple moments taken from a
METRANET Archive file.

"""
print(__doc__)

# Author: fvj
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

datapath = '/data/rad4alp/rawdata/15279/PHA15279/'
filename = 'PHA1527913307U.003'
iel=0
limx=(-150, 150)
limy=(-150, 150)

radar = pyart.aux_io.read_metranet(datapath+filename)
display = pyart.graph.RadarDisplay(radar)

# figure 1
fig = plt.figure(figsize=(20, 20))

ax = fig.add_subplot(4,3,1)
display.plot('reflectivity', iel, ax=ax, title='Horizontal Reflectivity', vmin=-20., vmax=60., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,2)
display.plot('reflectivity_vv', iel, ax=ax, title='Vertical Reflectivity', vmin=-20., vmax=60., colorbar_label='', axislabels=('', ''))
display.set_limits((-150, 150), (-150, 150), ax=ax)

ax = fig.add_subplot(4,3,3)
display.plot('differential_reflectivity', iel, ax=ax, title='Differential Reflectivity', vmin=-1., vmax=5., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,4)
display.plot('uncorrected_differential_phase', iel, ax=ax, title='Uncorrected Differential Phase', vmin=-180., vmax=180., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,5)
display.plot('cross_correlation_ratio', iel, ax=ax, title='Correlation Coefficient', vmin=0.7, vmax=1., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,6)
display.plot('mean_phase', iel, ax=ax, title='Mean Phase', vmin=-10., vmax=10., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,7)
display.plot('velocity', iel, ax=ax, title='Doppler Velocity', vmin=-10., vmax=10., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,8)
display.plot('spectrum_width', iel, ax=ax, title='Spectrum Width', vmin=0., vmax=4., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,9)
display.plot('wide_band_noise', iel, ax=ax, title='Wide Band Noise', vmin=0., vmax=10., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,10)
display.plot('stat_test_lag1', iel, ax=ax, title='Test lag 1', vmin=0., vmax=8., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,11)
display.plot('stat_test_lag2', iel, ax=ax, title='test lag 2', vmin=0., vmax=8., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

ax = fig.add_subplot(4,3,12)
display.plot('clutter_exit_code', iel, ax=ax, title='Clutter', vmin=0., vmax=100., colorbar_label='', axislabels=('', ''))
display.set_limits(limx, limy, ax=ax)

plt.show()