"""
======================================
Create a PPI plot from a RAINBOW file
======================================

An example which creates a PPI plot of a RAINBOW file using a RadarDisplay object.

"""
print(__doc__)

# Author: fvj
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

# # MCH DX50 file
# filepath='/data/DX50/rawdata/MALS_PAY_NO_FILTERING.vol/2016-05-23/'
# datetime = '2016052300065100'
# voltype='.vol'

filepath='/data/DX50/rawdata/MALStest02.vol/2012-09-25/'
datetime = '2012092523515000'
voltype='.vol'

# filepath='/data/DX50/rawdata/ThunderTracking_PPI.azi/2015-06-07/'
# datetime='2015060713431500'
# voltype='.azi'

# # ARPAP DX50 file
# filepath='/data/ARPAP/rawdata/vercelli_longpulse2.vol/2014-05-22/'
# datetime = '2014052215590900'
# voltype='.vol'

# wradlib example DX50 file
# filepath='/data/pyart_examples/'
# datetime = '2013070308340000'
# voltype='.azi'

iel=1
datatype='dBuZ'

if datatype == 'dBZ':
	field='reflectivity'
	minv=-3.
	maxv=60.
elif datatype == 'dBuZ':
	field='unfiltered_reflectivity'
	minv=-3.
	maxv=60.
elif datatype == 'ZDR':
	field='differential_reflectivity'
	minv=-1.
	maxv=5.
elif datatype =='PhiDP':
	field='differential_phase'
	minv=-10.
	maxv=130.
elif datatype =='uPhiDP':
	field='uncorrected_differential_phase'
	minv=-180.
	maxv=180.
elif datatype == 'RhoHV':
	field='cross_correlation_ratio'
	minv=0.7
	maxv=1.
elif datatype == 'KDP':
	field='specific_differential_phase'
	minv=-2.
	maxv=5.
elif datatype == 'V':
	field='velocity'
	minv=-30.
	maxv=30.
elif datatype == 'W':
	field='spectrum_width'
	minv=0.
	maxv=4.

filename=datetime+datatype+voltype

# create the plot using RadarDisplay
#radar = pyart.aux_io.read_rainbow(filepath+filename)
radar = pyart.aux_io.read_rainbow_wrl(filepath+filename)
display = pyart.graph.RadarDisplay(radar)

fig = plt.figure(figsize=[12, 9])
ax = fig.add_subplot(111)
display.plot(field, iel, vmin=minv, vmax=maxv)
display.set_limits(ylim=[-50, 50], xlim=[-50, 50])
display.plot_range_rings([10, 20, 30, 40])
display.plot_cross_hair(5.)
plt.show()
