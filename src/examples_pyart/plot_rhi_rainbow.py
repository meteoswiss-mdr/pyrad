"""
=====================================
Create a RHI plot from a RAINBOW file
=====================================

An example which creates a RHI plot of a RAINBOW file using a RadarDisplay object.

"""
print(__doc__)

# Author: fvj
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

# # MCH DX50 file
# filepath = '/data/DX50/rawdata/MALS_PAY_049_up_nopsr.ele/2016-05-23/'
# datetime = '2016052300025900'
# 
# filepath = '/data/DX50/rawdata/MALS_PAY_047_dw_nopsr.ele/2016-05-23/'
# datetime = '2016052300035100'

filepath = '/data/DX50/rawdata/ThunderTracking_00_up.ele/2015-06-07/'
datetime = '2015060713435700'

# filepath = '/data/cosmo/TEMP/MALS_PAY_049_up_nopsr.ele/2016-05-23/'
# datetime = '2016052300000000'
# runtime = '2016052300000000'

# filepath = '/data/cosmo/ISO0/MALS_PAY_049_up_nopsr.ele/2016-05-23/'
# datetime = '2016052300000000'
# runtime = '2016052300000000'

limy=[0, 12]
limx=[0, 60]

# ARPAP DX50 file
# filepath = '/data/ARPAP/rawdata/vercelli_az162.ele/2014-05-22/'
# datetime = '2014052215570300'
#
# limy=[0, 12]
# limx=[-80, 80]

# filepath = '/data/ARPAP/rawdata/vercelli_az241.ele/2014-05-22/'
# datetime = '2014052215580900'
# 
# limy=[0, 12]
# limx=[-80, 0]



voltype='.ele'
iaz=0
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
elif datatype == 'TEMP':
	field='temperature'
	minv=-60.
	maxv=20.
elif datatype == 'ISO0':
	field='iso0'
	minv=0.
	maxv=4.
    


if (datatype =='ISO0') or (datatype =='TEMP'):
    filename=datatype+'_RUN'+runtime+'_DX50'+datetime+voltype
else:
    filename=datetime+datatype+voltype

# create the plot using RadarDisplay
#radar = pyart.aux_io.read_rainbow(filepath+filename)
radar = pyart.aux_io.read_rainbow_wrl(filepath+filename)
display = pyart.graph.RadarDisplay(radar)
fig = plt.figure(figsize=[15, 5])
ax = fig.add_subplot(111)
display.plot(field, iaz, vmin=minv, vmax=maxv)
display.set_limits(ylim=limy, xlim=limx)
plt.show()
