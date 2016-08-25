"""
======================================================
Create a plot of multiple moments from a RAINBOW file
======================================================

An example which creates a plot containing multiple moments taken from a
RAINBOW Archive file.

"""
print(__doc__)

# Author: fvj
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

#datatypes=['dBuZ', 'dBZ', 'dBuZv', 'dBZv', 
#		 'Vu', 'V', 'Vvu', 'Vv',
#		 'Wu', 'W', 'Wvu', 'Wv',
#		 'PhiDP', 'uPhiDP', 'uPhiDPu',
#		 'KDP', 'KDPu', 'uKDPu',
#		 'RhoHV', 'RhoHVu',
#		 'ZDR', 'ZDRu',
#		 'SQI', 'SQIu', 'SQIvu', 'SQIv']

# # MCH DX50 file
# datapath='/data/DX50/rawdata/MALS_PAY_NO_FILTERING.vol/2016-05-23/'
# datetime='2016052300065100'
# filetype='.vol'
# iel=1
# 
# datatypes=['dBuZ', 'dBZ', 'dBuZv', 'dBZv', 
# 		 'Vu', 'V', 'Vvu', 'Vv',
# 		 'Wu', 'W', 'Wvu', 'Wv',
# 		 'PhiDP', 'uPhiDP', 'uPhiDPu',
# 		 'KDP',
# 		 'RhoHV', 'RhoHVu',
# 		 'ZDR', 'ZDRu']
		 
		 
# datapath='/data/DX50/rawdata/ThunderTracking_PPI.azi/2015-06-07/'
# datetime='2015060713431500'
# filetype='.azi'
# iel=0
# 
# datatypes=['dBuZ', 'dBZ',
# 		 'Vu', 'V',
# 		 'Wu', 'W',
# 		 'PhiDP', 'uPhiDP', 'uPhiDPu',
# 		 'KDP',
# 		 'RhoHV', 'RhoHVu',
# 		 'ZDR', 'ZDRu']

		 
# ARPAP DX50 file
datapath='/data/ARPAP/rawdata/vercelli_longpulse2.vol/2014-05-22/'
datetime = '2014052215590900'
filetype='.vol'
iel=1

datatypes=['dBuZ', 'dBZ',
		 'V',
		 'W',
		 'uPhiDP',
		 'RhoHV',
		 'ZDR']
		 
		 
		 


ndatatypes=len(datatypes)
limx=(-50, 50)
limy=(-50, 50) 

# read master data type
filename=datetime+datatypes[0]+filetype
print(filename)
radar = pyart.aux_io.read_rainbow_wrl(datapath+filename)

# add other fields
for i in range(1, ndatatypes):
	filename=datetime+datatypes[i]+filetype
	radar_aux = pyart.aux_io.read_rainbow_wrl(datapath+filename)
	for field_name in radar_aux.fields.keys():
		break	
	print(filename, field_name)	
	field_data = radar_aux.fields[field_name]['data']	
	field_metadata = pyart.config.get_metadata(field_name)
	field_metadata['data'] = field_data
	radar.add_field(field_name, field_metadata)

display = pyart.graph.RadarDisplay(radar)

# figure 1
fig = plt.figure(figsize=(15, 30))

if 'reflectivity' in radar.fields.keys():
	ax = fig.add_subplot(7,4,1)
	display.plot('reflectivity', iel, ax=ax, title='Horizontal Reflectivity', vmin=-20., vmax=60., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'reflectivity_vv' in radar.fields.keys():
	ax = fig.add_subplot(7,4,2)
	display.plot('reflectivity_vv', iel, ax=ax, title='Vertical Reflectivity', vmin=-20., vmax=60., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_reflectivity' in radar.fields.keys():
	ax = fig.add_subplot(7,4,3)
	display.plot('unfiltered_reflectivity', iel, ax=ax, title='Horizontal Unfiltered', vmin=-20., vmax=60., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_reflectivity_vv' in radar.fields.keys():
	ax = fig.add_subplot(7,4,4)
	display.plot('unfiltered_reflectivity_vv', iel, ax=ax, title='Vertical Unfiltered', vmin=-20., vmax=60., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)


if 'differential_reflectivity' in radar.fields.keys():
	ax = fig.add_subplot(7,4,5)
	display.plot('differential_reflectivity', iel, ax=ax, title='Differential Reflectivity', vmin=-1., vmax=5., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_differential_reflectivity' in radar.fields.keys():
	ax = fig.add_subplot(7,4,6)
	display.plot('unfiltered_differential_reflectivity', iel, ax=ax, title='Unfiltered ZDR', vmin=-1., vmax=5., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'cross_correlation_ratio' in radar.fields.keys():
	ax = fig.add_subplot(7,4,7)
	display.plot('cross_correlation_ratio', iel, ax=ax, title='Correlation Coefficient', vmin=0.7, vmax=1., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_cross_correlation_ratio' in radar.fields.keys():
	ax = fig.add_subplot(7,4,8)
	display.plot('unfiltered_cross_correlation_ratio', iel, ax=ax, title='Unfiltered RhoHV', vmin=0.7, vmax=1., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)


if 'differential_phase' in radar.fields.keys():
	ax = fig.add_subplot(7,4,9)
	display.plot('differential_phase', iel, ax=ax, title='Differential Phase', vmin=-10., vmax=130., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'uncorrected_differential_phase' in radar.fields.keys():
	ax = fig.add_subplot(7,4,10)
	display.plot('uncorrected_differential_phase', iel, ax=ax, title='Uncorrected', vmin=-180., vmax=180., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'uncorrected_unfiltered_differential_phase' in radar.fields.keys():
	ax = fig.add_subplot(7,4,11)
	display.plot('uncorrected_unfiltered_differential_phase', iel, ax=ax, title='Uncorrected Unfiltered', vmin=-180., vmax=180., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)


if 'specific_differential_phase' in radar.fields.keys():
	ax = fig.add_subplot(7,4,13)
	display.plot('specific_differential_phase', iel, ax=ax, title='Specific Differential Phase', vmin=-2., vmax=5., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'uncorrected_specific_differential_phase' in radar.fields.keys():
	ax = fig.add_subplot(7,4,14)
	display.plot('uncorrected_specific_differential_phase', iel, ax=ax, title='Uncorrected', vmin=-2., vmax=5., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'uncorrected_unfiltered_specific_differential_phase' in radar.fields.keys():
	ax = fig.add_subplot(7,4,15)
	display.plot('uncorrected_unfiltered_specific_differential_phase', iel, ax=ax, title='Uncorrected Unfiltered', vmin=-2., vmax=5., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)


if 'velocity' in radar.fields.keys():
	ax = fig.add_subplot(7,4,17)
	display.plot('velocity', iel, ax=ax, title='Doppler Velocity', vmin=-20., vmax=20., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'velocity_vv' in radar.fields.keys():
	ax = fig.add_subplot(7,4,18)
	display.plot('velocity_vv', iel, ax=ax, title='Vertical', vmin=-20., vmax=20., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_velocity' in radar.fields.keys():
	ax = fig.add_subplot(7,4,19)
	display.plot('unfiltered_velocity', iel, ax=ax, title='Unfiltered', vmin=-20., vmax=20., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_velocity_vv' in radar.fields.keys():
	ax = fig.add_subplot(7,4,20)
	display.plot('unfiltered_velocity_vv', iel, ax=ax, title='Unfiltered Vertical', vmin=-20., vmax=20., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)


if 'spectrum_width' in radar.fields.keys():
	ax = fig.add_subplot(7,4,21)
	display.plot('spectrum_width', iel, ax=ax, title='Spectrum Width', vmin=0., vmax=4., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'spectrum_width_vv' in radar.fields.keys():
	ax = fig.add_subplot(7,4,22)
	display.plot('spectrum_width_vv', iel, ax=ax, title='Vertical', vmin=0., vmax=4., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_spectrum_width' in radar.fields.keys():
	ax = fig.add_subplot(7,4,23)
	display.plot('unfiltered_spectrum_width', iel, ax=ax, title='Unfiltered', vmin=0., vmax=4., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_spectrum_width_vv' in radar.fields.keys():
	ax = fig.add_subplot(7,4,24)
	display.plot('unfiltered_spectrum_width_vv', iel, ax=ax, title='Unfiltered Vertical', vmin=0., vmax=4., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)
	

if 'signal_quality_index' in radar.fields.keys():
	ax = fig.add_subplot(7,4,25)
	display.plot('signal_quality_index', iel, ax=ax, title='Quality Index', vmin=0., vmax=4., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'signal_quality_index_vv' in radar.fields.keys():
	ax = fig.add_subplot(7,4,26)
	display.plot('signal_quality_index_vv', iel, ax=ax, title='Vertical', vmin=0., vmax=4., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_signal_quality_index' in radar.fields.keys():
	ax = fig.add_subplot(7,4,27)
	display.plot('unfiltered_signal_quality_index', iel, ax=ax, title='Unfiltered', vmin=0., vmax=4., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

if 'unfiltered_signal_quality_index_vv' in radar.fields.keys():
	ax = fig.add_subplot(7,4,28)
	display.plot('unfiltered_signal_quality_index_vv', iel, ax=ax, title='Unfiltered Vertical', vmin=0., vmax=4., colorbar_label='', axislabels=('', ''))
	display.set_limits(limx, limy, ax=ax)

plt.show()