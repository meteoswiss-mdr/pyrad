"""
================================================
Create a RHI plot from a group of metranet files
================================================

An example which creates a RHI plot from a group of metranet files using a RadarDisplay object.

"""
print(__doc__)

# Author: fvj
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

import copy

radar='A'
res='H'
year='15'
day='279'
time='1330'
elevs=['001', '002', '003', '004', '005', '006', '007', '008', '009', '010', '011', '012', '013', '014', '015', '016', '017', '018', '019', '020']
elevs=['001', '003', '005', '007', '009', '011']
nelevs=len(elevs)
azim=[270.5] # azimuth of the RHI to plot
limx=[0, 120]
limy=[0, 10]

filepath='/data/rad4alp/rawdata/'+year+day+'/P'+res+radar+year+day+'/'
filebase='P'+res+radar+year+day+time+'7U.'

# read master data type
filename=filebase+elevs[0]
print(filename)
radar = pyart.aux_io.read_metranet(filepath+filename)

# merge the elevations into a single radar instance
for i in range(1, nelevs):
	filename=filebase+elevs[i]
	print(filename)
	radar_aux = pyart.aux_io.read_metranet(filepath+filename)	
	radar=pyart.util.radar_utils.join_radar(radar, radar_aux)
	
# extract cross-section
xsect = pyart.util.cross_section_ppi(radar, azim)

# display figure
display = pyart.graph.RadarDisplay(xsect)
fig = plt.figure(figsize=(10, 10))

ax = fig.add_subplot(221)
display.plot('reflectivity', 0, vmin=-32, vmax=64.)
display.set_limits(xlim=limx, ylim=limy )

ax = fig.add_subplot(222)
display.plot('differential_reflectivity', 0, vmin=-1., vmax=5.)
display.set_limits(xlim=limx, ylim=limy )

ax = fig.add_subplot(223)
display.plot('uncorrected_differential_phase', 0, vmin=-10., vmax=130.)
display.set_limits(xlim=limx, ylim=limy )

ax = fig.add_subplot(224)
display.plot('cross_correlation_ratio', 0, vmin=0.7, vmax=1.)
display.set_limits(xlim=limx, ylim=limy )

plt.tight_layout()
plt.show()
