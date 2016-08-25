"""
=========================================
Create a b-scope plot from a RAINBOW file
=========================================

An example which creates a b-scope plot of a RAINBOW file

"""
print(__doc__)

# Author: fvj
# License: BSD 3 clause

import matplotlib.pyplot as plt
import pyart

# filepath='/data/DX50/rawdata/MALS_PAY_NO_FILTERING.vol/2016-05-23/'
# filename = '2016052300065100dBZ.vol'
# iel=1

filepath='/data/DX50/rawdata/MALS_PAY_zdrcal01.azi/2016-05-23/'
filename = '2016052300022500dBZ.azi'
iel=0

# read the data
# radar = pyart.aux_io.read_rainbow(filepath+filename)
radar = pyart.aux_io.read_rainbow_wrl(filepath+filename)

# plot data from the entire volume in a b-scope representation
fig = plt.figure()
plt.imshow(radar.fields['reflectivity']['data'], aspect=0.5, origin='bottom', vmin=-3, vmax=60.0)
plt.xlabel('range gate')
plt.ylabel('ray number')
plt.show()

# plot data from a single volume in a b-scope representation
refl_sweep_data = radar.get_field(sweep=iel, field_name='reflectivity')

fig = plt.figure()
plt.imshow(refl_sweep_data, aspect=0.5, origin='bottom', vmin=-3, vmax=60.0)
plt.xlabel('range gate')
plt.ylabel('ray number')
plt.show()

# plot iteratively
fig = plt.figure(figsize=[30, 5])
for i, refl_sweep_data in enumerate(radar.iter_field('reflectivity')):
    ax = fig.add_subplot(1, radar.nsweeps, i+1) 
    ax.imshow(refl_sweep_data, vmin=-3, vmax=60.0)
plt.show()