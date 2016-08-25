"""
================================
Example ams pyart course 6
================================

This is the 6th example in the AMS Radar conference 2015 pyart course.
Pyart data model

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.figsize'] = [12.0, 9.0]
import pyart


radarpath='/data/pyart_examples/'
radarfile='XSW110520113537.RAW7HHL'

radar = pyart.io.read(radarpath+radarfile)

# check data type
type(radar)
print('-------------------------------------------------------------\n\n\n')

# information on class radar
#help(radar)
#print('-------------------------------------------------------------\n\n\n')

# information on data
radar.info('compact')
print('-------------------------------------------------------------\n\n\n')

# information on specific attribute
print(radar.scan_type)
print('-------------------------------------------------------------\n\n\n')

# information on elevation field
print(type(radar.elevation))
print(radar.elevation.keys())
print(radar.elevation['long_name'])
print(radar.elevation['standard_name'])
print(radar.elevation['units'])

print(radar.elevation['data'])

print('-------------------------------------------------------------\n\n\n')

# plot elevation angles data
fig = plt.figure()
plt.plot(radar.elevation['data'])
plt.show()

# plot elevation angles data of a particular sweep
fig = plt.figure()
plt.plot(radar.elevation['data'][radar.get_start(0):radar.get_end(0)])
plt.show()

# plot elevation angles data for a particular sweep. Use get_slice to get the slice limits
sweep_0_slice = radar.get_slice(0)

fig = plt.figure()
plt.plot(radar.elevation['data'][sweep_0_slice])
plt.show()

sweep_1_slice = radar.get_slice(1)
fig = plt.figure()
plt.plot(radar.elevation['data'][sweep_1_slice])
plt.show()

# plot elevation angles data using the get_elevation method. There is also the get_azimuth method
fig = plt.figure()
plt.plot(radar.get_elevation(0))
plt.show()


# plot data from the entire volume in a b-scope representation
fig = plt.figure()
plt.imshow(radar.fields['reflectivity']['data'], aspect=0.5, origin='bottom')
plt.xlabel('range gate')
plt.ylabel('ray number')
plt.show()

# plot data from a single volume in a b-scope representation
refl_sweep_data = radar.get_field(sweep=0, field_name='reflectivity')

fig = plt.figure()
plt.imshow(refl_sweep_data, aspect=0.5, origin='bottom')
plt.xlabel('range gate')
plt.ylabel('ray number')
plt.show()

# create a new radar instance containing a single sweep
print(radar.nrays)
print(radar.nsweeps)

radar2 = radar.extract_sweeps([0])

print(radar2.nrays)
print(radar2.nsweeps)

fig = plt.figure()
plt.imshow(radar2.fields['reflectivity']['data'], aspect=0.5, origin='bottom')
plt.xlabel('range gate')
plt.ylabel('ray number')
plt.show()

print('-------------------------------------------------------------\n\n\n')

# iterate over sweeps
for start in radar.iter_start():
    print(start)
	
print('-------------------------------------------------------------\n\n\n')

for end in radar.iter_end():
    print(end)
	
print('-------------------------------------------------------------\n\n\n')

for start, end in radar.iter_start_end():
    print(start, end)

print('-------------------------------------------------------------\n\n\n')

fig = plt.figure()
ax = fig.add_subplot(111)
for i, sweep_slice in enumerate(radar.iter_slice()):
    ax.plot(radar.elevation['data'][sweep_slice] + i * 20)
plt.show()

fig = plt.figure()
for i, refl_sweep_data in enumerate(radar.iter_field('reflectivity')):
    ax = fig.add_subplot(1, radar.nsweeps, i+1) 
    ax.imshow(refl_sweep_data)
plt.show()
