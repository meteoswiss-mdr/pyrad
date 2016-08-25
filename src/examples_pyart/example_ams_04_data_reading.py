"""
================================
Example ams pyart course 4
================================

This is the 4th example in the AMS Radar conference 2015 pyart course.
Advance reading of file fields

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import pyart

radarpath='/data/pyart_examples/'
radarfile='XSW110520113537.RAW7HHL'

# read a file keeping the file field names
radar = pyart.io.read(radarpath+radarfile, file_field_names=False)
for field_name in radar.fields.keys():
    print(field_name)
print('-----')

# read a file and use the standard internal field names of pyart
radar = pyart.io.read(radarpath+radarfile, file_field_names=True)
for field_name in radar.fields.keys():
    print(field_name)
print('-----')

# read a file and costumize the field name
custom_mapping = {
    'DBZ2': 'refl',
    'VEL2': 'vel',
    'WIDTH2': 'width',
    'SQI2': 'sqi',
    'PHIDP2': None,
}

radar = pyart.io.read(radarpath+radarfile, field_names=custom_mapping)
for field_name in radar.fields.keys():
    print(field_name)
print('-----')

# exclude some fields from the reading
radar = pyart.io.read(radarpath+radarfile)
for field_name in radar.fields.keys():
    print(field_name)
print('-----')
	
radar = pyart.io.read(radarpath+radarfile, exclude_fields=['normalized_coherent_power'])
for field_name in radar.fields.keys():
    print(field_name)
print('-----')
	
custom_mapping = {
    'DBZ2': 'refl',
    'VEL2': 'vel',
    'WIDTH2': 'width',
    'SQI2': 'sqi',
    'PHIDP2': None,
}

radar = pyart.io.read(radarpath+radarfile, field_names=custom_mapping, exclude_fields=['sqi'])
for field_name in radar.fields.keys():
    print(field_name)
print('-----')


