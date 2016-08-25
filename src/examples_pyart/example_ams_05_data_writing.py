"""
================================
Example ams pyart course 5
================================

This is the 5th example in the AMS Radar conference 2015 pyart course.
Writing out files

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import pyart

radarpath='/data/pyart_examples/'
radarfile='XSW110520113537.RAW7HHL'

radar = pyart.io.read(radarpath+radarfile)
pyart.io.write_cfradial(radarpath+'example_converted_sigmet_file.nc', radar)
