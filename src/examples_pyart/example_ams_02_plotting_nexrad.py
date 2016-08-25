"""
================================
Example ams pyart course 2
================================

This is the second example in the AMS Radar conference 2015 pyart course.
Downloads the latest NEXRAD reflectivity scan from the US National Weather Service and plots it in a Cartesian map. 

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import matplotlib.pyplot as plt
import pyart
import urllib

# The NEXRAD site from which data will be 
nexrad_site = 'ktlx'


url = ('ftp://tgftp.nws.noaa.gov/SL.us008001/DF.of/'
       'DC.radar/DS.p19r0/SI.' + nexrad_site.lower() + '/sn.last')
handle = urllib.request.urlopen(url)
radar = pyart.io.read_nexrad_level3(handle)

fig = plt.figure()
display = pyart.graph.RadarMapDisplay(radar)
display.plot_ppi_map(
    'reflectivity', vmin=-32, vmax=80, cmap='pyart_NWSRef',
    resolution='c', embelish=False)
#display.basemap.drawcounties()
plt.show()








