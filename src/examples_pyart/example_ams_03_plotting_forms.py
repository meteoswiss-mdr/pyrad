"""
================================
Example ams pyart course 3
================================

This is the third example in the AMS Radar conference 2015 pyart course.
Plots data read in various formats

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [12.0, 9.0]
import pyart

radarpath='/data/pyart_examples/'

# reading and plotting MDV file from an ARM XSAPR radar
radarfile='110635.mdv'
radar = pyart.io.read_mdv(radarpath+radarfile)

fig = plt.figure()
display = pyart.graph.RadarDisplay(radar)
display.plot('reflectivity', vmin=-16, vmax=80, cmap='pyart_NWSRef')
plt.show()

# reading and plotting Sigmet/IRIS file from an ARM XSAPR radar
radarfile='XSW110520113537.RAW7HHL'
radar = pyart.io.read_sigmet(radarpath+radarfile)

fig = plt.figure()
display = pyart.graph.RadarDisplay(radar)
display.plot('reflectivity', vmin=-32, vmax=80, cmap='pyart_NWSRef')
plt.show()


# reading and plotting NEXRAD Level 2 file
radarfile='KOKX20110828_072107_V03.gz'
radar = pyart.io.read_nexrad_archive(radarpath+radarfile)

fig = plt.figure()
display = pyart.graph.RadarDisplay(radar)
display.plot('reflectivity', vmin=0, vmax=80, cmap='pyart_NWSRef')
plt.show()

# # reading and plotting ARM's KAZR radar (produces an error)
# radarfile='sgpkazrhiC1.a1.20110503.000001.cdf'
# radar = pyart.aux_io.read_kazr(radarpath+radarfile)
# 
# fig = plt.figure()
# display = pyart.graph.RadarDisplay(radar)
# display.plot('reflectivity_xpol', vmin=-70, vmax=0)
# plt.show()

# reading and plotting any file with the pyart.io.read
radarfile='110635.mdv'
radar = pyart.io.read(radarpath+radarfile)
# try the pyart.io.read function with othe other files used thus far..

fig = plt.figure()
display = pyart.graph.RadarDisplay(radar)
display.plot('reflectivity', vmin=-20, vmax=60, cmap='pyart_NWSRef')
plt.show()


