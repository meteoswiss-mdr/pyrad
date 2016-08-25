"""
================================
Example ams pyart course 12
================================

This is the 12th example in the AMS Radar conference 2015 pyart course.
Plot current nexrad using artview

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import os
import urllib
import pyart
#import numpy as np
import matplotlib.pyplot as plt
#import artview

class Get88D(object):
	"""
	A class method to retrieve and plot NEXRAD data.
	
	The metar retrieval is a modification of code found at:
	https://github.com/akrherz/iem/blob/master/scripts/asos/iem_scraper_example.py
	
	It is dependent upon the Iowa State Mesonet database.
	"""
	def __init__(self, radarID):
		'''Initialize the class'''
		
		# Set date formats to be used with datetime
		self.d_fmt = "%Y-%m-%d %H:%M"
		self.dout_fmt = "%Y-%m-%d_%H:%M"
		
		# Use passed arguments
		self.radarID = radarID
		
	def get_data(self):
		'''Function to return metar data and create output text file'''
		# Query the radarID directory to get file list
		SERVICE = "http://nomads.ncep.noaa.gov/pub/data/nccf/radar/nexrad_level2"
		nexurl = '%s/%s/'%(SERVICE, self.radarID)
		response = urllib.request.urlopen("%s%s"%(nexurl, "dir.list"))
		help(response.read())
		self.list88D = response.read().split("\n")
		
		# At this point you have a list of data files, BUT there are 2 columns
		# col 1 = file size, col2 = filename
		
		# Now grab the latest data file and save it locally to open
		data = urllib.request.urlopen("%s%s"%(nexurl,self.list88D[-2].split(" ")[1]))
		with open("latest88D.bz", "wb") as code:
			code.write(data.read())
		
	def plot_nexrad(self, vmin=None, vmax=None, xlims=None, ylims=None):
		'''Create a plot'''
		# Create a PyArt radar instance
		fig, ax = plt.subplots()
		#ax = plt.axes()
		self.r=pyart.io.read_nexrad_archive("latest88D.bz")
		d=pyart.graph.RadarDisplay(self.r)
		d.plot('reflectivity',0, vmin=vmin, vmax=vmax, cmap="pyart_Carbone42")
		if xlims is None:
			xlims = (-250., 250)
		if ylims is None:
			ylims = (-250., 250)
		d.set_limits(xlims, ylims)
		
	def remove_bz(self):
		os.remove("latest88D.bz")
		
radarID = "KABR"#KCYS"
Radar = Get88D(radarID)
Radar.get_data()
Radar.plot_nexrad(vmin=-10, vmax=60.)
Radar.remove_bz()

plt.show()

