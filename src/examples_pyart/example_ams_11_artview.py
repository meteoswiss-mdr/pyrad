"""
================================
Example ams pyart course 11
================================

This is the 11th example in the AMS Radar conference 2015 pyart course.
Using ARTview

"""
print(__doc__)

# Author: fvj
# 2016.05.10

import os
import pyart
import artview
#from artview.view import view
import matplotlib.pyplot as plt

fdir = '/data/pyart_examples/'
fname = 'KVNX20120501_040957_V06'
fpath = os.path.join(fdir, fname)
print(fpath)

r = pyart.io.read(fpath)
print(r.fields.keys())

artview.view.view(r)
