"""
================================
Example ams pyart course 13
================================

This is the 13th example in the AMS Radar conference 2015 pyart course.
Demonstrates the use of CSU_RadarTools: test DSD retrieval at C-band

"""

# Author: fvj
# 2016.05.12

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pyart
import glob
from pyart.io.common import radar_coords_to_cart
from skewt import SkewT
#from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain, csu_dsd, csu_kdp, csu_misc)
from csu_radartools import (csu_liquid_ice_mass, csu_blended_rain, csu_dsd, csu_misc)

def two_panel_plot(radar, sweep=0, var1='reflectivity', vmin1=0, vmax1=65, cmap1='RdYlBu_r', 
                   units1='dBZ', var2='differential_reflectivity', vmin2=-5, vmax2=5, 
                   cmap2='RdYlBu_r', units2='dB', return_flag=False, xlim=[-150,150],
                   ylim=[-150,150]):
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure(figsize=(13, 5))
    ax1 = fig.add_subplot(121)
    display.plot_ppi(var1, sweep=sweep, vmin=vmin1, vmax=vmax1, cmap=cmap1, 
                     colorbar_label=units1, mask_outside=True)
    display.set_limits(xlim=xlim, ylim=ylim)
    ax2 = fig.add_subplot(122)
    display.plot_ppi(var2, sweep=sweep, vmin=vmin2, vmax=vmax2, cmap=cmap2, 
                     colorbar_label=units2, mask_outside=True)
    display.set_limits(xlim=xlim, ylim=ylim)
    if return_flag:
        return fig, ax1, ax2, display

def add_field_to_radar_object(field, radar, field_name='FH', units='unitless', 
                              long_name='Hydrometeor ID', standard_name='Hydrometeor ID',
                              dz_field='ZC'):
    """
    Adds a newly created field to the Py-ART radar object. If reflectivity is a masked array,
    make the new field masked the same as reflectivity.
    """
    masked_field = np.ma.asanyarray(field)
    fill_value = -32768
    if hasattr(radar.fields[dz_field]['data'], 'mask'):
        setattr(masked_field, 'mask', radar.fields[dz_field]['data'].mask)
        fill_value = radar.fields[dz_field]['_FillValue']
    field_dict = {'data': masked_field,
                  'units': units,
                  'long_name': long_name,
                  'standard_name': standard_name,
                  '_FillValue': fill_value}
    radar.add_field(field_name, field_dict, replace_existing=True)
    return radar

def adjust_fhc_colorbar_for_pyart(cb):
    cb.set_ticks(np.arange(1.4, 10, 0.9))
    cb.ax.set_yticklabels(['Drizzle', 'Rain', 'Ice Crystals', 'Aggregates',
                           'Wet Snow', 'Vertical Ice', 'LD Graupel',
                           'HD Graupel', 'Hail', 'Big Drops'])
    cb.ax.set_ylabel('')
    cb.ax.tick_params(length=0)
    return cb

def adjust_meth_colorbar_for_pyart(cb):
    cb.set_ticks(np.arange(1.25, 5, 0.833))
    cb.ax.set_yticklabels(['R(Kdp, Zdr)', 'R(Kdp)', 'R(Z, Zdr)', 'R(Z)', 'R(Zrain)'])
    cb.ax.set_ylabel('')
    cb.ax.tick_params(length=0)
    return cb

	
# ----------------------------------------------------------------------
# test of fuzzy logic hydrometeor classification
# ----------------------------------------------------------------------

# Read in the data
filepath='/data/pyart_examples/'
radarfile = 'cfrad.20060119_170029.000_to_20060121_020810.000_CPOL_v1_PPI.nc'
radar = pyart.io.read(filepath+radarfile)
print(radar.fields.keys())

# Extract relevant radar fields	
dz = radar.fields['ZC']['data']
dr = radar.fields['ZD']['data']
kd = radar.fields['KD']['data']
rh = radar.fields['RH']['data']

# ---------------------------------------------------
# test of DSD retrieval
# ---------------------------------------------------
d0, Nw, mu = csu_dsd.calc_dsd(dz=dz, zdr=dr, kdp=kd, band='C')
radar = add_field_to_radar_object(d0, radar, field_name='D0', units='mm', 
                                  long_name='Median Volume Diameter',
                                  standard_name='Median Volume Diameter')
logNw = np.log10(Nw)
radar = add_field_to_radar_object(logNw, radar, field_name='NW', units='', 
                                  long_name='Normalized Intercept Parameter',
                                  standard_name='Normalized Intercept Parameter')
radar = add_field_to_radar_object(mu, radar, field_name='MU', units='', 
                                  long_name='Mu', standard_name='Mu')
								  
lim = [-80, 80]
two_panel_plot(radar, sweep=0, var1='D0', vmin1=0, vmax1=2, cmap1='GnBu', units1='mm',
               var2='NW', vmin2=0, vmax2=8, cmap2='cubehelix', units2='log10(Nw)', xlim=lim, 
               ylim=lim)
								  
plt.show()