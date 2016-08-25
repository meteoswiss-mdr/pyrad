"""
================================
Example ams pyart course 15
================================

This is the 15th example in the AMS Radar conference 2015 pyart course.
Demonstrates the use of CSU_RadarTools: fuzzy logic hydrometeor classification and blended rainfall estimation

"""

# Author: fvj
# 2016.05.12
# there are problems with modules csu_kdp and csu_fhc in csu_radartools

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pyart
import glob
from pyart.io.common import radar_coords_to_cart
from skewt import SkewT
from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain, csu_dsd, csu_kdp, csu_misc)

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

def get_z_from_radar(radar):
    """Input radar object, return z from radar (km, 2D)"""
    azimuth_1D = radar.azimuth['data']
    elevation_1D = radar.elevation['data']
    srange_1D = radar.range['data']
    sr_2d, az_2d = np.meshgrid(srange_1D, azimuth_1D)
    el_2d = np.meshgrid(srange_1D, elevation_1D)[1]
    xx, yy, zz = radar_coords_to_cart(sr_2d/1000.0, az_2d, el_2d)
    return zz + radar.altitude['data']


def check_sounding_for_montonic(sounding):
    """
    So the sounding interpolation doesn't fail, force the sounding to behave
    monotonically so that z always increases. This eliminates data from
    descending balloons.
    """
    snd_T = sounding.soundingdata['temp']  # In old SkewT, was sounding.data
    snd_z = sounding.soundingdata['hght']  # In old SkewT, was sounding.data
    dummy_z = []
    dummy_T = []
    if not snd_T.mask[0]:  # May cause issue for specific soundings
        dummy_z.append(snd_z[0])
        dummy_T.append(snd_T[0])
        for i, height in enumerate(snd_z):
            if i > 0:
                if snd_z[i] > snd_z[i-1] and not snd_T.mask[i]:
                    dummy_z.append(snd_z[i])
                    dummy_T.append(snd_T[i])
        snd_z = np.array(dummy_z)
        snd_T = np.array(dummy_T)
    return snd_T, snd_z


def interpolate_sounding_to_radar(sounding, radar):
    """Takes sounding data and interpolates it to every radar gate."""
    radar_z = get_z_from_radar(radar)
    radar_T = None
    snd_T, snd_z = check_sounding_for_montonic(sounding)
    shape = np.shape(radar_z)
    rad_z1d = radar_z.ravel()
    rad_T1d = np.interp(rad_z1d, snd_z, snd_T)
    return np.reshape(rad_T1d, shape), radar_z

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
sndfile = 'snd_Darwin.txt'
radarfile = 'cfrad.20060119_170029.000_to_20060121_020810.000_CPOL_v1_PPI.nc'
radar = pyart.io.read(filepath+radarfile)
print(radar.fields.keys())
sounding = SkewT.Sounding(filepath+sndfile)

# Extract relevant radar fields	
dz = radar.fields['ZC']['data']
dr = radar.fields['ZD']['data']
kd = radar.fields['KD']['data']
rh = radar.fields['RH']['data']

# interpolate the sounding to radar grid
radar_T, radar_z = interpolate_sounding_to_radar(sounding, radar)

# run the classification
scores = csu_fhc.csu_fhc_summer(dz=dz, zdr=dr, rho=rh, kdp=kd, use_temp=True, band='C', T=radar_T)
fh = np.argmax(scores, axis=0) + 1

# add classification field to radar object
radar = add_field_to_radar_object(fh, radar)

# plot the result
hid_colors = ['White', 'LightBlue', 'MediumBlue', 'DarkOrange', 'LightPink',
              'Cyan', 'DarkGray', 'Lime', 'Yellow', 'Red', 'Fuchsia']
cmaphid = colors.ListedColormap(hid_colors)
cmapmeth = colors.ListedColormap(hid_colors[0:6])

lim = [-80, 80]
fig, ax1, ax2, display = two_panel_plot(radar, sweep=5, var1='ZC', var2='ZD', vmin2=0, 
                                        vmax2=10, cmap2=cmaphid, units2='', return_flag=True, 
                                        xlim=lim, ylim=lim)
display.cbs[1] = adjust_fhc_colorbar_for_pyart(display.cbs[1])

plt.show()

# ---------------------------------------------------
# test of rainfall estimation
# ---------------------------------------------------
rain, method = csu_blended_rain.csu_hidro_rain(dz=dz, zdr=dr, kdp=kd, fhc=fh)
radar = add_field_to_radar_object(rain, radar, field_name='rain', units='mm h-1',
                                  long_name='HIDRO Rainfall Rate', 
                                  standard_name='Rainfall Rate')
radar = add_field_to_radar_object(method, radar, field_name='method', units='',
                                  long_name='HIDRO Rainfall Method', 
                                  standard_name='Rainfall Method')

lim = [-80, 80]
fig, ax1, ax2, display = two_panel_plot(radar, sweep=0, var1='rain', vmin1=0, vmax1=150,
                                        cmap1='GnBu', var2='method', vmin2=0, vmax2=5, 
                                        cmap2=cmapmeth, units2='', return_flag=True, 
                                        xlim=lim, ylim=lim, units1='mm h-1')
display.cbs[1] = adjust_meth_colorbar_for_pyart(display.cbs[1])

plt.show()

# ---------------------------------------------------
# test of blended rain rainfall estimation
# ---------------------------------------------------
rain, method, zdp, fi = csu_blended_rain.calc_blended_rain(dz=dz, zdr=dr, 
                                                           kdp=kd, ice_flag=True)
radar = add_field_to_radar_object(rain, radar, field_name='rain_blend', units='mm h-1',
                                  long_name='Blended Rainfall Rate', 
                                  standard_name='Rainfall Rate')
radar = add_field_to_radar_object(method, radar, field_name='method_blend', units='',
                                  long_name='Blended Rainfall Method', 
                                  standard_name='Rainfall Method')
radar = add_field_to_radar_object(zdp, radar, field_name='ZDP', units='dB',
                                  long_name='Difference Reflectivity',
                                  standard_name='Difference Reflectivity')
radar = add_field_to_radar_object(fi, radar, field_name='FI', units='', 
                                  long_name='Ice Fraction',
                                  standard_name='Ice Fraction')

lim = [-80, 80] 
fig, ax1, ax2, display = two_panel_plot(radar, sweep=0, var1='rain_blend', vmin1=0, vmax1=150,
                                        cmap1='GnBu', var2='method_blend', vmin2=0, vmax2=5, 
                                        cmap2=cmapmeth, units2='', return_flag=True, xlim=lim, 
                                        ylim=lim, units1='mm h-1')
display.cbs[1] = adjust_meth_colorbar_for_pyart(display.cbs[1])

plt.show()

# plot difference in reflectivity and ice fraction
two_panel_plot(radar, sweep=0, var1='ZDP', units1='dB', vmin1=0, vmax1=65,
               cmap1='cubehelix', var2='FI', vmin2=0, vmax2=1, 
               cmap2='YlOrRd_r', units2='', xlim=lim, ylim=lim)

plt.show()
