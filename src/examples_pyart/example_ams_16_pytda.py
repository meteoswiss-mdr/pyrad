"""
================================
Example ams pyart course 16
================================

This is the 16th example in the AMS Radar conference 2015 pyart course.
Demonstrates the use of PyTDA: a software to estimate turbulence

"""

# Author: fvj
# 2016.05.12

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import pyart
import pytda


def plot_list_of_fields(radar, sweep=0, fields=['reflectivity'], vmins=[0],
                        vmaxs=[65], units=['dBZ'], cmaps=['RdYlBu_r'],
                        return_flag=False, xlim=[-150, 150], ylim=[-150, 150],
                        mask_tuple=None):
    num_fields = len(fields)
    if mask_tuple is None:
        mask_tuple = []
        for i in np.arange(num_fields):
            mask_tuple.append(None)
    nrows = (num_fields + 1) // 2
    ncols = (num_fields + 1) % 2 + 1
    fig = plt.figure(figsize=(14.0, float(nrows)*5.5))
    display = pyart.graph.RadarDisplay(radar)
    for index, field in enumerate(fields):
        ax = fig.add_subplot(nrows, 2, index+1)
        display.plot_ppi(field, sweep=sweep, vmin=vmins[index],
                         vmax=vmaxs[index],
                         colorbar_label=units[index], cmap=cmaps[index],
                         mask_tuple=mask_tuple[index])
        display.set_limits(xlim=xlim, ylim=ylim)
    plt.tight_layout()
    if return_flag:
        return display

		
# read radar data
filepath='/data/pyart_examples/'
radarfile = 'cfrad.20101026_151323.000_to_20101026_151734.000_KGWX_v284_SUR.nc'
radar = pyart.io.read(filepath+radarfile)
print(radar.fields.keys())
print(radar.instrument_parameters['radar_beam_width_h']['data'][0])

# compute Eddy dissipation rate
pytda.calc_turb_vol(radar, name_sw='SW', name_dz='DZ', verbose=False,
                    gate_spacing=250.0/1000.0, use_ntda=False,
                    beamwidth=radar.instrument_parameters['radar_beam_width_h']['data'][0])

# plot result
plot_list_of_fields(radar, sweep=4, xlim=[-50, 50], ylim=[50, 150],
                    fields=['DZ', 'VR', 'SW', 'turbulence'],
                    vmins=[0, -30, 0, 0], vmaxs=[60, 30, 10, 1.0], 
                    units=['dBZ', 'm/s', 'm/s', 'EDR^0.33'],
                    cmaps=['pyart_LangRainbow12', 'seismic', 'YlOrRd', 'cubehelix'])

plt.show()

# write result
pyart.io.write_cfradial(filepath+'file_with_turbulence.nc', radar)


