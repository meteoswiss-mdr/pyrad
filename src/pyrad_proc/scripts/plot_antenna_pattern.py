#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
plot_antenna_pattern
================================================

This program plots an antenna pattern provided in a .csv file

"""

# Author: fvj
# License: BSD 3 clause

from pyrad.io import read_antenna_pattern
from pyrad.graph import plot_antenna_pattern

print(__doc__)


def main():
    """
    """

    file_path = '/home/lom/users/fvj/malsgit/config/antenna/'
    base_name_vec = ['ASR_HighBeamAzimuthPattern',
                     'ASR_HighBeamElevationPattern',
                     'ASR_LowBeamAzimuthPattern',
                     'ASR_LowBeamElevationPattern',
                     'PAR_AzAntenna_AzimuthPattern',
                     'PAR_AzAntenna_ElevationPattern']
    img_ext = '.png'
    linear = False
    twoway = False

    for base_name in base_name_vec:
        antpattern = read_antenna_pattern(
            file_path+base_name+'.csv', linear=linear, twoway=twoway)
        plot_antenna_pattern(
            antpattern, [file_path+base_name+img_ext], labelx='Angle [Deg]',
            linear=linear, twoway=twoway, title=base_name, ymin=None,
            ymax=None)

# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
