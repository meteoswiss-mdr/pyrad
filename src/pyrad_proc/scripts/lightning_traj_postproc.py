#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
lightning_traj_postproc
================================================

This program post-processes lightning trajectories

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import numpy as np

from pyrad.io import read_lightning_traj
from pyrad.graph import plot_histogram, get_colobar_label
from pyrad.io import get_fieldname_pyart
from pyrad.util import compute_histogram


print(__doc__)


def main():
    """
    """

    fname = '/data/pyrad_products/rad4alp_hydro_PHA/2017-06-29/hydroclass_traj/AT_FLASH/20170629130112_allflash_ts_trajlightning_hydro.csv'
    fname2 = '/data/pyrad_products/rad4alp_hydro_PHA/2017-06-29/hydroclass_traj/AT_FLASH/20170629130112_firstsource_ts_trajlightning_hydro.png'
    datatype = 'hydro'

    print("====== Lightning post-processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Lightning post-processing finished: ")

    (time_flash, flashnr, dBm, val_at_flash, val_mean, val_min, val_max,
     nval) = read_lightning_traj(fname)
     
    flashnr_first, unique_ind = np.unique(flashnr, return_index=True)

    val_first = val_at_flash[unique_ind]
    time_flash_first = time_flash[[unique_ind]
    
    bins, vals = compute_histogram(val_first, get_fieldname_pyart(datatype))
    
    labelx = get_colobar_label(dataset.fields[field_name], field_name)
    plot_histogram(bins, values, [fname2], labelx=labelx,
                   titl=("Trajectory Histogram First Source %s" %
                          time_flash_first[0].strftime("%Y-%m-%d")))
                          
    print("----- plot to '%s'" % fname2)


def _print_end_msg(text):
    """
    prints end message

    Parameters
    ----------
    text : str
        the text to be printed

    Returns
    -------
    Nothing

    """
    print(text + datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
