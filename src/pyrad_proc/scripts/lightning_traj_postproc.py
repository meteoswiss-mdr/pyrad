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
import glob
from warnings import warn
import os

import numpy as np

from pyrad.io import read_lightning_traj
from pyrad.graph import plot_histogram, get_colobar_label
from pyrad.io import get_fieldname_pyart
from pyrad.util import compute_histogram

from pyart.config import get_metadata



print(__doc__)


def main():
    """
    """
    # basepath = '/data/pyrad_products/rad4alp_hydro_PHA/'
    basepath = '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/'
    day_vec = [
        datetime.datetime(2017, 7, 19),
        datetime.datetime(2017, 7, 30),
        datetime.datetime(2017, 8, 1)]

    datatype_vec = [
        'hydro',
        'KDPc',
        'dBZc',
        'RhoHVc',
        'TEMP',
        'ZDRc']

    step_list = [
        None,
        0.05,
        0.5,
        0.001,
        1.,
        0.1]


    print("====== Lightning post-processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Lightning post-processing finished: ")

    for i, day in enumerate(day_vec):
        day_dir = day.strftime('%Y-%m-%d')
        day_str = day.strftime('%Y%m%d')

        for j, datatype in enumerate(datatype_vec):
            step = step_list[j]
            fname_test = (
                basepath+day_dir+'/*_traj/AT_FLASH/'+day_str +
                '*_allflash_ts_trajlightning_'+datatype+'.csv')
            fname_list = glob.glob(fname_test)
            if not fname_list:
                warn('No file found in '+fname_test)
                continue

            fname = fname_list[0]

            basepath_out = os.path.dirname(fname)
            fname_first_source = (
                basepath_out+'/'+day_str+'_firstsource_ts_trajlightning_' +
                datatype+'.png')
                
            fname_all_sources = (
                basepath_out+'/'+day_str+'_allsources_ts_trajlightning_' +
                datatype+'.png')

            print('\nReading file '+fname)
            (time_flash, flashnr, dBm, val_at_flash, val_mean, val_min,
             val_max, nval) = read_lightning_traj(fname)
             
            field_name = get_fieldname_pyart(datatype)
            field_dict = get_metadata(field_name)
            
            labelx = get_colobar_label(field_dict, field_name)
            
            # Plot all sources histogram
            bins, values = compute_histogram(val_at_flash, field_name, step=step)            
            
            plot_histogram(
                bins, values, [fname_all_sources], labelx=labelx,
                titl=("Trajectory Histogram %s" %
                      time_flash[0].strftime("%Y-%m-%d")))
                      
            print("----- plot to '%s'" % fname_all_sources)
            
            
            # Get and plot first sources histogram
            flashnr_first, unique_ind = np.unique(flashnr, return_index=True)

            val_first = val_at_flash[unique_ind]
            time_flash_first = time_flash[unique_ind]

            bins, values = compute_histogram(val_first, field_name, step=step)
            
            plot_histogram(
                bins, values, [fname_first_source], labelx=labelx,
                titl=("Trajectory Histogram First Source %s" %
                      time_flash_first[0].strftime("%Y-%m-%d")))

            print("----- plot to '%s'" % fname_first_source)


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
