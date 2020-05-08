#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_windmill_hist
================================================

This program compiles histograms out of time series of windmill radar returns
from a particular point in space

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os
import glob
from warnings import warn

import numpy as np

from pyrad.io import read_timeseries, read_proc_periods, get_fieldname_pyart
from pyrad.io import write_histogram
from pyrad.graph import plot_histogram, get_colobar_label
from pyrad.util import compute_histogram

from pyart.config import get_metadata

print(__doc__)


def main():
    """
    """

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # keyword arguments
    parser.add_argument(
        '--database', type=str,
        default='/store/msrad/radar/pyrad_products/mals_sha_windmills_point_WM1_20200321-20200325/',
        help='base path to the radar data')

    parser.add_argument(
        '--datatypes', type=str, default='dBuZ,dBuZv,rcs_h,rcs_v,uPhiDPu,RhoHVu,ZDRu,Vu,Wu',
        help='Data types. Coma separated')

    # RhoHVu: 0.003937 starting (0 to 1)
    # ZDRu:   0.078740 starting (-8 to 12)
    # Vu:     0.125    starting (-15.8 to 15.8)
    # Wu:     0-125    starting (0 to 15.8)
    parser.add_argument(
        '--steps', type=str, default='0.5,0.5,0.5,0.5,1.,0.003937,0.078740,0.125,0.125',
        help='Histogram steps. Coma separated')

    args = parser.parse_args()

    print("====== PYRAD windmill data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD windmill data processing finished: ")

    datatypes = args.datatypes.split(',')
    steps = np.asarray(args.steps.split(','), dtype=float)

    # Read periods of processing
    for step, datatype in zip(steps, datatypes):
        # Read data time series files
        flist = glob.glob(
            args.database+datatype+'_TS/TS/ts_POINT_MEASUREMENT_'+datatype +
            '_*.csv')
        if not flist:
            continue

        dt_ts, values = read_timeseries(flist[0])

        # Histogram plots
        field_name = get_fieldname_pyart(datatype)
        field_dict = get_metadata(field_name)

        fname = (
            args.database+datatype+'_TS/TS/ts_POINT_MEASUREMENT_hist_' +
            datatype+'.png')

        bin_edges, vals = compute_histogram(values, field_name, step=step)
        fname = plot_histogram(
            bin_edges, vals, [fname],
            labelx=get_colobar_label(field_dict, field_name), titl=datatype)
        print('Plotted '+' '.join(fname))

        fname = (
            args.database+datatype+'_TS/TS/ts_POINT_MEASUREMENT_hist_' +
            datatype+'.csv')
        hist, _ = np.histogram(vals, bins=bin_edges)
        fname = write_histogram(bin_edges, hist, fname)
        print('Written '+fname)


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
