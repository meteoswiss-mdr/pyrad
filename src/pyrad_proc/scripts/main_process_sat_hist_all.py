#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_sat_hist_all
================================================

This program puts together histograms of satellite data

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

from pyrad.io import read_histogram, get_fieldname_pyart
from pyrad.io import write_histogram
from pyrad.graph import plot_histogram2, get_colobar_label

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
        default='/store/msrad/radar/pyrad_products/MSG_ML_hist/',
        help='base path to the radar data')

    parser.add_argument(
        '--datadirs', type=str,
        default='HISTOGRAMS_ALL_DATA,HISTOGRAMS_NOPOH90,HISTOGRAMS_POH90',
        help='base path to the radar data')

    parser.add_argument(
        '--datatypes', type=str,
        default='HRV,HRV_norm,IR_108,WV_062-IR_108,HRV_norm_text,HRV_text,IR_108_text,WV_062-IR_108_text',
        help='Data types. Coma separated')

    args = parser.parse_args()

    print("====== PYRAD windmill data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD windmill data processing finished: ")

    datadirs = args.datadirs.split(',')
    datatypes = args.datatypes.split(',')

    # Read periods of processing
    for datadir in datadirs:
        for datatype in datatypes:
            # Read data time series files
            flist = sorted(glob.glob(
                args.database+'*/'+datadir+'/'+datatype +
                '/??????????????_histogram_*.csv'))

            if not flist:
                warn('No file in '+args.database+'*/'+datadir+'/'+datatype)
                continue
                
            first_read = False
            for fname in flist:
                print('Reading histogram file '+fname)
                hist_aux, bin_edges_aux = read_histogram(fname)
                if not first_read:
                    hist = hist_aux
                    bin_edges = bin_edges_aux
                    first_read = True
                    continue
                    
                if bin_edges_aux.size != bin_edges.size:
                    warn('Number of bins not identical: ' +
                         str(bin_edges_aux.size)+' '+str(bin_edges.size))
                    continue
                    
                if not np.allclose(bin_edges_aux, bin_edges):
                    warn('Histogram bins not identical')
                    continue
                hist += hist_aux

            fname = args.database+datadir+'_'+datatype+'.png'

            bin_centers = bin_edges[:-1]+((bin_edges[1]-bin_edges[0])/2.)
            fname = plot_histogram2(
                bin_centers, hist, [fname],
                labelx=datatype, titl=datadir+' '+datatype)
            print('Plotted '+' '.join(fname))

            fname = args.database+datadir+'_'+datatype+'.csv'
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
