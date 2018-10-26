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
from copy import deepcopy

import numpy as np

from pyrad.io import read_histogram
from pyrad.graph import plot_histogram2

print(__doc__)


def main():
    """
    """
    # basepath = '/data/pyrad_products/rad4alp_hydro_PHA/'
    basepath = '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/data_analysis/'

    print("====== Lightning post-processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Lightning post-processing finished: ")


    # sources = 'allsources'
    sources = 'firstsource'
    datatype = 'ZDRc'
    hist_all, bin_edges_all = read_histogram(
        basepath+'all_'+sources+'_ts_trajlightning_'+datatype+'.csv')
    hist_solid, bin_edges_solid = read_histogram(
        basepath+'solid_'+sources+'_ts_trajlightning_'+datatype+'.csv')
    hist_liquid, bin_edges_liquid = read_histogram(
        basepath+'liquid_'+sources+'_ts_trajlightning_'+datatype+'.csv')

    if (not np.array_equal(bin_edges_all, bin_edges_solid) or
            not np.array_equal(bin_edges_all, bin_edges_liquid)):
        warn('Bin edges should be identical to group histograms')
        return


    bin_res = bin_edges_all[1]-bin_edges_all[0]
    bin_centers = bin_edges_all[1:]-bin_res/2.
    fname = basepath+'group_Santis_hist_ZDRc_first_source.png'
    labelx = 'differential reflectivity [dB]'
    titl = 'ZDR at flash origin location'
    invert_xaxis=False

    fig, ax = plot_histogram2(
        bin_centers, hist_all, [fname], labelx=labelx, titl=titl, alpha=0.25,
        save_fig=False, color='b', invert_xaxis=invert_xaxis)

    fig, ax = plot_histogram2(
        bin_centers, hist_liquid, [fname], labelx=labelx, titl=titl,
        ax=ax, fig=fig, save_fig=False, color='g', alpha=0.25,
        invert_xaxis=invert_xaxis)

    fname_list = plot_histogram2(
        bin_centers, hist_solid, [fname], labelx=labelx, titl=titl,
        ax=ax, fig=fig, save_fig=True, color='r', alpha=0.25,
        invert_xaxis=invert_xaxis)

    print('plotted '+''.join(fname_list))





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
