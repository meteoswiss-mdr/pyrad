#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_extract_trt
================================================

This program extracts individual TRT cell data from the original files
and puts it in a separate file for each cell

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit

import numpy as np

from pyrad.io import read_trt_info_all, write_histogram, write_trt_info
from pyrad.graph import plot_histogram

print(__doc__)


def main():
    """
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # keyword arguments
    parser.add_argument(
        '--infobase', type=str,
        default='/store/msrad/radar/thundertracking/info/',
        help='name of folder containing the TRT cell data')

    args = parser.parse_args()

    print("====== TRT thundertracking info started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== TRT thundertracking info finished: ")

    trt_time, ids, rank, nscans, _, _, _, _, _, _, _, _, _, _ = (
        read_trt_info_all(args.infobase))

    bin_edges = [0, 12, 15, 25, 35, 40]

    fname_hist = args.infobase+'rank_hist.png'
    plot_histogram(
        bin_edges, rank, [fname_hist], labelx='Category',
        titl='TRT cell category histogram')

    print("----- plot to '%s'" % fname_hist)

    fname_hist = args.infobase+'rank_hist.csv'
    hist_values, _ = np.histogram(rank, bins=bin_edges)
    write_histogram(bin_edges, hist_values, fname_hist)

    print("----- written to '%s'" % fname_hist)

    id_uniques = np.unique(ids)
    print("Number of TRT scans : "+str(trt_time.size))
    print("Number of X-band radar scans: "+str(np.sum(nscans)))
    print("Number of unique cells: "+str(id_uniques.size))

    rank_max = np.array([])
    nscans_total = np.array([], dtype=int)
    trt_time_start = np.array([], dtype=datetime.datetime)
    trt_time_end = np.array([], dtype=datetime.datetime)
    for id_unique in id_uniques:
        rank_max_cell = np.max(rank[ids == id_unique])
        rank_max = np.append(rank_max, rank_max_cell)

        nscans_total = np.append(
            nscans_total, np.sum(nscans[ids == id_unique], dtype=int))

        trt_time_cell = sorted(trt_time[ids == id_unique])
        trt_time_start = np.append(trt_time_start, trt_time_cell[0])
        trt_time_end = np.append(trt_time_end, trt_time_cell[-1])

    max_scans = np.max(nscans_total)
    print("Max number of X-band radar scans per TRT cell: " +
          str(max_scans))

    print("Max rank per TRT cell: " + str(np.max(rank_max)))

    fname_hist = args.infobase+'max_rank_hist.png'
    plot_histogram(
        bin_edges, rank_max, [fname_hist], labelx='Category',
        titl='TRT cell max category histogram')

    print("----- plot to '%s'" % fname_hist)

    fname_hist = args.infobase+'max_rank_hist.csv'
    hist_values, _ = np.histogram(rank_max, bins=bin_edges)
    write_histogram(bin_edges, hist_values, fname_hist)

    print("----- written to '%s'" % fname_hist)

    bin_edges = np.append(np.arange(1, max_scans+1)-0.5, max_scans+0.5)

    fname_hist = args.infobase+'nscans_hist.png'
    plot_histogram(
        bin_edges, nscans_total, [fname_hist], labelx='Number of X-band scans',
        titl='Number of X-band scans per TRT cell')

    print("----- plot to '%s'" % fname_hist)

    fname_hist = args.infobase+'nscans_hist.csv'
    hist_values, _ = np.histogram(nscans_total, bins=bin_edges)
    write_histogram(bin_edges, hist_values, fname_hist)

    print("----- written to '%s'" % fname_hist)

    # list from more to less severe
    ind = np.flip(np.argsort(rank_max))

    id_uniques = id_uniques[ind]
    rank_max = rank_max[ind]
    nscans_total = nscans_total[ind]
    trt_time_start = trt_time_start[ind]
    trt_time_end = trt_time_end[ind]

    fname = args.infobase+'thundertracking_info.csv'
    write_trt_info(
        id_uniques, rank_max, nscans_total, trt_time_start, trt_time_end,
        fname)


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
