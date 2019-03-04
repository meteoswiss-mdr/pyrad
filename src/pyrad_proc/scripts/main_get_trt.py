#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_get_trt
================================================

This program reads a file containing dates where TRT data is needed and gets
the corresponding TRT files from the archive

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import subprocess
import os

import numpy as np

from pyrad.io import read_thundertracking_info

print(__doc__)


def main():
    """
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # keyword arguments
    parser.add_argument(
        '--info_file', type=str,
        default='/store/msrad/radar/thundertracking/info/thundertracking_info.csv',
        help='configuration file path')

    parser.add_argument(
        '--raw_trtbase', type=str,
        default='/store/msrad/radar/rad4alp/TRT/',
        help='name of folder containing the TRT cell data')

    args = parser.parse_args()

    print("====== TRT data fetching started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== TRT data fetching finished: ")

    _, _, _, trt_time_start, trt_time_end = (
        read_thundertracking_info(args.info_file))

    trt_times = np.append(trt_time_start, trt_time_end)

    trt_dates = np.array([], dtype=datetime.date)
    for trt_time in trt_times:
        trt_dates = np.append(trt_dates, trt_time.date())

    trt_dates = np.sort(np.unique(trt_dates))

    print(
        os.path.expanduser('~')+'/pyrad/tools/copyData/get_trt_data_cscs.sh')
    for trt_date in trt_dates:
        print(trt_date.strftime("%Y%m%d"))
        subprocess.run([
            os.path.expanduser('~') +
            '/pyrad/tools/copyData/get_trt_data_cscs.sh',
            '-d', trt_date.strftime("%Y%m%d"), '-p', args.raw_trtbase])


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
