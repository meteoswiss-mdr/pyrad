#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_data_thundertracking
================================================

This program activates the Pyrad processing to process
the polarimetric variables for all thundertracking data contained
in a csv file

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os

from pyrad.flow.flow_control import main as pyrad_main
from pyrad.io import read_thundertracking_info

print(__doc__)


def main():
    """
    """

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments
    parser.add_argument(
        'cfgfile_list', type=str, help='list of main configuration file. Coma separated')

    # keyword arguments
    parser.add_argument(
        '--cfgpath', type=str,
        default=os.path.expanduser('~')+'/pyrad/config/processing/',
        help='configuration file path')

    parser.add_argument(
        '--info_file', type=str,
        default='/store/msrad/radar/thundertracking/info/thundertracking_info.csv',
        help='File containing the dates to process')

    parser.add_argument(
        '--years', type=str,
        default=None,
        help='Years to process. If None all years in file will be processed')

    parser.add_argument(
        '--max_rank', type=float,
        default=None,
        help='Max rank to process')

    args = parser.parse_args()

    print("====== PYRAD thundertracking data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD thundertracking data processing finished: ")

    cfgfile_list = args.cfgfile_list.split(',')

    # get the years to process
    years = None
    if args.years is not None:
        years = list(map(int, args.years.split(',')))

    _, max_rank, _, trt_time_start, trt_time_end = read_thundertracking_info(
        args.info_file)

    for rank, time_start, time_end in zip(max_rank, trt_time_start, trt_time_end):
        # Do not process data if not in years
        if years is not None:
            if time_start.year not in years:
                continue

        if args.max_rank is not None:
            if rank > args.max_rank:
                continue

        for cfgfile in cfgfile_list:
            try:
                pyrad_main(
                    args.cfgpath+cfgfile, starttime=time_start,
                    endtime=time_end+datetime.timedelta(minutes=6))
            except ValueError:
                print(ValueError)


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
