#!/usr/bin/env python

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program processes data over a time span

To run the processing framework type:
    python main_process_data.py \
[config_file] [process_start_time] [process_end_time]
Example:
    python main_process_data.py \
'/home/lom/users/fvj/pyrad/config/processing/paradiso_fvj_vol.txt' \
'20140523000000' '20140523001000'

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse

from pyrad.flow import main

print(__doc__)


if __name__ == '__main__':
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments
    parser.add_argument(
        'proccfgfile', type=str, help='path to main configuration file')
    parser.add_argument(
        'starttime', type=str,
        help='starting time of the data to be processed')
    parser.add_argument(
        'endtime', type=str, help='end time of the data to be processed ')

    args = parser.parse_args()

    print('config file: '+args.proccfgfile)
    print('start time: '+args.starttime)
    print('end time: '+args.endtime)

    proc_starttime = datetime.datetime.strptime(
        args.starttime, '%Y%m%d%H%M%S')
    proc_endtime = datetime.datetime.strptime(args.endtime, '%Y%m%d%H%M%S')

    main(args.proccfgfile, proc_starttime, proc_endtime)
