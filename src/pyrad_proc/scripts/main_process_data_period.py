#!/usr/bin/env python

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program does the daily processes over a period of time

To run the processing framework type:
    python main_process_data.py \
[config_file] [process_start_date] [process_end_date]
Example:
    python main_process_data.py \
'/home/lom/users/fvj/pyrad/config/processing/paradiso_fvj_vol.txt' \
'20140523' '20140525'

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

    parser.add_argument(
        'proccfgfile', type=str, help='path to main configuration file')
    parser.add_argument(
        'starttime', type=str,
        help='starting date of the data to be processed')
    parser.add_argument(
        'endtime', type=str, help='end date of the data to be processed ')

    args = parser.parse_args()

    print('config file: '+args.proccfgfile)
    print('start date: '+args.starttime)
    print('end date: '+args.endtime)

    proc_startdate = datetime.datetime.strptime(
        args.starttime, '%Y%m%d')
    proc_enddate = datetime.datetime.strptime(args.endtime, '%Y%m%d')

    ndays = (proc_enddate - proc_startdate).days + 1
    print('Number of days to process: '+str(ndays)+'\n\n')

    for day in range(ndays):
        proc_starttime = proc_startdate + datetime.timedelta(days=day)
        #proc_endtime = proc_starttime + datetime.timedelta(days=1)
        proc_endtime = proc_starttime + datetime.timedelta(minutes=10)
        try:
            main(args.proccfgfile, proc_starttime, proc_endtime)
        except ValueError:
            print(ValueError)
