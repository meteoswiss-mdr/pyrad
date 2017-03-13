#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program processes and post-processes data over a time span

To run the processing framework type:
    python main_process_data.py \
[config_file] --starttime [process_start_time] --endtime [process_end_time] \
--postproc_cfgfile [postproc_config_file] --cfgpath [cfgpath]

If startime and endtime are not specified the program does realtime processing
postproc_cfgfile is an optional argument with default: None
cfgpath is an optional argument with default: \
'$HOME/pyrad/config/processing/'

Example:
    python main_process_data.py 'paradiso_fvj_vol.txt' --starttime \
'20140523000000' --endtime '20140523001000' --postproc_cfgfile \
'paradiso_fvj_vol_postproc.txt' --cfgpath '$HOME/pyrad/config/processing/'

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os

from pyrad.flow.flow_control import main as pyrad_main

print(__doc__)


def main():
    """
    """

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments
    parser.add_argument(
        'proc_cfgfile', type=str, help='name of main configuration file')

    # keyword arguments
    parser.add_argument(
        '--starttime', type=str, default=None,
        help=('starting time of the data to be processed. ' +
              'Format ''YYYYMMDDhhmmss'''))
    parser.add_argument(
        '--endtime', type=str, default=None,
        help='end time of the data to be processed. Format ''YYYYMMDDhhmmss''')
    parser.add_argument(
        '--postproc_cfgfile', type=str, default=None,
        help='name of main post-processing configuration file')
    parser.add_argument(
        '--cfgpath', type=str,
        default=os.path.expanduser('~')+'/pyrad/config/processing/',
        help='configuration file path')
    parser.add_argument(
        '--rt', type=int, default=0,
        help='if 1 the processing is real time.')

    args = parser.parse_args()

    print("====== PYRAD data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD data processing finished: ")

    print('config path: '+args.cfgpath)
    print('config file: '+args.proc_cfgfile)
    print('postproc config file: '+str(args.postproc_cfgfile))
    print('start time: '+args.starttime)
    print('end time: '+args.endtime)

    proc_starttime = None
    if args.starttime is not None:
        proc_starttime = datetime.datetime.strptime(
            args.starttime, '%Y%m%d%H%M%S')
    proc_endtime = None
    if args.endtime is not None:
        proc_endtime = datetime.datetime.strptime(args.endtime, '%Y%m%d%H%M%S')
    cfgfile_proc = args.cfgpath+args.proc_cfgfile

    pyrad_main(cfgfile_proc, proc_starttime, proc_endtime, rt=args.rt)

    # if not real time allow post-processing
    if args.postproc_cfgfile is not None and args.rt == 0:
        cfgfile_postproc = args.cfgpath+args.postproc_cfgfile
        pyrad_main(cfgfile_postproc, proc_starttime, proc_endtime, rt=args.rt)


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
