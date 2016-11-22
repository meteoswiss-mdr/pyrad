#!/usr/bin/env python

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program processes and post-processes data over a time span

To run the processing framework type:
    python main_process_data.py \
[config_file] [process_start_time] [process_end_time] \
--postproc_cfgfile [postproc_config_file] --cfgpath [cfgpath]

postproc_cfgfile is an optional argument with default: None
cfgpath is an optional argument with default: \
'/home/lom/users/fvj/pyrad/config/processing/'

Example:
    python main_process_data.py 'paradiso_fvj_vol.txt' '20140523000000' \
'20140523001000' --postproc_cfgfile 'paradiso_fvj_vol_postproc.txt' \
--cfgpath '/home/lom/users/fvj/pyrad/config/processing/'

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit

from pyrad.flow import main

print(__doc__)


if __name__ == '__main__':
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments
    parser.add_argument(
        'proc_cfgfile', type=str, help='name of main configuration file')
    parser.add_argument(
        'starttime', type=str,
        help='starting time of the data to be processed')
    parser.add_argument(
        'endtime', type=str, help='end time of the data to be processed ')

    # keyword arguments
    parser.add_argument(
        '--postproc_cfgfile', type=str, default=None,
        help='name of main post-processing configuration file')
    parser.add_argument(
        '--cfgpath', type=str,
        default='/home/lom/users/fvj/pyrad/config/processing/',
        help='configuration file path')

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

    proc_starttime = datetime.datetime.strptime(
        args.starttime, '%Y%m%d%H%M%S')
    proc_endtime = datetime.datetime.strptime(args.endtime, '%Y%m%d%H%M%S')
    cfgfile_proc = args.cfgpath+args.proc_cfgfile

    main(cfgfile_proc, proc_starttime, proc_endtime)
    if args.postproc_cfgfile is not None:
        cfgfile_postproc = args.cfgpath+args.postproc_cfgfile
        main(cfgfile_postproc, proc_starttime, proc_endtime)

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
    print(text + datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))

