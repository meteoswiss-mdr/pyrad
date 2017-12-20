#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Create PYRAD products using a plane trajectory.
"""

# Author: Andreas Leuenberger
# License: BSD 3 clause

# ---------------------------------------------------------
# Imports
# ---------------------------------------------------------
from __future__ import print_function
import sys
import argparse
from datetime import datetime
import atexit

from pyrad.flow import main as pyrad_main


def main():
    """
    """

    # Parse the arguments
    parser = argparse.ArgumentParser(
        description="Create PYRAD products using a plane trajectory",
        epilog="Example:\n"
        "  process_trajectory.py -c $HOME/pyrad/config/processing/mals_emm_rw22_traj.txt\n"
        "     --preproc_cfgfile $HOME/pyrad/config/processing/mals_emm_rw22_traj_preproc.txt\n"
        "     -i TS011 /data/mals_plane_traj/EMM/gnv_20161026_ts011_seat_emmen_flt01_ADS.txt",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Input arguments:
    parser.add_argument("-c", "--cfgfile", type=str,
                        help="Main configuration file. Defines the ",
                        default="")

    parser.add_argument(
        '--preproc_cfgfile', type=str, default=None,
        help='name of main pre-processing configuration file')

    parser.add_argument("trajfile", type=str,
                        help="Definition file of plane trajectory. "
                        "Configuration of scan sector, products, ...")

    parser.add_argument("starttime", nargs='?', type=str,
                        help="Starting time of the data to be processed. "
                        "Format: YYYYMMDDhhmm[ss]. If not given, the time "
                        "of the first sample is used.",
                        default="")
    parser.add_argument('endtime', nargs='?', type=str,
                        help="End time of the data to be processed. "
                        "Format: YYYYMMDDhhmm[ss]. If not given, the time "
                        "of the last sample is used.",
                        default="")

    parser.add_argument("-i", "--infostr",
                        help="Information string about the actual data "
                        "processing (e.g. 'RUN57'). This string is added "
                        "to the filenames of the product files.",
                        default="")

    args = parser.parse_args()
    dt_starttime = None
    dt_endtime = None

    if (len(args.starttime) > 0):
        try:
            if (len(args.starttime) == 14):
                dt_starttime = datetime.strptime(args.starttime,
                                                 '%Y%m%d%H%M%S')
            elif (len(args.starttime) == 12):
                dt_starttime = datetime.strptime(args.starttime, '%Y%m%d%H%M')
            else:
                raise
        except:
            print("process_trajectory.py: Format error: Argument 'starttime' "
                  "must be in format 'YYYYMMDDhhmm[ss]' (is '%s')" %
                  args.starttime,
                  file=sys.stderr)
            sys.exit(1)

    if (len(args.endtime) > 0):
        try:
            if (len(args.endtime) == 14):
                dt_endtime = datetime.strptime(args.endtime, '%Y%m%d%H%M%S')
            elif (len(args.endtime) == 12):
                dt_endtime = datetime.strptime(args.endtime, '%Y%m%d%H%M')
            else:
                raise
        except:
            print("process_trajectory.py: Format error: Argument 'endtime' "
                  "must be in format 'YYYYMMDDhhmm[ss]' (is '%s')" %
                  args.endtime,
                  file=sys.stderr)
            sys.exit(1)

    print("====== PYRAD trajectory processing started: %s" %
          datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))

    atexit.register(_print_end_msg,
                    "====== PYRAD trajectory processing finished: ")

    # try:  #XXX
    #if args.preproc_cfgfile is not None:
    #    pyrad_main(args.preproc_cfgfile, dt_starttime, dt_endtime,
    #               trajfile=args.trajfile, infostr=args.infostr)

    pyrad_main(args.cfgfile, dt_starttime, dt_endtime,
               trajfile=args.trajfile, infostr=args.infostr)
    # except Exception as ee:
    #    print(str(ee), file=sys.stderr)


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

# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
