#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_poh_hist
================================================

This program processes POH histogram data

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import glob
import os

import pandas as pd

from pyrad.io import read_histogram
from pyart.correct import sun_position_pysolar

print(__doc__)


def main():
    """
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments


    # keyword arguments
    parser.add_argument(
        '--starttime', type=str, default=None,
        help=('starting time of the data to be processed. ' +
              'Format ''YYYYMMDD'''))

    parser.add_argument(
        '--ndays', type=int, default=1,
        help='Number of days to process')

    parser.add_argument(
        '--basepath', type=str,
        default='/store/msrad/radar/pyrad_products/rad4alp_POH/',
        help='Base path of the data')

    args = parser.parse_args()

    print("====== POH histogram processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== POH histogram processing finished: ")

    starttime = datetime.datetime.strptime(
        args.starttime, '%Y%m%d')

    datetime_list = []
    POH90_list = []
    daynight_list = []
    for iday in range(args.ndays):
        time_dir = (starttime+datetime.timedelta(days=iday)).strftime(
            "%Y-%m-%d")
        data_input_path = args.basepath+time_dir+'/POH/HISTOGRAM/'

        flist = glob.glob(data_input_path+'*_POH_histogram_RAW_GRID_POH.csv')

        for fname in flist:
            print('Reading POH histogram file '+fname)
            hist, _ = read_histogram(fname)
            if hist[9] == 0:
                print('No POH above 90% in file')
                continue

            datetimestr = os.path.basename(fname)[0:14]
            if datetimestr[11] not in ('0', '5'):
                print('Time step has to be multiple of 5')
                continue
                
            dt = datetime.datetime.strptime(datetimestr, '%Y%m%d%H%M%S')

            # check if image in day light, night or transition
            el_sun_southwest, _ = sun_position_pysolar(dt, 45.5, 5.5)
            el_sun_northwest, _ = sun_position_pysolar(dt, 48., 5.5)
            el_sun_southeast, _ = sun_position_pysolar(dt, 45.5, 11.)
            el_sun_northeast, _ = sun_position_pysolar(dt, 48., 11.)

            # check sun position
            if (el_sun_southwest<0. and el_sun_northwest<0.
                    and el_sun_southeast<0. and el_sun_southeast<0.):
                # all points with elevation angle below 0°: night time
                daynight_list.append('night')
            elif (el_sun_southwest>10. and el_sun_northwest>10.
                    and el_sun_southeast>10. and el_sun_southeast>10.):
                # all points with elevation angle above 10°: day time
                daynight_list.append('day')
            else:
                daynight_list.append('transition')


            datetime_list.append(dt)
            POH90_list.append(hist[9])

    df = pd.DataFrame(data={
        'nominal time': datetime_list,
        'day/night': daynight_list,
        'npixels with POH > 90%': POH90_list})

    df.sort_values(by='nominal time', ignore_index=True, inplace=True)
    df.to_csv(args.basepath+'POH_above_90_percent.csv', index=True)


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
