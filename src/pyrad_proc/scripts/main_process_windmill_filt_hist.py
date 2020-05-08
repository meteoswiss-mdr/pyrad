#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_windmill_hist
================================================

This program creates histograms out of time series of windmill radar returns
from a particular point in space filtered according to the windmill
characteristics (orientation, rotor speed, etc.)

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os
import glob
from warnings import warn

import numpy as np

from pyrad.io import read_timeseries, read_proc_periods, get_fieldname_pyart
from pyrad.io import write_histogram
from pyrad.graph import plot_histogram, get_colobar_label
from pyrad.util import compute_histogram

from pyart.config import get_metadata

print(__doc__)


def main():
    """
    """

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # keyword arguments
    parser.add_argument(
        '--database', type=str,
        default='/store/msrad/radar/pyrad_products/mals_sha_windmills_point_WM1_20200321-20200325/',
        help='base path to the radar data')

    parser.add_argument(
        '--datatypes', type=str,
        default='dBuZ,dBuZv,rcs_h,rcs_v,uPhiDPu,RhoHVu,ZDRu,Vu,Wu',
        help='Data types. Coma separated')

    # RhoHVu: 0.003937 starting (0 to 1)
    # ZDRu:   0.078740 starting (-8 to 12)
    # Vu:     0.125    starting (-15.8 to 15.8)
    # Wu:     0-125    starting (0 to 15.8)
    parser.add_argument(
        '--steps', type=str,
        default='0.5,0.5,1.,0.5,0.5,0.003937,0.078740,0.125,0.125',
        help='Histogram steps. Coma separated')

    parser.add_argument(
        '--period_file_base', type=str,
        default='/users/jfigui/windmills_params/2020_Schaffhausen/rawdata/',
        help='name of folder containing the windmill data')

    parser.add_argument(
        '--wind_IDs', type=str, default='WM1',
        help='Windmill ID. Coma separated')

    parser.add_argument(
        '--orientations', type=str,
        default='0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350',
        help='Orientation respect to radar')

    parser.add_argument(
        '--span', type=float, default=10.,
        help='Span')

    parser.add_argument(
        '--vel_limit', type=float, default=0.,
        help='Velocity limit')

    args = parser.parse_args()

    print("====== PYRAD windmill data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD windmill data processing finished: ")

    datatypes = args.datatypes.split(',')
    steps = np.asarray(args.steps.split(','), dtype=float)

    wind_IDs = args.wind_IDs.split(',')
    orientations = np.asarray(args.orientations.split(','), dtype=float)
    speeds = ['speed_GT'+str(args.vel_limit), 'speed_LE'+str(args.vel_limit)]

    scan_type = 'ppi'

    # Read data time series
    for step, datatype in zip(steps, datatypes):
        # Read data time series files
        flist = glob.glob(
            args.database+datatype+'_TS/TS/ts_POINT_MEASUREMENT_'+datatype +
            '_*.csv')
        if not flist:
            continue

        dt_ts, values = read_timeseries(flist[0])
#        print('number of values read '+str(values.size))
#
#        # filter data extending to next day
#        val_ind = []
#        for ind, dt_ts_el in enumerate(dt_ts):
#            if dt_ts_el.day != 12:
#                val_ind.append(ind)
#        values = values[val_ind]
#        print('number of values kept '+str(values.size))

         # Read periods of processing
        for ori in orientations:
            for speed in speeds:
                for wind_ID in wind_IDs:
                    if wind_ID == 'WM1':
                        wind_ID2 = 'nx85215'
                    elif wind_ID == 'WM2':
                        wind_ID2 = 'nx85213'
                    elif wind_ID == 'WM3':
                        wind_ID2 = 'nx85214'
                    period_file = (
                        args.period_file_base+scan_type+'_'+wind_ID2+'_span' +
                        str(args.span)+'_ori'+str(ori)+'_'+speed+'.csv')

                    start_times, end_times = read_proc_periods(period_file)
                    if start_times is None:
                        continue

                    # keep only data within time periods
                    val_filt = []
                    for t_start, t_end in zip(start_times, end_times):
                        for dt, val in zip(dt_ts, values):
                            if dt >= t_start and dt <= t_end:
                                val_filt.append(val)
                    val_filt = np.ma.array(val_filt)

                    if val_filt.size == 0:
                        warn('No data for time periods in file '+period_file)
                        continue

                    # Histogram plots
                    field_name = get_fieldname_pyart(datatype)
                    field_dict = get_metadata(field_name)

                    fname = (
                        args.database+datatype+'_TS/TS/' +
                        datatype+'_span'+str(args.span)+'_ori'+str(ori)+'_' +
                        speed+'_hist.png')
                    titl = (
                        datatype+' span '+str(args.span) +
                        ' orientation '+str(ori)+' '+speed)

                    bin_edges, vals = compute_histogram(
                        val_filt, field_name, step=step)
                    fname = plot_histogram(
                        bin_edges, vals, [fname],
                        labelx=get_colobar_label(field_dict, field_name),
                        titl=titl)
                    print('Plotted '+' '.join(fname))

                    fname = (
                        args.database+datatype+'_TS/TS/' +
                        datatype+'_span'+str(args.span)+'_ori'+str(ori)+'_' +
                        speed+'_hist.csv')
                    hist, _ = np.histogram(vals, bins=bin_edges)
                    fname = write_histogram(bin_edges, hist, fname)
                    print('Written '+fname)


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
