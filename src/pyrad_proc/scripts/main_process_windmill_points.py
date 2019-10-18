#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_data_windmills
================================================

This program compiles histograms of windmill radar returns

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
        '--period_file_base', type=str,
        default='/users/jfigui/windmills_params/',
        help='name of folder containing the windmill data')

    parser.add_argument(
        '--database', type=str,
        default='/store/msrad/radar/pyrad_products/mals_sha_windmills_rhi/',
        help='base path to the radar data')

    parser.add_argument(
        '--wind_IDs', type=str, default='WM1,WM2,WM3',
        help='Windmill ID. Coma separated')

    parser.add_argument(
        '--points', type=str, default='Pmax',
        help='Point ID. Coma separated')

    parser.add_argument(
        '--datatypes', type=str, default='dBuZ,ZDRu,RhoHVu,uPhiDPu,Vu,Wu',
        help='Data types. Coma separated')

    parser.add_argument(
        '--steps', type=str, default='0.5,0.1,0.01,1,0.2,0.1',
        help='Histogram steps. Coma separated')

    parser.add_argument(
        '--orientations', type=str, default='0,45,90,135,180,225,270,315',
        help='Orientation respect to radar')

    parser.add_argument(
        '--span', type=float, default=45.,
        help='Span')

    parser.add_argument(
        '--vel_limit', type=float, default=0.,
        help='Velocity limit')

    parser.add_argument(
        '--dBuZ_threshold', type=float, default=None,
        help='dBuZ threshold')

    args = parser.parse_args()

    print("====== PYRAD windmill data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD windmill data processing finished: ")

    wind_IDs = args.wind_IDs.split(',')
    orientations = np.asarray(args.orientations.split(','), dtype=float)
    points = args.points.split(',')
    datatypes = args.datatypes.split(',')
    steps = np.asarray(args.steps.split(','), dtype=float)

    speeds = ['speed_GT'+str(args.vel_limit), 'speed_LE'+str(args.vel_limit)]
    if args.dBuZ_threshold is not None:
        wet_dry = [
            'dBuZ_LT'+str(args.dBuZ_threshold),
            'dBuZ_GE'+str(args.dBuZ_threshold)]

    if '_rhi/' in args.database:
        scan_type = 'rhi'
    else:
        scan_type = 'ppi'

    # Read periods of processing
    for ori in orientations:
        for speed in speeds:
            if args.dBuZ_threshold is not None:
                for wd in wet_dry:
                    for wind_ID in wind_IDs:
                        if wind_ID == 'WM1':
                            wind_ID2 = 'nx85215'
                        elif wind_ID == 'WM2':
                            wind_ID2 = 'nx85213'
                        elif wind_ID == 'WM3':
                            wind_ID2 = 'nx85214'
                        period_file = (
                            args.period_file_base+scan_type+'_'+wind_ID2+'_span' +
                            str(args.span)+'_ori'+str(ori)+'_'+speed+'_'+wd +
                            '.csv')

                        start_times, end_times = read_proc_periods(period_file)
                        if start_times is None:
                            continue
                        for point in points:
                            for step, datatype in zip(steps, datatypes):
                                # Read data time series files
                                flist = glob.glob(
                                    args.database+wind_ID+'_'+datatype+'_' +
                                    point+'_TS/TS/ts_POINT_MEASUREMENT_' +
                                    datatype+'*.csv')
                                if not flist:
                                    continue

                                dt_ts, values = read_timeseries(flist[0])

                                basepath = os.path.dirname(flist[0])+'/'

                                # keep only data within time periods
                                val_filt = np.ma.array([])
                                for t_start, t_end in zip(start_times, end_times):
                                    for dt, val in zip(dt_ts, values):
                                        if dt >= t_start and dt <= t_end:
                                            val_filt = np.ma.append(
                                                val_filt, val)

                                if val_filt.size == 0:
                                    warn('No data for point '+point +
                                         ' for time periods in file ' +
                                         period_file)
                                    continue

                                # Histogram plots
                                field_name = get_fieldname_pyart(datatype)
                                field_dict = get_metadata(field_name)

                                fname = (
                                    basepath+datatype+'_'+point+'_span' +
                                    str(args.span)+'_ori'+str(ori)+'_'+speed +
                                    '_'+wd+'_hist.png')
                                titl = (
                                    datatype+' '+point+' span ' +
                                    str(args.span)+' orientation '+str(ori) +
                                    ' '+speed+' '+wd)

                                bin_edges, vals = compute_histogram(
                                    val_filt, field_name, step=step)
                                fname = plot_histogram(
                                    bin_edges, vals, [fname],
                                    labelx=get_colobar_label(
                                        field_dict, field_name),
                                    titl=titl)
                                print('Plotted '+' '.join(fname))

                                fname = (
                                    basepath+datatype+'_'+point+'_span' +
                                    str(args.span)+'_ori'+str(ori)+'_'+speed +
                                    '_'+wd+'_hist.csv')
                                hist, _ = np.histogram(vals, bins=bin_edges)
                                fname = write_histogram(
                                    bin_edges, hist, fname)
                                print('Written '+fname)
            else:
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
                    for point in points:
                        for step, datatype in zip(steps, datatypes):
                            # Read data time series files
                            flist = glob.glob(
                                args.database+wind_ID+'_'+datatype+'_'+point +
                                '_TS/TS/ts_POINT_MEASUREMENT_'+datatype +
                                '*.csv')
                            if not flist:
                                continue

                            dt_ts, values = read_timeseries(flist[0])

                            basepath = os.path.dirname(flist[0])+'/'

                            # keep only data within time periods
                            val_filt = np.ma.array([])
                            for t_start, t_end in zip(start_times, end_times):
                                for dt, val in zip(dt_ts, values):
                                    if dt >= t_start and dt <= t_end:
                                        val_filt = np.ma.append(val_filt, val)

                            if val_filt.size == 0:
                                warn('No data for point '+point +
                                    ' for time periods in file '+period_file)
                                continue

                            # Histogram plots
                            field_name = get_fieldname_pyart(datatype)
                            field_dict = get_metadata(field_name)

                            fname = (
                                basepath+datatype+'_'+point+'_span' +
                                str(args.span)+'_ori'+str(ori)+'_'+speed +
                                '_hist.png')
                            titl = (
                                datatype+' '+point+' span '+str(args.span) +
                                ' orientation '+str(ori)+' '+speed)

                            bin_edges, vals = compute_histogram(
                                val_filt, field_name, step=step)
                            fname = plot_histogram(
                                bin_edges, vals, [fname],
                                labelx=get_colobar_label(
                                    field_dict, field_name),
                                titl=titl)
                            print('Plotted '+' '.join(fname))

                            fname = (
                                basepath+datatype+'_'+point+'_span' +
                                str(args.span)+'_ori'+str(ori)+'_'+speed +
                                '_hist.csv')
                            hist, _ = np.histogram(vals, bins=bin_edges)
                            fname = write_histogram(bin_edges, hist, fname)
                            print('Written '+fname)


    # All data
    for wind_ID in wind_IDs:
        for point in points:
            for step, datatype in zip(steps, datatypes):
                # Read data time series files
                flist = glob.glob(
                    args.database+wind_ID+'_'+datatype+'_'+point +
                    '_TS/TS/ts_POINT_MEASUREMENT_'+datatype+'*.csv')
                if not flist:
                    continue

                dt_ts, values = read_timeseries(flist[0])

                basepath = os.path.dirname(flist[0])+'/'

                # Histogram plots
                field_name = get_fieldname_pyart(datatype)
                field_dict = get_metadata(field_name)

                fname = basepath+datatype+'_'+point+'_all_data_hist.png'
                titl = datatype+' '+point+' all data '

                bin_edges, vals = compute_histogram(
                    values, field_name, step=step)
                fname = plot_histogram(
                    bin_edges, vals, [fname],
                    labelx=get_colobar_label(field_dict, field_name),
                    titl=titl)
                print('Plotted '+' '.join(fname))

                fname = basepath+datatype+'_'+point+'_all_data_hist.csv'
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
