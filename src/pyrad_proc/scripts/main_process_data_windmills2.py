#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_data_windmills2
================================================

This program gets time series of median values (obtained from histograms)
of windmill radar returns that occurred when the windmill had similar
characteristics of orientation respect to the radar or rotor speed, etc.

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import glob
import csv
from warnings import warn
import os

import numpy as np

from pyrad.io import read_histogram, get_fieldname_pyart
from pyrad.io import read_proc_periods
from pyrad.graph import plot_timeseries, get_field_name

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
        '--startdate', type=str, default=None,
        help='starting date of the data to be processed. Format ''YYYYMMDDhhmmss'' ')
    parser.add_argument(
        '--enddate', type=str, default=None,
        help='end date of the data to be processed. Format ''YYYYMMDDhhmmss'' ')

    parser.add_argument(
        '--radarbase', type=str,
        default='/store/msrad/radar/pyrad_products/mals_sha_windmills_rhi/',
        help='name of folder containing the radar data')

    parser.add_argument(
        '--windmill', type=str, default='WM1', help='Windmill ID')

    parser.add_argument(
        '--datatypes', type=str,
        default='dBuZ',
        help='Name of the polarimetric moments to process. Coma separated')


    parser.add_argument("-p", "--periodfile", type=str, default='',
                        help="Periods file")

    parser.add_argument("-i", "--infostr", type=str,
                        help="Information string about the actual data "
                        "processing (e.g. 'RUN57'). This string is added "
                        "to the filenames of the product files.",
                        default="")



    args = parser.parse_args()

    print("====== PYRAD windmill data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD windmill data processing finished: ")


    datatypes = args.datatypes.split(',')

    if args.infostr == "":
        infostr = ''
    else:
        infostr = '_'+args.infostr

    proc_startdates = None
    proc_enddates = None
    if args.startdate is not None and args.enddate is not None:
        proc_startdate = datetime.datetime.strptime(
            args.startdate, '%Y%m%d%H%M%S')
        proc_enddate = datetime.datetime.strptime(
            args.enddate, '%Y%m%d%H%M%S')
        ndays = (proc_enddate - proc_startdate).days + 1
        print('Number of days to process: '+str(ndays)+'\n\n')

        proc_startdates = [proc_startdate]
        proc_enddates = [proc_enddate]
    elif args.periodfile != '':
        proc_startdates, proc_enddates = read_proc_periods(args.periodfile)

    if proc_startdates is None or proc_enddates is None:
        warn('Time periods to process not defined')

    for datatype in datatypes:
        first_hist = True
        hist = None
        bin_edges = None

        flist = get_file_list(
            proc_startdates, proc_enddates, args.radarbase, args.windmill,
            datatype)

        if not flist:
            warn('No histogram files found')
            continue

        median_vals = np.ma.masked_all(np.size(flist))
        record_times = np.ma.masked_all(np.size(flist), dtype=datetime.datetime)
        for i, fname in enumerate(flist):
            datetimestr = os.path.basename(fname)[0:14]
            record_times[i] = datetime.datetime.strptime(
                datetimestr, '%Y%m%d%H%M%S')
            hist_aux, bin_edges_aux = read_histogram(fname)

            bin_centers = bin_edges_aux[:-1] + (bin_edges_aux[1:]-bin_edges_aux[:-1])/2.
            cum_sum = np.cumsum(hist_aux)
            median_val = bin_centers[cum_sum >= cum_sum[-1]/2.][0]
            if first_hist:
                hist = hist_aux
                bin_edges = bin_edges_aux
                first_hist = False
                median_vals[i] = median_val
            else:
                if not np.allclose(bin_edges, bin_edges_aux):
                    warn('Bin edges not identical')
                    continue
                hist += hist_aux
                median_vals[i] = median_val
        if hist is None:
            continue

        fname = (
            args.radarbase+proc_startdates[0].strftime("%Y%m%d-") +
            proc_enddates[-1].strftime("%Y%m%d")+'_'+args.windmill+'_median_ts_' +
            datatype+infostr+'.png')

        field_name = get_fieldname_pyart(datatype)
        field_dict = get_metadata(field_name)
        titl = (
            args.windmill+' '+proc_startdates[0].strftime("%Y%m%d-") +
            proc_enddates[-1].strftime("%Y%m%d")+infostr.replace('_', ' ')+'\n' +
            'median '+get_field_name(field_dict, field_name))

        plot_timeseries(
            record_times, [median_vals], [fname], labelx='Time UTC',
            labely='median '+datatype, title=titl)

        print("----- plot to '%s'" % fname)


        fname = (
            args.radarbase+proc_startdates[0].strftime("%Y%m%d-") +
            proc_enddates[-1].strftime("%Y%m%d")+'_'+args.windmill+'_median_ts_' +
            datatype+infostr+'.csv')
        write_ts(record_times, median_vals, fname)

        print("----- written to '%s'" % fname)



def get_file_list(starttimes, endtimes, radarbase, windmill, datatype):
    filelist = []
    for starttime, endtime in zip(starttimes, endtimes):
        startdate = starttime.replace(
            hour=0, minute=0, second=0, microsecond=0)
        enddate = endtime.replace(hour=0, minute=0, second=0, microsecond=0)
        ndays = int((enddate-startdate).days)+1
        t_filelist = []
        for i in range(ndays):
            daydir = (
                startdate+datetime.timedelta(days=i)).strftime('%Y-%m-%d')
            datapath = (
                radarbase+daydir+'/'+windmill+'_postproc/HISTOGRAM_' +
                datatype+'/')
            if not os.path.isdir(datapath):
                continue
            dayfilelist = glob.glob(datapath+'*_histogram_*'+datatype+'.csv')
            for filename in dayfilelist:
                t_filelist.append(filename)

        for filename in t_filelist:
            filenamestr = str(filename)
            datetimestr = os.path.basename(filenamestr)[0:14]
            fdatetime = datetime.datetime.strptime(
                datetimestr, '%Y%m%d%H%M%S')
            if (fdatetime >= starttime) and (fdatetime <= endtime):
                filelist.append(filenamestr)

    return sorted(filelist)


def write_ts(record_times, vals, fname):
    """
    writes time series of data

    Parameters
    ----------
    dataset : dict
        dictionary containing the time series parameters

    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'a', newline='') as csvfile:
        fieldnames = ['date', 'value']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()

        for val, dt in zip(vals, record_times):
            writer.writerow(
                {'date': dt.strftime('%Y%m%d%H%M%S'),
                'value': val})
        csvfile.close()

    return fname


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
