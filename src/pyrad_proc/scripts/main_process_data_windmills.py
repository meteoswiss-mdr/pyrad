#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_data_windmills
================================================

This program groups histograms of windmill radar returns that occurred when
the windmill had similar characteristics of orientation respect to the radar
or rotor speed

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import glob
from warnings import warn
import os

import numpy as np

from pyrad.io import read_histogram, write_histogram, get_fieldname_pyart
from pyrad.io import read_proc_periods
from pyrad.graph import plot_histogram2, get_colobar_label, get_field_name

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
        default='/store/msrad/radar/pyrad_products/mals_sha_windmills_hr/',
        help='name of folder containing the radar data')

    parser.add_argument(
        '--windmill', type=str, default='WM1', help='Windmill ID')

    parser.add_argument(
        '--datatypes', type=str,
        default='dBuZ,rcs_h,ZDRu,RhoHVu,uPhiDPu,Vu,Wu',
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

        for fname in flist:
            hist_aux, bin_edges_aux = read_histogram(fname)
            if first_hist:
                hist = hist_aux
                bin_edges = bin_edges_aux
                first_hist = False
            else:
                if not np.allclose(bin_edges, bin_edges_aux):
                    warn('Bin edges not identical')
                    continue
                hist += hist_aux

        if hist is None:
            continue

        fname = (
            args.radarbase+proc_startdates[0].strftime("%Y%m%d-") +
            proc_enddates[-1].strftime("%Y%m%d")+'_'+args.windmill+'_histogram_' +
            datatype+infostr+'.csv')
        write_histogram(bin_edges, hist, fname, datatype=datatype)

        bin_centers = bin_edges[:-1] + (bin_edges[1:]-bin_edges[:-1])/2.
        width = bin_edges[1:]-bin_edges[:-1]

        field_name = get_fieldname_pyart(datatype)
        field_dict = get_metadata(field_name)
        xlabel = get_colobar_label(field_dict, field_name)
        titl = (
            args.windmill+' '+proc_startdates[0].strftime("%Y%m%d-") +
            proc_enddates[-1].strftime("%Y%m%d")+infostr.replace('_', ' ')+'\n' +
            get_field_name(field_dict, field_name))
        fname = (
            args.radarbase+proc_startdates[0].strftime("%Y%m%d-") +
            proc_enddates[-1].strftime("%Y%m%d")+'_'+args.windmill+'_histogram_' +
            datatype+infostr+'.png')
        plot_histogram2(
            bin_centers, hist, [fname], width=width, labelx=xlabel, titl=titl)

        print("----- plot to '%s'" % fname)



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
