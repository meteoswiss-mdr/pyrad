#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
format_conversion_intercomp
================================================

This program converts radar intercomparison time series files in the old idl
output format into the new Pyrad format and plots the result

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import numpy as np

from pyrad.io import read_intercomp_scores_ts_old, write_intercomp_scores_ts
from pyrad.graph import plot_intercomp_scores_ts
from pyrad.io import get_fieldname_pyart


print(__doc__)


def main():
    """
    """

    input_path = (
        '/home/lom/users/fvj/tmp/')
    output_path = (
        '/home/lom/users/fvj/tmp/')
    rad1_vec = [
        'A', 'A', 'A', 'A', 'D', 'D', 'D', 'D', 'L', 'L', 'L', 'L', 'P', 'P',
        'P', 'P', 'W', 'W', 'W', 'W']
    rad2_vec = [
        'D', 'L', 'P', 'W', 'A', 'L', 'P', 'W', 'A', 'D', 'P', 'W', 'A', 'D',
        'L', 'W', 'A', 'D', 'L', 'P']
    var_input_vec = ['Zh', 'Zv']
    var_output_vec = ['dBZ', 'dBZv']
    np_min = 200.
    corr_min = 0.8
    year_vec = [
        datetime.datetime(2013, 1, 1), datetime.datetime(2014, 1, 1),
        datetime.datetime(2015, 1, 1), datetime.datetime(2016, 1, 1),
        datetime.datetime(2017, 1, 1)]

    print("====== Intercomparison format conversion started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Intercomparison format conversion finished: ")

    for i, rad1 in enumerate(rad1_vec):
        rad2 = rad2_vec[i]

        print('Processing Radars '+rad1+'-'+rad2)
        for j, var_input in enumerate(var_input_vec):
            var_output = var_output_vec[j]
            print('- Processing Variable '+var_output)
            for k, year in enumerate(year_vec):
                print('-- Processing Year '+year.strftime('%Y'))
                fname_input = (
                    input_path+'intercomp'+var_input+'PL'+rad1+'-'+rad2 +
                    year.strftime('%y')+'.all.txt')
                fname_output = (
                    output_path+year.strftime('%Y') +
                    '_ts_INTERCOMP_TIME_AVG_'+var_output+'_'+rad1+'-'+rad2 +
                    '.csv')
                figfname = [
                    output_path+year.strftime('%Y') +
                    '_ts_INTERCOMP_TIME_AVG_'+var_output+'_'+rad1+'-'+rad2 +
                    '.png']

                (date_vec, np_vec, meanbias_vec, medianbias_vec,
                 quant25bias_vec, quant75bias_vec, modebias_vec, corr_vec,
                 slope_vec, intercep_vec, intercep_slope1_vec) = (
                    read_intercomp_scores_ts_old(fname_input))

                if date_vec is None:
                    continue

                stats = {
                    'npoints': np_vec,
                    'meanbias': meanbias_vec,
                    'medianbias': medianbias_vec,
                    'quant25bias': quant25bias_vec,
                    'quant75bias': quant75bias_vec,
                    'modebias': modebias_vec,
                    'corr': corr_vec,
                    'slope': slope_vec,
                    'intercep': intercep_vec,
                    'intercep_slope_1': intercep_slope1_vec
                }
                field_name = get_fieldname_pyart(var_output)
                write_intercomp_scores_ts(
                    date_vec, stats, field_name, fname_output,
                    rad1_name=rad1, rad2_name=rad2)

                titldate = (date_vec[0].strftime('%Y%m%d')+'-' +
                            date_vec[-1].strftime('%Y%m%d'))
                titl = (
                    rad1+'-'+rad2+' '+field_name+' intercomparison '+titldate)

                plot_intercomp_scores_ts(
                    date_vec, np_vec, meanbias_vec, medianbias_vec,
                    quant25bias_vec, quant75bias_vec, modebias_vec, corr_vec,
                    slope_vec, intercep_vec, intercep_slope1_vec, figfname,
                    ref_value=0., np_min=np_min, corr_min=corr_min,
                    labelx='Time UTC', titl=titl)


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
