#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
format_conversion_monitoring
================================================

This program converts monitoring time series files in the old idl output
format into the new Pyrad format and plots the result

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import numpy as np

from pyrad.io import read_monitoring_ts_old, write_monitoring_ts
from pyrad.graph import plot_monitoring_ts
from pyrad.io import generate_field_name_str, get_fieldname_pyart


print(__doc__)


def main():
    """
    """

    input_path = (
        '/home/lom/users/fvj/tmp/')
    output_path = (
        '/home/lom/users/fvj/tmp/')
    rad_vec = ['A', 'D', 'L', 'P', 'W']
    var_input_vec = ['phidp0', 'rhoav', 'zdrbias', 'zhbias']
    var_output_vec = ['PhiDP0', 'RhoHV_rain', 'ZDR_prec', 'dBZ_bias']
    quantiles_vec = np.asarray([
        [25., 50., 75.], [65., 80., 95.], [25., 50., 75.], [25., 50., 75.]])
    ref_value_vec = [0., 0.99, 0.2, 0.]
    vmin_vec = [-20., 0.95, -2., -30.]
    vmax_vec = [20., 1.01, 2., 30.]
    np_min_vec = [500., 5000., 5000., 500.]
    year_vec = [
        datetime.datetime(2013, 1, 1), datetime.datetime(2014, 1, 1),
        datetime.datetime(2015, 1, 1), datetime.datetime(2016, 1, 1),
        datetime.datetime(2017, 1, 1)]

    print("====== Monitoring format conversion started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Monitoring format conversion finished: ")

    for i, rad in enumerate(rad_vec):
        print('Processing Radar '+rad)
        for j, var_input in enumerate(var_input_vec):
            var_output = var_output_vec[j]
            quantiles = quantiles_vec[j, :]
            ref_value = ref_value_vec[j]
            vmin = vmin_vec[j]
            vmax = vmax_vec[j]
            np_min = np_min_vec[j]
            print('- Processing Variable '+var_output)
            for k, year in enumerate(year_vec):
                print('-- Processing Year '+year.strftime('%Y'))
                fname_input = (
                    input_path+'PL'+rad+year.strftime('%y')+'_monitoring_' +
                    var_input+'.csv')
                fname_output = (
                    output_path+year.strftime('%Y')+'_'+rad +
                    '_ts_MONITORING_'+var_output+'.csv')
                figfname = [
                    output_path+year.strftime('%Y')+'_'+rad +
                    '_ts_MONITORING_'+var_output+'.png']

                date, np_t_vec, cquant_vec, lquant_vec, hquant_vec = (
                    read_monitoring_ts_old(fname_input))

                if date is None:
                    continue

                val_vec = np.ma.asarray(
                    [lquant_vec, cquant_vec, hquant_vec]).T
                print(np.shape(val_vec))
                write_monitoring_ts(
                    date, np_t_vec, val_vec, quantiles, var_output,
                    fname_output)

                titldate = (date[0].strftime('%Y%m%d')+'-' +
                            date[-1].strftime('%Y%m%d'))
                titl = rad+' Monitoring '+titldate

                labely = generate_field_name_str(var_output)

                plot_monitoring_ts(
                    date, np_t_vec, cquant_vec, lquant_vec, hquant_vec,
                    get_fieldname_pyart(var_output), figfname,
                    ref_value=ref_value, vmin=vmin, vmax=vmax, np_min=np_min,
                    labelx='Time UTC', labely=labely, titl=titl)


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
