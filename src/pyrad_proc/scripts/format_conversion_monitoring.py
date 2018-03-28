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
        '/store/msrad/radar/monitoring/polarimetry/monitoring_ts/')
    output_path = (
        '/store/msrad/radar/monitoring/polarimetry/')
    rad_vec = ['D']
    var_input_vec = ['phidp0', 'rhoav', 'zdrbias', 'zhbias']
    var_output_vec = ['PhiDP0', 'RhoHV_rain', 'ZDR_prec', 'dBZ_bias']
    year_vec = [
        datetime.datetime(2013, 1, 1), datetime.datetime(2014, 1, 1),
        datetime.datetime(2015, 1, 1), datetime.datetime(2016, 1, 1),
        datetime.datetime(2017, 1, 1)]

    plot_data = True

    print("====== Monitoring format conversion started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Monitoring format conversion finished: ")

    for i, rad in enumerate(rad_vec):
        print('Processing Radar '+rad)
        for j, var_input in enumerate(var_input_vec):
            var_output = var_output_vec[j]
            if var_output == 'dBZ':
                basedir = 'rad4alp_gc_PH'+rad
                dsdir = 'monitoring_clt_Zh'
                mon_type = 'GC_MONITORING'
                quantiles = [50., 95., 99.]
            elif var_output == 'dBZv':
                basedir = 'rad4alp_gc_PH'+rad
                dsdir = 'monitoring_clt_Zv'
                mon_type = 'GC_MONITORING'
                quantiles = [50., 95., 99.]
            elif var_output == 'RhoHV_rain':
                basedir = 'rad4alp_dataquality_PL'+rad
                dsdir = 'monitoring_RhoHV'
                mon_type = 'MONITORING'
                quantiles = [65., 80., 95.]
            elif var_output == 'PhiDP0':
                basedir = 'rad4alp_dataquality_PL'+rad
                dsdir = 'monitoring_PhiDP0'
                mon_type = 'MONITORING'
                quantiles = [25., 50., 75.]
            elif var_output == 'ZDR_prec':
                basedir = 'rad4alp_dataquality_PL'+rad
                dsdir = 'monitoring_ZDR'
                mon_type = 'MONITORING'
                quantiles = [25., 50., 75.]
            elif var_output == 'ZDR_snow':
                basedir = 'rad4alp_dataquality_PL'+rad
                dsdir = 'monitoring_ZDR_snow'
                mon_type = 'MONITORING'
                quantiles = [25., 50., 75.]
            elif var_output == 'dBZ_bias':
                basedir = 'rad4alp_dataquality_PL'+rad
                dsdir = 'monitoring_Zh_bias'
                mon_type = 'MONITORING'
                quantiles = [25., 50., 75.]

            print('- Processing Variable '+var_output)
            for k, year in enumerate(year_vec):
                print('-- Processing Year '+year.strftime('%Y'))
                fname_input = (
                    input_path+'PL'+rad+year.strftime('%y')+'_monitoring_' +
                    var_input+'.csv')
                fname_output = (
                    output_path+year.strftime('%Y')+'_'+rad +
                    '_ts_'+mon_type+'_'+var_output+'.csv')
                figfname = [
                    output_path+year.strftime('%Y')+'_'+rad +
                    '_ts_'+mon_type+'_'+var_output+'.png']

                date, np_t_vec, cquant_vec, lquant_vec, hquant_vec = (
                    read_monitoring_ts_old(fname_input))

                if date is None:
                    continue

                val_vec = np.ma.asarray(
                    [lquant_vec, cquant_vec, hquant_vec]).T
                print(np.shape(val_vec))
                fname = write_monitoring_ts(
                    date, np_t_vec, val_vec, quantiles, var_output,
                    fname_output)
                print('written file '+fname)

                if not plot_data:
                    continue

                titldate = (date[0].strftime('%Y%m%d')+'-' +
                            date[-1].strftime('%Y%m%d'))
                titl = rad+' Monitoring '+titldate

                labely = generate_field_name_str(var_output)

                if var_output == 'RhoHV_rain':
                    ref_value = 0.99
                    vmin = 0.95
                    vmax = 1.01
                    np_min = 5000
                elif var_output == 'PhiDP0':
                    ref_value = 0.
                    vmin = -20.
                    vmax = 20.
                    np_min = 5000
                elif var_output == 'ZDR_prec':
                    ref_value = 0.2
                    vmin = -2.
                    vmax = 2.
                    np_min = 5000
                elif var_output == 'ZDR_snow':
                    ref_value = 0.2
                    vmin = -2.
                    vmax = 2.
                    np_min = 5000
                elif var_output == 'dBZ_bias':
                    ref_value = 0.
                    vmin = -30.
                    vmax = 30.
                    np_min = 500

                fname = plot_monitoring_ts(
                    date, np_t_vec, cquant_vec, lquant_vec, hquant_vec,
                    get_fieldname_pyart(var_output), figfname,
                    ref_value=ref_value, vmin=vmin, vmax=vmax, np_min=np_min,
                    labelx='Time UTC', labely=labely, titl=titl)
                print('plotted file '+' '.join(fname))


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
