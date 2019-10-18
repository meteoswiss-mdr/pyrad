#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
rewrite_monitoring
================================================

This program rewrites a monitoring time series files into the correct
time order

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import numpy as np
import os

from pyrad.io import read_monitoring_ts, write_monitoring_ts
from pyrad.graph import plot_monitoring_ts
from pyrad.io import generate_field_name_str, get_fieldname_pyart


print(__doc__)


def main():
    """
    """

    input_base = (
        '/store/msrad/radar/pyrad_products/mals_tor_dataquality/')
    output_base = (
        '/store/msrad/radar/pyrad_products_test/mals_tor_dataquality/')

    var_vec = ['PhiDP0', 'RhoHV_rain', 'ZDR_prec', 'ZDR_snow', 'dBZ_bias']
    var_vec = ['ZDR_prec']

    plot_data = True

    print("====== Monitoring rewriting started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Monitoring rewriting finished: ")

    for j, var in enumerate(var_vec):
        if var == 'dBZ':
            dsdir = 'monitoring_clt_Zh'
            mon_type = 'GC_MONITORING'
            quantiles = [50., 95., 99.]
        elif var == 'dBZv':
            dsdir = 'monitoring_clt_Zv'
            mon_type = 'GC_MONITORING'
            quantiles = [50., 95., 99.]
        elif var == 'RhoHV_rain':
            dsdir = 'monitoring_RhoHV'
            mon_type = 'MONITORING'
            quantiles = [65., 80., 95.]
        elif var == 'PhiDP0':
            dsdir = 'monitoring_PhiDP0'
            mon_type = 'MONITORING'
            quantiles = [25., 50., 75.]
        elif var == 'ZDR_prec':
            dsdir = 'monitoring_ZDR_vert'
            # dsdir = 'monitoring_ZDR'
            mon_type = 'MONITORING'
            quantiles = [25., 50., 75.]
        elif var == 'ZDR_snow':
            dsdir = 'monitoring_ZDR_snow'
            mon_type = 'MONITORING'
            quantiles = [25., 50., 75.]
        elif var == 'dBZ_bias':
            dsdir = 'monitoring_Zh_bias'
            mon_type = 'MONITORING'
            quantiles = [25., 50., 75.]

        input_path = input_base+'/'+dsdir+'/VOL_TS/'
        output_path = output_base+'/'+dsdir+'/VOL_TS/'
        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        print('- Processing Variable '+var)
        fname_input = input_path+'ts_'+mon_type+'_'+var+'.csv'
        fname_output = output_path+'ts_'+mon_type+'_'+var+'.csv'
        figfname = [output_path+'ts_'+mon_type+'_'+var+'.png']

        date, np_t_vec, cquant_vec, lquant_vec, hquant_vec = (
            read_monitoring_ts(fname_input, sort_by_date=True))

        if date is None:
            continue

        val_vec = np.ma.asarray([lquant_vec, cquant_vec, hquant_vec]).T
        fname = write_monitoring_ts(
            date, np_t_vec, val_vec, quantiles, var, fname_output,
            rewrite=True)

        print('written file '+fname)

        if not plot_data:
            continue

        titldate = (date[0].strftime('%Y%m%d')+'-' +
                    date[-1].strftime('%Y%m%d'))
        titl = ' Monitoring '+titldate

        labely = generate_field_name_str(var)

        if var == 'dBZ':
            ref_value = 49.5
            vmin = 44.5
            vmax = 54.5
            np_min = 100000
        elif var == 'dBZv':
            ref_value = 51.5
            vmin = 46.5
            vmax = 56.5
            np_min = 100000
        elif var == 'RhoHV_rain':
            ref_value = 0.99
            vmin = 0.95
            vmax = 1.01
            np_min = 5000
        elif var == 'PhiDP0':
            ref_value = 0.
            vmin = -20.
            vmax = 20.
            np_min = 500000
        elif var == 'ZDR_prec':
            ref_value = 0.2
            vmin = -2.
            vmax = 2.
            np_min = 5000
        elif var == 'ZDR_snow':
            ref_value = 0.2
            vmin = -2.
            vmax = 2.
            np_min = 5000
        elif var == 'dBZ_bias':
            ref_value = 0.
            vmin = -30.
            vmax = 30.
            np_min = 100

        fname = plot_monitoring_ts(
            date, np_t_vec, cquant_vec, lquant_vec, hquant_vec,
            get_fieldname_pyart(var), figfname,
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
