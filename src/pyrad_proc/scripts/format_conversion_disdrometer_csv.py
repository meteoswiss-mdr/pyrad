#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
format_conversion_disdrometer_csv
================================================

This program converts the disdrometer csv time series files in the old
output format into the new format and saved the csv files.

"""

# Author: sue

from __future__ import print_function
import datetime
import glob
import os
import csv
import pathlib
import numpy as np

from pyrad.io import read_data_sensor


print(__doc__)


def main():
    """
    """

    basepath = '/data/disdrometer/mals_parsivel/amfortas/scattering/'
    os.chdir(basepath)
    file_list = glob.glob('*.txt')

    for fname in file_list:
        var = read_data_sensor.read_disdro_scattering(fname)
        write_csv_file(var)


def get_precip_type(precipcode):
    """
    maps the old precip type number into the new type name

    Parameters
    ----------
    datatype : float
        old precip type number

    Returns
    -------
    precip_type : str
        new precip type name

    """
    if precipcode == 0.0:
        field_name = 'RA'
    elif precipcode == 1.0:
        field_name = 'SN'
    elif np.isnan(precipcode):
        field_name = 'NP'
    else:
        raise ValueError('ERROR: Unknown data type '+precipcode)
    return field_name


def write_csv_file(data):
    """
    writes time series of data

    Parameters
    ----------
    var : tuple containing timeseries data of each variable

    """
    basepath = '/data/disdrometer/mals_parsivel/amfortas/scattering/'
    datapath = (data[0][0].strftime('%Y%m%d')[0:4]+'/' +
                data[0][0].strftime('%Y%m%d')[0:6]+'/')
    pathlib.Path(basepath+datapath).mkdir(parents=True, exist_ok=True)
    varnames = ['LWC', 'RR', 'dBZ', 'dBZv', 'ZDR', 'LDR', 'Ah',
                'Av', 'Adp', 'KDP', 'DeltaCo', 'RhoHV']
    i = 0
    for var in data[2:]:
        fname = (data[0][0].strftime('%Y%m%d')+'_amfortas_pay'+'_9.7GHz_' +
                 varnames[i]+'_el90.0.csv')
        with open(basepath+datapath+fname, 'w', newline='') as csvfile:
            csvfile.write('# Disdrometer timeseries data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Description: \n')
            csvfile.write('# Time series of '+varnames[i]+'\n')
            csvfile.write(
                '# Location [lat  lon]: ' +
                '46.84900  '+'6.92809\n')
            csvfile.write(
                '# Elevation: ' + '443.5m MSL\n')
            csvfile.write(
                '# Elevation angle: ' + '90\n')
            csvfile.write(
                '# Scattering Frequency: ' +
                '9.7 GHz\n')
            csvfile.write('# Data: ' + varnames[i] + '\n')
            csvfile.write(
                '# Start: ' +
                data[0][0].strftime(
                    '%Y-%m-%d %H:%M:%S UTC') + '\n')
            csvfile.write('#\n')

            fieldnames = ['date', 'Precip Code', varnames[i],
                          'Scattering Temp [deg C]']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for j in range(len(data[0])):
                writer.writerow(
                    {'date': data[0][j].strftime('%Y-%m-%d %H:%M:%S'),
                     'Precip Code': get_precip_type(data[1][j]),
                     varnames[i]: float(-9999) if np.isnan(var[j]) else var[j],
                     'Scattering Temp [deg C]': '10.'})
            i += 1
            csvfile.close()


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
