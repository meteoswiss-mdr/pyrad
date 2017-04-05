#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_retrieve_dwh
================================================

This program retrieves parameters from selected SwissMetNet stations, stores
them in a file and plots them and computes the average and standard deviation
over a period and stores it and plots it.

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
from subprocess import run

from pyrad.io import read_smn2, write_smn
from pyrad.graph import plot_timeseries
from pyrad.util import time_series_statistics, join_time_series

print(__doc__)


def main():
    """
    """

    file_path = '/data/FLORAKO/v2/'
    img_ext = 'png'
    avg_time = 3600
    base_time = 0

    smn_station_vec = ['PIL', 'WFJ', 'MTR', 'GUE', 'DIA']
#    smn_station_vec = ['VSTSN']
    tstart_vec = ['2007', '2008', '2009', '2010', '2011', '2012', '2013',
                  '2014', '2015', '2016']
    tend_vec = ['2008', '2009', '2010', '2011', '2012', '2013', '2014',
                '2015', '2016', '2017']
#    param_vec = ['tre200s0', 'tde200s0', 'ure200s0', 'gor000z0', 'fkl010z0']
#    ymin_vec = [-40., -40., 0., -40., 0.]
#    ymax_vec = [30., 30., 110., 1500., 20.]

    param_vec = ['prestas0']
    ymin_vec = [600.]
    ymax_vec = [900.]

    print("====== retrieval from DWH started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== retrieval from DWH finished: ")

    for station in smn_station_vec:
        i = 0
        for param in param_vec:
            ymin = ymin_vec[i]
            ymax = ymax_vec[i]
            i += 1
            for time in range(len(tstart_vec)):
                tstart = tstart_vec[time]+'0101000000'
                tend = tend_vec[time]+'0101000000'

                print('\n--- Processing station '+station+' parameter ' +
                      param+' between '+tstart+' and '+tend)

                base_name = station+'_'+param+'_'+tstart+'-'+tend
                file_name = file_path+base_name+'.csv'
                with open(file_name, 'wb') as out_file:
                    result = run(
                        ['jretrievedwh.sh',  '-d', 'COMMA', '-s', 'surface',
                         '-i', 'nat_abbr,'+station, '-t', tstart+','+tend,
                         '-n', param], stdout=out_file)
                id, date, value = read_smn2(file_name)

                if date is None:
                    continue

                plot_timeseries(
                    date, [value], [file_path+base_name+'.'+img_ext],
                    labelx='Time [UTC]', labely=param, labels=[station],
                    title=station+' '+tstart+'-'+tend, period=0,
                    timeformat=None, colors=None, linestyles=None, ymin=ymin,
                    ymax=ymax)

                date_avg, value_avg = time_series_statistics(
                    date, value, avg_time=avg_time, base_time=base_time,
                    method='mean', dropnan=True)
                date_std, value_std = time_series_statistics(
                    date, value, avg_time=avg_time, base_time=base_time,
                    method='std', dropnan=True)

                date_series, value_avg, value_std = join_time_series(
                    date_avg, value_avg, date_std, value_std, dropnan=True)

                plot_timeseries(
                    date_series,
                    [value_avg, value_avg+value_std, value_avg-value_std],
                    [file_path+base_name+'_avg'+str(avg_time)+'s.'+img_ext],
                    labelx='Time [UTC]', labely=param, labels=None,
                    title=station+' '+tstart+'-'+tend+' avg '+str(
                        avg_time)+' s', period=0, timeformat=None,
                    colors=['b', 'r', 'r'], linestyles=['-', '--', '--'],
                    ymin=ymin, ymax=ymax)

                write_smn(
                    date_series, value_avg, value_std,
                    file_path+base_name+'_avg'+str(avg_time)+'s.csv')


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
