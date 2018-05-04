#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program processes and post-processes data over a time span

To run the processing framework type:
    python main_process_data.py \
[config_file] --starttime [process_start_time] --endtime [process_end_time] \
--postproc_cfgfile [postproc_config_file] --cfgpath [cfgpath]

If startime and endtime are not specified the program determines them from
the trajectory file or the last processed volume.
postproc_cfgfile is an optional argument with default: None
cfgpath is an optional argument with default: \
'$HOME/pyrad/config/processing/'
The trajectory file can be of type plane or type lightning. If it is of type \
lightning the flash number can be specified

Example:
    python main_process_data.py 'paradiso_fvj_vol.txt' --starttime \
'20140523000000' --endtime '20140523001000' --postproc_cfgfile \
'paradiso_fvj_vol_postproc.txt' --cfgpath '$HOME/pyrad/config/processing/'

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import os
import glob

import numpy as np

from pyrad.io import read_rhi_profile, get_fieldname_pyart
from pyrad.graph import _plot_time_range, get_field_name

from pyart.config import get_metadata

print(__doc__)


def main():
    """
    """
    file_path = '/data/pyrad_products/rad4alp_hydro_PHA/2017-06-29/reflectivity_trt_traj/PROFILE/'
    trt_cell_id = '2017062913000182'
    datatype = 'dBZc'
    hres = 250

    labels = ['50.0-percentile', '25.0-percentile', '75.0-percentile']
    if datatype == 'RhoHVc':
        labels = ['80.0-percentile', '65.0-percentile', '95.0-percentile']
    elif datatype == 'hydro':
        labels = ['Mode ', '% points mode']

    print("====== Plot time-height started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Plot time-height finished: ")

    flist = glob.glob(
        file_path+'*_'+trt_cell_id+'_rhi_profile_*_'+datatype +
        '_hres'+str(hres)+'.csv')

    data_ma = []
    datetime_arr = np.ma.array([], dtype=datetime.datetime)
    for fname in flist:
        bfile = os.path.basename(fname)
        datetimestr = bfile[0:14]
        fdatetime = datetime.datetime.strptime(datetimestr, '%Y%m%d%H%M%S')
        datetime_arr = np.append(datetime_arr, fdatetime)

        print('\nReading file '+fname)
        height, np_t, vals = read_rhi_profile(fname, labels)
        val = vals[:, 0]
        print(np.shape(val))
        data_ma.append(val)
    data_ma = np.ma.asarray(data_ma)

    # sort data as a function of time
    ind = np.argsort(datetime_arr)
    datetime_arr = datetime_arr[ind]
    data_ma = data_ma[ind, :]

    basepath_out = os.path.dirname(flist[0])
    fname = (
        basepath_out+'/'+trt_cell_id+'_trt_TIME_HEIGHT_'+datatype +
        '_hres'+str(hres)+'.png')
    field_name = get_fieldname_pyart(datatype)
    field_dict = get_metadata(field_name)
    titl = 'TRT cell '+trt_cell_id+'\n'+get_field_name(field_dict, field_name)

    # put date time array as seconds from start of TRT cell
    dt_s = np.empty(datetime_arr.size, dtype=float)
    for i, dt in enumerate(datetime_arr):
        dt_s[i] = (dt-datetime_arr[0]).total_seconds()

    _plot_time_range(
        dt_s, height, data_ma, field_name, [fname], titl=titl,
        figsize=[10, 8], dpi=72)

    print("----- plot to '%s'" % fname)


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
