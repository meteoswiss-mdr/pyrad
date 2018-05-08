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
from warnings import warn

from pyrad.io import read_profile_ts, get_fieldname_pyart
from pyrad.graph import _plot_time_range, get_field_name

from pyart.config import get_metadata

print(__doc__)


def main():
    """
    """
    file_base = '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/'
    time_dir_list = ['2017-06-29']
    trt_cell_id = '2017062913000174'
    hres = 250

    datatype_list = ['dBZc', 'ZDRc', 'RhoHVc', 'KDPc', 'TEMP', 'hydro']
    dataset_list = ['reflectivity', 'ZDRc', 'RhoHVc', 'KDPc', 'temperature', 'hydroclass']

    print("====== Plot time-height started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Plot time-height finished: ")

    for time_dir in time_dir_list:
        for i, datatype in enumerate(datatype_list):

            labels = ['50.0-percentile', '25.0-percentile', '75.0-percentile']
            if datatype == 'RhoHVc':
                labels = ['80.0-percentile', '65.0-percentile', '95.0-percentile']
            elif datatype == 'hydro':
                labels = [
                    'Mode', '2nd most common', '3rd most common',
                    '% points mode', '% points 2nd most common',
                    '% points 3nd most common']

            dataset = dataset_list[i]
            file_path = file_base+time_dir+'/'+dataset+'_trt_traj/PROFILE/'
            flist = glob.glob(
                file_path+'*_'+trt_cell_id+'_rhi_profile_*_'+datatype +
                '_hres'+str(hres)+'.csv')

            if not flist:
                warn('No profile files found in '+file_path+' for TRT cell ' +
                     trt_cell_id+' with resolution '+str(hres))
                continue

            tbin_edges, hbin_edges, data_ma = read_profile_ts(
                flist, labels, hres=hres)

            basepath_out = os.path.dirname(flist[0])
            fname = (
                basepath_out+'/'+trt_cell_id+'_trt_TIME_HEIGHT_'+datatype +
                '_hres'+str(hres)+'.png')
            field_name = get_fieldname_pyart(datatype)
            field_dict = get_metadata(field_name)
            titl = 'TRT cell '+trt_cell_id+'\n'+get_field_name(field_dict, field_name)

            vmin = vmax = None
            if datatype == 'RhoHVc':
                vmin = 0.95
                vmax = 1.00
            _plot_time_range(
                tbin_edges, hbin_edges, data_ma, field_name, [fname],
                titl=titl, figsize=[10, 8], vmin=vmin, vmax=vmax, dpi=72)

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
