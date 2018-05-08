#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

This program processes TRT data

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
import argparse
import atexit
import os
import glob
from warnings import warn

import numpy as np

from pyrad.flow.flow_control import main as pyrad_main
from pyrad.io import get_fieldname_pyart
from pyrad.io import read_profile_ts, read_histogram_ts, read_quantiles_ts
from pyrad.graph import get_field_name, get_colobar_label
from pyrad.graph import _plot_time_range

from pyart.config import get_metadata

print(__doc__)


def main():
    """
    """

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments
    parser.add_argument(
        'proc_cfgfile', type=str, help='name of main configuration file')

    parser.add_argument(
        'trtpath', type=str,
        help='name of folder containing the TRT cell data')

    # keyword arguments
    parser.add_argument(
        '--cfgpath', type=str,
        default=os.path.expanduser('~')+'/pyrad/config/processing/',
        help='configuration file path')

    args = parser.parse_args()

    print("====== PYRAD data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD data processing finished: ")

    print('config path: '+args.cfgpath)
    print('config file: '+args.proc_cfgfile)
    print('trt path: '+args.trtpath)

    cfgfile_proc = args.cfgpath+args.proc_cfgfile
    trajtype = 'trt'

    trt_list = glob.glob(args.trtpath+'*.trt')
    trt_cell_id_list = []
    for fname in trt_list:
        print('processing TRT cell file '+fname)
        try:
            infostr = os.path.basename(fname).split('.')[0]
            pyrad_main(
                cfgfile_proc, trajfile=fname, infostr=infostr,
                trajtype=trajtype)
            trt_cell_id_list.append(infostr)
        except ValueError:
            print(ValueError)

    # plot time series
    datatype_list = ['dBZc', 'ZDRc', 'RhoHVc', 'KDPc', 'TEMP', 'hydro']
    dataset_list = [
        'reflectivity', 'ZDRc', 'RhoHVc', 'KDPc', 'temperature', 'hydroclass']
    file_base = '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/'
    hres = 250

    for i, trt_cell_id in enumerate(trt_cell_id_list):
        dt_str = trt_cell_id[0:12]
        dt = datetime.datetime.strptime(dt_str, "%Y%m%d%H%M")
        time_dir = dt.strftime("%Y-%m-%d")
        for j, datatype in enumerate(datatype_list):
            dataset = dataset_list[j]
            file_base2 = file_base+time_dir+'/'+dataset+'_trt_traj/'

            field_name = get_fieldname_pyart(datatype)
            field_dict = get_metadata(field_name)
            titl = 'TRT cell '+trt_cell_id+'\n'+get_field_name(
                field_dict, field_name)

            # plot time-height
            flist = glob.glob(
                file_base2+'PROFILE/*_'+trt_cell_id+'_rhi_profile_*_' +
                datatype+'_hres'+str(hres)+'.csv')
            if not flist:
                warn('No profile files found in '+file_base2 +
                     'PROFILE/ for TRT cell ' +
                     trt_cell_id+' with resolution '+str(hres))
            else:
                labels = [
                    '50.0-percentile', '25.0-percentile', '75.0-percentile']
                if datatype == 'RhoHVc':
                    labels = [
                        '80.0-percentile', '65.0-percentile',
                        '95.0-percentile']
                elif datatype == 'hydro':
                    labels = [
                        'Mode', '2nd most common', '3rd most common',
                        '% points mode', '% points 2nd most common',
                        '% points 3nd most common']
                tbin_edges, hbin_edges, data_ma = read_profile_ts(
                    flist, labels, hres=hres)

                basepath_out = os.path.dirname(flist[0])
                fname = (
                    basepath_out+'/'+trt_cell_id+'_trt_TIME_HEIGHT_' +
                    datatype+'_hres'+str(hres)+'.png')

                vmin = vmax = None
                if datatype == 'RhoHVc':
                    vmin = 0.95
                    vmax = 1.00
                _plot_time_range(
                    tbin_edges, hbin_edges, data_ma, field_name, [fname],
                    titl=titl, figsize=[10, 8], vmin=vmin, vmax=vmax, dpi=72)

                print("----- plot to '%s'" % fname)

            # plot time-hist
            flist = glob.glob(
                file_base2+'HISTOGRAM/*_'+trt_cell_id+'_histogram_*_' +
                datatype+'.csv')

            if not flist:
                warn('No histogram files found in '+file_base2 +
                     'HISTOGRAM/ for TRT cell '+trt_cell_id)
            else:
                tbin_edges, bin_edges, data_ma = read_histogram_ts(
                    flist, datatype)

                basepath_out = os.path.dirname(flist[0])
                fname = (
                    basepath_out+'/'+trt_cell_id+'_trt_HISTOGRAM_'+datatype +
                    '.png')

                _plot_time_range(
                    tbin_edges, bin_edges, data_ma, 'frequency_of_occurrence',
                    [fname], titl=titl,
                    ylabel=get_colobar_label(field_dict, field_name),
                    vmin=0., vmax=np.max(data_ma), figsize=[10, 8], dpi=72)

                print("----- plot to '%s'" % fname)

            # plot quantiles
            flist = glob.glob(
                file_base2+'QUANTILES/*_'+trt_cell_id+'_quantiles_*_' +
                datatype+'.csv')

            if not flist:
                warn('No quantiles files found in '+file_base2 +
                     'QUANTILES/ for TRT cell '+trt_cell_id)
                continue

            tbin_edges, qbin_edges, data_ma = read_quantiles_ts(
                flist, step=5., qmin=0., qmax=100.)

            basepath_out = os.path.dirname(flist[0])
            fname = (
                basepath_out+'/'+trt_cell_id+'_trt_QUANTILES_'+datatype +
                '.png')

            vmin = vmax = None
            if datatype == 'RhoHVc':
                vmin = 0.95
                vmax = 1.00
            _plot_time_range(
                tbin_edges, qbin_edges, data_ma, field_name, [fname],
                titl=titl, ylabel='Quantile', vmin=vmin, vmax=vmax,
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
