#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_extract_trt
================================================

This program extracts individual TRT cell data from the original files
and puts it in a separate file for each cell

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import glob
import os
from warnings import warn

import numpy as np

from pyrad.io import read_trt_info2, read_trt_traj_data
from pyrad.io import write_trt_thundertracking_data

print(__doc__)


def main():
    """
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # keyword arguments
    parser.add_argument(
        '--infobase', type=str,
        default='/store/msrad/radar/thundertracking/info/',
        help='name of folder containing the TRT cell data')

    parser.add_argument(
        '--radarbase', type=str,
        default='/store/msrad/radar/DX50/rawdata/ThunderTracking_00_up.ele/',
        help='name of folder containing the TRT cell data')

    parser.add_argument(
        '--trtbase', type=str,
        default='/store/msrad/radar/trt/',
        help='name of folder containing the TRT cell data')

    args = parser.parse_args()

    print("====== TRT thundertracking info started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== TRT thundertracking info finished: ")

    file_list = glob.glob(args.infobase+'*.txt')
    if not file_list:
        warn('No info files in '+args.infobase)
        return None

    nrecords_total = 0
    nrecords_scan = 0
    nrecords_trt = 0
    for file in file_list:
        print('Reading TRT info file '+file)
        (trt_time_tt, ids_tt, rank_tt, scan_time_tt, azi_tt, rng_tt, lat_tt,
         lon_tt, ell_l_tt, ell_s_tt, ell_or_tt, vel_x_tt, vel_y_tt,
         det_tt) = read_trt_info2(file)

        if trt_time_tt is None:
            continue

        nrecords_total += trt_time_tt.size

        # look for scan
        scan_start_time_tt = np.ma.masked_all(
            trt_time_tt.size, dtype=datetime.datetime)
        nrecords_scan_local = 0
        for i, tscan in enumerate(scan_time_tt):
            scan_dir = tscan.strftime("%Y-%m-%d/")
            t_file_list = glob.glob(args.radarbase+scan_dir+'*dBZ.ele')

            if t_file_list is None:
                warn('No DX50 radar file matching scan time')
                continue

            t_file_list2 = []
            for filename in t_file_list:
                datetimestr = os.path.basename(filename)[0:14]
                fdatetime = datetime.datetime.strptime(
                    datetimestr, '%Y%m%d%H%M%S')
                if fdatetime is not None:
                    if ((fdatetime >= tscan) and
                            (fdatetime <= tscan+datetime.timedelta(seconds=60.))):
                        t_file_list2.append(fdatetime)

            if t_file_list2:
                scan_start_time_tt[i] = t_file_list2[0]
                nrecords_scan += 1
                nrecords_scan_local += 1

        # look for complete TRT cell info
        trt_files = glob.glob(
            args.trtbase+'*/TRTC_cell/'+str(ids_tt[0])+'.trt')
        if trt_files is None:
            warn('No trt file found for TRT ID '+str(ids_tt[0]))
            continue

        trt_file = trt_files[0]

        print('Reading TRT file '+trt_file)
        (ids_trt, trt_time_trt, lon_trt, lat_trt, _, _, _, area_trt, _, _, _,
         _, CG_n_trt, CG_p_trt, CG_trt, CG_percent_p_trt, ET45_trt, ET45m_trt,
         ET15_trt, ET15m_trt, VIL_trt, maxH_trt, maxHm_trt, POH_trt,
         rank2_trt, Dvel_x_trt, Dvel_y_trt, cell_contour_trt) = (
             read_trt_traj_data(trt_file))

        if ids_trt is None:
            warn('No TRT file found matching ID '+str(ids_tt[0]))
            continue

        nscans = trt_time_tt.size

        nrecords_trt += trt_time_trt.size
        print('TRT steps: '+str(trt_time_trt.size))
        print('Nominal X-band scans: '+str(nscans))
        print('Actual X-band scans: '+str(nrecords_scan_local))

        lon_tt = np.ma.masked_all(nscans)
        lat_tt = np.ma.masked_all(nscans)
        area_tt = np.ma.masked_all(nscans)
        CG_n_tt = np.ma.masked_all(nscans, dtype=int)
        CG_p_tt = np.ma.masked_all(nscans, dtype=int)
        CG_tt = np.ma.masked_all(nscans, dtype=int)
        CG_percent_p_tt = np.ma.masked_all(nscans)
        ET45_tt = np.ma.masked_all(nscans)
        ET45m_tt = np.ma.masked_all(nscans)
        ET15_tt = np.ma.masked_all(nscans)
        ET15m_tt = np.ma.masked_all(nscans)
        VIL_tt = np.ma.masked_all(nscans)
        maxH_tt = np.ma.masked_all(nscans)
        maxHm_tt = np.ma.masked_all(nscans)
        POH_tt = np.ma.masked_all(nscans)
        rank2_tt = np.ma.masked_all(nscans)
        Dvel_x_tt = np.ma.masked_all(nscans)
        Dvel_y_tt = np.ma.masked_all(nscans)
        cell_contour_tt = np.ma.masked_all(nscans, dtype=dict)

        for i, dt_trt in enumerate(trt_time_trt):
            lon_tt[trt_time_tt == dt_trt] = lon_trt[i]
            lat_tt[trt_time_tt == dt_trt] = lat_trt[i]
            area_tt[trt_time_tt == dt_trt] = area_trt[i]
            CG_n_tt[trt_time_tt == dt_trt] = CG_n_trt[i]
            CG_p_tt[trt_time_tt == dt_trt] = CG_p_trt[i]
            CG_tt[trt_time_tt == dt_trt] = CG_trt[i]
            CG_percent_p_tt[trt_time_tt == dt_trt] = CG_percent_p_trt[i]
            ET45_tt[trt_time_tt == dt_trt] = ET45_trt[i]
            ET45m_tt[trt_time_tt == dt_trt] = ET45m_trt[i]
            ET15_tt[trt_time_tt == dt_trt] = ET15_trt[i]
            ET15m_tt[trt_time_tt == dt_trt] = ET15m_trt[i]
            VIL_tt[trt_time_tt == dt_trt] = VIL_trt[i]
            maxH_tt[trt_time_tt == dt_trt] = maxH_trt[i]
            maxHm_tt[trt_time_tt == dt_trt] = maxHm_trt[i]
            POH_tt[trt_time_tt == dt_trt] = POH_trt[i]
            rank2_tt[trt_time_tt == dt_trt] = rank2_trt[i]
            Dvel_x_tt[trt_time_tt == dt_trt] = Dvel_x_trt[i]
            Dvel_y_tt[trt_time_tt == dt_trt] = Dvel_y_trt[i]
            cell_contour_tt[trt_time_tt == dt_trt] = cell_contour_trt[i]

        fname = trt_file.replace('.trt', '_tt.trt')
        write_trt_thundertracking_data(
            ids_tt, scan_time_tt, scan_start_time_tt, azi_tt, rng_tt,
            trt_time_tt, lon_tt, lat_tt, ell_l_tt, ell_s_tt, ell_or_tt,
            area_tt, vel_x_tt, vel_y_tt, det_tt, rank_tt, CG_n_tt, CG_p_tt,
            CG_tt, CG_percent_p_tt, ET45_tt, ET45m_tt, ET15_tt, ET15m_tt,
            VIL_tt, maxH_tt, maxHm_tt, POH_tt, rank2_tt, Dvel_x_tt, Dvel_y_tt,
            cell_contour_tt, fname)

    print('\nTotal TRT steps '+str(nrecords_trt))
    print('Total X-band records '+str(nrecords_total))
    print('Total X-band scans '+str(nrecords_scan))

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
