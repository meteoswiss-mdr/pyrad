#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
get_rpc
================================================

This program computes the rimed particles column

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

from pyrad.io import write_trt_rpc, read_profile_ts
from pyrad.io import read_trt_thundertracking_traj_data
from pyrad.io import read_thundertracking_info

print(__doc__)


def main():
    """
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # keyword arguments
    parser.add_argument(
        '--days', type=str, default=None,
        help='Dates to process. Format YYYY-MM-DD. Coma separated')

    parser.add_argument(
        '--info_file', type=str,
        default='/store/msrad/radar/thundertracking/info/thundertracking_info.csv',
        help='configuration file path')

    parser.add_argument(
        '--years', type=str,
        default=None,
        help='Years to process. If None all years in file will be processed')

    parser.add_argument(
        '--max_rank', type=float,
        default=None,
        help='Max rank to process')

    parser.add_argument(
        '--trtbase', type=str,
        default='/store/msrad/radar/trt/',
        help='name of folder containing the TRT cell data')

    parser.add_argument(
        '--radarbase', type=str,
        default='/store/msrad/radar/pyrad_products/thundertracking/',
        help='name of folder containing the radar data')

    parser.add_argument(
        '--dataset', type=str,
        default='hydro',
        help='Name of the directory containing the datasets')

    parser.add_argument(
        '--datatype', type=str,
        default='hydro',
        help='Name of the directory containing the datasets')

    parser.add_argument(
        '--hres', type=float, default=250., help='Height resolution')

    parser.add_argument(
        '--path_structure', type=int, default=2,
        help='If true the data is at the cell center')

    args = parser.parse_args()

    print("====== PYRAD RPC processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD RPC processing finished: ")

    datatype = args.datatype
    labels = [
        'Mode', '2nd most common', '3rd most common', '% points mode',
        '% points 2nd most common', '% points 3rd most common']

    print('trt path: '+args.trtbase)
    print('radar data path: '+args.radarbase)

    if args.days is not None:
        time_dir_list = args.days.split(',')
    else:
        # get the years to process
        years = None
        if args.years is not None:
            years = list(map(int, args.years.split(',')))

        _, max_rank, _, trt_time_start, trt_time_end = (
            read_thundertracking_info(args.info_file))
        trt_times = np.append(trt_time_start, trt_time_end)

        trt_dates = np.array([], dtype=datetime.date)
        for trt_time in trt_times:
            trt_dates = np.append(trt_dates, trt_time.date())
        trt_dates = np.sort(np.unique(trt_dates))
        time_dir_list = []
        for rank, trt_date in zip(max_rank, trt_dates):
            if years is not None:
                if trt_date.year not in years:
                    continue

            if args.max_rank is not None:
                if rank > args.max_rank:
                    continue

            time_dir_list.append(trt_date.strftime("%Y-%m-%d"))

    # Find all TRT files in directory
    trt_list = []
    if args.days is not None:
        for time_dir in time_dir_list:
            trt_list.extend(glob.glob(
                args.trtbase+time_dir+'/TRTC_cell_plots/All/*.trt'))
            trt_list.extend(glob.glob(
                args.trtbase+time_dir+'/TRTC_cell_plots/Some/*.trt'))
    else:
        for time_dir in time_dir_list:
            trt_list.extend(glob.glob(
                args.trtbase+time_dir+'/TRTC_cell/*_tt.trt'))

    if len(trt_list) == 0:
        warn('No valid TRT files found in '+args.trtbase)
        return

    # get altitude of graupel column
    for trt_fname in trt_list:
        trt_cell_id = os.path.basename(trt_fname).split('.')[0]
        trt_cell_id = trt_cell_id.replace('_tt', '')

        dt_str = trt_cell_id[0:12]
        dt_cell1 = datetime.datetime.strptime(dt_str, "%Y%m%d%H%M")
        dt_cell2 = dt_cell1+datetime.timedelta(days=1)
        time_dir = [
            dt_cell1.strftime("%Y-%m-%d"), dt_cell2.strftime("%Y-%m-%d")]

        if args.path_structure == 1:
            file_base2 = [
                args.radarbase+time_dir[0]+'/'+args.dataset+'_trt_center_traj/',
                args.radarbase+time_dir[1]+'/'+args.dataset+'_trt_center_traj/']
        elif args.path_structure == 0:
            file_base2 = [
                args.radarbase+time_dir[0]+'/'+args.dataset+'_trt_traj/',
                args.radarbase+time_dir[1]+'/'+args.dataset+'_trt_traj/']
        elif args.path_structure == 2:
            file_base2 = [
                args.radarbase+time_dir[0]+'/trt_traj_tt/',
                args.radarbase+time_dir[1]+'/trt_traj_tt/']
        elif args.path_structure == 3:
            file_base2 = [
                args.radarbase+time_dir[0]+'/trt_traj_vis_filt/',
                args.radarbase+time_dir[1]+'/trt_traj_vis_filt/']

        # plot time-height
        if args.path_structure in (2, 3):
            flist = []
            for fbase in file_base2:
                flist.extend(glob.glob(
                    fbase+'PROFILE_'+args.dataset+'/*_'+trt_cell_id +
                    '_rhi_profile_*_'+datatype+'_hres'+str(int(args.hres)) +
                    '.csv'))
        else:
            flist = []
            for fbase in file_base2:
                flist.extend(glob.glob(
                    fbase+'PROFILE/*_'+trt_cell_id+'_rhi_profile_*_' +
                    datatype+'_hres'+str(int(args.hres))+'.csv'))

        if not flist:
            warn('No profile files found in '+file_base2[0] +
                 'PROFILE/ for TRT cell ' +
                 trt_cell_id+' with resolution '+str(args.hres))
            continue

        tbin_edges, hbin_edges, _, data_ma, start_time = (
            read_profile_ts(flist, labels, hres=args.hres, t_res=None))

        # Get min and max altitude of graupel/hail area
        (traj_ID, _, scan_time, _, _, _, lon, lat, _, _, _, area, _, _, _,
         RANKr, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _) = (
             read_trt_thundertracking_traj_data(trt_fname))

        # Mask data where there is no scan
        ind = np.where(np.ma.getmaskarray(scan_time) == 0)
        traj_ID = traj_ID[ind]
        lon = lon[ind]
        lat = lat[ind]
        area = area[ind]
        RANKr = RANKr[ind]
        scan_time = scan_time[ind]

        hmin, hmax, freq = get_graupel_column(
            tbin_edges, hbin_edges, data_ma, start_time, scan_time)

        fname = (
            args.trtbase+'rpc/'+trt_cell_id +
            '_cell_rimed_particles_column.csv')
        write_trt_rpc(
            traj_ID, scan_time, lon, lat, area, RANKr, hmin, hmax, freq,
            fname, timeformat='%Y%m%d%H%M%S')

        print("----- written to '%s'" % fname)



def get_graupel_column(tbin_edges, hbin_edges, data_ma, start_time,
                       yyyymmddHHMM):
    """
    Gets the minimum and maximum heigth of the graupel column. This is defined
    as the maximum extend where rimed particles (code 5) or solid hail (code
    10) in the hydrometeor classification are dominant

    Parameters
    ----------
    tbin_edges : 1D array of floats
        The time bin edges [s]
    hbin_edges : 1D array of floats
        The height bin edges [m MSL]
    data_ma : 2D array of ints
        Matrix containing the time-height hydrometeor classification
        information
    start_time : datetime object
        start time of the radar scan
    yyyymmddHHMM : 1D array of datetime objects
        time steps of the TRT cell

    Returns
    -------
    hmin, hmax : 1D float arrays
        the minimum and maximum altitude of the rimed particles column

    """
    tbin_rights = tbin_edges[1:]
    hbin_lefts = hbin_edges[:-1]
    hbin_rights = hbin_edges[1:]
    H_lefts, _ = np.meshgrid(hbin_lefts, tbin_rights)
    H_lefts = np.ma.asarray(H_lefts)
    H_rights, _ = np.meshgrid(hbin_rights, tbin_rights)
    H_rights = np.ma.asarray(H_rights)
    valid = np.logical_or(data_ma == 5, data_ma == 10)
    mask = np.logical_not(valid)
    H_lefts[mask] = np.ma.masked
    H_rights[mask] = np.ma.masked
    hmin = np.ma.min(H_lefts, axis=1)
    hmax = np.ma.max(H_rights, axis=1)
    ind_hmin = np.ma.argmin(H_lefts, axis=1)
    ind_hmax = np.ma.argmax(H_rights, axis=1)
    n_rp = (ind_hmax-ind_hmin+1).astype(float)

    freq = np.ma.zeros((hmin.size))
    for tt in range(tbin_rights.size):
        ind = np.where(valid[tt, :] == 1)[0]
        if ind.size == 0:
            continue
        freq[tt] = 100.*float(ind.size)/n_rp[tt]

    # All TRT cell time steps have a corresponding rimed column height value
    # Return the values
    if tbin_rights.size == yyyymmddHHMM.size:
        return hmin, hmax, freq

    print('Missing '+str(yyyymmddHHMM.size-tbin_rights.size)+' time steps')

    # Missing rimed height values. Determine those missing and put the data
    # in the right time step
    hmin_aux = np.ma.masked_all(yyyymmddHHMM.size, dtype=float)
    hmax_aux = np.ma.masked_all(yyyymmddHHMM.size, dtype=float)
    freq_aux = np.ma.masked_all(yyyymmddHHMM.size, dtype=float)
    for k, dt_trt in enumerate(yyyymmddHHMM):
        for l, dt_rad in enumerate(start_time):
            if dt_trt == dt_rad:
                hmin_aux[k] = hmin[l]
                hmax_aux[k] = hmax[l]
                freq_aux[k] = freq[l]
                break
    return hmin_aux, hmax_aux, freq_aux


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
