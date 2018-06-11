#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_lightning_trt
================================================

This program processess lightning information on TRT cells

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import glob
from warnings import warn
import os

import numpy as np

from pyrad.io import read_trt_traj_data, read_lightning, write_histogram
from pyrad.io import read_histogram_ts
from pyrad.graph import plot_histogram, plot_pos, plot_timeseries
from pyrad.graph import _plot_time_range
from pyrad.util import belongs_roi_indices

print(__doc__)


def main():
    """
    """
    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    parser.add_argument(
        'trtpath', type=str, help='name of folder containing the TRT cell data')

    # keyword arguments
    parser.add_argument(
        '--flashpath', type=str,
        default='/store/msrad/lightning/LMA/Santis/',
        help='name of the folder containing the Santis LMA data')

    args = parser.parse_args()

    print("====== Lightning TRT cells started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Lightning TRT cells finished: ")

    print('trt path: '+args.trtpath)

    # Get bins altitude
    alt_min = 0.
    alt_max = 20000.
    step = 100.
    bin_edges = np.linspace(
        alt_min-step/2., alt_max+step/2, num=int((alt_max-alt_min)/step)+2)

    trt_list = glob.glob(args.trtpath+'*.trt')
    for trt_fname in trt_list:
        print('processing TRT cell file '+trt_fname)

        # reading TRT cell file
        (traj_ID, yyyymmddHHMM, lon, lat, ell_L, ell_S, ell_or, area,
         vel_x, vel_y, det, RANKr, CG_n, CG_p, CG, CG_percent_p, ET45,
         ET45m, ET15, ET15m, VIL, maxH, maxHm, POH, RANK, Dvel_x,
         Dvel_y, cell_contour) = read_trt_traj_data(trt_fname)

        # reading lightning file
        infostr = os.path.basename(trt_fname).split('.')[0]
        dt_str = infostr[0:12]
        dt = datetime.datetime.strptime(dt_str, '%Y%m%d%H%M')

        flashnr, time_data, time_in_flash, lat, lon, alt, dBm = (
            read_lightning(args.flashpath+dt.strftime("%y%m%d")+'.txt'))

        if flashnr is None:
            continue

        # Get first sources data
        flashnr_first, unique_ind = np.unique(flashnr, return_index=True)
        lat_first = lat[unique_ind]
        lon_first = lon[unique_ind]
        alt_first = alt[unique_ind]
        dBm_first = dBm[unique_ind]
        time_first = time_data[unique_ind]

        # Find cell period
        cell_dt_s = np.empty(yyyymmddHHMM.size, dtype=float)
        for i, cell_time in enumerate(yyyymmddHHMM):
            cell_dt_s[i] = (cell_time-yyyymmddHHMM[0]).total_seconds()
        t_res = np.mean(cell_dt_s[1:]-cell_dt_s[:-1])

        # analyze first sources
        print('\n\n--- Processing first sources ----')
        nflashes = np.zeros(yyyymmddHHMM.size, dtype=int)
        for i, cell_time in enumerate(yyyymmddHHMM):
            cell_time_str = cell_time.strftime("%Y%m%d%H%M%S")

            # Find flashes within time step of cell
            tstart_cell_step = cell_time-datetime.timedelta(seconds=t_res)
            inds = np.where(np.logical_and(
                time_first > tstart_cell_step, time_first <= cell_time))[0]
            if inds.size == 0:
                warn('No flashes within time step')
                fname_hist = (
                    args.trtpath+cell_time_str+'_'+infostr +
                    '_hist_alt_first_source.csv')
                fname_hist = write_histogram(
                    bin_edges, np.zeros(bin_edges.size-1, dtype=int),
                    fname_hist, step=step)
                print('----- save to '+fname_hist)
                continue

            lat_cell = lat_first[inds]
            lon_cell = lon_first[inds]
            alt_cell = alt_first[inds]

            # Find flashes within cell contour
            inds, is_roi = belongs_roi_indices(
                lat_cell, lon_cell, cell_contour[i])
            if is_roi == 'None':
                warn('No flashes within cell contour')
                fname_hist = (
                    args.trtpath+cell_time_str+'_'+infostr +
                    '_hist_alt_first_source.csv')
                fname_hist = write_histogram(
                    bin_edges, np.zeros(bin_edges.size-1, dtype=int),
                    fname_hist, step=step)
                print('----- save to '+fname_hist)
                continue
            elif is_roi == 'All':
                inds = inds[0]

            lat_cell = lat_cell[inds]
            lon_cell = lon_cell[inds]
            alt_cell = alt_cell[inds]

            # compute number of flashes
            nflashes[i] = lat_cell.size

            # Plot altitude histogram
            fname_hist = (
                args.trtpath+cell_time_str+'_'+infostr +
                '_hist_alt_first_source.png')
            titl = (
                cell_time_str+' '+infostr +
                ' TRT cell. Flash first source altitude')
            fname_hist = plot_histogram(
                bin_edges, alt_cell, [fname_hist], labelx='Altitude [m MSL]',
                titl=titl)
            print('Plotted '+' '.join(fname_hist))

            fname_hist = (
                args.trtpath+cell_time_str+'_'+infostr +
                '_hist_alt_first_source.csv')
            hist, bin_edges = np.histogram(alt_cell, bins=bin_edges)
            fname_hist = write_histogram(
                bin_edges, hist, fname_hist, step=step)
            print('----- save to '+fname_hist)

            # plot position first source
            figfname = (
                args.trtpath+cell_time_str+'_'+infostr +
                '_first_source_pos_max_height_on_top.png')
            titl = (
                cell_time_str+' '+infostr +
                ' first flash source position. Highest on top')
            figfname = plot_pos(
                lat_cell, lon_cell, alt_cell, [figfname],
                sort_altitude='Highest_on_top', cb_label='Source height [m MSL]',
                titl=titl)
            print('Plotted '+' '.join(figfname))

        # plot time series lightning
        figfname = (
            args.trtpath+cell_time_str+'_'+infostr+'_N_first_sources.png')
        titl = infostr+' Number of first sources in cell'
        figfname = plot_timeseries(
            yyyymmddHHMM, [nflashes], [figfname], labelx='Time [UTC]',
            labely='N flashes', labels=['N flashes'], title=titl,
            period=0, timeformat=None, colors=None, linestyles=None,
            markers=None, ymin=None, ymax=None, dpi=72)
        print('Plotted '+' '.join(figfname))

        # plot time-hist_height
        flist = glob.glob(
            args.trtpath+'*_'+infostr+'_hist_alt_first_source.csv')

        if not flist:
            warn('No histogram files found in '+args.trtpath +
                 ' for TRT cell '+infostr)
        else:
            tbin_edges, bin_edges, data_ma = read_histogram_ts(
                flist, 'flash_altitude')

            vmax = np.max(data_ma)
            if vmax == 0.:
                warn('Unable to plot histogram. No valid data')
            else:
                fname_hist = (
                    args.trtpath+infostr+'_trt_HISTOGRAM_alt_first_source.png')
                titl = 'TRT cell '+infostr+'\n'+'Altitude of first flash source'
                _plot_time_range(
                    tbin_edges, bin_edges, data_ma, 'frequency_of_occurrence',
                    [fname_hist], titl=titl,
                    ylabel='Altitude [m MSL]',
                    vmin=0., vmax=vmax, figsize=[10, 8], dpi=72)

                print("----- plot to '%s'" % fname_hist)

        # analyze all flashes
        print('\n\n--- Processing all sources ----')
        nflashes = np.zeros(yyyymmddHHMM.size, dtype=int)
        for i, cell_time in enumerate(yyyymmddHHMM):
            cell_time_str = cell_time.strftime("%Y%m%d%H%M%S")

            # Find flashes within time step of cell
            tstart_cell_step = cell_time-datetime.timedelta(seconds=t_res)
            inds = np.where(np.logical_and(
                time_data > tstart_cell_step, time_data <= cell_time))[0]
            if inds.size == 0:
                warn('No flashes within time step')
                fname_hist = (
                    args.trtpath+cell_time_str+'_'+infostr +
                    '_hist_alt_all_sources.csv')
                fname_hist = write_histogram(
                    bin_edges, np.zeros(bin_edges.size-1, dtype=int),
                    fname_hist, step=step)
                print('----- save to '+fname_hist)
                continue

            lat_cell = lat[inds]
            lon_cell = lon[inds]
            alt_cell = alt[inds]

            # Find flashes within cell contour
            inds, is_roi = belongs_roi_indices(
                lat_cell, lon_cell, cell_contour[i])
            if is_roi == 'None':
                warn('No flashes within cell contour')
                fname_hist = (
                    args.trtpath+cell_time_str+'_'+infostr +
                    '_hist_alt_all_sources.csv')
                fname_hist = write_histogram(
                    bin_edges, np.zeros(bin_edges.size-1, dtype=int),
                    fname_hist, step=step)
                print('----- save to '+fname_hist)
                continue
            elif is_roi == 'All':
                inds = inds[0]

            lat_cell = lat_cell[inds]
            lon_cell = lon_cell[inds]
            alt_cell = alt_cell[inds]

            # compute number of flashes
            nflashes[i] = lat_cell.size

            # Plot altitude histogram
            fname_hist = (
                args.trtpath+cell_time_str+'_'+infostr +
                '_hist_alt_all_sources.png')
            titl = (
                cell_time_str+' '+infostr +
                ' TRT cell. Flash all sources altitude')
            fname_hist = plot_histogram(
                bin_edges, alt_cell, [fname_hist], labelx='Altitude [m MSL]',
                titl=titl)
            print('Plotted '+' '.join(fname_hist))

            fname_hist = (
                args.trtpath+cell_time_str+'_'+infostr +
                '_hist_alt_all_sources.csv')
            hist, bin_edges = np.histogram(alt_cell, bins=bin_edges)
            fname_hist = write_histogram(
                bin_edges, hist, fname_hist, step=step)
            print('----- save to '+fname_hist)

            # plot position first source
            figfname = (
                args.trtpath+cell_time_str+'_'+infostr +
                '_all_sources_pos_max_height_on_top.png')
            titl = (
                cell_time_str+' '+infostr +
                ' all flash sources position. Highest on top')
            figfname = plot_pos(
                lat_cell, lon_cell, alt_cell, [figfname],
                sort_altitude='Highest_on_top', cb_label='Source height [m MSL]',
                titl=titl)
            print('Plotted '+' '.join(figfname))

        # plot time series lightning
        figfname = (
            args.trtpath+cell_time_str+'_'+infostr+'_N_all_sources.png')
        titl = infostr+' Number of sources in cell'
        figfname = plot_timeseries(
            yyyymmddHHMM, [nflashes], [figfname], labelx='Time [UTC]',
            labely='N flashes', labels=['N flashes'], title=titl,
            period=0, timeformat=None, colors=None, linestyles=None,
            markers=None, ymin=None, ymax=None, dpi=72)
        print('Plotted '+' '.join(figfname))

        # plot time-hist_height
        flist = glob.glob(
            args.trtpath+'*_'+infostr+'_hist_alt_all_sources.csv')

        if not flist:
            warn('No histogram files found in '+args.trtpath +
                 ' for TRT cell '+infostr)
        else:
            tbin_edges, bin_edges, data_ma = read_histogram_ts(
                flist, 'flash_altitude')

            vmax = np.max(data_ma)
            if vmax == 0.:
                warn('Unable to plot histogram. No valid data')
            else:
                fname_hist = (
                    args.trtpath+'/'+infostr+'_trt_HISTOGRAM_alt_all_source.png')
                titl = 'TRT cell '+infostr+'\n'+'Altitude of all flash sources'
                _plot_time_range(
                    tbin_edges, bin_edges, data_ma, 'frequency_of_occurrence',
                    [fname_hist], titl=titl,
                    ylabel='Altitude [m MSL]',
                    vmin=0., vmax=np.max(data_ma), figsize=[10, 8], dpi=72)

                print("----- plot to '%s'" % fname_hist)


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
