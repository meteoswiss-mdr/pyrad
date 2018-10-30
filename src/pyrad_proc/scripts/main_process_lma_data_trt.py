#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_lma_data_trt
================================================

This program gets the LMA lightning data generated within a TRT cell and
provides outputs such as time-height plot of flashes/sources in each cell
and files containing the number of flashes/sources within the cell at each
time step

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
from pyrad.io import read_histogram_ts, write_trt_cell_lightning
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

    # positional arguments
    parser.add_argument(
        'days', nargs='+', type=str,
        help='Dates to process. Format YYYY-MM-DD')

    # keyword arguments
    parser.add_argument(
        '--nsources_min', type=int, default=10,
        help='Minimum number of sources to consider the LMA flash valid')

    parser.add_argument(
        '--trtbase', type=str,
        default='/store/msrad/radar/trt/',
        help='name of folder containing the TRT cell data')

    parser.add_argument(
        '--flashpath', type=str,
        default='/store/msrad/lightning/LMA/Santis/',
        help='name of the folder containing the Santis LMA data')

    args = parser.parse_args()

    print("====== LMA data TRT processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== LMA data TRT processing finished: ")

    print('trt path: '+args.trtbase)

    time_dir_list = args.days

    # Get bins altitude
    alt_min = 0.
    alt_max = 14000.
    step = 100.
    bin_edges = np.linspace(
        alt_min-step/2., alt_max+step/2, num=int((alt_max-alt_min)/step)+2)

    trt_list = []
    for time_dir in time_dir_list:
        trt_list.extend(glob.glob(
            args.trtbase+time_dir+'/TRTC_cell_plots/All/*.trt'))
        trt_list.extend(glob.glob(
            args.trtbase+time_dir+'/TRTC_cell_plots/Some/*.trt'))

    cell_ID_list = np.ma.asarray([], dtype=int)
    time_list = np.ma.asarray([], dtype=datetime.datetime)
    lon_list = np.ma.asarray([], dtype=float)
    lat_list = np.ma.asarray([], dtype=float)
    flash_density_list = np.ma.asarray([], dtype=float)
    sources_density_list = np.ma.asarray([], dtype=float)
    rank_flash_density_list = np.ma.asarray([], dtype=float)
    area_list = np.ma.asarray([], dtype=float)
    nflash_list = np.ma.asarray([], dtype=int)
    nsources_list = np.ma.asarray([], dtype=int)
    for trt_fname in trt_list:
        print('processing TRT cell file '+trt_fname)
        trtpath = os.path.dirname(trt_fname)+'/'

        # reading TRT cell file
        (traj_ID, yyyymmddHHMM, lon_trt, lat_trt, _, _, _, area, _, _, _,
         RANKr, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, cell_contour) = (
             read_trt_traj_data(trt_fname))

        # reading lightning file
        infostr = os.path.basename(trt_fname).split('.')[0]
        dt_str = infostr[0:12]
        dt = datetime.datetime.strptime(dt_str, '%Y%m%d%H%M')

        flashnr, time_data, _, lat, lon, alt, _ = read_lightning(
            args.flashpath+dt.strftime("%y%m%d")+'.txt')

        if flashnr is None:
            continue

        # Filter data with less than nsources_min sources
        unique_flashnr = np.unique(flashnr, return_index=False)

        ind = []
        for flash in unique_flashnr:
            ind_flash = np.where(flashnr == flash)[0]
            if ind_flash.size < args.nsources_min:
                continue
            ind.extend(ind_flash)

        if np.size(ind) == 0:
            continue

        flashnr = flashnr[ind]
        time_data = time_data[ind]
        lat = lat[ind]
        lon = lon[ind]
        alt = alt[ind]

        # Get first sources data
        _, unique_ind = np.unique(flashnr, return_index=True)
        lat_first = lat[unique_ind]
        lon_first = lon[unique_ind]
        alt_first = alt[unique_ind]
        # dBm_first = dBm[unique_ind]
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
                    trtpath+cell_time_str+'_'+infostr +
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
                    trtpath+cell_time_str+'_'+infostr +
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
                trtpath+cell_time_str+'_'+infostr +
                '_hist_alt_first_source.png')
            titl = (
                cell_time_str+' '+infostr +
                ' TRT cell. Flash first source altitude')
            fname_hist = plot_histogram(
                bin_edges, alt_cell, [fname_hist], labelx='Altitude [m MSL]',
                titl=titl)
            print('Plotted '+' '.join(fname_hist))

            fname_hist = (
                trtpath+cell_time_str+'_'+infostr +
                '_hist_alt_first_source.csv')
            hist, bin_edges = np.histogram(alt_cell, bins=bin_edges)
            fname_hist = write_histogram(
                bin_edges, hist, fname_hist, step=step)
            print('----- save to '+fname_hist)

            # plot position first source
            figfname = (
                trtpath+cell_time_str+'_'+infostr +
                '_first_source_pos_max_height_on_top.png')
            titl = (
                cell_time_str+' '+infostr +
                ' first flash source position. Highest on top')
            figfname = plot_pos(
                lat_cell, lon_cell, alt_cell, [figfname],
                sort_altitude='Highest_on_top',
                cb_label='Source height [m MSL]', titl=titl)
            print('Plotted '+' '.join(figfname))

        # plot time series lightning density
        figfname = (
            trtpath+cell_time_str+'_'+infostr+'_dens_first_sources.png')
        titl = infostr+' First sources density'
        figfname = plot_timeseries(
            yyyymmddHHMM, [nflashes/area], [figfname], labelx='Time [UTC]',
            labely='Flash dens [Flashes/Km2]', labels=['Flash density'],
            title=titl, period=0, timeformat=None, colors=None,
            linestyles=None, markers=None, ymin=None, ymax=None, dpi=72)
        print('Plotted '+' '.join(figfname))

        # plot time series lightning
        figfname = (
            trtpath+cell_time_str+'_'+infostr+'_N_first_sources.png')
        titl = infostr+' Number of first sources in cell'
        figfname = plot_timeseries(
            yyyymmddHHMM, [nflashes], [figfname], labelx='Time [UTC]',
            labely='N flashes', labels=['N flashes'], title=titl,
            period=0, timeformat=None, colors=None, linestyles=None,
            markers=None, ymin=None, ymax=None, dpi=72)
        print('Plotted '+' '.join(figfname))

        # plot time-hist_height
        flist = glob.glob(
            trtpath+'*_'+infostr+'_hist_alt_first_source.csv')

        if not flist:
            warn('No histogram files found in '+trtpath +
                 ' for TRT cell '+infostr)
        else:
            tbin_edges, bin_edges, data_ma, _ = read_histogram_ts(
                flist, 'flash_altitude')

            vmax = np.max(data_ma)
            if vmax == 0.:
                warn('Unable to plot histogram. No valid data')
            else:
                data_ma[data_ma == 0.] = np.ma.masked
                fname_hist = (
                    trtpath+infostr+'_trt_HISTOGRAM_alt_first_source.png')
                titl = ('TRT cell '+infostr+'\n' +
                        'Altitude of first flash source')
                _plot_time_range(
                    tbin_edges, bin_edges, data_ma, 'frequency_of_occurrence',
                    [fname_hist], titl=titl,
                    ylabel='Altitude [m MSL]',
                    vmin=0., vmax=vmax, figsize=[10, 8], dpi=72)

                print("----- plot to '%s'" % fname_hist)

        # Append flash data
        cell_ID_list = np.ma.append(cell_ID_list, traj_ID)
        time_list = np.ma.append(time_list, yyyymmddHHMM)
        lon_list = np.ma.append(lon_list, lon_trt)
        lat_list = np.ma.append(lat_list, lat_trt)
        flash_density_list = np.ma.append(flash_density_list, nflashes/area)
        rank_flash_density_list = np.ma.append(
            rank_flash_density_list, RANKr)
        area_list = np.ma.append(area_list, area)
        nflash_list = np.ma.append(nflash_list, nflashes)

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
                    trtpath+cell_time_str+'_'+infostr +
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
                    trtpath+cell_time_str+'_'+infostr +
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
                trtpath+cell_time_str+'_'+infostr +
                '_hist_alt_all_sources.png')
            titl = (
                cell_time_str+' '+infostr +
                ' TRT cell. Flash all sources altitude')
            fname_hist = plot_histogram(
                bin_edges, alt_cell, [fname_hist], labelx='Altitude [m MSL]',
                titl=titl)
            print('Plotted '+' '.join(fname_hist))

            fname_hist = (
                trtpath+cell_time_str+'_'+infostr +
                '_hist_alt_all_sources.csv')
            hist, bin_edges = np.histogram(alt_cell, bins=bin_edges)
            fname_hist = write_histogram(
                bin_edges, hist, fname_hist, step=step)
            print('----- save to '+fname_hist)

            # plot position first source
            figfname = (
                trtpath+cell_time_str+'_'+infostr +
                '_all_sources_pos_max_height_on_top.png')
            titl = (
                cell_time_str+' '+infostr +
                ' all flash sources position. Highest on top')
            figfname = plot_pos(
                lat_cell, lon_cell, alt_cell, [figfname],
                sort_altitude='Highest_on_top',
                cb_label='Source height [m MSL]', titl=titl)
            print('Plotted '+' '.join(figfname))

        # plot time series lightning
        figfname = (
            trtpath+cell_time_str+'_'+infostr+'_N_all_sources.png')
        titl = infostr+' Number of sources in cell'
        figfname = plot_timeseries(
            yyyymmddHHMM, [nflashes], [figfname], labelx='Time [UTC]',
            labely='N flashes', labels=['N flashes'], title=titl,
            period=0, timeformat=None, colors=None, linestyles=None,
            markers=None, ymin=None, ymax=None, dpi=72)
        print('Plotted '+' '.join(figfname))

        # plot time series lightning density
        figfname = (
            trtpath+cell_time_str+'_'+infostr+'_dens_all_sources.png')
        titl = infostr+' Sources density'
        figfname = plot_timeseries(
            yyyymmddHHMM, [nflashes/area], [figfname], labelx='Time [UTC]',
            labely='Source dens [Sources/Km2]', labels=['Source density'],
            title=titl, period=0, timeformat=None, colors=None,
            linestyles=None, markers=None, ymin=None, ymax=None, dpi=72)
        print('Plotted '+' '.join(figfname))

        # plot time-hist_height
        flist = glob.glob(
            trtpath+'*_'+infostr+'_hist_alt_all_sources.csv')

        if not flist:
            warn('No histogram files found in '+trtpath +
                 ' for TRT cell '+infostr)
        else:
            tbin_edges, bin_edges, data_ma, _ = read_histogram_ts(
                flist, 'flash_altitude')

            vmax = np.max(data_ma)
            if vmax == 0.:
                warn('Unable to plot histogram. No valid data')
            else:
                data_ma[data_ma == 0.] = np.ma.masked
                fname_hist = (
                    trtpath+'/'+infostr+'_trt_HISTOGRAM_alt_all_source.png')
                titl = (
                    'TRT cell '+infostr+'\n'+'Altitude of all flash sources')
                _plot_time_range(
                    tbin_edges, bin_edges, data_ma, 'frequency_of_occurrence',
                    [fname_hist], titl=titl,
                    ylabel='Altitude [m MSL]',
                    vmin=0., vmax=np.max(data_ma), figsize=[10, 8], dpi=72)

                print("----- plot to '%s'" % fname_hist)

        # Append sources data
        sources_density_list = np.ma.append(sources_density_list, nflashes/area)
        nsources_list = np.ma.append(nsources_list, nflashes)

    fname = args.trtbase+'cell_LMA_flashes.csv'
    write_trt_cell_lightning(
        cell_ID_list, time_list, lon_list, lat_list, area_list,
        rank_flash_density_list, nflash_list, flash_density_list, fname)

    fname = args.trtbase+'cell_LMA_sources.csv'
    write_trt_cell_lightning(
        cell_ID_list, time_list, lon_list, lat_list, area_list,
        rank_flash_density_list, nsources_list, sources_density_list, fname)


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
