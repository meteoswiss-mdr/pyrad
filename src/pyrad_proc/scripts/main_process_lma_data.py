#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_lma_data
================================================

This program reads LMA raw data and plots several features such as sources
position, histograms, etc.

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit

import numpy as np

from pyrad.io import read_lightning, write_histogram
from pyrad.graph import plot_histogram, plot_pos, _plot_time_range

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
        help='Dates to process. Format YYYYMMDD')

    # keyword arguments
    parser.add_argument(
        '--basepath', type=str,
        default='/store/msrad/lightning/LMA/Santis/',
        help='name of folder containing the LMA lightning data')

    parser.add_argument(
        '--nsources_min', type=int, default=10,
        help='Minimum number of sources to consider the LMA flash valid')

    args = parser.parse_args()

    print("====== LMA data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== LMA data processing finished: ")

    day_vec = []
    for day in args.days:
        day_vec.append(datetime.datetime.strptime(day, '%Y%m%d'))

    flashnr = np.asarray([], dtype=int)
    # time_data = np.asarray([], dtype=datetime.datetime)
    # time_in_flash = np.asarray([], dtype=float)
    lat = np.asarray([], dtype=float)
    lon = np.asarray([], dtype=float)
    alt = np.asarray([], dtype=float)
    dBm = np.asarray([], dtype=float)

    # define histogram bin edges
    bin_edges_alt = np.arange(-50., 14150., 100.)
    bin_edges_dBm = np.arange(-17., 47., 1.)

    flash_cnt = 0
    for day in day_vec:
        day_str = day.strftime('%y%m%d')
        fname = args.basepath+day_str+'.txt'

        print('Reading LMA data file '+fname)
        (flashnr_aux, _, _, lat_aux, lon_aux,
         alt_aux, dBm_aux) = read_lightning(fname)

        # Filter data with less than nsources_min sources
        unique_flashnr = np.unique(flashnr_aux, return_index=False)

        ind = []
        for flash in unique_flashnr:
            ind_flash = np.where(flashnr_aux == flash)[0]
            if ind_flash.size < args.nsources_min:
                continue
            ind.extend(ind_flash)

        flashnr_aux = flashnr_aux[ind]
        lat_aux = lat_aux[ind]
        lon_aux = lon_aux[ind]
        alt_aux = alt_aux[ind]
        dBm_aux = dBm_aux[ind]

        # get first flash
        _, unique_ind = np.unique(flashnr_aux, return_index=True)
        lat_first = lat_aux[unique_ind]
        lon_first = lon_aux[unique_ind]
        alt_first = alt_aux[unique_ind]
        dBm_first = dBm_aux[unique_ind]

        print('N sources: '+str(flashnr_aux.size))
        print('Min alt: '+str(alt_aux.min()))
        print('Max alt: '+str(alt_aux.max()))
        print('Min power: '+str(dBm_aux.min()))
        print('Max power: '+str(dBm_aux.max()))
        print('N flashes: '+str(alt_first.size))

        # put data in global list
        flashnr = np.append(flashnr, flashnr_aux+flash_cnt)
        flash_cnt += flashnr_aux.max()
        # time_data = np.append(time_data, time_data_aux)
        # time_in_flash = np.append(time_in_flash, time_in_flash_aux)
        lat = np.append(lat, lat_aux)
        lon = np.append(lon, lon_aux)
        alt = np.append(alt, alt_aux)
        dBm = np.append(dBm, dBm_aux)

        # Plot histogram altitude all sources
        fname_hist = args.basepath+day_str+'_hist_alt.png'
        fname_hist = plot_histogram(
            bin_edges_alt, alt_aux, [fname_hist], labelx='Altitude [m MSL]',
            titl=day_str+' Flash sources altitude')
        print('Plotted '+' '.join(fname_hist))

        # Plot histogram altitude first sources
        fname_hist = args.basepath+day_str+'_hist_alt_first_source.png'
        fname_hist = plot_histogram(
            bin_edges_alt, alt_first, [fname_hist], labelx='Altitude [m MSL]',
            titl=day_str+' Flash first source altitude')
        print('Plotted '+' '.join(fname_hist))

        # Plot histogram power all sources
        fname_hist = args.basepath+day_str+'_hist_dBm.png'
        fname_hist = plot_histogram(
            bin_edges_dBm, dBm_aux, [fname_hist], labelx='Power [dBm]',
            titl=day_str+' Flash sources power')
        print('Plotted '+' '.join(fname_hist))

        # Plot histogram power first sources
        fname_hist = args.basepath+day_str+'_hist_dBm_first_source.png'
        fname_hist = plot_histogram(
            bin_edges_dBm, dBm_first, [fname_hist], labelx='Power [dBm]',
            titl=day_str+' Flash first source power')
        print('Plotted '+' '.join(fname_hist))

        # plot position all sources
        figfname = args.basepath+day_str+'_sources_pos_max_height_on_top.png'
        figfname = plot_pos(
            lat_aux, lon_aux, alt_aux, [figfname],
            sort_altitude='Highest_on_top', cb_label='Source height [m MSL]',
            titl=day_str+' flash sources position. Highest on top')
        print('Plotted '+' '.join(figfname))

        figfname = args.basepath+day_str+'_sources_pos_min_height_on_top.png'
        figfname = plot_pos(
            lat_aux, lon_aux, alt_aux, [figfname],
            sort_altitude='Lowest_on_top', cb_label='Source height [m MSL]',
            titl=day_str+' flash sources position. Lowest on top')
        print('Plotted '+' '.join(figfname))

        # plot position first source
        figfname = (
            args.basepath+day_str+'_first_source_pos_max_height_on_top.png')
        figfname = plot_pos(
            lat_first, lon_first, alt_first, [figfname],
            sort_altitude='Highest_on_top', cb_label='Source height [m MSL]',
            titl=day_str+' first flash source position. Highest on top')
        print('Plotted '+' '.join(figfname))

        figfname = (
            args.basepath+day_str+'_first_source_pos_min_height_on_top.png')
        figfname = plot_pos(
            lat_first, lon_first, alt_first, [figfname],
            sort_altitude='Lowest_on_top', cb_label='Source height [m MSL]',
            titl=day_str+' first flash source position. Lowest on top')
        print('Plotted '+' '.join(figfname))

        # Plot 2D histogram all sources
        H, _, _ = np.histogram2d(
            alt_aux, dBm_aux, bins=[bin_edges_alt, bin_edges_dBm])

        # set 0 values to blank
        H = np.ma.asarray(H)
        H[H == 0] = np.ma.masked

        fname_hist = args.basepath+day_str+'_Santis_2Dhist_alt_dBm.png'
        fname_hist = _plot_time_range(
            bin_edges_alt, bin_edges_dBm, H, None, [fname_hist],
            titl='LMA sources Altitude-Power histogram',
            xlabel='Altitude [m MSL]', ylabel='Power [dBm]',
            clabel='Occurrence',
            vmin=0, vmax=None, figsize=[10, 8], dpi=72)
        print('Plotted '+' '.join(fname_hist))

        # Plot 2D histogram first sources
        H, _, _ = np.histogram2d(
            alt_first, dBm_first, bins=[bin_edges_alt, bin_edges_dBm])

        # set 0 values to blank
        H = np.ma.asarray(H)
        H[H == 0] = np.ma.masked

        fname_hist = (
            args.basepath+day_str+'_Santis_2Dhist_alt_dBm_first_source.png')
        fname_hist = _plot_time_range(
            bin_edges_alt, bin_edges_dBm, H, None, [fname_hist],
            titl='LMA first sources Altitude-Power histogram',
            xlabel='Altitude [m MSL]', ylabel='Power [dBm]',
            clabel='Occurrence',
            vmin=0, vmax=None, figsize=[10, 8], dpi=72)
        print('Plotted '+' '.join(fname_hist))

    print('N sources total: '+str(alt.size))
    print('Min alt: '+str(alt.min()))
    print('Max alt: '+str(alt.max()))
    print('Min power: '+str(dBm.min()))
    print('Max power: '+str(dBm.max()))

    # Get first sources
    _, unique_ind = np.unique(flashnr, return_index=True)
    lat_first = lat[unique_ind]
    lon_first = lon[unique_ind]
    alt_first = alt[unique_ind]
    dBm_first = dBm[unique_ind]

    # plot position all sources
    figfname = args.basepath+'Santis_LMA_sources_pos_max_height_on_top.png'
    figfname = plot_pos(
        lat, lon, alt, [figfname], sort_altitude='Highest_on_top',
        cb_label='Source height [m MSL]',
        titl='Flash sources position. Highest on top')
    print('Plotted '+' '.join(figfname))

    figfname = args.basepath+'Santis_LMA_sources_pos_min_height_on_top.png'
    figfname = plot_pos(
        lat, lon, alt, [figfname], sort_altitude='Lowest_on_top',
        cb_label='Source height [m MSL]',
        titl='Flash sources position. Lowest on top')
    print('Plotted '+' '.join(figfname))

    # plot position first source
    print('N flashes total: '+str(alt_first.size))

    figfname = (
        args.basepath+'Santis_LMA_first_source_pos_max_height_on_top.png')
    figfname = plot_pos(
        lat_first, lon_first, alt_first, [figfname],
        sort_altitude='Highest_on_top', cb_label='Source height [m MSL]',
        titl='First flash source position. Highest on top')
    print('Plotted '+' '.join(figfname))

    figfname = (
        args.basepath+'Santis_LMA_first_source_pos_min_height_on_top.png')
    figfname = plot_pos(
        lat_first, lon_first, alt_first, [figfname],
        sort_altitude='Lowest_on_top', cb_label='Source height [m MSL]',
        titl='First flash source position. Lowest on top')
    print('Plotted '+' '.join(figfname))

    # Plot histogram altitude all sources
    fname_hist = args.basepath+'Santis_hist_alt.png'
    fname_hist = plot_histogram(
        bin_edges_alt, alt, [fname_hist], labelx='Altitude [m MSL]',
        titl='Flash sources altitude')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.basepath+'Santis_hist_alt.csv'
    hist_alt, _ = np.histogram(alt, bins=bin_edges_alt)
    fname_hist = write_histogram(bin_edges_alt, hist_alt, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram altitude first sources
    fname_hist = args.basepath+'Santis_hist_alt_first_source.png'
    fname_hist = plot_histogram(
        bin_edges_alt, alt_first, [fname_hist], labelx='Altitude [m MSL]',
        titl='Flash first source altitude')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.basepath+'Santis_hist_alt_first_source.csv'
    hist_alt_first, _ = np.histogram(alt_first, bins=bin_edges_alt)
    fname_hist = write_histogram(bin_edges_alt, hist_alt_first, fname_hist)
    print('Written '+fname_hist)

    fname_hist = args.basepath+'Santis_hist_dBm.png'
    fname_hist = plot_histogram(
        bin_edges_dBm, dBm, [fname_hist], labelx='Power [dBm]',
        titl='Flash sources power')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.basepath+'Santis_hist_dBm.csv'
    hist_dBm, _ = np.histogram(dBm, bins=bin_edges_dBm)
    fname_hist = write_histogram(bin_edges_dBm, hist_dBm, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram power first sources
    fname_hist = args.basepath+'Santis_hist_dBm_first_source.png'
    fname_hist = plot_histogram(
        bin_edges_dBm, dBm_first, [fname_hist], labelx='Power [dBm]',
        titl='Flash first source power')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.basepath+'Santis_hist_dBm_first_source.csv'
    hist_dBm_first, _ = np.histogram(dBm_first, bins=bin_edges_dBm)
    fname_hist = write_histogram(bin_edges_dBm, hist_dBm_first, fname_hist)
    print('Written '+fname_hist)

    # Plot 2D histogram all sources
    H, _, _ = np.histogram2d(alt, dBm, bins=[bin_edges_alt, bin_edges_dBm])

    # set 0 values to blank
    H = np.ma.asarray(H)
    H[H == 0] = np.ma.masked

    fname_hist = args.basepath+'Santis_2Dhist_alt_dBm.png'
    fname_hist = _plot_time_range(
        bin_edges_alt, bin_edges_dBm, H, None, [fname_hist],
        titl='LMA sources Altitude-Power histogram',
        xlabel='Altitude [m MSL]', ylabel='Power [dBm]',
        clabel='Occurrence',
        vmin=0, vmax=None, figsize=[10, 8], dpi=72)
    print('Plotted '+' '.join(fname_hist))

    # Plot 2D histogram first sources
    H, _, _ = np.histogram2d(
        alt_first, dBm_first, bins=[bin_edges_alt, bin_edges_dBm])

    # set 0 values to blank
    H = np.ma.asarray(H)
    H[H == 0] = np.ma.masked

    fname_hist = args.basepath+'Santis_2Dhist_alt_dBm_first_source.png'
    fname_hist = _plot_time_range(
        bin_edges_alt, bin_edges_dBm, H, None, [fname_hist],
        titl='LMA first sources Altitude-Power histogram',
        xlabel='Altitude [m MSL]', ylabel='Power [dBm]',
        clabel='Occurrence',
        vmin=0, vmax=None, figsize=[10, 8], dpi=72)
    print('Plotted '+' '.join(fname_hist))


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
