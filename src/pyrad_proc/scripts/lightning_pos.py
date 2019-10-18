#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
lightning_pos
================================================

This program plots the position of LMA lightning data
and several data histograms

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit

import numpy as np

from pyrad.io import read_lightning, write_histogram
from pyrad.graph import plot_histogram, plot_pos

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

    args = parser.parse_args()

    print("====== Lightning Position plotting started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Lightning Position plotting finished: ")

    day_vec = []
    for day in enumerate(args.days):
        day_vec.append(datetime.datetime.strptime(day, '%Y%m%d'))

    flashnr = np.asarray([], dtype=int)
    # time_data = np.asarray([], dtype=datetime.datetime)
    # time_in_flash = np.asarray([], dtype=float)
    lat = np.asarray([], dtype=float)
    lon = np.asarray([], dtype=float)
    alt = np.asarray([], dtype=float)
    dBm = np.asarray([], dtype=float)

    flash_cnt = 0
    for day in day_vec:
        day_str = day.strftime('%y%m%d')
        fname = args.basepath+day_str+'.txt'

        print('Reading LMA data file '+fname)
        (flashnr_aux, _, _, lat_aux, lon_aux,
         alt_aux, dBm_aux) = read_lightning(fname)

        print('N sources: '+str(flashnr_aux.size))

        flashnr = np.append(flashnr, flashnr_aux+flash_cnt)
        flash_cnt += flashnr_aux.max()
        # time_data = np.append(time_data, time_data_aux)
        # time_in_flash = np.append(time_in_flash, time_in_flash_aux)
        lat = np.append(lat, lat_aux)
        lon = np.append(lon, lon_aux)
        alt = np.append(alt, alt_aux)
        dBm = np.append(dBm, dBm_aux)

        # Get first sources data
        _, unique_ind = np.unique(flashnr_aux, return_index=True)
        lat_first = lat_aux[unique_ind]
        lon_first = lon_aux[unique_ind]
        alt_first = alt_aux[unique_ind]
        dBm_first = dBm_aux[unique_ind]

        print('N flashes: '+str(alt_first.size))

        # Get bins altitude
        alt_min = alt_aux.min()
        alt_max = alt_aux.max()
        step = 100.
        bins = np.linspace(
            alt_min-step/2., alt_max+step/2.,
            num=int((alt_max-alt_min)/step)+2)

        # Plot histogram altitude all sources
        fname_hist = args.basepath+day_str+'_hist_alt.png'
        fname_hist = plot_histogram(
            bins, alt_aux, [fname_hist], labelx='Altitude [m MSL]',
            titl=day_str+' Flash sources altitude')
        print('Plotted '+' '.join(fname_hist))

        # Plot histogram altitude first sources
        fname_hist = args.basepath+day_str+'_hist_alt_first_source.png'
        fname_hist = plot_histogram(
            bins, alt_first, [fname_hist], labelx='Altitude [m MSL]',
            titl=day_str+' Flash first source altitude')
        print('Plotted '+' '.join(fname_hist))

        # Get bins dBm
        dBm_min = dBm_aux.min()
        dBm_max = dBm_aux.max()
        step = 1.
        bins = np.linspace(
            dBm_min-step/2., dBm_max+step/2.,
            num=int((dBm_max-dBm_min)/step)+2)

        # Plot histogram power all sources
        fname_hist = args.basepath+day_str+'_hist_dBm.png'
        fname_hist = plot_histogram(
            bins, dBm_aux, [fname_hist], labelx='Power [dBm]',
            titl=day_str+' Flash sources power')
        print('Plotted '+' '.join(fname_hist))

        # Plot histogram power first sources
        fname_hist = args.basepath+day_str+'_hist_dBm_first_source.png'
        fname_hist = plot_histogram(
            bins, dBm_first, [fname_hist], labelx='Power [dBm]',
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

    print('N sources total: '+str(alt.size))

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

    # Get bins altitude
    alt_min = alt.min()
    alt_max = alt.max()
    step = 100.
    bins = np.linspace(
        alt_min-step/2., alt_max+step/2., num=int((alt_max-alt_min)/step)+2)

    # Plot histogram altitude all sources
    fname_hist = args.basepath+'Santis_hist_alt.png'
    fname_hist = plot_histogram(
        bins, alt, [fname_hist], labelx='Altitude [m MSL]',
        titl='Flash sources altitude')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.basepath+'Santis_hist_alt.csv'
    hist_alt, _ = np.histogram(alt, bins=bins)
    fname_hist = write_histogram(bins, hist_alt, fname_hist)
    print('Written '+' '.join(fname_hist))

    # Plot histogram altitude first sources
    fname_hist = args.basepath+'Santis_hist_alt_first_source.png'
    fname_hist = plot_histogram(
        bins, alt_first, [fname_hist], labelx='Altitude [m MSL]',
        titl='Flash first source altitude')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.basepath+'Santis_hist_alt_fist_source.csv'
    hist_alt_fist, _ = np.histogram(alt_first, bins=bins)
    fname_hist = write_histogram(bins, hist_alt_fist, fname_hist)
    print('Written '+' '.join(fname_hist))

    # Get bins dBm
    dBm_min = dBm.min()
    dBm_max = dBm.max()
    step = 1.
    bins = np.linspace(
        dBm_min-step/2., dBm_max+step/2., num=int((dBm_max-dBm_min)/step)+2)

    fname_hist = args.basepath+'Santis_hist_dBm.png'
    fname_hist = plot_histogram(
        bins, dBm, [fname_hist], labelx='Power [dBm]',
        titl='Flash sources power')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.basepath+'Santis_hist_dBm.csv'
    hist_dBm, _ = np.histogram(dBm, bins=bins)
    fname_hist = write_histogram(bins, hist_dBm, fname_hist)
    print('Written '+' '.join(fname_hist))

    # Plot histogram power first sources
    fname_hist = args.basepath+'Santis_hist_dBm_first_source.png'
    fname_hist = plot_histogram(
        bins, dBm_first, [fname_hist], labelx='Power [dBm]',
        titl='Flash first source power')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.basepath+'Santis_hist_dBm_first_source.csv'
    hist_dBm_first, _ = np.histogram(dBm_first, bins=bins)
    fname_hist = write_histogram(bins, hist_dBm_first, fname_hist)
    print('Written '+' '.join(fname_hist))


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
