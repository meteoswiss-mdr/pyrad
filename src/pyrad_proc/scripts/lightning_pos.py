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
import atexit

import numpy as np

from pyrad.io import read_lightning
from pyrad.graph import plot_histogram, plot_pos

print(__doc__)


def main():
    """
    """
    # basepath = '/data/lightning/Santis/'
    basepath = '/store/msrad/lightning/LMA/Santis/'
    day_vec = [
        datetime.datetime(2017, 7, 19),
        datetime.datetime(2017, 7, 30),
        datetime.datetime(2017, 8, 1)]


    print("====== Lightning Position plotting started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Lightning Position plotting finished: ")

    flashnr = np.asarray([], dtype=int)
    time_data = np.asarray([], dtype=datetime.datetime)
    time_in_flash = np.asarray([], dtype=float)
    lat = np.asarray([], dtype=float)
    lon = np.asarray([], dtype=float)
    alt = np.asarray([], dtype=float)
    dBm = np.asarray([], dtype=float)

    flash_cnt = 0
    for i, day in enumerate(day_vec):
        day_str = day.strftime('%y%m%d')
        fname = basepath+day_str+'.txt'

        print('Reading LMA data file '+fname)
        (flashnr_aux, time_data_aux, time_in_flash_aux, lat_aux, lon_aux,
         alt_aux, dBm_aux) = read_lightning(fname)

        flashnr = np.append(flashnr, flashnr_aux+flash_cnt)
        flash_cnt += flashnr_aux.max()
        # time_data = np.append(time_data, time_data_aux)
        # time_in_flash = np.append(time_in_flash, time_in_flash_aux)
        lat = np.append(lat, lat_aux)
        lon = np.append(lon, lon_aux)
        alt = np.append(alt, alt_aux)
        dBm = np.append(dBm, dBm_aux)

        # Get first sources data
        flashnr_first, unique_ind = np.unique(flashnr_aux, return_index=True)
        lat_first = lat_aux[unique_ind]
        lon_first = lon_aux[unique_ind]
        alt_first = alt_aux[unique_ind]
        dBm_first = dBm_aux[unique_ind]

        # Get bins altitude
        alt_min = alt_aux.min()
        alt_max = alt_aux.max()
        step = 100.
        bins = np.linspace(alt_min-step/2., alt_max+step/2., num=int((alt_max-alt_min)/step)+2)

        fname_hist = basepath+day_str+'_hist_alt.png'
        fname_hist = plot_histogram(
            bins, alt_aux, [fname_hist], labelx='Altitude [m MSL]',
            titl=day_str+' Flash sources altitude')
        print('Plotted '+' '.join(fname_hist))

        fname_hist = basepath+day_str+'_hist_alt_first_source.png'
        fname_hist = plot_histogram(
            bins, alt_first, [fname_hist], labelx='Altitude [m MSL]',
            titl=day_str+' Flash first source altitude')
        print('Plotted '+' '.join(fname_hist))

        # Get bins dBm
        dBm_min = dBm_aux.min()
        dBm_max = dBm_aux.max()
        step = 1.
        bins = np.linspace(dBm_min-step/2., dBm_max+step/2., num=int((dBm_max-dBm_min)/step)+2)

        fname_hist = basepath+day_str+'_hist_dBm.png'
        fname_hist = plot_histogram(
            bins, dBm_aux, [fname_hist], labelx='Power [dBm]',
            titl=day_str+' Flash sources power')
        print('Plotted '+' '.join(fname_hist))


        fname_hist = basepath+day_str+'_hist_dBm_first_source.png'
        fname_hist = plot_histogram(
            bins, dBm_first, [fname_hist], labelx='Power [dBm]',
            titl=day_str+' Flash first source power')
        print('Plotted '+' '.join(fname_hist))

        print('N sources: '+str(alt_aux.size))

        # plot position all sources
        figfname = basepath+day_str+'_sources_pos_max_height_on_top.png'
        figfname = plot_pos(
            lat_aux, lon_aux, alt_aux, [figfname],
            sort_altitude='Highest_on_top', cb_label='Source height [m MSL]',
            titl=day_str+' flash sources position. Highest on top')
        print('Plotted '+' '.join(figfname))

        figfname = basepath+day_str+'_sources_pos_min_height_on_top.png'
        figfname = plot_pos(
            lat_aux, lon_aux, alt_aux, [figfname],
            sort_altitude='Lowest_on_top', cb_label='Source height [m MSL]',
            titl=day_str+' flash sources position. Lowest on top')
        print('Plotted '+' '.join(figfname))

        # plot position first source
        print('N flashes: '+str(alt_first.size))

        # plot position first source
        figfname = basepath+day_str+'_first_source_pos_max_height_on_top.png'
        figfname = plot_pos(
            lat_first, lon_first, alt_first, [figfname],
            sort_altitude='Highest_on_top', cb_label='Source height [m MSL]',
            titl=day_str+' first flash source position. Highest on top')
        print('Plotted '+' '.join(figfname))

        figfname = basepath+day_str+'_first_source_pos_min_height_on_top.png'
        figfname = plot_pos(
            lat_first, lon_first, alt_first, [figfname],
            sort_altitude='Lowest_on_top', cb_label='Source height [m MSL]',
            titl=day_str+' first flash source position. Lowest on top')
        print('Plotted '+' '.join(figfname))

    print('N sources total: '+str(alt.size))

    # plot position all sources
    figfname = basepath+'Santis_LMA_sources_pos_max_height_on_top.png'
    figfname = plot_pos(
        lat, lon, alt, [figfname], sort_altitude='Highest_on_top',
        cb_label='Source height [m MSL]',
        titl='Flash sources position. Highest on top')
    print('Plotted '+' '.join(figfname))

    figfname = basepath+'Santis_LMA_sources_pos_min_height_on_top.png'
    figfname = plot_pos(
        lat, lon, alt, [figfname], sort_altitude='Lowest_on_top',
        cb_label='Source height [m MSL]',
        titl='Flash sources position. Lowest on top')
    print('Plotted '+' '.join(figfname))

    # plot position first source
    flashnr_first, unique_ind = np.unique(flashnr, return_index=True)
    lat_first = lat[unique_ind]
    lon_first = lon[unique_ind]
    alt_first = alt[unique_ind]

    print('N flashes total: '+str(alt_first.size))

    # plot position all sources
    figfname = basepath+'Santis_LMA_first_source_pos_max_height_on_top.png'
    figfname = plot_pos(
        lat_first, lon_first, alt_first, [figfname],
        sort_altitude='Highest_on_top', cb_label='Source height [m MSL]',
        titl='First flash source position. Highest on top')
    print('Plotted '+' '.join(figfname))

    figfname = basepath+'Santis_LMA_first_source_pos_min_height_on_top.png'
    figfname = plot_pos(
        lat_first, lon_first, alt_first, [figfname],
        sort_altitude='Lowest_on_top', cb_label='Source height [m MSL]',
        titl='First flash source position. Lowest on top')
    print('Plotted '+' '.join(figfname))


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
