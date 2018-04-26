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

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

from pyrad.io import read_lightning
from pyrad.graph import plot_histogram

print(__doc__)


def main():
    """
    """
    basepath = '/data/lightning/Santis/'
    # basepath = '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/'
    day_vec = [
        datetime.datetime(2017, 6, 29),
        datetime.datetime(2017, 6, 30),
        datetime.datetime(2017, 7, 10),
        datetime.datetime(2017, 7, 14),
        datetime.datetime(2017, 7, 18)]


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
    for i, day in enumerate(day_vec):
        day_str = day.strftime('%y%m%d')
        fname = basepath+day_str+'.txt'

        print('Reading LMA data file '+fname)
        (flashnr_aux, time_data_aux, time_in_flash_aux, lat_aux, lon_aux,
         alt_aux, dBm_aux) = read_lightning(fname)

        # flashnr = np.append(flashnr, flashnr_aux)
        # time_data = np.append(time_data, time_data_aux)
        # time_in_flash = np.append(time_in_flash, time_in_flash_aux)
        lat = np.append(lat, lat_aux)
        lon = np.append(lon, lon_aux)
        alt = np.append(alt, alt_aux)
        dBm = np.append(dBm, dBm_aux)

        # Get first sources indices
        flashnr_first, unique_ind = np.unique(flashnr_aux, return_index=True)

        # Get bins altitude
        alt_min = alt_aux.min()
        alt_max = alt_aux.max()
        step = 100.
        bins = np.linspace(alt_min, alt_max, num=int((alt_max-alt_min)/step))

        fname_hist = basepath+day_str+'_hist_alt.png'
        plot_histogram(bins, alt_aux, [fname_hist], labelx='Altitude [m MSL]',
                       titl=day_str+' Flash sources altitude')
        print('Plotted '+fname_hist)

        alt_first = alt_aux[unique_ind]
        fname_hist = basepath+day_str+'_hist_alt_first_source.png'
        plot_histogram(bins, alt_first, [fname_hist], labelx='Altitude [m MSL]',
                       titl=day_str+' Flash first source altitude')
        print('Plotted '+fname_hist)

        # Get bins dBm
        dBm_min = dBm_aux.min()
        dBm_max = dBm_aux.max()
        step = 1.
        bins = np.linspace(dBm_min, dBm_max, num=int((dBm_max-dBm_min)/step))

        fname_hist = basepath+day_str+'_hist_dBm.png'
        plot_histogram(bins, dBm_aux, [fname_hist], labelx='Power [dBm]',
                       titl=day_str+' Flash sources power')
        print('Plotted '+fname_hist)

        dBm_first = dBm_aux[unique_ind]
        fname_hist = basepath+day_str+'_hist_dBm_first_source.png'
        plot_histogram(bins, dBm_first, [fname_hist], labelx='Power [dBm]',
                       titl=day_str+' Flash first source power')
        print('Plotted '+fname_hist)

    # plot position
    figfname = basepath+'Santis_LMA_flashes_pos.png'

    # sort data by altitude
    ind = np.argsort(alt)
    # inverse indices
    # ind = ind[::-1]
    lat = lat[ind]
    lon = lon[ind]
    alt = alt[ind]

    print('N flashes: '+str(alt.size))

    marker = 'x'
    col = alt
    cmap = 'viridis'
    norm = plt.Normalize(alt.min(), alt.max())
    cb_label = 'Flash source height [m MSL]'

    fig = plt.figure(figsize=[10, 8], dpi=72)
    ax = fig.add_subplot(111, aspect='equal')

    cax = ax.scatter(
        lon, lat, c=col, marker=marker, alpha=1, cmap=cmap, norm=norm)
    cb = fig.colorbar(cax, orientation='horizontal')
    cb.set_label(cb_label)

    plt.title('Flash sources position')
    plt.xlabel('Lon [Deg]')
    plt.ylabel('Lat [Deg]')

    # Turn on the grid
    ax.grid()

    fig.savefig(figfname, dpi=72)
    plt.close()
    print('Plotted '+figfname)


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
