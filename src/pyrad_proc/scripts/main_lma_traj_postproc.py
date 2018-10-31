#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
lightning_traj_postproc
================================================

This program post-processes lightning trajectories

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
from copy import deepcopy
from warnings import warn

import numpy as np

from pyrad.io import read_lightning_all, get_fieldname_pyart, write_histogram
from pyrad.io import write_ts_lightning
from pyrad.graph import plot_histogram, plot_pos, get_colobar_label
from pyrad.graph import _plot_time_range, plot_histogram2
from pyrad.util import compute_histogram

from pyart.config import get_metadata


print(__doc__)


def main():
    """
    """
    # basepath = '/data/pyrad_products/rad4alp_hydro_PHA/'
    basepath = '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/'
    day_vec = [
        datetime.datetime(2017, 6, 29),
        datetime.datetime(2017, 6, 30),
        datetime.datetime(2017, 7, 10),
        datetime.datetime(2017, 7, 14),
        datetime.datetime(2017, 7, 18),
        datetime.datetime(2017, 7, 19),
        datetime.datetime(2017, 7, 30),
        datetime.datetime(2017, 8, 1)]

#    day_vec = [
#        datetime.datetime(2017, 7, 14)]
#
#    pol_vals_labels = [
#        'hydro [-]', 'KDPc [deg/Km]', 'dBZc [dBZ]', 'RhoHVc [-]',
#        'TEMP [deg C]', 'ZDRc [dB]']
#
#    datatype_vec = [
#        'hydro',
#        'KDPc',
#        'dBZc',
#        'RhoHVc',
#        'TEMP',
#        'ZDRc']
#
#    step_list = [
#        None,
#        0.05,
#        0.5,
#        0.001,
#        1.,
#        0.1]

    pol_vals_labels = [
        'hydro', 'entropy', 'propAG', 'propCR', 'propIH', 'propLR', 'propMH', 'propRN',
        'propRP', 'propVI', 'propWS']

    datatype_vec = [
        'hydro', 'entropy', 'propAG', 'propCR', 'propIH', 'propLR', 'propMH', 'propRN',
        'propRP', 'propVI', 'propWS']

    step_list = [None, 0.1, 1., 1., 1., 1., 1., 1., 1., 1., 1.]

    basename = 'Santis_data_entropy_IC'

    filt_type = 'keep_all'
    nsources_min = 10

    for label in pol_vals_labels:
        if 'hydro' in label:
            hydro_label = label
            break

    print("====== Lightning post-processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Lightning post-processing finished: ")

    # read all the data to analyze
    flashnr, time_data, time_in_flash, lat, lon, alt, dBm, pol_vals_dict = (
        read_data(
            basepath, day_vec, basename=basename,
            pol_vals_labels=pol_vals_labels))

    # Get indices of data to keep
    if filt_type == 'keep_all':
        ind, data_ID, subtitl = get_indices_all_data(
            flashnr, nsources_min=nsources_min)
    elif filt_type == 'keep_solid':
        ind, data_ID, subtitl = get_indices_solid_phase(
            flashnr, pol_vals_dict[hydro_label], nsources_min=nsources_min)
    elif filt_type == 'keep_liquid':
        ind, data_ID, subtitl = get_indices_liquid_phase(
            flashnr, pol_vals_dict[hydro_label], nsources_min=nsources_min)
    elif filt_type == 'keep_liquid_origin':
        ind, data_ID, subtitl = get_indices_liquid_phase_origin(
            flashnr, pol_vals_dict[hydro_label], nsources_min=nsources_min)
    else:
        warn('Unknown filter type '+filt_type)
        return

    flashnr_filt = flashnr[ind]
    time_data_filt = time_data[ind]
    time_in_flash_filt = time_in_flash[ind]
    lat_filt = lat[ind]
    lon_filt = lon[ind]
    alt_filt = alt[ind]
    dBm_filt = dBm[ind]
    pol_vals_dict_filt = deepcopy(pol_vals_dict)
    for key in pol_vals_dict.keys():
        pol_vals_dict_filt[key] = pol_vals_dict[key][ind]

#    # write the filtered data in a file
#    vals_list = []
#    for label in pol_vals_labels:
#        vals_list.append(pol_vals_dict_filt[label])
#
#    fname = basepath+basename'_'+data_ID+'.csv'
#    write_ts_lightning(
#        flashnr_filt, time_data_filt, time_in_flash_filt, lat_filt, lon_filt,
#        alt_filt, dBm_filt, vals_list, fname, pol_vals_labels)
#    print('written to '+fname)

    # get flashes origin of filtered data
    flashnr_first, ind_first = np.unique(
        flashnr_filt, return_index=True)
    time_data_first = time_data_filt[ind_first]
    time_in_flash_first = time_in_flash_filt[ind_first]
    lat_first = lat_filt[ind_first]
    lon_first = lon_filt[ind_first]
    alt_first = alt_filt[ind_first]
    dBm_first = dBm_filt[ind_first]
    pol_vals_dict_first = deepcopy(pol_vals_dict_filt)
    for key in pol_vals_dict_filt.keys():
        pol_vals_dict_first[key] = pol_vals_dict_filt[key][ind_first]

    print('N flashes: '+str(flashnr_first.size))
    print('N sources: '+str(flashnr_filt.size))

    # Analyse the data

    # create histogram of hydrometeor proportions
    if 'propAG' in pol_vals_dict_filt:
        hydro_hist = np.zeros(10)
        bins_centers = np.arange(0, 10, 1)
        bins_edges = np.arange(-0.5, 10.5, 1)

        # Plot all sources histogram
        nradar_bins = pol_vals_dict_filt['propAG'].size

        # look for empty gates
        ind = np.ma.where(pol_vals_dict_filt['hydro'] == 0)[0]

        hydro_hist[0] = ind.size/nradar_bins*100.
        hydro_hist[1] = np.ma.sum(pol_vals_dict_filt['propAG'])/nradar_bins
        hydro_hist[2] = np.ma.sum(pol_vals_dict_filt['propCR'])/nradar_bins
        hydro_hist[3] = np.ma.sum(pol_vals_dict_filt['propLR'])/nradar_bins
        hydro_hist[4] = np.ma.sum(pol_vals_dict_filt['propRP'])/nradar_bins
        hydro_hist[5] = np.ma.sum(pol_vals_dict_filt['propRN'])/nradar_bins
        hydro_hist[6] = np.ma.sum(pol_vals_dict_filt['propVI'])/nradar_bins
        hydro_hist[7] = np.ma.sum(pol_vals_dict_filt['propWS'])/nradar_bins
        hydro_hist[8] = np.ma.sum(pol_vals_dict_filt['propMH'])/nradar_bins
        hydro_hist[9] = np.ma.sum(pol_vals_dict_filt['propIH'])/nradar_bins

        fname = (
            basepath+data_ID+'_allsources_ts_trajlightning_hydro_prop.png')
        plot_histogram2(
            bins_centers, hydro_hist, [fname],
            labelx='radar echo classification (-)', labely='percentage',
            titl='Trajectory Histogram All Sources'+subtitl)

        print("----- plot to '%s'" % fname)

        # store histogram
        fname = (
            basepath+data_ID+'_allsources_ts_trajlightning_hydro_prop.csv')
        write_histogram(bins_edges, hydro_hist, fname)
        print('Written '+fname)

        # Plot first sources histogram
        nradar_bins = pol_vals_dict_first['propAG'].size

        # look for empty gates
        ind = np.ma.where(pol_vals_dict_first['hydro'] == 0)[0]

        hydro_hist[0] = ind.size/nradar_bins*100.
        hydro_hist[1] = np.ma.sum(pol_vals_dict_first['propAG'])/nradar_bins
        hydro_hist[2] = np.ma.sum(pol_vals_dict_first['propCR'])/nradar_bins
        hydro_hist[3] = np.ma.sum(pol_vals_dict_first['propLR'])/nradar_bins
        hydro_hist[4] = np.ma.sum(pol_vals_dict_first['propRP'])/nradar_bins
        hydro_hist[5] = np.ma.sum(pol_vals_dict_first['propRN'])/nradar_bins
        hydro_hist[6] = np.ma.sum(pol_vals_dict_first['propVI'])/nradar_bins
        hydro_hist[7] = np.ma.sum(pol_vals_dict_first['propWS'])/nradar_bins
        hydro_hist[8] = np.ma.sum(pol_vals_dict_first['propMH'])/nradar_bins
        hydro_hist[9] = np.ma.sum(pol_vals_dict_first['propIH'])/nradar_bins

        fname = (
            basepath+data_ID+'_firstsource_ts_trajlightning_hydro_prop.png')
        plot_histogram2(
            bins_centers, hydro_hist, [fname],
            labelx='radar echo classification (-)', labely='percentage',
            titl='Trajectory Histogram First Sources'+subtitl)

        print("----- plot to '%s'" % fname)

        # store histogram
        fname = (
            basepath+data_ID+'_firstsource_ts_trajlightning_hydro_prop.csv')
        write_histogram(bins_edges, hydro_hist, fname)
        print('Written '+fname)

    for i, key in enumerate(pol_vals_labels):
        step = step_list[i]
        datatype = datatype_vec[i]

        field_name = get_fieldname_pyart(datatype)
        field_dict = get_metadata(field_name)

        labelx = get_colobar_label(field_dict, field_name)

        vals = pol_vals_dict_filt[key]
        bins, values = compute_histogram(vals, field_name, step=step)

        print(datatype+' min: '+str(vals.min()))
        print(datatype+' max: '+str(vals.max()))

        # Plot all sources histogram
        fname_first_source = (
            basepath+data_ID+'_allsources_ts_trajlightning_' +
            datatype+'.png')
        plot_histogram(
            bins, values, [fname_first_source], labelx=labelx,
            titl='Trajectory Histogram All Sources'+subtitl)

        print("----- plot to '%s'" % fname_first_source)

        # store histogram
        fname_first_source = (
            basepath+data_ID+'_allsources_ts_trajlightning_' +
            datatype+'.csv')
        hist_values, _ = np.histogram(values, bins=bins)
        write_histogram(bins, hist_values, fname_first_source)
        print('Written '+fname_first_source)

        # First sources
        vals = pol_vals_dict_first[key]
        bins, values = compute_histogram(vals, field_name, step=step)

        # Plot first source histogram
        fname_first_source = (
            basepath+data_ID+'_firstsource_ts_trajlightning_' +
            datatype+'.png')
        plot_histogram(
            bins, values, [fname_first_source], labelx=labelx,
            titl='Trajectory Histogram First Source'+subtitl)

        print("----- plot to '%s'" % fname_first_source)

        # store histogram
        fname_first_source = (
            basepath+data_ID+'_firstsource_ts_trajlightning_' +
            datatype+'.csv')
        hist_values_first, _ = np.histogram(values, bins=bins)
        write_histogram(bins, hist_values_first, fname_first_source)
        print('Written '+fname_first_source)

    # Get histograms of sources altitude and power

    # define histogram bin edges
    bin_edges_alt = np.arange(-50., 14150., 100.)
    bin_edges_dBm = np.arange(-17., 47., 1.)
    bin_edges_time = np.arange(0, 25, 1)

    # Plot histogram time of occurrence
    time_hour_first = np.empty(time_data_first.size)
    for i, dt_first in enumerate(time_data_first):
        time_first = dt_first.time()
        time_hour_first[i] = (
            time_first.hour+time_first.minute/60.+time_first.second/3600. +
            time_first.microsecond/3600e6)

    fname_hist = basepath+data_ID+'_Santis_hist_time.png'
    fname_hist = plot_histogram(
        bin_edges_time, time_hour_first, [fname_hist], labelx='Hour [UTC]',
        titl='Flash occurrence time'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_time.csv'
    hist_time, _ = np.histogram(time_hour_first, bins=bin_edges_time)
    fname_hist = write_histogram(bin_edges_time, hist_time, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram altitude all sources
    fname_hist = basepath+data_ID+'_Santis_hist_alt.png'
    fname_hist = plot_histogram(
        bin_edges_alt, alt_filt, [fname_hist], labelx='Altitude [m MSL]',
        titl='Flash sources altitude'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_alt.csv'
    hist_alt, _ = np.histogram(alt_filt, bins=bin_edges_alt)
    fname_hist = write_histogram(bin_edges_alt, hist_alt, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram altitude first sources
    fname_hist = basepath+data_ID+'_Santis_hist_alt_first_source.png'
    fname_hist = plot_histogram(
        bin_edges_alt, alt_first, [fname_hist], labelx='Altitude [m MSL]',
        titl='Flash first source altitude'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_alt_first_source.csv'
    hist_alt_fist, _ = np.histogram(alt_first, bins=bin_edges_alt)
    fname_hist = write_histogram(bin_edges_alt, hist_alt_fist, fname_hist)
    print('Written '+fname_hist)

    fname_hist = basepath+data_ID+'_Santis_hist_dBm.png'
    fname_hist = plot_histogram(
        bin_edges_dBm, dBm_filt, [fname_hist], labelx='Power [dBm]',
        titl='Flash sources power'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_dBm.csv'
    hist_dBm, _ = np.histogram(dBm_filt, bins=bin_edges_dBm)
    fname_hist = write_histogram(bin_edges_dBm, hist_dBm, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram power first sources
    fname_hist = basepath+data_ID+'_Santis_hist_dBm_first_source.png'
    fname_hist = plot_histogram(
        bin_edges_dBm, dBm_first, [fname_hist], labelx='Power [dBm]',
        titl='Flash first source power'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_dBm_first_source.csv'
    hist_dBm_first, _ = np.histogram(dBm_first, bins=bin_edges_dBm)
    fname_hist = write_histogram(bin_edges_dBm, hist_dBm_first, fname_hist)
    print('Written '+fname_hist)

    # Plot 2D histogram all sources
    H, _, _ = np.histogram2d(alt_filt, dBm_filt, bins=[bin_edges_alt, bin_edges_dBm])

    # set 0 values to blank
    H = np.ma.asarray(H)
    H[H == 0] = np.ma.masked

    fname_hist = basepath+data_ID+'_Santis_2Dhist_alt_dBm.png'
    fname_hist = _plot_time_range(
        bin_edges_alt, bin_edges_dBm, H, None, [fname_hist],
        titl='LMA sources Altitude-Power histogram'+subtitl,
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

    fname_hist = basepath+data_ID+'_Santis_2Dhist_alt_dBm_first_source.png'
    fname_hist = _plot_time_range(
        bin_edges_alt, bin_edges_dBm, H, None, [fname_hist],
        titl='LMA first sources Altitude-Power histogram'+subtitl,
        xlabel='Altitude [m MSL]', ylabel='Power [dBm]',
        clabel='Occurrence',
        vmin=0, vmax=None, figsize=[10, 8], dpi=72)
    print('Plotted '+' '.join(fname_hist))

    # plot position all sources
    figfname = basepath+data_ID+'_Santis_LMA_sources_pos_max_height_on_top.png'
    figfname = plot_pos(
        lat_filt, lon_filt, alt_filt, [figfname], sort_altitude='Highest_on_top',
        cb_label='Source height [m MSL]',
        titl='Flash sources position. Highest on top'+subtitl)
    print('Plotted '+' '.join(figfname))

    figfname = basepath+data_ID+'_Santis_LMA_sources_pos_min_height_on_top.png'
    figfname = plot_pos(
        lat_filt, lon_filt, alt_filt, [figfname], sort_altitude='Lowest_on_top',
        cb_label='Source height [m MSL]',
        titl='Flash sources position. Lowest on top'+subtitl)
    print('Plotted '+' '.join(figfname))

    # plot position first source
    figfname = (
        basepath+data_ID+'_Santis_LMA_first_source_pos_max_height_on_top.png')
    figfname = plot_pos(
        lat_first, lon_first, alt_first, [figfname],
        sort_altitude='Highest_on_top', cb_label='Source height [m MSL]',
        titl='First flash source position. Highest on top'+subtitl)
    print('Plotted '+' '.join(figfname))

    figfname = (
        basepath+data_ID+'_Santis_LMA_first_source_pos_min_height_on_top.png')
    figfname = plot_pos(
        lat_first, lon_first, alt_first, [figfname],
        sort_altitude='Lowest_on_top', cb_label='Source height [m MSL]',
        titl='First flash source position. Lowest on top'+subtitl)
    print('Plotted '+' '.join(figfname))


def read_data(basepath, day_vec, basename='Santis_data',
              pol_vals_labels=['hydro [-]', 'KDPc [deg/Km]', 'dBZc [dBZ]',
                               'RhoHVc [-]', 'TEMP [deg C]', 'ZDRc [dB]']):
    """
    reads data to analyze

    Parameters
    ----------
    day_vec : array of datetime objects
        The dates to analyze
    basename : str
        the base name of the files containing the data
    pol_vals_labels : list of str
        the labels corresponding to the polarimetric radar variables contained
        in the file to be analyzed

    Returns
    -------
    flashnr, time_data, time_in_flash, lat, lon, alt, dBm,
    pol_vals_dict : tupple
        tupple containing all the data read

    """
    flashnr = np.asarray([], dtype=int)
    time_data = np.asarray([], dtype=datetime.datetime)
    time_in_flash = np.asarray([], dtype=float)
    lat = np.asarray([], dtype=float)
    lon = np.asarray([], dtype=float)
    alt = np.asarray([], dtype=float)
    dBm = np.asarray([], dtype=float)

    pol_vals_dict = dict()
    for label in pol_vals_labels:
        pol_vals_dict.update({label: np.ma.asarray([])})

    # Read all data
    for day in day_vec:
        day_str = day.strftime('%Y%m%d')
        fname = basepath+day_str+'_'+basename+'.csv'
        print('Reading data file '+fname)
        (flashnr_aux, time_data_aux, time_in_flash_aux, lat_aux, lon_aux,
         alt_aux, dBm_aux, pol_vals_dict_aux) = read_lightning_all(
             fname, labels=pol_vals_labels)

        flashnr = np.append(flashnr, int(day_str)*1000000+flashnr_aux)
        time_data = np.append(time_data, time_data_aux)
        time_in_flash = np.append(time_in_flash, time_in_flash_aux)
        lat = np.append(lat, lat_aux)
        lon = np.append(lon, lon_aux)
        alt = np.append(alt, alt_aux)
        dBm = np.append(dBm, dBm_aux)
        for key in pol_vals_dict.keys():
            pol_vals_dict[key] = np.ma.append(
                pol_vals_dict[key], pol_vals_dict_aux[key])

    return (
        flashnr, time_data, time_in_flash, lat, lon, alt, dBm, pol_vals_dict)


def read_data_two_sources(basepath, day_vec, basename1='Santis_data',
                          basename2='Santis_data', basename_out='Santis_data',
                          pol_vals_labels=['hydro [-]', 'KDPc [deg/Km]',
                                           'dBZc [dBZ]', 'RhoHVc [-]',
                                           'TEMP [deg C]', 'ZDRc [dB]']):
    """
    reads data to analyze

    Parameters
    ----------
    day_vec : array of datetime objects
        The dates to analyze
    basename : str
        the base name of the files containing the data
    pol_vals_labels : list of str
        the labels corresponding to the polarimetric radar variables contained
        in the file to be analyzed

    Returns
    -------
    flashnr, time_data, time_in_flash, lat, lon, alt, dBm,
    pol_vals_dict : tupple
        tupple containing all the data read

    """
    flashnr = np.asarray([], dtype=int)
    time_data = np.asarray([], dtype=datetime.datetime)
    time_in_flash = np.asarray([], dtype=float)
    lat = np.asarray([], dtype=float)
    lon = np.asarray([], dtype=float)
    alt = np.asarray([], dtype=float)
    dBm = np.asarray([], dtype=float)

    pol_vals_dict = dict()
    for label in pol_vals_labels:
        pol_vals_dict.update({label: np.ma.asarray([])})

    # Read all data
    flash_cnt = 0
    for day in day_vec:
        day_str = day.strftime('%Y%m%d')
        fname = basepath+day_str+'_'+basename1+'.csv'
        print('Reading data file '+fname)
        (flashnr_aux, time_data_aux, time_in_flash_aux, lat_aux, lon_aux,
         alt_aux, dBm_aux, pol_vals_dict_aux) = read_lightning_all(
             fname, labels=pol_vals_labels)

        fname = basepath+day_str+'_'+basename2+'.csv'
        print('Reading data file '+fname)
        (flashnr_aux2, time_data_aux2, time_in_flash_aux2, lat_aux2, lon_aux2,
         alt_aux2, dBm_aux2, pol_vals_dict_aux2) = read_lightning_all(
             fname, labels=pol_vals_labels)

        # get unique flashes in both files
        flashnr_first1 = np.unique(flashnr_aux, return_index=False)
        flashnr_first2 = np.unique(flashnr_aux2, return_index=True)

        # keep flashes common in both files
        flashnr_both = flashnr_first1[
            np.isin(flashnr_first1, flashnr_first2, assume_unique=True,
                    invert=False)]

        ind = []
        for flash in flashnr_both:
            ind.extend(np.where(flashnr_aux == flash)[0])

        flashnr_filt = flashnr_aux[ind]
        time_data_filt = time_data_aux[ind]
        time_in_flash_filt = time_in_flash_aux[ind]
        lat_filt = lat_aux[ind]
        lon_filt = lon_aux[ind]
        alt_filt = alt_aux[ind]
        dBm_filt = dBm_aux[ind]
        pol_vals_dict_filt = deepcopy(pol_vals_dict_aux)
        for key in pol_vals_dict_aux.keys():
            pol_vals_dict_filt[key] = pol_vals_dict_aux[key][ind]

        flashnr = np.append(flashnr, flashnr_aux+flash_cnt)
        flash_cnt += flashnr_aux.max()
        time_data = np.append(time_data, time_data_aux)
        time_in_flash = np.append(time_in_flash, time_in_flash_aux)
        lat = np.append(lat, lat_aux)
        lon = np.append(lon, lon_aux)
        alt = np.append(alt, alt_aux)
        dBm = np.append(dBm, dBm_aux)
        for key in pol_vals_dict.keys():
            pol_vals_dict[key] = np.ma.append(
                pol_vals_dict[key], pol_vals_dict_aux[key])

        # write the results in a file
        if basename_out is not None:
            vals_list = []
            for label in pol_vals_labels:
                vals_list.append(pol_vals_dict_filt[label])

            fname = basepath+day_str+'_'+basename_out+'.csv'
            write_ts_lightning(
                flashnr_filt, time_data_filt, time_in_flash_filt, lat_filt,
                lon_filt, alt_filt, dBm_filt, vals_list, fname,
                pol_vals_labels)
            print('written to '+fname)

    return flashnr, time_data, time_in_flash, lat, lon, alt, dBm, pol_vals_dict


def get_indices_all_data(flashnr, nsources_min=0):
    """
    Get indices of all flash sources

    Parameters
    ----------
    flashnr: 1D array
        The flash number of each source

    Returns
    -------
    ind : 1D array
        the indices corresponding to the data to keep
    data_ID : str
        an identifier of the type of data kept
    subtitl : str
        the subtitle to add to the plots generated

    """
#    data_ID = 'All'
#    subtitl = ''

    data_ID = 'IC'
    subtitl = '\nLMA flashes with associated '+data_ID+' EUCLID flashes'

    # Get unique flashes
    unique_flashnr = np.unique(flashnr, return_index=False)

    # get sources of each flashes
    ind = []
    for flash in unique_flashnr:
        ind_flash = np.where(flashnr == flash)[0]
        if ind_flash.size < nsources_min:
            continue
        ind.extend(ind_flash)

    return ind, data_ID, subtitl


def get_indices_solid_phase(flashnr, hydro, nsources_min=0):
    """
    Get indices of flash sources of flashes propagating exclusively within the
    solid phase

    Parameters
    ----------
    flashnr: 1D array
        The flash number of each source
    hydro : 1D array
        The dominant hydrometeor class in the region of each source
    nsources_min : float
        Minimum number of sources to accept the flash

    Returns
    -------
    ind : 1D array
        the indices corresponding to the data to keep
    data_ID : str
        an identifier of the type of data kept
    subtitl : str
        the subtitle to add to the plots generated

    """
    data_ID = 'solid'
    subtitl = '\nFlashes propagating into the solid phase only'

    ind = np.ma.where(np.logical_or.reduce((
        np.ma.getmaskarray(hydro) == 1, hydro == 0., hydro == 3., hydro == 5.,
        hydro == 7., hydro == 8.)))[0]

    # Get unique flash that contain these sources
    flashnr_filter = np.unique(flashnr[ind], return_index=False)

    # get Flash ID of all flashes
    flashnr_ID = np.unique(flashnr, return_index=False)

    # get flashes not in filtered flashes
    unique_flashnr_filt = flashnr_ID[
        np.isin(flashnr_ID, flashnr_filter, assume_unique=True, invert=True)]

    ind = []
    for flash in unique_flashnr_filt:
        ind_flash = np.where(flashnr == flash)[0]
        if ind_flash.size < nsources_min:
            continue
        ind.extend(ind_flash)

    return ind, data_ID, subtitl


def get_indices_liquid_phase(flashnr, hydro, nsources_min=0):
    """
    Get indices of flash sources of flashes propagating into the liquid phase

    Parameters
    ----------
    flashnr: 1D array
        The flash number of each source
    hydro : 1D array
        The dominant hydrometeor class in the region of each source
    nsources_min : float
        Minimum number of sources to accept the flash

    Returns
    -------
    ind : 1D array
        the indices corresponding to the data to keep
    data_ID : str
        an identifier of the type of data kept
    subtitl : str
        the subtitle to add to the plots generated

    """
    data_ID = 'liquid'
    subtitl = '\nFlashes propagating into the liquid phase'
    ind = np.ma.where(np.logical_or.reduce((
        hydro == 3., hydro == 5., hydro == 8.)))[0]

    # Get unique flashes that contain these sources
    unique_flashnr_filt = np.unique(flashnr[ind], return_index=False)

    # get sources of those flashes
    ind = []
    for flash in unique_flashnr_filt:
        ind_flash = np.where(flashnr == flash)[0]
        if ind_flash.size < nsources_min:
            continue
        ind.extend(ind_flash)

    return ind, data_ID, subtitl


def get_indices_liquid_phase_origin(flashnr, hydro, nsources_min=0):
    """
    Get indices of flash sources of flashes generated in the liquid phase

    Parameters
    ----------
    flashnr: 1D array
        The flash number of each source
    hydro : 1D array
        The dominant hydrometeor class in the region of each source
    nsources_min : float
        Minimum number of sources to accept the flash

    Returns
    -------
    ind : 1D array
        the indices corresponding to the data to keep
    data_ID : str
        an identifier of the type of data kept
    subtitl : str
        the subtitle to add to the plots generated

    """
    data_ID = 'liquid_origin'
    subtitl = '\nFlashes with origin in the liquid phase'

    # Get unique flashes origin
    unique_flashnr, unique_ind = np.unique(
        flashnr, return_index=True)

    # get sources with origin in liquid phase
    ind = []
    for i, flash in enumerate(unique_flashnr):
        orig_hydro = hydro[unique_ind[i]]
        if orig_hydro != 3. and orig_hydro != 5. and orig_hydro != 8.:
            continue

        ind_flash = np.where(flashnr == flash)[0]
        if ind_flash.size < nsources_min:
            continue
        ind.extend(ind_flash)

    return ind, data_ID, subtitl


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
