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
import shapely

from pyrad.io import read_lightning_all, get_fieldname_pyart, write_histogram
from pyrad.io import write_ts_lightning
from pyrad.graph import plot_histogram, plot_pos, get_colobar_label
from pyrad.graph import _plot_time_range, plot_histogram2
from pyrad.util import compute_histogram

from pyart.config import get_metadata
from pyart.core import wgs84_to_swissCH1903


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

    basename = 'Santis_data_entropy'
    filt_type = 'keep_non_solid_phase_origin'
    nsources_min = 10

    if 'entropy' in basename:
        pol_vals_labels = [
        'hydro', 'entropy', 'propAG', 'propCR', 'propIH', 'propLR', 'propMH',
        'propRN', 'propRP', 'propVI', 'propWS']

        datatype_vec = [
            'hydro', 'entropy', 'propAG', 'propCR', 'propIH', 'propLR',
            'propMH', 'propRN', 'propRP', 'propVI', 'propWS']

        step_list = [None, 0.1, 1., 1., 1., 1., 1., 1., 1., 1., 1.]
    else:
        pol_vals_labels = [
            'hydro [-]', 'KDPc [deg/Km]', 'dBZc [dBZ]', 'RhoHVc [-]',
            'TEMP [deg C]', 'ZDRc [dB]']

        datatype_vec = ['hydro', 'KDPc', 'dBZc', 'RhoHVc', 'TEMP', 'ZDRc']

        step_list = [None, 0.05, 0.5, 0.001, 1., 0.1]


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

#    flashnr, time_data, time_in_flash, lat, lon, alt, dBm, pol_vals_dict = (
#        read_data_two_sources(
#            basepath, day_vec, basename1='Santis_data_entropy',
#            basename2='Santis_data_entropy_CGt', basename_out='Santis_data_entropy_no_CG',
#            keep_common=False,
#            pol_vals_labels=pol_vals_labels))

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
    elif filt_type == 'keep_mixed_phase_origin':
        ind, data_ID, subtitl = get_indices_mixed_phase_origin(
            flashnr, pol_vals_dict[hydro_label], nsources_min=nsources_min)
    elif filt_type == 'keep_non_solid_phase_origin':
        ind, data_ID, subtitl = get_indices_non_solid_phase_origin(
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

    # get duration and area of flash
    duration_filt = np.ma.masked_all(flashnr_first.size)
    area_filt = np.ma.masked_all(flashnr_first.size)

    chy_filt, chx_filt, _ = wgs84_to_swissCH1903(
        lon_filt, lat_filt, alt_filt, no_altitude_transform=True)

    for i, flash_ID in enumerate(flashnr_first):
        time_data_flash = time_data_filt[flashnr_filt == flash_ID]
        duration_filt[i] = (
            1e3*(time_data_flash[-1]-time_data_flash[0]).total_seconds())

        chy_flash = chy_filt[flashnr_filt == flash_ID]
        chx_flash = chx_filt[flashnr_filt == flash_ID]

        points_flash = shapely.geometry.MultiPoint(
            list(zip(chy_flash, chx_flash)))
        area_filt[i] = points_flash.minimum_rotated_rectangle.area*1e-6

    print('N flashes: '+str(flashnr_first.size))
    print('N sources: '+str(flashnr_filt.size))

    # Analyse the data

    # create histograms of hydrometeor proportions
    if 'propAG' in pol_vals_dict_filt:
        bins_centers = np.arange(0, 10, 1)
        bins_edges = np.arange(-0.5, 10.5, 1)

        # Create histogram of number of different hydrometeors types in each
        # radar range gate. All sources
        nhydros_hist = hist_nhydros_gate(pol_vals_dict_filt, percent_min=10.)

        fname = (
            basepath+data_ID+'_allsources_ts_trajlightning_nhydro.png')
        plot_histogram2(
            bins_centers, nhydros_hist, [fname],
            labelx='Number of hydrometeors in radar range gate', labely='occurrence',
            titl='Trajectory Histogram All Sources'+subtitl)

        print("----- plot to '%s'" % fname)

        # store histogram
        fname = (
            basepath+data_ID+'_allsources_ts_trajlightning_nhydro.csv')
        write_histogram(bins_edges, nhydros_hist, fname)
        print('Written '+fname)


        # Create histogram of number of different hydrometeors types in each
        # radar range gate. First source
        nhydros_hist = hist_nhydros_gate(pol_vals_dict_first, percent_min=10.)

        fname = (
            basepath+data_ID+'_firstsource_ts_trajlightning_nhydro.png')
        plot_histogram2(
            bins_centers, nhydros_hist, [fname],
            labelx='Number of hydrometeors in radar range gate', labely='occurrence',
            titl='Trajectory Histogram First Sources'+subtitl)

        print("----- plot to '%s'" % fname)

        # store histogram
        fname = (
            basepath+data_ID+'_firstsource_ts_trajlightning_nhydro.csv')
        write_histogram(bins_edges, nhydros_hist, fname)
        print('Written '+fname)

        # Create histograms of dominant hydrometeors all sources
        hydro_hist2 = hist_dominant_hydrometeors(
            pol_vals_dict_filt, percent_min=10.)

        fname_hist = basepath+data_ID+'_allsources_ts_trajlightning_hydro_dominant.png'
        fname_hist = _plot_time_range(
            bins_edges, bins_edges, hydro_hist2, None, [fname_hist],
            titl='Trajectory Histogram All Sources'+subtitl,
            xlabel='Dominant hydrometeor', ylabel='2nd most dominant hydrometeor',
            vmin=0, clabel='Occurrence', figsize=[10, 8], dpi=72)
        print('Plotted '+' '.join(fname_hist))

        # Create histogram of dominant hydrometeors first sources
        hydro_hist2 = hist_dominant_hydrometeors(
            pol_vals_dict_first, percent_min=10.)

        fname_hist = basepath+data_ID+'_firstsource_ts_trajlightning_hydro_dominant.png'
        fname_hist = _plot_time_range(
            bins_edges, bins_edges, hydro_hist2, None, [fname_hist],
            titl='Trajectory Histogram First Sources'+subtitl,
            xlabel='Dominant hydrometeor', ylabel='2nd most dominant hydrometeor',
            vmin=0, clabel='Occurrence', figsize=[10, 8], dpi=72)
        print('Plotted '+' '.join(fname_hist))

        # create histogram of percentage of dominant hydrometeor all sources
        hydro_hist = hist_hydrometeor_mixtures(pol_vals_dict_filt)

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

        # create histogram of percentage of dominant hydrometeor first sources
        hydro_hist = hist_hydrometeor_mixtures(pol_vals_dict_first)

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
    bin_edges_area = np.arange(0., 2100., 100.)
    bin_edges_duration = np.arange(0., 1100., 100.)

    # Plot histogram of LMA flash area
    _, area_filt_values = compute_histogram(
        area_filt, None, bin_edges=bin_edges_area)
    fname_hist = basepath+data_ID+'_Santis_hist_area.png'
    fname_hist = plot_histogram(
        bin_edges_area, area_filt_values, [fname_hist], labelx='Area [km2]',
        titl='Flash area'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_area.csv'
    hist_area, _ = np.histogram(area_filt_values, bins=bin_edges_area)
    fname_hist = write_histogram(bin_edges_area, hist_area, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram of LMA flash duration
    _, duration_filt_values = compute_histogram(
        duration_filt, None, bin_edges=bin_edges_duration)
    fname_hist = basepath+data_ID+'_Santis_hist_duration.png'
    fname_hist = plot_histogram(
        bin_edges_duration, duration_filt_values, [fname_hist],
        labelx='Duration [ms]',
        titl='Flash duration'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_duration.csv'
    hist_duration, _ = np.histogram(duration_filt_values, bins=bin_edges_duration)
    fname_hist = write_histogram(bin_edges_duration, hist_duration, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram time of occurrence
    time_hour_first = occurrence_time(time_data_first)

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
    _, alt_filt_values = compute_histogram(
        alt_filt, None, bin_edges=bin_edges_alt)
    fname_hist = basepath+data_ID+'_Santis_hist_alt.png'
    fname_hist = plot_histogram(
        bin_edges_alt, alt_filt_values, [fname_hist], labelx='Altitude [m MSL]',
        titl='Flash sources altitude'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_alt.csv'
    hist_alt, _ = np.histogram(alt_filt_values, bins=bin_edges_alt)
    fname_hist = write_histogram(bin_edges_alt, hist_alt, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram altitude first sources
    _, alt_first_values = compute_histogram(
        alt_first, None, bin_edges=bin_edges_alt)
    fname_hist = basepath+data_ID+'_Santis_hist_alt_first_source.png'
    fname_hist = plot_histogram(
        bin_edges_alt, alt_first_values, [fname_hist], labelx='Altitude [m MSL]',
        titl='Flash first source altitude'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_alt_first_source.csv'
    hist_alt_fist, _ = np.histogram(alt_first_values, bins=bin_edges_alt)
    fname_hist = write_histogram(bin_edges_alt, hist_alt_fist, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram power all sources
    _, dBm_filt_values = compute_histogram(
        dBm_filt, None, bin_edges=bin_edges_dBm)
    fname_hist = basepath+data_ID+'_Santis_hist_dBm.png'
    fname_hist = plot_histogram(
        bin_edges_dBm, dBm_filt_values, [fname_hist], labelx='Power [dBm]',
        titl='Flash sources power'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_dBm.csv'
    hist_dBm, _ = np.histogram(dBm_filt_values, bins=bin_edges_dBm)
    fname_hist = write_histogram(bin_edges_dBm, hist_dBm, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram power first sources
    _, dBm_first_values = compute_histogram(
        dBm_first, None, bin_edges=bin_edges_dBm)
    fname_hist = basepath+data_ID+'_Santis_hist_dBm_first_source.png'
    fname_hist = plot_histogram(
        bin_edges_dBm, dBm_first_values, [fname_hist], labelx='Power [dBm]',
        titl='Flash first source power'+subtitl)
    print('Plotted '+' '.join(fname_hist))

    fname_hist = basepath+data_ID+'_Santis_hist_dBm_first_source.csv'
    hist_dBm_first, _ = np.histogram(dBm_first_values, bins=bin_edges_dBm)
    fname_hist = write_histogram(bin_edges_dBm, hist_dBm_first, fname_hist)
    print('Written '+fname_hist)

    # Plot 2D histogram all sources
    H, _, _ = np.histogram2d(
        alt_filt_values, dBm_filt_values, bins=[bin_edges_alt, bin_edges_dBm])

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
        alt_first_values, dBm_first_values, bins=[bin_edges_alt, bin_edges_dBm])

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
                          keep_common=True,
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

        if keep_common:
            # get unique flashes in both files
            flashnr_first1 = np.unique(flashnr_aux, return_index=False)
            flashnr_first2 = np.unique(flashnr_aux2, return_index=True)

            # get common flashes
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
        else:
            (flashnr_filt, time_data_filt, time_in_flash_filt, lat_filt,
             lon_filt, alt_filt, dBm_filt, pol_vals_dict_filt,
             flashnr_unique1) = get_unique_flashes_data(
                 flashnr_aux, time_data_aux, time_in_flash_aux, lat_aux,
                 lon_aux, alt_aux, dBm_aux, pol_vals_dict_aux, flashnr_aux2)

            (flashnr_filt2, time_data_filt2, time_in_flash_filt2, lat_filt2,
             lon_filt2, alt_filt2, dBm_filt2, pol_vals_dict_filt2,
             flashnr_unique2) = get_unique_flashes_data(
                 flashnr_aux2, time_data_aux2, time_in_flash_aux2, lat_aux2,
                 lon_aux2, alt_aux2, dBm_aux2, pol_vals_dict_aux2, flashnr_aux)

            # join unique flashes
            flashnr_filt = np.ma.append(flashnr_filt, flashnr_filt2)
            time_data_filt = np.ma.append(time_data_filt, time_data_filt2)
            time_in_flash_filt = np.ma.append(
                time_in_flash_filt, time_in_flash_filt2)
            lat_filt = np.ma.append(lat_filt, lat_filt2)
            lon_filt = np.ma.append(lon_filt, lon_filt2)
            alt_filt = np.ma.append(alt_filt, alt_filt2)
            dBm_filt = np.ma.append(dBm_filt, dBm_filt2)
            for key in pol_vals_dict_filt.keys():
                pol_vals_dict_filt[key] = np.ma.append(
                    pol_vals_dict_filt[key],
                    pol_vals_dict_filt2[key])

            # reorder flashes according to flash ID
            flashnr_unique = np.ma.append(flashnr_unique1, flashnr_unique2)
            flashnr_unique = np.ma.sort(flashnr_unique)

            ind = []
            for flash in flashnr_unique:
                ind.extend(np.where(flashnr_filt == flash)[0])

            flashnr_filt = flashnr_filt[ind]
            time_data_filt = time_data_filt[ind]
            time_in_flash_filt = time_in_flash_filt[ind]
            lat_filt = lat_filt[ind]
            lon_filt = lon_filt[ind]
            alt_filt = alt_filt[ind]
            dBm_filt = dBm_filt[ind]
            for key in pol_vals_dict_filt.keys():
                pol_vals_dict_filt[key] = pol_vals_dict_filt[key][ind]

        flashnr = np.append(flashnr, flashnr_filt+flash_cnt)
        flash_cnt += flashnr_aux.max()
        time_data = np.append(time_data, time_data_filt)
        time_in_flash = np.append(time_in_flash, time_in_flash_filt)
        lat = np.append(lat, lat_filt)
        lon = np.append(lon, lon_filt)
        alt = np.append(alt, alt_filt)
        dBm = np.append(dBm, dBm_filt)
        for key in pol_vals_dict.keys():
            pol_vals_dict[key] = np.ma.append(
                pol_vals_dict[key], pol_vals_dict_filt[key])

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


def get_unique_flashes_data(flashnr1, time_data1, time_in_flash1, lat1, lon1,
                            alt1, dBm1, pol_vals_dict1, flashnr2):
    """
    Get flash data that is in one dataset but not in the other

    Parameters
    ----------
    flashnr1, time_data1, time_in_flash1, lat1, lon1,
    alt1, dBm1, pol_vals_dict1,: 1D arrays
        The data of dataset one
    flashnr2 : 1D array
        the identifier of flashes in dataset 2

    Returns
    -------
    flashnr_filt, time_data_filt, time_in_flash_filt, lat_filt, lon_filt,
    alt_filt, dBm_filt, pol_vals_dict_filt : 1D arrays
        The filtered data

    """
    # get unique flashes in both datasets
    flashnr_first1 = np.unique(flashnr1, return_index=False)
    flashnr_first2 = np.unique(flashnr2, return_index=True)

    # get flashes in dataset 1 not in dataset2
    flashnr_unique1 = flashnr_first1[
        np.isin(flashnr_first1, flashnr_first2, assume_unique=True,
                invert=True)]

    ind = []
    for flash in flashnr_unique1:
        ind.extend(np.where(flashnr1 == flash)[0])

    # flashes unique to first dataset
    flashnr_filt = flashnr1[ind]
    time_data_filt = time_data1[ind]
    time_in_flash_filt = time_in_flash1[ind]
    lat_filt = lat1[ind]
    lon_filt = lon1[ind]
    alt_filt = alt1[ind]
    dBm_filt = dBm1[ind]
    pol_vals_dict_filt = deepcopy(pol_vals_dict1)
    for key in pol_vals_dict_filt.keys():
        pol_vals_dict_filt[key] = pol_vals_dict1[key][ind]

    return (
        flashnr_filt, time_data_filt, time_in_flash_filt, lat_filt, lon_filt,
        alt_filt, dBm_filt, pol_vals_dict_filt, flashnr_unique1)


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
    data_ID = 'All'
    subtitl = ''

    data_ID = 'CGp'
    subtitl = '\nLMA flashes with associated '+data_ID+' EUCLID flashes'

    data_ID = 'no_CG'
    subtitl = '\nLMA flashes without associated CG EUCLID flashes'

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


def get_indices_non_solid_phase_origin(flashnr, hydro, nsources_min=0):
    """
    Get indices of flash sources of flashes generated in the liquid and mixed
    phases

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
    data_ID = 'non_solid_origin'
    subtitl = '\nFlashes with origin in the liquid and mixed phases'

    # Get unique flashes origin
    unique_flashnr, unique_ind = np.unique(
        flashnr, return_index=True)

    # Get indices of flashes with origin in liquid and mixed phase layer
    ind = np.ma.where(np.logical_or.reduce((
        hydro[unique_ind] == 3., hydro[unique_ind] == 5.,
        hydro[unique_ind] == 7., hydro[unique_ind] == 8.)))[0]

    # Get unique flashes that contain these sources
    unique_flashnr_filt = np.unique(unique_flashnr[ind], return_index=False)

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

    # Get indices of flashes with origin in liquid layer
    ind = np.ma.where(np.logical_or.reduce((
        hydro[unique_ind] == 3., hydro[unique_ind] == 5.,
        hydro[unique_ind] == 8.)))[0]

    # Get unique flashes that contain these sources
    unique_flashnr_filt = np.unique(unique_flashnr[ind], return_index=False)

    # get sources of those flashes
    ind = []
    for flash in unique_flashnr_filt:
        ind_flash = np.where(flashnr == flash)[0]
        if ind_flash.size < nsources_min:
            continue
        ind.extend(ind_flash)

    return ind, data_ID, subtitl


def get_indices_mixed_phase_origin(flashnr, hydro, nsources_min=0):
    """
    Get indices of flash sources of flashes generated in the mixed phase

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
    data_ID = 'mixed_phase_origin'
    subtitl = '\nFlashes with origin in the mixed phase'

    # Get unique flashes origin
    unique_flashnr, unique_ind = np.unique(
        flashnr, return_index=True)

    # Get indices of flashes with origin in wet snow
    ind = np.ma.where(hydro[unique_ind] == 7.)[0]

    # Get unique flashes that contain these sources
    unique_flashnr_filt = np.unique(unique_flashnr[ind], return_index=False)

    # get sources of those flashes
    ind = []
    for flash in unique_flashnr_filt:
        ind_flash = np.where(flashnr == flash)[0]
        if ind_flash.size < nsources_min:
            continue
        ind.extend(ind_flash)

    return ind, data_ID, subtitl


def hist_dominant_hydrometeors(pol_vals_dict, percent_min=10.):
    """
    Computes a 2D histogram with the dominant and second dominant
    hydrometeor

    Parameters
    ----------
    pol_vals_dict : dict
        dictionary containing the proportions of hydrometeors and the
        dominant hydrometeor type for each radar range gate where and
        LMA source was located
    percent_min : float
        Minimum percentage of the proportion of hydrometeors to consider
        the hydrometeor relevant

    Returns
    -------
    hydro_hist : 2D-array
        the 2D histogram

    """
    # keep only data that has been classified by the hydrometeor classification
    valid = np.logical_not(np.ma.getmaskarray(pol_vals_dict['hydro']))
    propAG = pol_vals_dict['propAG'][valid]
    propCR = pol_vals_dict['propCR'][valid]
    propLR = pol_vals_dict['propLR'][valid]
    propRP = pol_vals_dict['propRP'][valid]
    propRN = pol_vals_dict['propRN'][valid]
    propVI = pol_vals_dict['propVI'][valid]
    propWS = pol_vals_dict['propWS'][valid]
    propMH = pol_vals_dict['propMH'][valid]
    propIH = pol_vals_dict['propIH'][valid]
    hydro = pol_vals_dict['hydro'][valid]

    nradar_bins = propAG.size
    ind_nc = np.ma.where(hydro == 0)[0]
    ind_c = np.ma.where(hydro != 0)[0]

    values = np.ma.masked_all((10, nradar_bins))
    values[0, ind_nc] = 100.
    values[1, ind_c] = propAG[ind_c]
    values[2, ind_c] = propCR[ind_c]
    values[3, ind_c] = propLR[ind_c]
    values[4, ind_c] = propRP[ind_c]
    values[5, ind_c] = propRN[ind_c]
    values[6, ind_c] = propVI[ind_c]
    values[7, ind_c] = propWS[ind_c]
    values[8, ind_c] = propMH[ind_c]
    values[9, ind_c] = propIH[ind_c]

    inds_sorted = np.ma.argsort(values, endwith=False, axis=0)
    values_sorted = np.ma.sort(values, endwith=False, axis=0)

    dominant_hydros = inds_sorted[9, :]
    second_hydros = inds_sorted[8, :]
    prop_seconds = values_sorted[8, :]

    # if less than percent_min the second dominant is negligible
    second_hydros[prop_seconds < percent_min] = (
        dominant_hydros[prop_seconds < percent_min])

    hydro_hist = np.ma.zeros((10, 10))
    for i in range(nradar_bins):
        hydro_hist[dominant_hydros[i], second_hydros[i]] += 1

    hydro_hist[hydro_hist == 0] = np.ma.masked

    return hydro_hist


def hist_hydrometeor_mixtures(pol_vals_dict):
    """
    Computes the histogram of percentage of dominant hydrometeors considering
    proportions of hydrometeors in a range bin

    Parameters
    ----------
    pol_vals_dict : dict
        dictionary containing the proportions of hydrometeors and the
        dominant hydrometeor type for each radar range gate where and
        LMA source was located

    Returns
    -------
    hydro_hist : 1D-array
        the histogram

    """
    # keep only data that has been classified by the hydrometeor classification
    valid = np.logical_not(np.ma.getmaskarray(pol_vals_dict['hydro']))
    propAG = pol_vals_dict['propAG'][valid]
    propCR = pol_vals_dict['propCR'][valid]
    propLR = pol_vals_dict['propLR'][valid]
    propRP = pol_vals_dict['propRP'][valid]
    propRN = pol_vals_dict['propRN'][valid]
    propVI = pol_vals_dict['propVI'][valid]
    propWS = pol_vals_dict['propWS'][valid]
    propMH = pol_vals_dict['propMH'][valid]
    propIH = pol_vals_dict['propIH'][valid]
    hydro = pol_vals_dict['hydro'][valid]

    nradar_bins = propAG.size
    ind_nc = np.ma.where(hydro == 0)[0]
    ind_c = np.ma.where(hydro != 0)[0]

    hydro_hist = np.ma.zeros(10)
    hydro_hist[0] = ind_nc.size/nradar_bins*100.
    hydro_hist[1] = np.ma.sum(propAG[ind_c])/nradar_bins
    hydro_hist[2] = np.ma.sum(propCR[ind_c])/nradar_bins
    hydro_hist[3] = np.ma.sum(propLR[ind_c])/nradar_bins
    hydro_hist[4] = np.ma.sum(propRP[ind_c])/nradar_bins
    hydro_hist[5] = np.ma.sum(propRN[ind_c])/nradar_bins
    hydro_hist[6] = np.ma.sum(propVI[ind_c])/nradar_bins
    hydro_hist[7] = np.ma.sum(propWS[ind_c])/nradar_bins
    hydro_hist[8] = np.ma.sum(propMH[ind_c])/nradar_bins
    hydro_hist[9] = np.ma.sum(propIH[ind_c])/nradar_bins

    return hydro_hist


def hist_nhydros_gate(pol_vals_dict, percent_min=10.):
    """
    Computes the histogram of the number of hydrometeor types in a radar
    range gate

    Parameters
    ----------
    pol_vals_dict : dict
        dictionary containing the proportions of hydrometeors and the
        dominant hydrometeor type for each radar range gate where and
        LMA source was located

    Returns
    -------
    nhydro_hist : 1D-array
        the histogram

    """
    # keep only data that has been classified by the hydrometeor classification
    valid = np.logical_not(np.ma.getmaskarray(pol_vals_dict['hydro']))
    propAG = pol_vals_dict['propAG'][valid]
    propCR = pol_vals_dict['propCR'][valid]
    propLR = pol_vals_dict['propLR'][valid]
    propRP = pol_vals_dict['propRP'][valid]
    propRN = pol_vals_dict['propRN'][valid]
    propVI = pol_vals_dict['propVI'][valid]
    propWS = pol_vals_dict['propWS'][valid]
    propMH = pol_vals_dict['propMH'][valid]
    propIH = pol_vals_dict['propIH'][valid]
    hydro = pol_vals_dict['hydro'][valid]

    nradar_bins = propAG.size
    ind_nc = np.ma.where(hydro == 0)[0]
    ind_c = np.ma.where(hydro != 0)[0]

    print(ind_c.size)
    print(ind_nc.size)

    values = np.ma.zeros((10, nradar_bins))
    values[0, ind_nc] = 100.
    values[1, ind_c] = propAG[ind_c]
    values[2, ind_c] = propCR[ind_c]
    values[3, ind_c] = propLR[ind_c]
    values[4, ind_c] = propRP[ind_c]
    values[5, ind_c] = propRN[ind_c]
    values[6, ind_c] = propVI[ind_c]
    values[7, ind_c] = propWS[ind_c]
    values[8, ind_c] = propMH[ind_c]
    values[9, ind_c] = propIH[ind_c]

    # Mask with presence of each hydrometeor in each gates
    mask = np.ma.ones((9, nradar_bins), dtype=int)

    # If proportion more than percent_min consider hydrometeor as present
    mask[values[1:, :] < percent_min] = 0

    # Number of hydrometeors present in each gate
    nhydros = np.ma.sum(mask, axis=0)

    nhydros_hist = np.ma.zeros(10)
    nhydros_hist[0] = ind_nc.size

    for i in range(nradar_bins):
        if nhydros[i] > 0:
            nhydros_hist[nhydros[i]] += 1

    print(np.ma.sum(nhydros_hist))

    return nhydros_hist


def occurrence_time(dt_data):
    """
    Computes the time of occurrence [h] within the day of the LMA sources

    Parameters
    ----------
    dt_data : 1D-array
        array containing the date and time of occurence of the LMA sources

    Returns
    -------
    hydro_hist : 1D-array
        the histogram

    """
    time_hour = np.empty(dt_data.size)
    for i, dt in enumerate(dt_data):
        t_data = dt.time()
        time_hour[i] = (
            t_data.hour+t_data.minute/60.+t_data.second/3600. +
            t_data.microsecond/3600e6)

    return time_hour


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
