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
import glob
from warnings import warn
import os
from copy import deepcopy

import numpy as np

from pyrad.io import read_lightning_all, get_fieldname_pyart, write_histogram
from pyrad.io import write_ts_lightning
from pyrad.graph import plot_histogram, plot_pos, get_colobar_label
from pyrad.graph import plot_scatter_comp, _plot_time_range
from pyrad.util import compute_histogram

from pyart.config import get_metadata

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.pyplot as plt



print(__doc__)


def main():
    """
    """
    # basepath = '/data/pyrad_products/rad4alp_hydro_PHA/'
    basepath = '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/data_analysis/'
    day_vec = [
        datetime.datetime(2017, 6, 29),
        datetime.datetime(2017, 6, 30),
        datetime.datetime(2017, 7, 10),
        datetime.datetime(2017, 7, 14),
        datetime.datetime(2017, 7, 18),
        datetime.datetime(2017, 7, 19),
        datetime.datetime(2017, 7, 30),
        datetime.datetime(2017, 8, 1)]

    pol_vals_labels = [
        'hydro [-]', 'KDPc [deg/Km]', 'dBZc [dBZ]', 'RhoHVc [-]',
        'TEMP [deg C]', 'ZDRc [dB]']

    print("====== Lightning post-processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Lightning post-processing finished: ")

    flashnr = np.asarray([], dtype=int)
    time_data = np.asarray([], dtype=datetime.datetime)
    time_in_flash = np.asarray([], dtype=float)
    lat = np.asarray([], dtype=float)
    lon = np.asarray([], dtype=float)
    alt = np.asarray([], dtype=float)
    dBm = np.asarray([], dtype=float)
    pol_vals_dict = {
        'hydro [-]': np.ma.asarray([], dtype=np.int16),
        'KDPc [deg/Km]': np.ma.asarray([], dtype=float),
        'dBZc [dBZ]': np.ma.asarray([], dtype=float),
        'RhoHVc [-]': np.ma.asarray([], dtype=float),
        'TEMP [deg C]': np.ma.asarray([], dtype=float),
        'ZDRc [dB]': np.ma.asarray([], dtype=float)
    }

    # Read all data
    flash_cnt = 0
    for i, day in enumerate(day_vec):
        day_str = day.strftime('%Y%m%d')
        fname = basepath+'CG_data/'+day_str+'_Santis_data_CGt.csv'
        print('Reading data file '+fname)
        flashnr_aux, time_data_aux, time_in_flash_aux, lat_aux, lon_aux, alt_aux, dBm_aux, pol_vals_dict_aux = (
            read_lightning_all(
                fname, labels=pol_vals_labels))

        fname = basepath+'IC_data/'+day_str+'_Santis_data_IC.csv'
        print('Reading data file '+fname)
        flashnr_aux2, time_data_aux2, time_in_flash_aux2, lat_aux2, lon_aux2, alt_aux2, dBm_aux2, pol_vals_dict_aux2 = (
            read_lightning_all(
                fname, labels=pol_vals_labels))

        # get flashes in both files
        flashnr_first1, ind_first1 = np.unique(
            flashnr_aux, return_index=True)
        flashnr_first2, ind_first2 = np.unique(
            flashnr_aux2, return_index=True)

        flashnr_unique1 = flashnr_first1[np.isin(flashnr_first1, flashnr_first2, assume_unique=True, invert=True)]
        flashnr_unique2 = flashnr_first2[np.isin(flashnr_first2, flashnr_first1, assume_unique=True, invert=True)]

        print('Number of unique flashes: '+str(flashnr_unique1.size+flashnr_unique2.size))
        
        continue
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


        # write the results in a file
        vals_list = []
        for label in pol_vals_labels:
            vals_list.append(pol_vals_dict_filt[label])

        fname = basepath+day_str+'_Santis_data_CGpn.csv'
        write_ts_lightning(
            flashnr_filt, time_data_filt, time_in_flash_filt, lat_filt, lon_filt,
            alt_filt, dBm_filt, vals_list, fname, pol_vals_labels)
        print('written to '+fname)

        flashnr = np.append(flashnr, flashnr_filt+flash_cnt)
        flash_cnt += flashnr_aux.max()
        time_data = np.append(time_data, time_data_filt)
        time_in_flash = np.append(time_in_flash, time_in_flash_filt)
        lat = np.append(lat, lat_filt)
        lon = np.append(lon, lon_filt)
        alt = np.append(alt, alt_filt)
        dBm = np.append(dBm, dBm_filt)
        for key in pol_vals_dict.keys():
            pol_vals_dict[key] = np.ma.append(pol_vals_dict[key], pol_vals_dict_filt[key])

    # Filter data according to predominant hydrometeor class

#    # get indices of flashes with origin in the liquid region
#    data_ID = 'liquid_origin'
#    subtitl = '\nFlashes with origin in the liquid phase'
#
#    # Get unique flashes origin
#    unique_flashnr, unique_ind = np.unique(
#        flashnr, return_index=True)
#
#    # get sources with origin in liquid phase
#    ind = []
#    for i, flash in enumerate(unique_flashnr):
#        orig_hydro = pol_vals_dict['hydro [-]'][unique_ind[i]]
#        if orig_hydro == 3. or orig_hydro == 5. or orig_hydro == 8.:
#            ind.extend(np.where(flashnr == flash)[0])
#    # End of getting flashes with origin in the liquid phase

    # get all data
    data_ID = 'CGpn'
    subtitl = '\nLMA flashes associated to both CGn and CGp EUCLID strokes'
    ind = np.arange(flashnr.size)
    # End of getting all data
#
#    # get indices of sources in liquid phase region
#    data_ID = 'liquid'
#    subtitl = '\nFlashes propagating into the liquid phase'
#    ind = np.ma.where(np.logical_or.reduce((
#        pol_vals_dict['hydro [-]'] == 3.,
#        pol_vals_dict['hydro [-]'] == 5.,
#        pol_vals_dict['hydro [-]'] == 8.)))[0]
#
#    # Get unique flashes that contain these sources
#    unique_flashnr_filt, unique_ind_filt = np.unique(
#        flashnr[ind], return_index=True)
#
#    # get sources of those flashes
#    ind = []
#    for flash in unique_flashnr_filt:
#        ind.extend(np.where(flashnr == flash)[0])
#    # End of getting sources in the liquid phase region
#
#
#    # Get indices of sources propagating exclusively through the solid phase
#    # get sources in liquid, mixed phase, unclassified and no data regions
#    data_ID = 'solid'
#    subtitl = '\nFlashes propagating into the solid phase only'
#    ind = np.ma.where(np.logical_or.reduce((
#        np.ma.getmaskarray(pol_vals_dict['hydro [-]']) == 1,
#        pol_vals_dict['hydro [-]'] == 0.,
#        pol_vals_dict['hydro [-]'] == 3.,
#        pol_vals_dict['hydro [-]'] == 5.,
#        pol_vals_dict['hydro [-]'] == 7.,
#        pol_vals_dict['hydro [-]'] == 8.)))[0]
#
#    # Get unique flash that contain these sources
#    flashnr_filter = np.unique(flashnr[ind], return_index=False)
#
#    # get Flash ID of all flashes
#    flashnr_ID = np.unique(flashnr, return_index=False)
#
#    # get flashes not in filtered flashes
#    unique_flashnr_filt = flashnr_ID[
#        np.isin(flashnr_ID, flashnr_filter, assume_unique=True, invert=True)]
#
#    ind = []
#    for flash in unique_flashnr_filt:
#        ind.extend(np.where(flashnr == flash)[0])
#    # End of getting sources in the solid phase region

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

#    # write first sources in a file
#    vals_list = []
#    for key in pol_vals_dict.keys():
#        vals_list.append(pol_vals_dict_first[key])
#
#    fname = basepath+data_ID+'first_source_data.csv'
#    write_ts_lightning(
#        flashnr_first, time_data_first, time_in_flash_first, lat_first, lon_first,
#        alt_first, dBm_first, vals_list, fname, pol_vals_labels)
#    print('written to '+fname)
#
#    return
#
#    # write the results in a file
#    vals_list = []
#    for key in pol_vals_dict.keys():
#        vals_list.append(pol_vals_dict_filt[key])
#
#    fname = basepath+data_ID+'_data.csv'
#    write_ts_lightning(
#        flashnr_filt, time_data_filt, time_in_flash_filt, lat_filt, lon_filt,
#        alt_filt, dBm_filt, vals_list, fname, pol_vals_labels)
#    print('written to '+fname)
#
#    return

    # Get histogram of polarimetric variables at flash location
    datatype_vec = [
        'hydro',
        'KDPc',
        'dBZc',
        'RhoHVc',
        'TEMP',
        'ZDRc']

    step_list = [
        None,
        0.05,
        0.5,
        0.001,
        1.,
        0.1]

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


    return











    time_hour_first = np.empty(datetimes_first.size)
    for i, datetime_first in enumerate(datetimes_first):
        time_first = datetime_first.time()
        time_hour_first[i] = time_first.hour+time_first.minute/60.+time_first.second/3600.+time_first.microsecond/3600e6

#    fname = basepath+'Santis_alt_vs_time.png'
#    plot_scatter_comp(time_hour_first, alt_first, [fname], labelx='Hour [UTC]',
#                      labely='Flash origin altitude [m MSL]', titl='Flash origin altitude respect to time', axis=None,
#                      metadata=None, dpi=72)

    # Compute 2D histogram and
    bin_edges_time = np.arange(25.)
    bin_edges_alt = np.arange(-100., 14100, 200.)
    H, _, _ = np.histogram2d(time_hour_first, alt_first, bins=[bin_edges_time, bin_edges_alt])

    # Compute time histogram
    hist_time, _ = np.histogram(time_hour_first, bins=bin_edges_time)

    # divide 2D histogram by total number of samples at each time
    for i, val_time in enumerate(hist_time):
        H[i, :] = H[i, :]/val_time*100.

    H = H.T

    # Set 0 values to blank
    H = np.ma.asarray(H)
    H[H == 0] = np.ma.masked

    fig = plt.figure(figsize=(10, 6), dpi=72)
    ax = fig.add_subplot(111, title='Flash origin altitude respect to time')
    X, Y = np.meshgrid(bin_edges_time, bin_edges_alt)
    cax = ax.pcolormesh(X, Y, H, vmin=0., vmax=5.)

    plt.xlabel('Hour [UTC]')
    plt.ylabel('Flash origin altitude [m MSL]')

    cb = fig.colorbar(cax)
    cb.set_label('Frequency of Occurrence [per cent]')

    # Make a tight layout
    fig.tight_layout()
    fig.savefig(basepath+'Santis_alt_vs_time_density_plot.png', dpi=72)
    plt.close()

    print('Plotted '+basepath+'Santis_alt_vs_time_density_plot.png')


#    # Compute and plot histograms
#    # Get bins altitude
#    bin_edges = np.arange(25.)
#
#    # Plot histogram time of occurrence
#    fname_hist = basepath+'Santis_hist_time.png'
#    fname_hist = plot_histogram(
#        bin_edges, time_hour_first, [fname_hist], labelx='hour [UTC]',
#        titl='Flashes time of occurrence')
#    print('Plotted '+' '.join(fname_hist))
#
#    fname_hist = basepath+'Santis_hist_time.csv'
#    hist_time, _ = np.histogram(time_hour_first, bins=bin_edges)
#    fname_hist = write_histogram(bin_edges, hist_time, fname_hist)
#    print('Written '+' '.join(fname_hist))




#    # get sources in wet snow region
#    ind = np.ma.where(np.logical_or.reduce((
#        pol_vals_dict['hydro [-]'] == 3.,
#        pol_vals_dict['hydro [-]'] == 5.,
#        pol_vals_dict['hydro [-]'] == 7.,
#        pol_vals_dict['hydro [-]'] == 8.)))[0]
#
#    # Get unique flash that contain these sources
#    flashnr_liquid, unique_ind_liquid = np.unique(flashnr[ind], return_index=True)
#
#    ind = []
#    for flash in flashnr_liquid:
#        ind.extend(np.where(flashnr == flash)[0])
#
#    # Get sources of flashes propagating into liquid phase region
#    time_data_liquid = time_data[ind]
#    lat_liquid = lat[ind]
#    lon_liquid = lon[ind]
#    alt_liquid = alt[ind]
#    dBm_liquid = dBm[ind]
#    flashnr_liquid = flashnr[ind]
#    hydro_liquid = pol_vals_dict['hydro [-]'][ind]
#
#    print('non solid flash sources: '+str(lat_liquid.size))
#
#    # flash origin only
#    flashnr_first, unique_ind = np.unique(flashnr_liquid, return_index=True)
#    time_data_first = time_data_liquid[unique_ind]
#    lat_first = lat_liquid[unique_ind]
#    lon_first = lon_liquid[unique_ind]
#    alt_first = alt_liquid[unique_ind]
#    dBm_first = dBm_liquid[unique_ind]
#    hydro_first = hydro_liquid[unique_ind]
#
#    print('non solid phase flashes: '+str(lat_first.size))
#
#    # flashes with origin in rain only
#    ind_rain = np.ma.where(hydro_first == 8)[0]
#
#    print('melting hail origin flashes: '+str(ind_rain.size))
#
#    # Plots
#    figfname = basepath+'Santis_LMA_firstsources_mh_pos_min_height_on_top.png'
#    figfname = plot_pos(
#        lat_first[ind_rain], lon_first[ind_rain], alt_first[ind_rain], [figfname], sort_altitude='Lowest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='First flash sources position. Lowest on top\nFlashes with origin in melting hail only')
#    print('Plotted '+' '.join(figfname))
#
#    figfname = basepath+'Santis_LMA_firstsources_mh_pos_max_height_on_top.png'
#    figfname = plot_pos(
#        lat_first[ind_rain], lon_first[ind_rain], alt_first[ind_rain], [figfname], sort_altitude='Highest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='First flash sources position. Highest on top\nFlashes with origin in melting hail only')
#    print('Plotted '+' '.join(figfname))
#
#
#
#
#
#
#
#    # get sources in liquid, mixed phase, unclassified and no data regions
#    ind = np.ma.where(np.logical_or.reduce((
#        np.ma.getmaskarray(pol_vals_dict['hydro [-]']) == 1,
#        pol_vals_dict['hydro [-]'] == 0.,
#        pol_vals_dict['hydro [-]'] == 3.,
#        pol_vals_dict['hydro [-]'] == 5.,
#        pol_vals_dict['hydro [-]'] == 7.,
#        pol_vals_dict['hydro [-]'] == 8.)))[0]
#
#    # Get unique flash that contain these sources
#    flashnr_filter = np.unique(flashnr[ind], return_index=False)
#
#    # get Flash ID
#    flashnr_ID = np.unique(flashnr, return_index=False)
#
#    # get flashes not in filtered flashes
#    flashnr_solid = flashnr_ID[np.isin(flashnr_ID, flashnr_filter, assume_unique=True, invert=True)]
#
#    ind = []
#    for flash in flashnr_solid:
#        ind.extend(np.where(flashnr == flash)[0])
#
#    # Get sources of flashes propagating into liquid phase region
#    time_data_solid = time_data[ind]
#    lat_solid = lat[ind]
#    lon_solid = lon[ind]
#    alt_solid = alt[ind]
#    dBm_solid = dBm[ind]
#    flashnr_solid = flashnr[ind]
#
#    print('Solid phase flash sources: '+str(lat_solid.size))
#
#    # flash origin only
#    flashnr_first, unique_ind = np.unique(flashnr_solid, return_index=True)
#    time_data_first = time_data_solid[unique_ind]
#    lat_first = lat_solid[unique_ind]
#    lon_first = lon_solid[unique_ind]
#    alt_first = alt_solid[unique_ind]
#    dBm_first = dBm_solid[unique_ind]
#
#    print('Solid phase flashes: '+str(lat_first.size))
#
#    # Plots
#    figfname = basepath+'Santis_LMA_sources_solid_data_pos_min_height_on_top.png'
#    figfname = plot_pos(
#        lat_solid, lon_solid, alt_solid, [figfname], sort_altitude='Lowest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='Flash sources position. Lowest on top\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(figfname))
#
#    figfname = basepath+'Santis_LMA_firstsources_solid_data_pos_min_height_on_top.png'
#    figfname = plot_pos(
#        lat_first, lon_first, alt_first, [figfname], sort_altitude='Lowest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='First flash sources position. Lowest on top\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(figfname))
#
#    figfname = basepath+'Santis_LMA_firstsources_solid_data_pos_max_height_on_top.png'
#    figfname = plot_pos(
#        lat_first, lon_first, alt_first, [figfname], sort_altitude='Highest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='First flash sources position. Highest on top\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(figfname))
#
#
#    # Compute and plot histograms
#    # Get bins altitude
#    alt_min = alt_solid.min()
#    alt_max = alt_solid.max()
#    step = 100.
#    bins = np.linspace(alt_min-step/2., alt_max+step/2., num=int((alt_max-alt_min)/step)+2)
#
#    # Plot histogram altitude all sources
#    fname_hist = basepath+'solid_Santis_hist_alt.png'
#    fname_hist = plot_histogram(
#        bins, alt_solid, [fname_hist], labelx='Altitude [m MSL]',
#        titl='Flash sources altitude\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(fname_hist))
#
#    fname_hist = basepath+'solid_Santis_hist_alt.csv'
#    hist_alt, _ = np.histogram(alt_solid, bins=bins)
#    fname_hist = write_histogram(bins, hist_alt, fname_hist)
#    print('Written '+' '.join(fname_hist))
#
#    # Plot histogram altitude first sources
#    fname_hist = basepath+'solid_Santis_hist_alt_first_source.png'
#    fname_hist = plot_histogram(
#        bins, alt_first, [fname_hist], labelx='Altitude [m MSL]',
#        titl='Flash first source altitude\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(fname_hist))
#
#    fname_hist = basepath+'solid_Santis_hist_alt_fist_source.csv'
#    hist_alt_fist, _ = np.histogram(alt_first, bins=bins)
#    fname_hist = write_histogram(bins, hist_alt_fist, fname_hist)
#    print('Written '+' '.join(fname_hist))
#
#    # Get bins dBm
#    dBm_min = dBm_solid.min()
#    dBm_max = dBm_solid.max()
#    step = 1.
#    bins = np.linspace(dBm_min-step/2., dBm_max+step/2., num=int((dBm_max-dBm_min)/step)+2)
#
#    fname_hist = basepath+'solid_Santis_hist_dBm.png'
#    fname_hist = plot_histogram(
#        bins, dBm_solid, [fname_hist], labelx='Power [dBm]',
#        titl='Flash sources power\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(fname_hist))
#
#    fname_hist = basepath+'solid_Santis_hist_dBm.csv'
#    hist_dBm, _ = np.histogram(dBm_solid, bins=bins)
#    fname_hist = write_histogram(bins, hist_dBm, fname_hist)
#    print('Written '+' '.join(fname_hist))
#
#    # Plot histogram power first sources
#    fname_hist = basepath+'solid_Santis_hist_dBm_first_source.png'
#    fname_hist = plot_histogram(
#        bins, dBm_first, [fname_hist], labelx='Power [dBm]',
#        titl='Flash first source power\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(fname_hist))
#
#    fname_hist = basepath+'solid_Santis_hist_dBm_first_source.csv'
#    hist_dBm_first, _ = np.histogram(dBm_first, bins=bins)
#    fname_hist = write_histogram(bins, hist_dBm_first, fname_hist)
#    print('Written '+' '.join(fname_hist))
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
#
#    for i, key in enumerate(pol_vals_labels):
#        step = step_list[i]
#        datatype = datatype_vec[i]
#
#        field_name = get_fieldname_pyart(datatype)
#        field_dict = get_metadata(field_name)
#
#        labelx = get_colobar_label(field_dict, field_name)
#
#        vals = pol_vals_dict[key][ind]
#        bins, values = compute_histogram(vals, field_name, step=step)
#
#        # Plot all sources histogram
#        fname_first_source = (
#            basepath+'solid_allsources_ts_trajlightning_' +
#            datatype+'.png')
#        plot_histogram(
#            bins, values, [fname_first_source], labelx=labelx,
#            titl="Trajectory Histogram All Sources\nFlashes propagating into the solid phase only")
#
#        print("----- plot to '%s'" % fname_first_source)
#
#        # store histogram
#        fname_first_source = (
#            basepath+'solid_allsources_ts_trajlightning_' +
#            datatype+'.csv')
#        hist_values, _ = np.histogram(values, bins=bins)
#        write_histogram(bins, hist_values, fname_first_source)
#        print('Written '+fname_first_source)
#
#        # First sources
#        vals_first = vals[unique_ind]
#        bins, values = compute_histogram(vals_first, field_name, step=step)
#
#        # Plot first source histogram
#        fname_first_source = (
#            basepath+'solid_firstsource_ts_trajlightning_' +
#            datatype+'.png')
#        plot_histogram(
#            bins, values, [fname_first_source], labelx=labelx,
#            titl="Trajectory Histogram First Source\nFlashes propagating into the solid phase only")
#
#        print("----- plot to '%s'" % fname_first_source)
#
#        # store histogram
#        fname_first_source = (
#            basepath+'solid_firstsource_ts_trajlightning_' +
#            datatype+'.csv')
#        hist_values_first, _ = np.histogram(values, bins=bins)
#        write_histogram(bins, hist_values_first, fname_first_source)
#        print('Written '+fname_first_source)






















#    # get sources in liquid regions
#    ind = np.ma.where(np.logical_or.reduce((
#        pol_vals_dict['hydro [-]'] == 3.,
#        pol_vals_dict['hydro [-]'] == 5.,
#        pol_vals_dict['hydro [-]'] == 7.,
#        pol_vals_dict['hydro [-]'] == 8.)))[0]
#
#    # Get unique flash that contain these sources
#    flashnr_solid, unique_ind_solid = np.unique(flashnr[ind], return_index=True)
#
#    ind = []
#    for flash in flashnr_solid:
#        ind.extend(np.where(flashnr == flash)[0])
#
#    # Get sources of flashes propagating into liquid phase region
#    time_data_solid = time_data[ind]
#    lat_solid = lat[ind]
#    lon_solid = lon[ind]
#    alt_solid = alt[ind]
#    dBm_solid = dBm[ind]
#    flashnr_solid = flashnr[ind]
#
#    print('Liquid and mixed phase flash sources: '+str(lat_solid.size))
#
#    # flash origin only
#    flashnr_first, unique_ind = np.unique(flashnr_solid, return_index=True)
#    time_data_first = time_data_solid[unique_ind]
#    lat_first = lat_solid[unique_ind]
#    lon_first = lon_solid[unique_ind]
#    alt_first = alt_solid[unique_ind]
#    dBm_first = dBm_solid[unique_ind]
#
#    print('Liquid and mixed phase flashes: '+str(lat_first.size))
#
#    # Plots
#    figfname = basepath+'Santis_LMA_sources_solid_data_pos_min_height_on_top.png'
#    figfname = plot_pos(
#        lat_solid, lon_solid, alt_solid, [figfname], sort_altitude='Lowest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='Flash sources position. Lowest on top\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(figfname))
#
#    figfname = basepath+'Santis_LMA_firstsources_solid_data_pos_min_height_on_top.png'
#    figfname = plot_pos(
#        lat_first, lon_first, alt_first, [figfname], sort_altitude='Lowest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='First flash sources position. Lowest on top\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(figfname))
#
#    figfname = basepath+'Santis_LMA_firstsources_solid_data_pos_max_height_on_top.png'
#    figfname = plot_pos(
#        lat_first, lon_first, alt_first, [figfname], sort_altitude='Highest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='First flash sources position. Highest on top\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(figfname))
#
#
#    # Compute and plot histograms
#    # Get bins altitude
#    alt_min = alt_solid.min()
#    alt_max = alt_solid.max()
#    step = 100.
#    bins = np.linspace(alt_min-step/2., alt_max+step/2., num=int((alt_max-alt_min)/step)+2)
#
#    # Plot histogram altitude all sources
#    fname_hist = basepath+'solid_Santis_hist_alt.png'
#    fname_hist = plot_histogram(
#        bins, alt_solid, [fname_hist], labelx='Altitude [m MSL]',
#        titl='Flash sources altitude\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(fname_hist))
#
#    fname_hist = basepath+'solid_Santis_hist_alt.csv'
#    hist_alt, _ = np.histogram(alt_solid, bins=bins)
#    fname_hist = write_histogram(bins, hist_alt, fname_hist)
#    print('Written '+' '.join(fname_hist))
#
#    # Plot histogram altitude first sources
#    fname_hist = basepath+'solid_Santis_hist_alt_first_source.png'
#    fname_hist = plot_histogram(
#        bins, alt_first, [fname_hist], labelx='Altitude [m MSL]',
#        titl='Flash first source altitude\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(fname_hist))
#
#    fname_hist = basepath+'solid_Santis_hist_alt_fist_source.csv'
#    hist_alt_fist, _ = np.histogram(alt_first, bins=bins)
#    fname_hist = write_histogram(bins, hist_alt_fist, fname_hist)
#    print('Written '+' '.join(fname_hist))
#
#    # Get bins dBm
#    dBm_min = dBm_solid.min()
#    dBm_max = dBm_solid.max()
#    step = 1.
#    bins = np.linspace(dBm_min-step/2., dBm_max+step/2., num=int((dBm_max-dBm_min)/step)+2)
#
#    fname_hist = basepath+'solid_Santis_hist_dBm.png'
#    fname_hist = plot_histogram(
#        bins, dBm_solid, [fname_hist], labelx='Power [dBm]',
#        titl='Flash sources power\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(fname_hist))
#
#    fname_hist = basepath+'solid_Santis_hist_dBm.csv'
#    hist_dBm, _ = np.histogram(dBm_solid, bins=bins)
#    fname_hist = write_histogram(bins, hist_dBm, fname_hist)
#    print('Written '+' '.join(fname_hist))
#
#    # Plot histogram power first sources
#    fname_hist = basepath+'solid_Santis_hist_dBm_first_source.png'
#    fname_hist = plot_histogram(
#        bins, dBm_first, [fname_hist], labelx='Power [dBm]',
#        titl='Flash first source power\nFlashes propagating into the solid phase only')
#    print('Plotted '+' '.join(fname_hist))
#
#    fname_hist = basepath+'solid_Santis_hist_dBm_first_source.csv'
#    hist_dBm_first, _ = np.histogram(dBm_first, bins=bins)
#    fname_hist = write_histogram(bins, hist_dBm_first, fname_hist)
#    print('Written '+' '.join(fname_hist))
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
#
#    for i, key in enumerate(pol_vals_labels):
#        step = step_list[i]
#        datatype = datatype_vec[i]
#
#        field_name = get_fieldname_pyart(datatype)
#        field_dict = get_metadata(field_name)
#
#        labelx = get_colobar_label(field_dict, field_name)
#
#        vals = pol_vals_dict[key][ind]
#        bins, values = compute_histogram(vals, field_name, step=step)
#
#        # Plot all sources histogram
#        fname_first_source = (
#            basepath+'solid_allsources_ts_trajlightning_' +
#            datatype+'.png')
#        plot_histogram(
#            bins, values, [fname_first_source], labelx=labelx,
#            titl="Trajectory Histogram All Sources\nFlashes propagating into the solid phase only")
#
#        print("----- plot to '%s'" % fname_first_source)
#
#        # store histogram
#        fname_first_source = (
#            basepath+'solid_allsources_ts_trajlightning_' +
#            datatype+'.csv')
#        hist_values, _ = np.histogram(values, bins=bins)
#        write_histogram(bins, hist_values, fname_first_source)
#        print('Written '+fname_first_source)
#
#        # First sources
#        vals_first = vals[unique_ind]
#        bins, values = compute_histogram(vals_first, field_name, step=step)
#
#        # Plot first source histogram
#        fname_first_source = (
#            basepath+'solid_firstsource_ts_trajlightning_' +
#            datatype+'.png')
#        plot_histogram(
#            bins, values, [fname_first_source], labelx=labelx,
#            titl="Trajectory Histogram First Source\nFlashes propagating into the solid phase only")
#
#        print("----- plot to '%s'" % fname_first_source)
#
#        # store histogram
#        fname_first_source = (
#            basepath+'solid_firstsource_ts_trajlightning_' +
#            datatype+'.csv')
#        hist_values_first, _ = np.histogram(values, bins=bins)
#        write_histogram(bins, hist_values_first, fname_first_source)
#        print('Written '+fname_first_source)
#
#
#
#    # get sources without radar data
#    mask_hydro = np.ma.getmaskarray(pol_vals_dict['hydro [-]'])
#
#    figfname = basepath+'Santis_LMA_sources_no_radar_data_pos_min_height_on_top.png'
#    figfname = plot_pos(
#        lat[mask_hydro], lon[mask_hydro], alt[mask_hydro], [figfname], sort_altitude='Lowest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='Flash sources with no radar data position. Lowest on top')
#    print('Plotted '+' '.join(figfname))
#
#    # get sources with hydroclass no classification
#    ind_nc = np.ma.where(pol_vals_dict['hydro [-]'] == 0.)[0]
#
#    print('Non classified flash sources: '+str(ind_nc.size))
#
#    figfname = basepath+'Santis_LMA_sources_no_class_data_pos_min_height_on_top.png'
#    figfname = plot_pos(
#        lat[ind_nc], lon[ind_nc], alt[ind_nc], [figfname], sort_altitude='Lowest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='Flash sources classified as NC position. Lowest on top')
#    print('Plotted '+' '.join(figfname))
#
#    # flash origin only
#    flashnr_first, unique_ind = np.unique(flashnr, return_index=True)
#    lat_first = lat[unique_ind]
#    lon_first = lon[unique_ind]
#    alt_first = alt[unique_ind]
#    hydro_first = pol_vals_dict['hydro [-]'][unique_ind]
#
#    ind_nc = np.ma.where(hydro_first == 0.)[0]
#
#    print('Non classified flash origin: '+str(ind_nc.size))
#
#    figfname = basepath+'Santis_LMA_firstsources_no_class_data_pos_min_height_on_top.png'
#    figfname = plot_pos(
#        lat_first[ind_nc], lon_first[ind_nc], alt_first[ind_nc], [figfname], sort_altitude='Lowest_on_top',
#        cb_label='Source height [m MSL]',
#        titl='First flash sources classified as NC position. Lowest on top')
#    print('Plotted '+' '.join(figfname))
#
#    # Plot data histograms
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
#
#    for i, key in enumerate(pol_vals_labels):
#        vals_first = pol_vals_dict[key][unique_ind]
#
#        step = step_list[i]
#        datatype = datatype_vec[i]
#
#        field_name = get_fieldname_pyart(datatype)
#        field_dict = get_metadata(field_name)
#
#        labelx = get_colobar_label(field_dict, field_name)
#
#        bins, values = compute_histogram(vals_first, field_name, step=step)
#
#        # Plot first source histogram
#        fname_first_source = (
#            basepath+'firstsource_ts_trajlightning_' +
#            datatype+'.png')
#        plot_histogram(
#            bins, values, [fname_first_source], labelx=labelx,
#            titl="Trajectory Histogram First Source")
#
#        print("----- plot to '%s'" % fname_first_source)
#
#        # store histogram
#        fname_first_source = (
#            basepath+'firstsource_ts_trajlightning_' +
#            datatype+'.csv')
#        hist_values_first, _ = np.histogram(values, bins=bins)
#        write_histogram(bins, hist_values_first, fname_first_source)
#        print('Written '+fname_first_source)



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
