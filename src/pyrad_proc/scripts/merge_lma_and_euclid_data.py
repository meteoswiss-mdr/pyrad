#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_euclid_data
================================================

This program reads EUCLID raw data and plots several features such as sources
position, histograms, etc.

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit

import numpy as np
import shapely

from pyrad.io import read_meteorage, read_lightning_all, write_ts_lightning
from pyrad.io import write_histogram
from pyrad.util import belongs_roi_indices, compute_histogram
from pyrad.graph import plot_pos, plot_histogram

from pyart.core import wgs84_to_swissCH1903

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
        '--euclid_basepath', type=str,
        default='/store/msrad/lightning/meteorage/',
        help='name of folder containing the EUCLID lightning data')

    parser.add_argument(
        '--lma_basepath', type=str,
        default='/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/',
        help='name of folder containing the LMA lightning data')

    parser.add_argument(
        '--lma_basename', type=str,
        default='Santis_data',
        help='base name of the files containing the LMA lightning data')

    parser.add_argument(
        '--datatypes', type=str,
        default='hydro,KDPc,dBZc,RhoHVc,TEMP,ZDRc',
        help='Name of the polarimetric moments to process. Coma separated')

    parser.add_argument(
        '--labels', type=str,
        default=('hydro [-],KDPc [deg/Km],dBZc [dBZ],RhoHVc [-],' +
                 'TEMP [deg C],ZDRc [dB]'),
        help='Labels in the csv file for each polarimetric variable')

    parser.add_argument(
        '--nsources_min', type=int, default=10,
        help='Minimum number of sources to consider the LMA flash valid')

    parser.add_argument(
        '--scale_factor', type=float, default=1.2,
        help='Factor by which the area covered by the LMA flash has to be ' +
             'enlarged to find EUCLID strokes')

    parser.add_argument(
        '--delay', type=float, default=100000.,
        help='delay after end of LMA flash where to look for EUCLID strokes [micros]')

    parser.add_argument(
        '--anticipation', type=float, default=100000.,
        help='anticipation of the start of an LMA flash where to look for EUCLID strokes [micros]')

    parser.add_argument(
        '--min_area', type=float, default=25.,
        help='Minimum size of the area where to look for an EUCLID stroke [km]')

    parser.add_argument(
        '--euclidtype', type=str, default='CGt',
        help='Type of Euclid stroke. Can be: CGt, CGp, CGn, IC')

    args = parser.parse_args()

    print("====== EUCLID data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== EUCLID data processing finished: ")

    day_vec = []
    for day in args.days:
        day_vec.append(datetime.datetime.strptime(day, '%Y%m%d'))

    datatype_vec = args.datatypes.split(',')
    pol_vals_labels = args.labels.split(',')

    if np.size(datatype_vec) != np.size(pol_vals_labels):
        warn(
            str(np.size(datatype_vec))+' datatypes but ' +
            str(np.size(pol_vals_labels)) +
            ' labels. Their number must be equal')
        return

    flashnr_sel_all = np.ma.asarray([], dtype=int)
    t_start_sel_all = np.ma.asarray([], dtype=datetime.datetime)
    t_end_sel_all = np.ma.asarray([], dtype=datetime.datetime)
    area_sel_all = np.ma.asarray([])
    neuclid_sel_all = np.ma.asarray([], dtype=int)
    intens_EU_sel_all = np.ma.asarray([])

    for day in day_vec:
        day_str = day.strftime('%Y%m%d')
        fname_lma = args.lma_basepath+day_str+'_'+args.lma_basename+'.csv'
        fname_euclid = args.euclid_basepath+'THX/THX'+day.strftime('%y%j0000')+'.prd'

        print('Reading EUCLID data file '+fname_euclid)
        (stroke_time, lon_euclid, lat_euclid, intens, ns, mode, intra, ax, ki2, ecc, incl,
         sind) = read_meteorage(fname_euclid)

        print('Reading LMA data file '+fname_lma)
        flashnr, time_data, time_in_flash, lat_lma, lon_lma, alt_lma, dBm, pol_vals_dict = (
            read_lightning_all(fname_lma, labels=pol_vals_labels))

        flashnr_first = np.unique(flashnr, return_index=False)

        # keep strokes of interest
        if args.euclidtype == 'CGt' or args.euclidtype == 'CGp' or args.euclidtype == 'CGn':
            time_EU_filt = stroke_time[intra == 0]
            lon_EU_filt = lon_euclid[intra == 0]
            lat_EU_filt = lat_euclid[intra == 0]
            intens_EU_filt = intens[intra == 0]

            if args.euclidtype == 'CGp':
                time_EU_filt = time_EU_filt[intens_EU_filt > 0.]
                lon_EU_filt = lon_EU_filt[intens_EU_filt > 0.]
                lat_EU_filt = lat_EU_filt[intens_EU_filt > 0.]
                intens_EU_filt = intens_EU_filt[intens_EU_filt > 0.]
            elif args.euclidtype == 'CGn':
                time_EU_filt = time_EU_filt[intens_EU_filt < 0.]
                lon_EU_filt = lon_EU_filt[intens_EU_filt < 0.]
                lat_EU_filt = lat_EU_filt[intens_EU_filt < 0.]
                intens_EU_filt = intens_EU_filt[intens_EU_filt < 0.]
        else:
            intens_EU_filt = intens[intra == 1]
            time_EU_filt = stroke_time[intra == 1]
            lon_EU_filt = lon_euclid[intra == 1]
            lat_EU_filt = lat_euclid[intra == 1]

        # get Swiss coordinates
        chy_EU, chx_EU, _ = wgs84_to_swissCH1903(
            lon_EU_filt, lat_EU_filt, np.zeros(lon_EU_filt.size), no_altitude_transform=True)

        chy_lma, chx_lma, _ = wgs84_to_swissCH1903(
            lon_lma, lat_lma, alt_lma, no_altitude_transform=True)

        flashnr_filt = np.asarray([], dtype=int)
        time_data_filt = np.asarray([], dtype=datetime.datetime)
        time_in_flash_filt = np.asarray([], dtype=float)
        lat_filt = np.asarray([], dtype=float)
        lon_filt = np.asarray([], dtype=float)
        alt_filt = np.asarray([], dtype=float)
        dBm_filt = np.asarray([], dtype=float)

        pol_vals_dict_filt = dict()
        for label in pol_vals_labels:
            pol_vals_dict_filt.update({label: np.ma.asarray([])})
        nflashes_time_rejected = 0
        nflashes_area_rejected = 0
        nflashes_insufficient_sources = 0
        nflashes_small_area_rejected = 0
        nstrokes_accepted = 0

        flashnr_sel = np.ma.asarray([], dtype=int)
        t_start_sel = np.ma.asarray([], dtype=datetime.datetime)
        t_end_sel = np.ma.asarray([], dtype=datetime.datetime)
        area_sel = np.ma.asarray([])
        neuclid_sel = np.ma.asarray([], dtype=int)
        intens_EU_sel = np.ma.asarray([])

        for flash_ID in flashnr_first:
            # get LMA data of flash
            flashnr_flash = flashnr[flashnr == flash_ID]
            time_data_flash = time_data[flashnr == flash_ID]
            time_in_flash_flash = time_in_flash[flashnr == flash_ID]
            lat_flash = lat_lma[flashnr == flash_ID]
            lon_flash = lon_lma[flashnr == flash_ID]
            alt_flash = alt_lma[flashnr == flash_ID]
            dBm_flash = dBm[flashnr == flash_ID]
            pol_vals_dict_flash = dict()
            for key in pol_vals_dict.keys():
                pol_vals_dict_flash.update({key: pol_vals_dict[key][flashnr == flash_ID]})

            if flashnr_flash.size < args.nsources_min:
                # print('Not enough sources for flash '+str(flash_ID))
                nflashes_insufficient_sources +=1
                continue

            chy_flash = chy_lma[flashnr == flash_ID]
            chx_flash = chx_lma[flashnr == flash_ID]

            # check if there are EUCLID strokes within LMA time
            t_start = time_data_flash[0]-datetime.timedelta(microseconds=args.anticipation)
            t_end = time_data_flash[-1]+datetime.timedelta(microseconds=args.delay)
            ind = np.where(np.logical_and(
                time_EU_filt >= t_start, time_EU_filt <=t_end))[0]
            if ind.size == 0:
                # print('No EUCLID '+args.euclidtype+' flashes within time of LMA flash '+str(flash_ID))
                nflashes_time_rejected += 1
                continue

            # check if there are EUCLID strokes within LMA area
            lon_EU_flash = lon_EU_filt[ind]
            lat_EU_flash = lat_EU_filt[ind]
            intens_EU_flash = intens_EU_filt[ind]

            chy_EU_flash = chy_EU[ind]
            chx_EU_flash = chx_EU[ind]

            points_flash_lma = shapely.geometry.MultiPoint(
                list(zip(chy_flash, chx_flash)))
            points_EU = shapely.geometry.MultiPoint(
                list(zip(chy_EU_flash, chx_EU_flash)))

            rectangle_lma = points_flash_lma.minimum_rotated_rectangle
            # print('LMA area before scaling : '+str(rectangle_lma.area*1e-6))

            if rectangle_lma.area*1e-6 == 0:
                # print('LMA area too small for flash '+str(flash_ID))
                nflashes_small_area_rejected += 1

                # Plot position of LMA sources AND EUCLID stroke
                types = np.zeros(lat_flash.size)

                figfname = args.lma_basepath+'rejected_'+day_str+'_'+str(flash_ID)+'_LMA_EUCLID_'+args.euclidtype+'_pos.png'
                figfname = plot_pos(
                    lat_flash, lon_flash, types, [figfname],
                    cb_label='Type of detection 0: LMA, 1: EUCLID',
                    titl=day_str+' '+str(flash_ID)+' rejected LMA flash\nLMA and EUCLID '+args.euclidtype+' positions')
                print('Plotted '+' '.join(figfname))

                continue

            roi_lma = shapely.affinity.scale(
                rectangle_lma, xfact=args.scale_factor, yfact=args.scale_factor)

            area_roi = roi_lma.area*1e-6
            if area_roi < args.min_area:
                scale_factor = args.scale_factor
                while area_roi < args.min_area:
                    scale_factor += 0.1
                    roi_lma = shapely.affinity.scale(
                        rectangle_lma, xfact=scale_factor, yfact=scale_factor)
                    area_roi = roi_lma.area*1e-6
                # print('scale_factor: '+str(scale_factor))
            # print('LMA area after scaling : '+str(roi_lma.area*1e-6))

            if roi_lma.disjoint(points_EU):
                # print('No EUCLID '+args.euclidtype+' flashes within area of LMA flash '+str(flash_ID))
                nflashes_area_rejected += 1

#                # Plot position of LMA sources AND EUCLID stroke
#                lats = np.append(lat_flash, lat_EU_flash)
#                lons = np.append(lon_flash, lon_EU_flash)
#                types = np.append(np.zeros(lat_flash.size), np.ones(lat_EU_flash.size))
#
#                figfname = args.lma_basepath+'rejected_'+day_str+'_'+str(flash_ID)+'_LMA_EUCLID_'+args.euclidtype+'_pos.png'
#                figfname = plot_pos(
#                    lats, lons, types, [figfname],
#                    cb_label='Type of detection 0: LMA, 1: EUCLID',
#                    titl=day_str+' '+str(flash_ID)+' rejected LMA flash\nLMA and EUCLID '+args.euclidtype+' positions')
#                print('Plotted '+' '.join(figfname))

                continue

            flashnr_filt = np.append(flashnr_filt, flashnr_flash)
            time_data_filt = np.append(time_data_filt, time_data_flash)
            time_in_flash_filt = np.append(time_in_flash_filt, time_in_flash_flash)
            lat_filt = np.append(lat_filt, lat_flash)
            lon_filt = np.append(lon_filt, lon_flash)
            alt_filt = np.append(alt_filt, alt_flash)
            dBm_filt = np.append(dBm_filt, dBm_flash)
            for key in pol_vals_dict_filt.keys():
                pol_vals_dict_filt[key] = np.ma.append(
                    pol_vals_dict_filt[key], pol_vals_dict_flash[key])

            if not roi_lma.contains(points_EU):
                points_EU = points_EU.intersection(roi_lma)

                inds = []
                if points_EU.geom_type == 'Point':
                    ind = np.where(np.logical_and(
                        chy_EU_flash == points_EU.x, chx_EU_flash == points_EU.y))
                    if len(ind) == 1:
                        ind = ind[0]
                    inds.extend(ind)
                else:
                    points_EU_list = list(points_EU)
                    for point in points_EU_list:
                        ind = np.where(np.logical_and(
                            chy_EU_flash == point.x, chx_EU_flash == point.y))
                        if len(ind) == 1:
                            ind = ind[0]
                        inds.extend(ind)
                lon_EU_flash = lon_EU_flash[inds]
                lat_EU_flash = lat_EU_flash[inds]
                intens_EU_flash = intens_EU_flash[inds]

            print(str(lon_EU_flash.size)+' EUCLID '+args.euclidtype+' flashes for LMA flash '+str(flash_ID))
            nstrokes_accepted += lon_EU_flash.size

#            # Plot position of LMA sources AND EUCLID stroke
#            lats = np.append(lat_flash, lat_EU_flash)
#            lons = np.append(lon_flash, lon_EU_flash)
#            types = np.append(np.zeros(lat_flash.size), np.ones(lat_EU_flash.size))
#
#            figfname = args.lma_basepath+day_str+'_'+str(flash_ID)+'_LMA_EUCLID_'+args.euclidtype+'_pos.png'
#            figfname = plot_pos(
#                lats, lons, types, [figfname],
#                cb_label='Type of detection 0: LMA, 1: EUCLID',
#                titl=day_str+' '+str(flash_ID)+' LMA flash\nLMA and EUCLID '+args.euclidtype+' positions')
#            print('Plotted '+' '.join(figfname))

            flashnr_sel = np.append(flashnr_sel, flash_ID)
            t_start_sel = np.append(t_start_sel, time_data_flash[0])
            t_end_sel = np.append(t_end_sel, time_data_flash[-1])
            area_sel = np.append(area_sel, rectangle_lma.area*1e-6)
            neuclid_sel = np.append(neuclid_sel, lon_EU_flash.size)
            intens_EU_sel = np.append(intens_EU_sel, intens_EU_flash)

        # Plot histogram of number of EUCLID strokes
        bins_edges = np.arange(-0.5, 21.5, 1)
        fname_hist = args.lma_basepath+day_str+'_'+args.lma_basename+'_n'+args.euclidtype+'_hist.png'
        fname_hist = plot_histogram(
            bins_edges, neuclid_sel, [fname_hist], labelx='Number of EUCLID strokes per LMA flash',
            titl=day_str+' EUCLID '+args.euclidtype+' strokes per LMA flash')
        print('Plotted '+' '.join(fname_hist))

        fname_hist = args.lma_basepath+day_str+'_'+args.lma_basename+'_n'+args.euclidtype+'_hist.csv'
        _, hist_neuclid = compute_histogram(neuclid_sel, None, bin_edges=bins_edges)
        hist_neuclid, _ = np.histogram(hist_neuclid, bins=bins_edges)
        fname_hist = write_histogram(bins_edges, hist_neuclid, fname_hist)
        print('Written '+fname_hist)

        # Plot histogram of LMA flash area
        bins_edges = np.arange(0., 2010., 10.)
        fname_hist = args.lma_basepath+day_str+'_'+args.lma_basename+'_'+args.euclidtype+'_LMA_area_hist.png'
        fname_hist = plot_histogram(
            bins_edges, area_sel, [fname_hist], labelx='Area of LMA flash [km2]',
            titl=day_str+' area of LMA flash with associated EUCLID '+args.euclidtype+' strokes')
        print('Plotted '+' '.join(fname_hist))

        fname_hist = args.lma_basepath+day_str+'_'+args.lma_basename+'_'+args.euclidtype+'_LMA_area_hist.csv'
        _, hist_area = compute_histogram(area_sel, None, bin_edges=bins_edges)
        hist_area, _ = np.histogram(hist_area, bins=bins_edges)
        fname_hist = write_histogram(bins_edges, hist_area, fname_hist)
        print('Written '+fname_hist)

        # Plot histogram of LMA flash duration [milliseconds]
        duration = np.ma.zeros(t_end_sel.size)
        for i in range(t_end_sel.size):
            duration[i] = 1e3*(t_end_sel[i]-t_start_sel[i]).total_seconds()

        bins_edges = np.arange(0., 1010., 10.)
        fname_hist = args.lma_basepath+day_str+'_'+args.lma_basename+'_'+args.euclidtype+'_LMA_duration_hist.png'
        fname_hist = plot_histogram(
            bins_edges, duration, [fname_hist], labelx='Duration of LMA flash [ms]',
            titl=day_str+' Duration of LMA flash with associated EUCLID '+args.euclidtype+' strokes')
        print('Plotted '+' '.join(fname_hist))

        fname_hist = args.lma_basepath+day_str+'_'+args.lma_basename+'_'+args.euclidtype+'_LMA_duration_hist.csv'
        _, hist_duration = compute_histogram(duration, None, bin_edges=bins_edges)
        hist_duration, _ = np.histogram(hist_duration, bins=bins_edges)
        fname_hist = write_histogram(bins_edges, hist_duration, fname_hist)
        print('Written '+fname_hist)

        # Plot histogram of intensity of strokes [kA]
        bins_edges = np.arange(-100., 101., 1.)
        fname_hist = args.lma_basepath+day_str+'_'+args.lma_basename+'_'+args.euclidtype+'_EUCLID_intensity_hist.png'
        fname_hist = plot_histogram(
            bins_edges, intens_EU_sel, [fname_hist], labelx='Intensity of EUCLID strokes [kA]',
            titl=day_str+' Intensity of EUCLID '+args.euclidtype+' strokes')
        print('Plotted '+' '.join(fname_hist))

        fname_hist = args.lma_basepath+day_str+'_'+args.lma_basename+'_'+args.euclidtype+'_EUCLID_intensity_hist.csv'
        _, hist_intensity = compute_histogram(intens_EU_sel, None, bin_edges=bins_edges)
        hist_intensity, _ = np.histogram(hist_intensity, bins=bins_edges)
        fname_hist = write_histogram(bins_edges, hist_intensity, fname_hist)
        print('Written '+fname_hist)


        flashnr_filt_first = np.unique(flashnr_filt, return_index=False)
        print('N EUCLID strokes: '+str(time_EU_filt.size))
        print('N EUCLID strokes in LMA flashes accepted: '+str(nstrokes_accepted)+'\n')

        print('N LMA flashes: '+str(flashnr_first.size))
        print('N LMA flashes insufficient sources: '+str(nflashes_insufficient_sources))
        print('N LMA flashes time rejected: '+str(nflashes_time_rejected))
        print('N LMA flashes small area rejected: '+str(nflashes_small_area_rejected))
        print('N LMA flashes area rejected: '+str(nflashes_area_rejected))
        print('N LMA flashes accepted: '+str(flashnr_filt_first.size))
        print('N LMA sources accepted: '+str(flashnr_filt.size)+'\n\n\n')


    #    # write the results in a file
    #    vals_list = []
    #    for label in pol_vals_labels:
    #        vals_list.append(pol_vals_dict_filt[label])
    #
    #    fname = args.lma_basepath+day_str+'_'+args.lma_basename+'_'+args.euclidtype+'.csv'
    #    write_ts_lightning(
    #        flashnr_filt, time_data_filt, time_in_flash_filt, lat_filt, lon_filt,
    #        alt_filt, dBm_filt, vals_list, fname, pol_vals_labels)
    #    print('written to '+fname)

        flashnr_sel_all = np.append(flashnr_sel_all, flashnr_sel)
        t_start_sel_all = np.append(t_start_sel_all, t_start_sel)
        t_end_sel_all = np.append(t_end_sel_all, t_end_sel)
        area_sel_all = np.append(area_sel_all, area_sel)
        neuclid_sel_all = np.append(neuclid_sel_all, neuclid_sel)
        intens_EU_sel_all = np.append(intens_EU_sel_all, intens_EU_sel)

    # Plot histogram of number of EUCLID strokes
    bins_edges = np.arange(-0.5, 21.5, 1)
    fname_hist = args.lma_basepath+args.lma_basename+'_n'+args.euclidtype+'_hist.png'
    fname_hist = plot_histogram(
        bins_edges, neuclid_sel_all, [fname_hist], labelx='Number of EUCLID strokes per LMA flash',
        titl='EUCLID '+args.euclidtype+' strokes per LMA flash')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.lma_basepath+args.lma_basename+'_n'+args.euclidtype+'_hist.csv'
    _, hist_neuclid = compute_histogram(neuclid_sel_all, None, bin_edges=bins_edges)
    hist_neuclid, _ = np.histogram(hist_neuclid, bins=bins_edges)
    fname_hist = write_histogram(bins_edges, hist_neuclid, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram of LMA flash area
    bins_edges = np.arange(0., 2010., 10.)
    fname_hist = args.lma_basepath+args.lma_basename+'_'+args.euclidtype+'_LMA_area_hist.png'
    fname_hist = plot_histogram(
        bins_edges, area_sel_all, [fname_hist], labelx='Area of LMA flash [km2]',
        titl='area of LMA flash with associated EUCLID '+args.euclidtype+' strokes')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.lma_basepath+args.lma_basename+'_'+args.euclidtype+'_LMA_area_hist.csv'
    _, hist_area = compute_histogram(area_sel_all, None, bin_edges=bins_edges)
    hist_area, _ = np.histogram(hist_area, bins=bins_edges)
    fname_hist = write_histogram(bins_edges, hist_area, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram of LMA flash duration [milliseconds]
    duration = np.ma.zeros(t_end_sel.size)
    for i in range(t_end_sel.size):
        duration[i] = 1e3*(t_end_sel_all[i]-t_start_sel_all[i]).total_seconds()
    print(duration[i])
    bins_edges = np.arange(0., 1010., 10.)
    fname_hist = args.lma_basepath+args.lma_basename+'_'+args.euclidtype+'_LMA_duration_hist.png'
    fname_hist = plot_histogram(
        bins_edges, duration, [fname_hist], labelx='Duration of LMA flash [ms]',
        titl='Duration of LMA flash with associated EUCLID '+args.euclidtype+' strokes')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.lma_basepath+args.lma_basename+'_'+args.euclidtype+'_LMA_duration_hist.csv'
    _, hist_duration = compute_histogram(duration, None, bin_edges=bins_edges)
    hist_duration, _ = np.histogram(hist_duration, bins=bins_edges)
    fname_hist = write_histogram(bins_edges, hist_duration, fname_hist)
    print('Written '+fname_hist)

    # Plot histogram of intensity of strokes [kA]
    bins_edges = np.arange(-100., 101., 1.)
    fname_hist = args.lma_basepath+args.lma_basename+'_'+args.euclidtype+'_EUCLID_intensity_hist.png'
    fname_hist = plot_histogram(
        bins_edges, intens_EU_sel_all, [fname_hist], labelx='Intensity of EUCLID strokes [kA]',
        titl='Intensity of EUCLID '+args.euclidtype+' strokes')
    print('Plotted '+' '.join(fname_hist))

    fname_hist = args.lma_basepath+args.lma_basename+'_'+args.euclidtype+'_EUCLID_intensity_hist.csv'
    _, hist_intensity = compute_histogram(intens_EU_sel_all, None, bin_edges=bins_edges)
    hist_intensity, _ = np.histogram(hist_intensity, bins=bins_edges)
    fname_hist = write_histogram(bins_edges, hist_intensity, fname_hist)
    print('Written '+fname_hist)


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
