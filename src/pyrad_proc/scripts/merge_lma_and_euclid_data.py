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
from pyrad.util import belongs_roi_indices
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
        default='/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/data_analysis/',
        help='name of folder containing the EUCLID lightning data')

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
        '--scale_factor', type=float, default=1.2,
        help='Factor by which the area covered by the LMA flash has to be ' +
             'enlarged to find EUCLID strokes')

    parser.add_argument(
        '--delay', type=float, default=1000.,
        help='delay after end of LMA flash where to look for EUCLID strokes [micros]')

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

    for day in day_vec:
        day_str = day.strftime('%Y%m%d')
        fname_lma = args.lma_basepath+day_str+'_Santis_data.csv'
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

            if args.euclidtype == 'CGp':
                intens_EU_filt = intens[intra == 0]

                time_EU_filt = time_EU_filt[intens_EU_filt > 0.]
                lon_EU_filt = lon_EU_filt[intens_EU_filt > 0.]
                lat_EU_filt = lat_EU_filt[intens_EU_filt > 0.]
            elif args.euclidtype == 'CGn':
                intens_EU_filt = intens[intra == 0]

                time_EU_filt = time_EU_filt[intens_EU_filt < 0.]
                lon_EU_filt = lon_EU_filt[intens_EU_filt < 0.]
                lat_EU_filt = lat_EU_filt[intens_EU_filt < 0.]
        else:
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
        pol_vals_dict_filt = {
            'hydro [-]': np.ma.asarray([], dtype=np.int16),
            'KDPc [deg/Km]': np.ma.asarray([], dtype=float),
            'dBZc [dBZ]': np.ma.asarray([], dtype=float),
            'RhoHVc [-]': np.ma.asarray([], dtype=float),
            'TEMP [deg C]': np.ma.asarray([], dtype=float),
            'ZDRc [dB]': np.ma.asarray([], dtype=float)
        }
        nflashes_time_rejected = 0
        nflashes_area_rejected = 0
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

            chy_flash = chy_lma[flashnr == flash_ID]
            chx_flash = chx_lma[flashnr == flash_ID]

            # check if there are EUCLID strokes within LMA time
            t_start = time_data_flash[0]
            t_end = time_data_flash[-1]+datetime.timedelta(microseconds=args.delay)
            ind = np.where(np.logical_and(
                time_EU_filt >= t_start, time_EU_filt <=t_end))[0]
            if ind.size == 0:
                print('No EUCLID '+args.euclidtype+' flashes within time of LMA flash '+str(flash_ID))
                nflashes_time_rejected += 1
                continue

            # check if there are EUCLID strokes within LMA area
            lon_EU_flash = lon_EU_filt[ind]
            lat_EU_flash = lat_EU_filt[ind]

            chy_EU_flash = chy_EU[ind]
            chx_EU_flash = chx_EU[ind]

            points_flash_lma = shapely.geometry.MultiPoint(
                list(zip(chy_flash, chx_flash)))
            points_EU = shapely.geometry.MultiPoint(
                list(zip(chy_EU_flash, chx_EU_flash)))
            roi_lma = shapely.affinity.scale(
                points_flash_lma.minimum_rotated_rectangle,
                xfact=args.scale_factor, yfact=args.scale_factor)

            if roi_lma.disjoint(points_EU):
                print('No EUCLID '+args.euclidtype+' flashes within area of LMA flash '+str(flash_ID))
                nflashes_area_rejected += 1
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

            if roi_lma.contains(points_EU):
                print(str(lon_EU_flash.size)+' EUCLID '+args.euclidtype+' flashes for LMA flash '+str(flash_ID))
                continue

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

            print(str(lon_EU_flash.size)+' EUCLID '+args.euclidtype+' flashes for LMA flash '+str(flash_ID))

        flashnr_first = np.unique(flashnr_filt, return_index=False)
        print('N EUCLID strokes: '+str(time_EU_filt.size))
        print('N LMA flashes time rejected: '+str(nflashes_time_rejected))
        print('N LMA flashes area rejected: '+str(nflashes_area_rejected))
        print('N LMA flashes: '+str(flashnr_first.size))
        print('N LMA sources: '+str(flashnr_filt.size))


        # write the results in a file
        vals_list = []
        for label in pol_vals_labels:
            vals_list.append(pol_vals_dict_filt[label])

        fname = args.lma_basepath+day_str+'_Santis_data_'+args.euclidtype+'.csv'
        write_ts_lightning(
            flashnr_filt, time_data_filt, time_in_flash_filt, lat_filt, lon_filt,
            alt_filt, dBm_filt, vals_list, fname, pol_vals_labels)
        print('written to '+fname)
        

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
