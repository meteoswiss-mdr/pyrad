#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
merge_lma_and_radar_data
================================================

This program merges the original LMA data with the
co-located radar data and puts it in a new file

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import glob
from warnings import warn

import numpy as np

from pyrad.io import read_lightning, read_lightning_traj
from pyrad.io import write_ts_lightning

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
        default='/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/',
        help='name of folder containing the radar data')

    parser.add_argument(
        '--basepath_lma', type=str,
        default='/store/msrad/lightning/LMA/Santis/',
        help='name of folder containing the original LMA data')

    parser.add_argument(
        '--datatypes', type=str,
        default='hydro,KDPc,dBZc,RhoHVc,TEMP,ZDRc',
        help='Name of the polarimetric moments to process. Coma separated')

    parser.add_argument(
        '--datasets', type=str,
        default='hydroclass,KDPc,reflectivity,RhoHVc,temperature,ZDRc',
        help='Name of the datasets (directories where the data is stored. Coma separated')

    parser.add_argument(
        '--labels', type=str,
        default=('hydro [-],KDPc [deg/Km],dBZc [dBZ],RhoHVc [-],' +
                 'TEMP [deg C],ZDRc [dB]'),
        help='Labels in the csv file for each polarimetric variable')

    args = parser.parse_args()

    print("====== Merging of LMA and radar data started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Merging of LMA and radar data started finished: ")

    day_vec = []
    for day in args.days:
        day_vec.append(datetime.datetime.strptime(day, '%Y%m%d'))

    datatype_vec = args.datatypes.split(',')
    dataset_vec = args.datasets.split(',')
    pol_vals_labels = args.labels.split(',')

    if (np.size(datatype_vec) != np.size(pol_vals_labels) or
            np.size(datatype_vec) != np.size(dataset_vec)):
        warn(
            str(np.size(datatype_vec))+' datatypes but ' +
            str(np.size(pol_vals_labels))+' labels and '+
            str(np.size(dataset_vec))+' datasets. Their number must be equal')
        return

    for day in day_vec:
        day_str = day.strftime('%y%m%d')
        fname = args.basepath_lma+day_str+'.txt'

        print('\nReading LMA data file '+fname)
        flashnrs, times, times_in_flash, lats, lons, alts, dBms = (
            read_lightning(fname))

        nsources_file = flashnrs.size

        # find unique flash numbers
        flashnrs_first = np.unique(flashnrs)

        vals_list = []
        for datatype, dataset in zip(datatype_vec, dataset_vec):
            day_dir = day.strftime('%Y-%m-%d')
            day_str = day.strftime('%Y%m%d')

            fname_test = (
                args.basepath+day_dir+'/'+dataset+'_traj/AT_FLASH/'+day_str +
                '*_allflash_ts_trajlightning_'+datatype+'.csv')
            fname_list = glob.glob(fname_test)
            if not fname_list:
                warn('No file found in '+fname_test)
                continue

            fname = fname_list[0]

            print('\nReading file '+fname)
            times_rad, flashnrs_rad, dBms_rad, vals_rad, _, _, _, _ = (
                read_lightning_traj(fname))

            vals_rad_aux = np.ma.masked_all(nsources_file)
            nsources_missing = 0
            nflashes_missing = 0
            for flashnr_first in flashnrs_first:
                vals_rad_flash = vals_rad[flashnrs_rad == flashnr_first]
                nsources_flash_rad = vals_rad_flash.size
                nsources_flash = np.size(flashnrs[flashnrs == flashnr_first])
                if nsources_flash_rad == 0:
                    # no radar data for flash
                    print('No radar data for flash nr '+str(flashnr_first))
                    nsources_missing += nsources_flash
                    nflashes_missing += 1
                    continue

                if nsources_flash == nsources_flash_rad:
                    # radar data contained all flash sources
                    vals_rad_aux[flashnrs == flashnr_first] = vals_rad_flash
                else:
                    # radar data contain some of the flash
                    # print(
                    #    'Partial radar coverage of flash nr ' +
                    #    str(flashnr_first))
                    times_flash = times[flashnrs == flashnr_first]
                    times_rad_flash = times_rad[flashnrs_rad == flashnr_first]
                    dBms_rad_flash = dBms_rad[flashnrs_rad == flashnr_first]
                    for time_flash in times_flash:
                        found = False
                        for time_rad_flash in  times_rad_flash:
                            delta_time = np.abs(
                                (time_flash-time_rad_flash).total_seconds())
                            if delta_time < 0.00006:
                                vals_rad_src = vals_rad_flash[
                                    times_rad_flash == time_rad_flash]
                                dBm_rad_src = dBms_rad_flash[
                                    times_rad_flash == time_rad_flash]
                                dBm_src = dBms[times == time_flash]
                                if (np.any(np.round(dBm_rad_src*10.)
                                           == np.round(dBm_src*10))):
                                    vals_rad_aux[times == time_flash] = (
                                        vals_rad_src[
                                            np.round(dBm_rad_src*10.) ==
                                            np.round(dBm_src*10)])
                                    found = True
                                    break
                        if not found:
                            nsources_missing += 1
                            # print(
                            #    'Missing radar data for flash source time ' +
                            #    str(time_in_single_flash[k]))

            vals_list.append(vals_rad_aux)
            print('Missing sources in radar file '+str(nsources_missing))
            print('Missing flashes in radar file '+str(nflashes_missing))

        fname = args.basepath+day_str+'_Santis_data.csv'
        fname = write_ts_lightning(
            flashnrs, times, times_in_flash, lats, lons, alts, dBms,
            vals_list, fname, pol_vals_labels)


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
