#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_windmill_filtering
================================================

This program finds periods in time where the windmills nacelle had a certain
orientation with respect to the radar (e.g. 0°+-5°) and the windmill rotor had
a certain velocity (typically either 0 or >0)

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
from warnings import warn

import numpy as np

from pyrad.io import read_windmills_data, write_proc_periods

print(__doc__)


def main():
    """
    """

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    parser.add_argument(
        'startdate', type=str,
        help='starting date of the data to be processed. Format ''YYYYMMDDhhmmss'' ')
    parser.add_argument(
        'enddate', type=str,
        help='end date of the data to be processed. Format ''YYYYMMDDhhmmss'' ')

    # keyword arguments
    parser.add_argument(
        '--basename', type=str,
        default='/users/jfigui/windmills_params/',
        help='name of folder containing the windmill data')

    parser.add_argument(
        '--wind_IDs', type=str, default='nx85213,nx85214,nx85215',
        help='Windmill ID. Coma separated')
        
    parser.add_argument(
        '--wind_periods', type=str, default='20200227-20200326',
        help='Windmill periods. Coma separated')

    parser.add_argument(
        '--orientations', type=str, default='0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350',
        help='Orientation respect to radar')

    parser.add_argument(
        '--span', type=float, default=10.,
        help='Span')

    parser.add_argument(
        '--vel_limit', type=float, default=0.,
        help='Velocity limit')

    args = parser.parse_args()

    print("====== PYRAD windmill data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD windmill data processing finished: ")


    start_date = datetime.datetime.strptime(args.startdate, '%Y%m%d%H%M%S')
    end_date = datetime.datetime.strptime(args.enddate, '%Y%m%d%H%M%S')
    wind_IDs = args.wind_IDs.split(',')
    wind_periods = args.wind_periods.split(',')
    orientations = np.asarray(args.orientations.split(','), dtype=float)

    for wind_ID in wind_IDs:
        print(wind_ID)

        fname = args.basename+wind_ID+'_rawdata_10m_'+wind_periods[0]+'.csv'
        windmill_dict = read_windmills_data(fname)

        if len(wind_periods) > 1:
            for wind_period in wind_periods[1:]:
                fname = (
                    args.basename+wind_ID+'_rawdata_10m_'+wind_period+'.csv')

                windmill_dict_aux = read_windmills_data(fname)
                windmill_dict['dt_remote'] = np.ma.append(
                    windmill_dict['dt_remote'],
                    windmill_dict_aux['dt_remote'])
                windmill_dict['dt_server'] = np.ma.append(
                    windmill_dict['dt_server'],
                    windmill_dict_aux['dt_server'])
                windmill_dict['rotor_speed_avg'] = np.ma.append(
                    windmill_dict['rotor_speed_avg'],
                    windmill_dict_aux['rotor_speed_avg'])
                windmill_dict['rotor_speed_min'] = np.ma.append(
                    windmill_dict['rotor_speed_min'],
                    windmill_dict_aux['rotor_speed_min'])
                windmill_dict['rotor_speed_max'] = np.ma.append(
                    windmill_dict['rotor_speed_max'],
                    windmill_dict_aux['rotor_speed_max'])
                windmill_dict['nacelle_pos'] = np.ma.append(
                    windmill_dict['nacelle_pos'],
                    windmill_dict_aux['nacelle_pos'])
                windmill_dict['blade_angle_1'] = np.ma.append(
                    windmill_dict['blade_angle_1'],
                    windmill_dict_aux['blade_angle_1'])
                windmill_dict['blade_angle_2'] = np.ma.append(
                    windmill_dict['blade_angle_2'],
                    windmill_dict_aux['blade_angle_2'])
                windmill_dict['blade_angle_3'] = np.ma.append(
                    windmill_dict['blade_angle_3'],
                    windmill_dict_aux['blade_angle_3'])
                windmill_dict['t_outside'] = np.ma.append(
                    windmill_dict['t_outside'],
                    windmill_dict_aux['t_outside'])
                windmill_dict['ice_level'] = np.ma.append(
                    windmill_dict['ice_level'],
                    windmill_dict_aux['ice_level'])

        ind = np.ma.where(np.logical_and(
            windmill_dict['dt_remote'] >= start_date,
            windmill_dict['dt_remote'] <= end_date))

        windmill_dict['dt_remote'] = windmill_dict['dt_remote'][ind]
        windmill_dict['dt_server'] = windmill_dict['dt_server'][ind]
        windmill_dict['rotor_speed_avg'] = (
            windmill_dict['rotor_speed_avg'][ind])
        windmill_dict['rotor_speed_min'] = (
            windmill_dict['rotor_speed_min'][ind])
        windmill_dict['rotor_speed_max'] = (
            windmill_dict['rotor_speed_max'][ind])
        windmill_dict['nacelle_pos'] = windmill_dict['nacelle_pos'][ind]
        windmill_dict['blade_angle_1'] = windmill_dict['blade_angle_1'][ind]
        windmill_dict['blade_angle_2'] = windmill_dict['blade_angle_2'][ind]
        windmill_dict['blade_angle_3'] = windmill_dict['blade_angle_3'][ind]
        windmill_dict['t_outside'] = windmill_dict['t_outside'][ind]
        windmill_dict['ice_level'] = windmill_dict['ice_level'][ind]

        # change angle to orientation respect to radar
        if wind_ID == 'nx85215':
            azi = 337.84
        elif wind_ID == 'nx85213':
            azi = 342.27
        elif wind_ID == 'nx85214':
            azi = 340.24

        nacelle_radar_ori = np.mod(
            windmill_dict['nacelle_pos']+(360.-azi+180.), 360.)

        print('total records: '+str(np.size(nacelle_radar_ori)))

        # filter according to orientation respect to radar
        for orientation in orientations:
            ind = get_indices_orientation(
                nacelle_radar_ori, orient=orientation, span=args.span)

            if ind is None:
                continue

            print('records with orientation '+str(orientation)+' deg +- ' +
                  str(args.span/2.)+': '+str(np.size(ind)))

            dt_remote = windmill_dict['dt_remote'][ind]
            rotor_speed_avg = windmill_dict['rotor_speed_avg'][ind]
            # nacelle_pos = nacelle_radar_ori[ind]
            # blade_angle_1 = windmill_dict['blade_angle_1'][ind]
            # blade_angle_2 = windmill_dict['blade_angle_2'][ind]
            # blade_angle_3 = windmill_dict['blade_angle_3'][ind]


            # Filter according to rotor speed
            ind_fast = np.where(rotor_speed_avg > args.vel_limit)[0]
            ind_slow = np.where(rotor_speed_avg <= args.vel_limit)[0]

            print('records with vel <= '+str(args.vel_limit)+': '+
                  str(np.size(ind_slow)))

            if ind_fast.size > 0:
                dt_remote_fast = dt_remote[ind_fast]
                # rotor_speed_avg_fast = rotor_speed_avg[ind_fast]
                # nacelle_pos_fast = nacelle_pos[ind_fast]
                # blade_angle_1_fast = blade_angle_1[ind_fast]
                # blade_angle_2_fast = blade_angle_2[ind_fast]
                # blade_angle_3_fast = blade_angle_3[ind_fast]

                start_times_fast, end_times_fast = find_contiguous_times(
                    dt_remote_fast)

                fname = (
                    args.basename+wind_ID+'_span'+str(args.span)+'_ori' +
                    str(orientation)+'_speedGT'+str(args.vel_limit)+'.csv')
                write_proc_periods(start_times_fast, end_times_fast, fname)
            else:
                warn('No data for rotor speed above '+str(args.vel_limit) +
                     ' rpm')

            if ind_slow.size > 0:
                dt_remote_slow = dt_remote[ind_slow]
                # rotor_speed_avg_slow = rotor_speed_avg[ind_slow]
                # nacelle_pos_slow = nacelle_pos[ind_slow]
                # blade_angle_1_slow = blade_angle_1[ind_slow]
                # blade_angle_2_slow = blade_angle_2[ind_slow]
                # blade_angle_3_slow = blade_angle_3[ind_slow]

                start_times_slow, end_times_slow = find_contiguous_times(
                    dt_remote_slow)

                fname = (
                    args.basename+wind_ID+'_span'+str(args.span)+'_ori' +
                    str(orientation)+'_speedLE'+str(args.vel_limit)+'.csv')
                write_proc_periods(start_times_slow, end_times_slow, fname)
            else:
                warn('No data for rotor speed below '+str(args.vel_limit) +
                     ' rpm')


def get_indices_orientation(orientations, orient=0., span=45.):
    """
    Get indices of data with the orientation within
    ]orient-span/2, orient+span/2]

    Parameters
    ----------
    orientations : array of floats
        The orientations array
    orient, span : float
        The orientation and the span

    Returns
    -------
    ind : array of ints or None
        The indices of data with the prescribed orientation. None
        otherwise

    """
    left_limit = orient-span/2.
    right_limit = orient+span/2.

    if left_limit < 0. or right_limit > 360.:
        if left_limit < 0.:
            left_limit = 360.-left_limit
        if right_limit > 360.:
            right_limit -= 360.
        ind = np.where(np.logical_or(
            orientations > left_limit, orientations <= right_limit))[0]
    else:
        ind = np.where(np.logical_and(
            orientations > left_limit, orientations <= right_limit))[0]

    if ind.size == 0:
        warn('No data for nacelle orientation '+str(orient)+' deg respect to radar')
        return None

    return ind


def find_contiguous_times(times, step=600):
    """
    Given and array of ordered times, find those contiguous according to
    a maximum time step

    Parameters
    ----------
    times : array of datetimes
        The array of times
    step : float
        The time step [s]

    Returns
    -------
    start_times, end_times : array of date times
        The start and end of each consecutive time period

    """
    run = []
    periods = []
    expect = None
    for time in times:
        if time == expect or expect is None:
            run.append(time)
        else:
            run = [time]
            periods.append(run)
        expect = time+datetime.timedelta(seconds=step)

    print('number of consecutive periods: '+str(len(periods)))

    start_times = np.array([], dtype=datetime.datetime)
    end_times = np.array([], dtype=datetime.datetime)
    for period in periods:
        start_times = np.append(
            start_times, period[0]-datetime.timedelta(seconds=step))
        end_times = np.append(end_times, period[-1])

    return start_times, end_times


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
