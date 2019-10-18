#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_process_data_windmills
================================================

This program compiles histograms of windmill radar returns

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import atexit
import os
import glob
from warnings import warn

import numpy as np

from pyrad.io import read_windmills_data, write_histogram
from pyrad.io import read_histogram, write_proc_periods
from pyrad.graph import plot_histogram
from pyrad.util import compute_histogram, find_contiguous_times

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
        '--loadbasepath', type=str,
        default='/store/msrad/radar/pyrad_products/',
        help='base path to the radar data')

    parser.add_argument(
        '--loadname', type=str,
        default='mals_sha_windmills_hr',
        help='Folder containing the radar data')

    parser.add_argument(
        '--wind_IDs', type=str, default='nx85213,nx85214,nx85215',
        help='Windmill ID. Coma separated')

    parser.add_argument(
        '--orientations', type=str, default='0,45,90,135,180,225,270,315',
        help='Orientation respect to radar')

    parser.add_argument(
        '--span', type=float, default=45.,
        help='Span')

    parser.add_argument(
        '--vel_limit', type=float, default=0.,
        help='Velocity limit')

    parser.add_argument(
        '--dBuZ_threshold', type=float, default=None,
        help='dBuZ threshold')

    args = parser.parse_args()

    print("====== PYRAD windmill data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD windmill data processing finished: ")


    start_date = datetime.datetime.strptime(args.startdate, '%Y%m%d%H%M%S')
    end_date = datetime.datetime.strptime(args.enddate, '%Y%m%d%H%M%S')
    wind_IDs = args.wind_IDs.split(',')
    orientations = np.asarray(args.orientations.split(','), dtype=float)

    cfg = {
        'loadbasepath': [args.loadbasepath],
        'loadname': [args.loadname]}

    if 'rhi' in args.loadname:
        scan_type = 'rhi'
    else:
        scan_type = 'ppi'

    speeds = ['speed_GT'+str(args.vel_limit), 'speed_LE'+str(args.vel_limit)]
    if args.dBuZ_threshold is not None:
        wet_dry = [
            'dBuZ_LT'+str(args.dBuZ_threshold),
            'dBuZ_GE'+str(args.dBuZ_threshold)]

    for wind_ID in wind_IDs:
        print('\n'+wind_ID)

        fname = args.basename+wind_ID+'_rawdata_10m20190227_20190326.csv'
        fname2 = args.basename+wind_ID+'_rawdata_10m20190327_20190426.csv'

        windmill_dict = read_windmills_data(fname)
        windmill_dict_aux = read_windmills_data(fname2)
        windmill_dict['dt_remote'] = np.ma.append(
            windmill_dict['dt_remote'], windmill_dict_aux['dt_remote'])
        windmill_dict['dt_server'] = np.ma.append(
            windmill_dict['dt_server'], windmill_dict_aux['dt_server'])
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
            windmill_dict['nacelle_pos'], windmill_dict_aux['nacelle_pos'])
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
            windmill_dict['t_outside'], windmill_dict_aux['t_outside'])
        windmill_dict['ice_level'] = np.ma.append(
            windmill_dict['ice_level'], windmill_dict_aux['ice_level'])

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
            wind_ID2 = 'WM1'
        elif wind_ID == 'nx85213':
            azi = 342.27
            wind_ID2 = 'WM2'
        elif wind_ID == 'nx85214':
            azi = 340.24
            wind_ID2 = 'WM3'

        nacelle_radar_ori = np.mod(
            windmill_dict['nacelle_pos']+(360.-azi+180.), 360.)

        print('total windmill records: '+str(np.size(nacelle_radar_ori)))

        # Compute median dBuZ if stratified according to wet/dry
        if args.dBuZ_threshold is not None:
            dt_radar, ori_radar, speed_radar, val_radar = compute_median(
                windmill_dict['dt_remote'], nacelle_radar_ori,
                windmill_dict['rotor_speed_avg'],
                cfg['loadbasepath'][0]+cfg['loadname'][0]+'/', wind_ID2)
        else:
            dt_radar = windmill_dict['dt_remote']
            ori_radar = nacelle_radar_ori
            speed_radar = windmill_dict['rotor_speed_avg']
            val_radar = None

        # Compute histograms and time periods
        for speed in speeds:
            dt_radar_aux, ori_aux, val_radar_aux = filter_speed(
                dt_radar, ori_radar, speed_radar, speed, val_radar=val_radar,
                vel_limit=args.vel_limit)

            if args.dBuZ_threshold is not None:
                for wd in wet_dry:
                    dt_radar_aux2, ori_aux2 = filter_wet_dry(
                        dt_radar_aux, ori_aux, val_radar_aux, wd,
                        dBuZ_threshold=args.dBuZ_threshold)

                    plot_and_write_histogram(
                        ori_aux2, args.basename, wind_ID, args.startdate,
                        args.enddate, scan_type, speed, wd=wd)

                    # filter according to orientation respect to radar
                    filter_orientation(
                        dt_radar_aux2, ori_aux2, orientations, args.basename,
                        wind_ID, scan_type, args.span, speed, wd=wd)
            else:
                plot_and_write_histogram(
                    ori_aux, args.basename, wind_ID, args.startdate,
                    args.enddate, scan_type, speed, wd=None)

                # filter according to orientation respect to radar
                filter_orientation(
                    dt_radar_aux, ori_aux, orientations, args.basename,
                    wind_ID, scan_type, args.span, speed, wd=None)


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


def get_file_list(starttimes, endtimes, radarbase, windmill, datatype):
    filelist = []
    for starttime, endtime in zip(starttimes, endtimes):
        startdate = starttime.replace(
            hour=0, minute=0, second=0, microsecond=0)
        enddate = endtime.replace(hour=0, minute=0, second=0, microsecond=0)
        ndays = int((enddate-startdate).days)+1
        t_filelist = []
        for i in range(ndays):
            daydir = (
                startdate+datetime.timedelta(days=i)).strftime('%Y-%m-%d')
            datapath = (
                radarbase+daydir+'/'+windmill+'_postproc/HISTOGRAM_' +
                datatype+'/')
            if not os.path.isdir(datapath):
                continue
            dayfilelist = glob.glob(datapath+'*_histogram_*'+datatype+'.csv')
            for filename in dayfilelist:
                t_filelist.append(filename)

        for filename in t_filelist:
            filenamestr = str(filename)
            datetimestr = os.path.basename(filenamestr)[0:14]
            fdatetime = datetime.datetime.strptime(
                datetimestr, '%Y%m%d%H%M%S')
            if (fdatetime >= starttime) and (fdatetime <= endtime):
                filelist.append(filenamestr)

    return sorted(filelist)


def compute_median(dt_remotes, origins, rotor_speeds, basepath, wind_ID):
    dt_radar = np.ma.array([], dtype=datetime.datetime)
    ori_radar = np.ma.array([])
    speed_radar = np.ma.array([])
    val_radar = np.ma.array([])
    for dt_remote, ori, rotor_speed in zip(dt_remotes, origins, rotor_speeds):
        file_list = get_file_list(
            [dt_remote-datetime.timedelta(minutes=10)], [dt_remote],
            basepath, wind_ID, 'dBuZ')

        if file_list is not None and file_list:
            dt_radar = np.ma.append(dt_radar, dt_remote)
            ori_radar = np.ma.append(ori_radar, ori)
            speed_radar = np.ma.append(speed_radar, rotor_speed)

            # Get median
            hist_aux, bin_edges_aux = read_histogram(file_list[0])
            bin_centers = (
                bin_edges_aux[:-1]+(bin_edges_aux[1:]-bin_edges_aux[:-1])/2.)
            cum_sum = np.cumsum(hist_aux)
            median_val = bin_centers[cum_sum >= cum_sum[-1]/2.][0]
            val_radar = np.ma.append(val_radar, median_val)

    print('total radar records: '+str(np.size(dt_radar)))

    return dt_radar, ori_radar, speed_radar, val_radar


def filter_speed(dt_remotes, origins, rotor_speeds, speed, val_radar=None,
                 vel_limit=0.):
    val_radar_filt = None
    if 'speed_GT' in speed:
        dt_filt = dt_remotes[rotor_speeds > vel_limit]
        ori_filt = origins[rotor_speeds > vel_limit]
        if val_radar is not None:
            val_radar_filt = val_radar[rotor_speeds > vel_limit]
    else:
        dt_filt = dt_remotes[rotor_speeds <= vel_limit]
        ori_filt = origins[rotor_speeds <= vel_limit]
        if val_radar is not None:
            val_radar_filt = val_radar[rotor_speeds <= vel_limit]

    return dt_filt, ori_filt, val_radar_filt


def filter_wet_dry(dt_remotes, origins, val_radar, wd, dBuZ_threshold=10.):
    if 'dBuZ_LT' in wd:
        dt_filt = dt_remotes[val_radar < dBuZ_threshold]
        ori_filt = origins[val_radar < dBuZ_threshold]
    else:
        dt_filt = dt_remotes[val_radar >= dBuZ_threshold]
        ori_filt = origins[val_radar >= dBuZ_threshold]

    return dt_filt, ori_filt


def plot_and_write_histogram(values, basename, wind_ID, startdate, enddate,
                             scan_type, speed, wd=None):
    bin_edges = np.arange(-0.5, 361.5, 1.)
    if wd is None:
        wd_str = ''
    else:
        wd_str = '_'+wd

    fname = (
        basename+wind_ID+'_'+startdate+'-' +enddate+'_'+scan_type +
        '_nacelle_ori_rotor_'+speed+wd_str)
    titl = (
        wind_ID+' '+startdate+'-'+enddate+'\n'+scan_type +
        ' Nacelle orientation when rotor '+speed.replace('_', ' ') +
        ' rpm\n '+wd_str)

    _, vals = compute_histogram(values, None, bin_edges=bin_edges)
    fig_fname = plot_histogram(
        bin_edges, vals, [fname+'.png'],
        labelx='Orientation respect to radar [deg]', titl=titl)
    print('Plotted '+' '.join(fig_fname))

    hist, _ = np.histogram(vals, bins=bin_edges)
    csv_fname = write_histogram(bin_edges, hist, fname+'.csv')
    print('Written '+csv_fname)


def filter_orientation(wm_dt, wm_oris, orientations, basename, wind_ID,
                       scan_type, span, speed, wd=None):
    if wd is None:
        wd_str = ''
    else:
        wd_str = '_'+wd
    for ori in orientations:
        ind = get_indices_orientation(wm_oris, orient=ori, span=span)

        if ind is None:
            print(
                'records with orientation '+str(ori) +
                ' deg +- '+str(span/2.)+', '+speed+', '+wd_str+': 0')
            continue

        print(
            'records with orientation '+str(ori)+' deg +- '+str(span/2.) +
            ', '+speed+', '+wd_str+': '+str(np.size(ind)))

        start_times, end_times = find_contiguous_times(wm_dt[ind])

        fname = (
            basename+scan_type+'_'+wind_ID+'_span'+str(span)+'_ori'+str(ori) +
            '_'+speed+wd_str+'.csv')
        write_proc_periods(start_times, end_times, fname)


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
