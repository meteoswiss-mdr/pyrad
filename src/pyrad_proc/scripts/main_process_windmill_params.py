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

import numpy as np

from pyrad.io import read_windmills_data, write_histogram
from pyrad.graph import plot_timeseries, plot_histogram, _plot_time_range
from pyrad.util import compute_histogram

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

    args = parser.parse_args()

    print("====== PYRAD windmill data processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== PYRAD windmill data processing finished: ")


    start_date = datetime.datetime.strptime(args.startdate, '%Y%m%d%H%M%S')
    end_date = datetime.datetime.strptime(args.enddate, '%Y%m%d%H%M%S')
    wind_IDs = args.wind_IDs.split(',')

    for wind_ID in wind_IDs:
        print(wind_ID)
        
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

        # Time series plots
        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_rotor_speed.png')
        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate+'\nRotor Speed')

        plot_timeseries(
            windmill_dict['dt_remote'],
            [windmill_dict['rotor_speed_avg'],
             windmill_dict['rotor_speed_min'],
             windmill_dict['rotor_speed_max']], [fname],
            labely='Speed [rpm]', labels=['avg', 'min', 'max'],
            title=titl, colors=['b', 'r', 'r'])

        print("----- plot to '%s'" % fname)

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_nacelle_pos.png')

        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate +
            '\nNacelle position')
        plot_timeseries(
            windmill_dict['dt_remote'], [windmill_dict['nacelle_pos']],
            [fname], labely='Position [deg from North]', labels=['pos'],
            title=titl)

        print("----- plot to '%s'" % fname)

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_blade_ang.png')

        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate +
            '\nBlade angle')
        plot_timeseries(
            windmill_dict['dt_remote'],
            [windmill_dict['blade_angle_1'], windmill_dict['blade_angle_2'],
             windmill_dict['blade_angle_3']], [fname],
            labely='Angle [deg]',
            labels=['angle 1', 'angle 2', 'angle 3'], title=titl)

        print("----- plot to '%s'" % fname)

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_t_outside.png')

        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate +
            '\nOutside temperature')
        plot_timeseries(
            windmill_dict['dt_remote'],
            [windmill_dict['t_outside']], [fname],
            labely='Temperature [deg Celsius]', labels=['T'], title=titl)

        print("----- plot to '%s'" % fname)

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_ice_level.png')

        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate +
            '\nIce level')
        plot_timeseries(
            windmill_dict['dt_remote'],
            [windmill_dict['ice_level']], [fname],
            labely='Ice level', labels=['Ice level'], title=titl)

        print("----- plot to '%s'" % fname)


        # Histogram plots
        bin_edges = np.arange(-0.5, 13.5, 1.)
        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_rotor_speed_avg_hist.png')
        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate+'\nRotor Speed avg')

        _, vals = compute_histogram(
            windmill_dict['rotor_speed_avg'], None, bin_edges=bin_edges)
        fname = plot_histogram(
            bin_edges, vals, [fname], labelx='Speed [rpm]', titl=titl)
        print('Plotted '+' '.join(fname))

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_rotor_speed_avg_hist.csv')
        hist, _ = np.histogram(vals, bins=bin_edges)
        fname = write_histogram(bin_edges, hist, fname)
        print('Written '+fname)


        bin_edges = np.arange(-0.5, 361.5, 1.)
        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_nacelle_pos_hist.png')
        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate+'\nNacelle position')

        _, vals = compute_histogram(
            windmill_dict['nacelle_pos'], None, bin_edges=bin_edges)
        fname = plot_histogram(
            bin_edges, vals, [fname], labelx='Position [deg from North]',
            titl=titl)
        print('Plotted '+' '.join(fname))

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_nacelle_pos_hist.csv')
        hist, _ = np.histogram(vals, bins=bin_edges)
        fname = write_histogram(bin_edges, hist, fname)
        print('Written '+fname)


        bin_edges = np.arange(-0.5, 91.5, 1.)
        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_blade_ang1_hist.png')
        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate+'\nBlade angle 1')

        _, vals = compute_histogram(
            windmill_dict['blade_angle_1'], None, bin_edges=bin_edges)
        fname = plot_histogram(
            bin_edges, vals, [fname], labelx='Angle [deg]',
            titl=titl)
        print('Plotted '+' '.join(fname))

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_blade_ang1_hist.csv')
        hist, _ = np.histogram(vals, bins=bin_edges)
        fname = write_histogram(bin_edges, hist, fname)
        print('Written '+fname)


        bin_edges = np.arange(-0.5, 91.5, 1.)
        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_blade_ang2_hist.png')
        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate+'\nBlade angle 2')

        _, vals = compute_histogram(
            windmill_dict['blade_angle_2'], None, bin_edges=bin_edges)
        fname = plot_histogram(
            bin_edges, vals, [fname], labelx='Angle [deg]',
            titl=titl)
        print('Plotted '+' '.join(fname))

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_blade_ang2_hist.csv')
        hist, _ = np.histogram(vals, bins=bin_edges)
        fname = write_histogram(bin_edges, hist, fname)
        print('Written '+fname)


        bin_edges = np.arange(-0.5, 91.5, 1.)
        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_blade_ang3_hist.png')
        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate+'\nBlade angle 3')

        _, vals = compute_histogram(
            windmill_dict['blade_angle_3'], None, bin_edges=bin_edges)
        fname = plot_histogram(
            bin_edges, vals, [fname], labelx='Angle [deg]',
            titl=titl)
        print('Plotted '+' '.join(fname))

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_blade_ang3_hist.csv')
        hist, _ = np.histogram(vals, bins=bin_edges)
        fname = write_histogram(bin_edges, hist, fname)
        print('Written '+fname)


        # 2D histogram plots
        bin_edges_rotor = np.arange(-0.5, 13.5, 1.)
        _, vals_rotor = compute_histogram(
            windmill_dict['rotor_speed_avg'], None, bin_edges=bin_edges_rotor)

        bin_edges_nacelle = np.arange(-0.5, 361.5, 1.)
        _, vals_nacelle = compute_histogram(
            windmill_dict['nacelle_pos'], None, bin_edges=bin_edges_nacelle)

        H, _, _ = np.histogram2d(
            vals_nacelle, vals_rotor,
            bins=[bin_edges_nacelle, bin_edges_rotor])

        # set 0 values to blank
        H = np.ma.asarray(H)
        H[H == 0] = np.ma.masked

        fname = (
            args.basename+wind_ID+'_'+args.startdate+'-'+args.enddate +
            '_nacelle-rotor_hist.png')
        titl = (
            wind_ID+' '+args.startdate+'-'+args.enddate +
            '\nNacelle position-Rotor speed')
        fname = _plot_time_range(
            bin_edges_nacelle, bin_edges_rotor, H, None, [fname],
            titl=titl,
            xlabel='Nacelle position [deg from north]',
            ylabel='Rotor speed [rpm]',
            clabel='Occurrence',
            vmin=0, vmax=None, figsize=[10, 8], dpi=72)
        print('Plotted '+' '.join(fname))



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
