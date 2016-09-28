"""
pyrad.io.write_data
====================

Functions for writing pyrad output data

.. autosummary::
    :toctree: generated/

    write_timeseries
    write_sun_hits
    write_sun_retrieval
    generate_field_name_str

"""

import glob
import csv

from pyart.config import get_fillvalue, get_metadata

from ..io.read_data import get_fieldname_rainbow


def write_timeseries(dataset, fname):
    """
    writes time series of data

    Parameters
    ----------
    dataset : dict
        dictionary containing the time series parameters

    fname : str
        file name where to store the data

    Returns
    -------
    No return

    """
    filelist = glob.glob(fname)
    if len(filelist) == 0:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar timeseries data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Description: \n')
            csvfile.write('# Time series of a weather radar data over a ' +
                          'fixed location.\n')
            csvfile.write(
                '# Location [lon, lat, alt]: ' +
                str(dataset['point_coordinates_WGS84_lon_lat_alt']) + '\n')
            csvfile.write(
                '# Nominal antenna coordinates used [az, el, r]: ' +
                str(dataset['antenna_coordinates_az_el_r'])+'\n')
            csvfile.write(
                '# Data: '+generate_field_name_str(dataset['datatype'])+'\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write(
                '# Start: ' +
                dataset['time'].strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
            csvfile.write('#\n')

            fieldnames = ['date', 'az', 'el', 'r', 'value']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            writer.writerow(
                {'date': dataset['time'],
                 'az': dataset['used_antenna_coordinates_az_el_r'][0],
                 'el': dataset['used_antenna_coordinates_az_el_r'][1],
                 'r': dataset['used_antenna_coordinates_az_el_r'][2],
                 'value': dataset['value']})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = ['date', 'az', 'el', 'r', 'value']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow(
                {'date': dataset['time'],
                 'az': dataset['used_antenna_coordinates_az_el_r'][0],
                 'el': dataset['used_antenna_coordinates_az_el_r'][1],
                 'r': dataset['used_antenna_coordinates_az_el_r'][2],
                 'value': dataset['value']})
            csvfile.close()


def write_sun_hits(sun_hits, fname):
    """
    Writes sun hits data.

    Parameters
    ----------
    sun_hits : dict
        dictionary containing the sun hits parameters

    fname : str
        file name where to store the data

    Returns
    -------
    No return

    """
    filelist = glob.glob(fname)
    nhits = len(sun_hits['time'])
    if len(filelist) == 0:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar sun hits data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write('#\n')

            fieldnames = ['time', 'ray', 'rad_el', 'rad_az', 'sun_el',
                          'sun_az',
                          'Ph', 'std(Ph)', 'NPh',
                          'Pv', 'std(Pv)', 'NPv',
                          'ZDR', 'std(ZDR)', 'NPzdr',
                          'NPval', 'NPrng']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i in range(nhits):
                writer.writerow(
                    {'time': sun_hits['time'][i],
                     'ray': sun_hits['ray'][i],
                     'rad_el': sun_hits['rad_el'][i],
                     'rad_az': sun_hits['rad_az'][i],
                     'sun_el': sun_hits['sun_el'][i],
                     'sun_az': sun_hits['sun_az'][i],
                     'Ph': sun_hits['pwrh'][i],
                     'std(Ph)': sun_hits['pwrh_std'][i],
                     'NPh': sun_hits['npointsh'][i],
                     'Pv': sun_hits['pwrv'][i],
                     'std(Pv)': sun_hits['pwrv_std'][i],
                     'NPv': sun_hits['npointsv'][i],
                     'ZDR': sun_hits['zdr'][i],
                     'std(ZDR)': sun_hits['zdr_std'][i],
                     'NPzdr': sun_hits['npointszdr'][i],
                     'NPval': sun_hits['nvalid'][i],
                     'NPrng': sun_hits['nrange'][i]})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = ['time', 'ray', 'rad_el', 'rad_az',
                          'sun_el', 'sun_az',
                          'Ph', 'std(Ph)', 'NPh',
                          'Pv', 'std(Pv)', 'NPv',
                          'ZDR', 'std(ZDR)', 'NPzdr',
                          'NPval', 'NPrng']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i in range(nhits):
                writer.writerow(
                    {'time': sun_hits['time'][i],
                     'ray': sun_hits['ray'][i],
                     'rad_el': sun_hits['rad_el'][i],
                     'rad_az': sun_hits['rad_az'][i],
                     'sun_el': sun_hits['sun_el'][i],
                     'sun_az': sun_hits['sun_az'][i],
                     'Ph': sun_hits['pwrh'][i],
                     'std(Ph)': sun_hits['pwrh_std'][i],
                     'NPh': sun_hits['npointsh'][i],
                     'Pv': sun_hits['pwrv'][i],
                     'std(Pv)': sun_hits['pwrv_std'][i],
                     'NPv': sun_hits['npointsv'][i],
                     'ZDR': sun_hits['zdr'][i],
                     'std(ZDR)': sun_hits['zdr_std'][i],
                     'NPzdr': sun_hits['npointszdr'][i],
                     'NPval': sun_hits['nvalid'][i],
                     'NPrng': sun_hits['nrange'][i]})
            csvfile.close()


def write_sun_retrieval(sun_retrieval, fname):
    """
    Writes sun retrieval data.

    Parameters
    ----------
    sun_retrieval : dict
        dictionary containing the sun retrieval parameters

    fname : str
        file name where to store the data

    Returns
    -------
    No return

    """
    filelist = glob.glob(fname)
    if len(filelist) == 0:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar sun retrievals data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write('#\n')

            fieldnames = [
                'time',
                'nhits_h', 'el_width_h', 'az_width_h',
                'el_bias_h', 'az_bias_h', 'Ph', 'std(Ph)',
                'nhits_v', 'el_width_v', 'az_width_v',
                'el_bias_v', 'az_bias_v', 'Pv', 'std(Pv)',
                'nhits_zdr', 'ZDR', 'std(ZDR)']

            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            writer.writerow(
                {'time': sun_retrieval['time'].strftime('%Y%m%d'),
                 'nhits_h': sun_retrieval['nhits_h'],
                 'el_width_h': sun_retrieval['el_width_h'],
                 'az_width_h': sun_retrieval['az_width_h'],
                 'el_bias_h': sun_retrieval['el_bias_h'],
                 'az_bias_h': sun_retrieval['az_bias_h'],
                 'Ph': sun_retrieval['Ph'],
                 'std(Ph)': sun_retrieval['std(Ph)'],
                 'nhits_v': sun_retrieval['nhits_v'],
                 'el_width_v': sun_retrieval['el_width_v'],
                 'az_width_v': sun_retrieval['az_width_v'],
                 'el_bias_v': sun_retrieval['el_bias_v'],
                 'az_bias_v': sun_retrieval['az_bias_v'],
                 'Pv': sun_retrieval['Pv'],
                 'std(Pv)': sun_retrieval['std(Pv)'],
                 'nhits_zdr': sun_retrieval['nhits_zdr'],
                 'ZDR': sun_retrieval['ZDR'],
                 'std(ZDR)': sun_retrieval['std(ZDR)']})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'time',
                'nhits_h', 'el_width_h', 'az_width_h',
                'el_bias_h', 'az_bias_h', 'Ph', 'std(Ph)',
                'nhits_v', 'el_width_v', 'az_width_v',
                'el_bias_v', 'az_bias_v', 'Pv', 'std(Pv)',
                'nhits_zdr', 'ZDR', 'std(ZDR)']

            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow(
                {'time': sun_retrieval['time'].strftime('%Y%m%d'),
                 'nhits_h': sun_retrieval['nhits_h'],
                 'el_width_h': sun_retrieval['el_width_h'],
                 'az_width_h': sun_retrieval['az_width_h'],
                 'el_bias_h': sun_retrieval['el_bias_h'],
                 'az_bias_h': sun_retrieval['az_bias_h'],
                 'Ph': sun_retrieval['Ph'],
                 'std(Ph)': sun_retrieval['std(Ph)'],
                 'nhits_v': sun_retrieval['nhits_v'],
                 'el_width_v': sun_retrieval['el_width_v'],
                 'az_width_v': sun_retrieval['az_width_v'],
                 'el_bias_v': sun_retrieval['el_bias_v'],
                 'az_bias_v': sun_retrieval['az_bias_v'],
                 'Pv': sun_retrieval['Pv'],
                 'std(Pv)': sun_retrieval['std(Pv)'],
                 'nhits_zdr': sun_retrieval['nhits_zdr'],
                 'ZDR': sun_retrieval['ZDR'],
                 'std(ZDR)': sun_retrieval['std(ZDR)']})
            csvfile.close()


def generate_field_name_str(datatype):
    """
    Generates a field name in a nice to read format.

    Parameters
    ----------
    datatype : str
        The data type

    Returns
    -------
    field_str : str
        The field name

    """
    field_name = get_fieldname_rainbow(datatype)
    field_dic = get_metadata(field_name)
    field_str = field_dic['standard_name'].replace('_', ' ')
    field_str = field_str[0].upper() + field_str[1:]
    field_str += ' ('+field_dic['units']+')'

    return field_str
