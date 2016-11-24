"""
pyrad.io.write_data
====================

Functions for writing pyrad output data

.. autosummary::
    :toctree: generated/

    write_timeseries
    write_ts_polar_data
    write_monitoring_ts
    write_sun_hits
    write_sun_retrieval
    generate_field_name_str

"""

from __future__ import print_function
import glob
import csv

from pyart.config import get_fillvalue, get_metadata

from .io_aux import generate_field_name_str


def write_timeseries(ts, fname):
    """
    """

    print("----- write to '%s'" % fname)

    try:
        tsfile = open(fname, "w")
    except:
        raise Exception("ERROR: Could not create file '%s'" % fname)

    print("# Weather radar timeseries data file", file=tsfile)
    print("# Project: MALSplus", file=tsfile)
    print("# Start : %s UTC" % ts.time_vector[0].strftime("%Y-%m-%d %H:%M:%S"),
          file=tsfile)
    print("# End   : %s UTC" % ts.time_vector[-1].strftime("%Y-%m-%d %H:%M:%S"),
          file=tsfile)
    print("# Header lines with comments are preceded by '#'", file=tsfile)
    for line in ts.description:
        print("# %s" % line, file=tsfile)
    print("#", file=tsfile)

    # Make raw header
    if (ts.timeformat is None):
        print("# Date, UTC [seconds since midnight]", end="", file=tsfile)
    else:
        print("# Date [%s]" % ts.timeformat, end="", file=tsfile)

    for ds in ts.dataseries:
        print(", %s [%s]" % (ds.label, ds.unit), end="", file=tsfile)
    print("", file=tsfile)

    # Store the data
    nsample = len(ts.time_vector)
    for kk in range(nsample):
        if (ts.timeformat is None):
            dt = ts.time_vector[kk]
            daystr = dt.strftime("%d-%b-%Y")
            secs = dt.hour*3600. + dt.minute*60. + dt.second + \
                dt.microsecond/1000000.
            print("%s, %14.4f" % (daystr, secs), end="", file=tsfile)
        else:
            print(ts.time_vector[kk].strftime(ts.timeformat), end="",
                  file=tsfile)

        for ds in ts.dataseries:
            print(", %14.4f" % (ds.data[kk]), end="", file=tsfile)

        print("", file=tsfile)

    tsfile.close()


def write_ts_polar_data(dataset, fname):
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
    fname : str
        the name of the file where data has written

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

    return fname


def write_monitoring_ts(start_time, np_t, values, quantiles, datatype, fname):
    """
    writes time series of data

    Parameters
    ----------
    start_time : datetime object
        the time of the monitoring
    np_t : int
        the total number of points
    values: float array
        the values at certain quantiles
    quantiles: float array
        the quantiles computed
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    values_aux = values.filled(fill_value=get_fillvalue())
    filelist = glob.glob(fname)
    if len(filelist) == 0:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar monitoring timeseries data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Description: \n')
            csvfile.write('# Time series of a monitoring of weather radar' +
                          ' data.\n')
            csvfile.write(
                '# Quantiles: '+str(quantiles[1])+', '+str(quantiles[0])+', ' +
                str(quantiles[2])+' percent.\n')
            csvfile.write(
                '# Data: '+generate_field_name_str(datatype)+'\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write(
                '# Start: ' +
                start_time.strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
            csvfile.write('#\n')

            fieldnames = ['date', 'NP', 'central_quantile', 'low_quantile',
                          'high_quantile']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            writer.writerow(
                {'date': start_time.strftime('%Y%m%d%H%M%S'),
                 'NP': np_t,
                 'central_quantile': values_aux[1],
                 'low_quantile': values_aux[0],
                 'high_quantile': values_aux[2]})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = ['date', 'NP', 'central_quantile', 'low_quantile',
                          'high_quantile']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow(
                {'date': start_time.strftime('%Y%m%d%H%M%S'),
                 'NP': np_t,
                 'central_quantile': values_aux[1],
                 'low_quantile': values_aux[0],
                 'high_quantile': values_aux[2]})
            csvfile.close()

    return fname


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
    fname : str
        the name of the file where data has written

    """
    dBm_sun_hit = sun_hits['dBm_sun_hit'].filled(fill_value=get_fillvalue())
    std_dBm_sun_hit = sun_hits['std(dBm_sun_hit)'].filled(
        fill_value=get_fillvalue())
    dBmv_sun_hit = sun_hits['dBmv_sun_hit'].filled(fill_value=get_fillvalue())
    std_dBmv_sun_hit = sun_hits['std(dBmv_sun_hit)'].filled(
        fill_value=get_fillvalue())
    zdr_sun_hit = sun_hits['ZDR_sun_hit'].filled(fill_value=get_fillvalue())
    std_zdr_sun_hit = sun_hits['std(ZDR_sun_hit)'].filled(
        fill_value=get_fillvalue())

    filelist = glob.glob(fname)
    nhits = len(sun_hits['time'])
    if len(filelist) == 0:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar sun hits data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write('#\n')

            fieldnames = [
                'time', 'ray', 'NPrng',
                'rad_el', 'rad_az', 'sun_el', 'sun_az',
                'dBm_sun_hit', 'std(dBm_sun_hit)', 'NPh', 'NPhval',
                'dBmv_sun_hit', 'std(dBmv_sun_hit)', 'NPv', 'NPvval',
                'ZDR_sun_hit', 'std(ZDR_sun_hit)', 'NPzdr', 'NPzdrval']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i in range(nhits):
                writer.writerow(
                    {'time': sun_hits['time'][i],
                     'ray': sun_hits['ray'][i],
                     'NPrng': sun_hits['NPrng'][i],
                     'rad_el': sun_hits['rad_el'][i],
                     'rad_az': sun_hits['rad_az'][i],
                     'sun_el': sun_hits['sun_el'][i],
                     'sun_az': sun_hits['sun_az'][i],
                     'dBm_sun_hit': dBm_sun_hit[i],
                     'std(dBm_sun_hit)': std_dBm_sun_hit[i],
                     'NPh': sun_hits['NPh'][i],
                     'NPhval': sun_hits['NPhval'][i],
                     'dBmv_sun_hit': dBmv_sun_hit[i],
                     'std(dBmv_sun_hit)': std_dBmv_sun_hit[i],
                     'NPv': sun_hits['NPv'][i],
                     'NPvval': sun_hits['NPvval'][i],
                     'ZDR_sun_hit': zdr_sun_hit[i],
                     'std(ZDR_sun_hit)': std_zdr_sun_hit[i],
                     'NPzdr': sun_hits['NPzdr'][i],
                     'NPzdrval': sun_hits['NPzdrval'][i]})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'time', 'ray', 'NPrng',
                'rad_el', 'rad_az', 'sun_el', 'sun_az',
                'dBm_sun_hit', 'std(dBm_sun_hit)', 'NPh', 'NPhval',
                'dBmv_sun_hit', 'std(dBmv_sun_hit)', 'NPv', 'NPvval',
                'ZDR_sun_hit', 'std(ZDR_sun_hit)', 'NPzdr', 'NPzdrval']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i in range(nhits):
                writer.writerow(
                    {'time': sun_hits['time'][i],
                     'ray': sun_hits['ray'][i],
                     'NPrng': sun_hits['NPrng'][i],
                     'rad_el': sun_hits['rad_el'][i],
                     'rad_az': sun_hits['rad_az'][i],
                     'sun_el': sun_hits['sun_el'][i],
                     'sun_az': sun_hits['sun_az'][i],
                     'dBm_sun_hit': dBm_sun_hit[i],
                     'std(dBm_sun_hit)': std_dBm_sun_hit[i],
                     'NPh': sun_hits['NPh'][i],
                     'NPhval': sun_hits['NPhval'][i],
                     'dBmv_sun_hit': dBmv_sun_hit[i],
                     'std(dBmv_sun_hit)': std_dBmv_sun_hit[i],
                     'NPv': sun_hits['NPv'][i],
                     'NPvval': sun_hits['NPvval'][i],
                     'ZDR_sun_hit': zdr_sun_hit[i],
                     'std(ZDR_sun_hit)': std_zdr_sun_hit[i],
                     'NPzdr': sun_hits['NPzdr'][i],
                     'NPzdrval': sun_hits['NPzdrval'][i]})
            csvfile.close()

    return fname


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
    fname : str
        the name of the file where data has written

    """
    first_hit_time = sun_retrieval['first_hit_time'].strftime('%Y%m%d%H%M%S')
    last_hit_time = sun_retrieval['last_hit_time'].strftime('%Y%m%d%H%M%S')
    el_width_h = sun_retrieval['el_width_h'].filled(fill_value=get_fillvalue())
    az_width_h = sun_retrieval['az_width_h'].filled(fill_value=get_fillvalue())
    el_bias_h = sun_retrieval['el_bias_h'].filled(fill_value=get_fillvalue())
    az_bias_h = sun_retrieval['az_bias_h'].filled(fill_value=get_fillvalue())
    dBm_sun_est = sun_retrieval['dBm_sun_est'].filled(
        fill_value=get_fillvalue())
    std_dBm_sun_est = sun_retrieval['std(dBm_sun_est)'].filled(
        fill_value=get_fillvalue())

    el_width_v = sun_retrieval['el_width_v'].filled(fill_value=get_fillvalue())
    az_width_v = sun_retrieval['az_width_v'].filled(fill_value=get_fillvalue())
    el_bias_v = sun_retrieval['el_bias_v'].filled(fill_value=get_fillvalue())
    az_bias_v = sun_retrieval['az_bias_v'].filled(fill_value=get_fillvalue())
    dBmv_sun_est = sun_retrieval['dBmv_sun_est'].filled(
        fill_value=get_fillvalue())
    std_dBmv_sun_est = sun_retrieval['std(dBmv_sun_est)'].filled(
        fill_value=get_fillvalue())

    zdr_sun_est = sun_retrieval['ZDR_sun_est'].filled(
        fill_value=get_fillvalue())
    std_zdr_sun_est = sun_retrieval['std(ZDR_sun_est)'].filled(
        fill_value=get_fillvalue())
    dBm_sun_ref = sun_retrieval['dBm_sun_ref'].filled(
        fill_value=get_fillvalue())
    ref_time = 'None'
    if sun_retrieval['ref_time'] is not None:
        ref_time = sun_retrieval['ref_time'].strftime('%Y%m%d%H%M%S')

    filelist = glob.glob(fname)
    if len(filelist) == 0:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar sun retrievals data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write('#\n')

            fieldnames = [
                'first_hit_time', 'last_hit_time',
                'nhits_h', 'el_width_h', 'az_width_h',
                'el_bias_h', 'az_bias_h', 'dBm_sun_est', 'std(dBm_sun_est)',
                'nhits_v', 'el_width_v', 'az_width_v',
                'el_bias_v', 'az_bias_v', 'dBmv_sun_est', 'std(dBmv_sun_est)',
                'nhits_zdr', 'ZDR_sun_est', 'std(ZDR_sun_est)', 'dBm_sun_ref',
                'ref_time']

            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            writer.writerow(
                {'first_hit_time': first_hit_time,
                 'last_hit_time': last_hit_time,
                 'nhits_h': sun_retrieval['nhits_h'],
                 'el_width_h': el_width_h,
                 'az_width_h': az_width_h,
                 'el_bias_h': el_bias_h,
                 'az_bias_h': az_bias_h,
                 'dBm_sun_est': dBm_sun_est,
                 'std(dBm_sun_est)': std_dBm_sun_est,
                 'nhits_v': sun_retrieval['nhits_v'],
                 'el_width_v': el_width_v,
                 'az_width_v': az_width_v,
                 'el_bias_v': el_bias_v,
                 'az_bias_v': az_bias_v,
                 'dBmv_sun_est': dBmv_sun_est,
                 'std(dBmv_sun_est)': std_dBmv_sun_est,
                 'nhits_zdr': sun_retrieval['nhits_zdr'],
                 'ZDR_sun_est': zdr_sun_est,
                 'std(ZDR_sun_est)': std_zdr_sun_est,
                 'dBm_sun_ref': dBm_sun_ref,
                 'ref_time': ref_time})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'first_hit_time', 'last_hit_time',
                'nhits_h', 'el_width_h', 'az_width_h',
                'el_bias_h', 'az_bias_h', 'dBm_sun_est', 'std(dBm_sun_est)',
                'nhits_v', 'el_width_v', 'az_width_v',
                'el_bias_v', 'az_bias_v', 'dBmv_sun_est', 'std(dBmv_sun_est)',
                'nhits_zdr', 'ZDR_sun_est', 'std(ZDR_sun_est)', 'dBm_sun_ref',
                'ref_time']

            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow(
                {'first_hit_time': first_hit_time,
                 'last_hit_time': last_hit_time,
                 'nhits_h': sun_retrieval['nhits_h'],
                 'el_width_h': el_width_h,
                 'az_width_h': az_width_h,
                 'el_bias_h': el_bias_h,
                 'az_bias_h': az_bias_h,
                 'dBm_sun_est': dBm_sun_est,
                 'std(dBm_sun_est)': std_dBm_sun_est,
                 'nhits_v': sun_retrieval['nhits_v'],
                 'el_width_v': el_width_v,
                 'az_width_v': az_width_v,
                 'el_bias_v': el_bias_v,
                 'az_bias_v': az_bias_v,
                 'dBmv_sun_est': dBmv_sun_est,
                 'std(dBmv_sun_est)': std_dBmv_sun_est,
                 'nhits_zdr': sun_retrieval['nhits_zdr'],
                 'ZDR_sun_est': zdr_sun_est,
                 'std(ZDR_sun_est)': std_zdr_sun_est,
                 'dBm_sun_ref': dBm_sun_ref,
                 'ref_time': ref_time})
            csvfile.close()

    return fname
