"""
pyrad.io.write_data
====================

Functions for writing pyrad output data

.. autosummary::
    :toctree: generated/

    send_msg
    write_alarm_msg
    write_last_state
    write_smn
    write_rhi_profile
    write_field_coverage
    write_cdf
    write_ts_polar_data
    write_ts_cum
    write_monitoring_ts
    write_excess_gates
    write_intercomp_scores_ts
    write_colocated_gates
    write_colocated_data
    write_colocated_data_time_avg
    write_sun_hits
    write_sun_retrieval
    generate_field_name_str

"""

from __future__ import print_function
import glob
import csv
from warnings import warn
import smtplib
from email.message import EmailMessage
import fcntl
import time

import numpy as np

from pyart.config import get_fillvalue

from .io_aux import generate_field_name_str


def send_msg(sender, receiver_list, subject, fname):
    """
    sends the content of a text file by email

    Parameters
    ----------
    sender : str
        the email address of the sender
    receiver_list : list of string
        list with the email addresses of the receiver
    subject : str
        the subject of the email
    fname : str
        name of the file containing the content of the email message

    Returns
    -------
    fname : str
        the name of the file containing the content

    """
    # Create the container email message.
    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = sender
    msg['To'] = ', '.join(receiver_list)

    # Open the plain text file whose name is in fname for reading.
    with open(fname) as fp:
        # Create a text/plain message
        msg.set_content(fp.read())

    # Send the message via our own SMTP server.
    with smtplib.SMTP('localhost') as s:
        s.send_message(msg)

    return fname


def write_alarm_msg(radar_name, param_name_unit, date_last, target, tol_abs,
                    np_trend, value_trend, tol_trend, nevents, np_last,
                    value_last, fname):
    """
    writes an alarm file

    Parameters
    ----------
    radar_name : str
        Name of the radar being controlled
    param_name_unit : str
        Parameter and units
    date_last : datetime object
        date of the current event
    target, tol_abs : float
        Target value and tolerance
    np_trend : int
        Total number of points in trend
    value_trend, tol_trend : float
        Trend value and tolerance
    nevents: int
        Number of events in trend
    np_last : int
        Number of points in the current event
    value_last : float
        Value of the current event
    fname : str
        Name of file where to store the alarm information

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'w', newline='') as txtfile:
        txtfile.write(
            'Weather radar polarimetric parameters monitoring alarm\n')
        if radar_name is not None:
            txtfile.write('Radar name: '+radar_name+'\n')
        txtfile.write('Parameter [Unit]: '+param_name_unit+'\n')
        txtfile.write('Date : '+date_last.strftime('%Y%m%d')+'\n')
        txtfile.write('Target value: '+str(target)+' +/- '+str(tol_abs)+'\n')
        if np_trend > 0:
            txtfile.write(
                'Number of points in trend: '+str(np_trend)+' in ' +
                str(nevents)+' events\n')
            txtfile.write(
                'Trend value: '+str(value_trend)+' (tolerance: +/- ' +
                str(tol_trend)+')\n')
        else:
            txtfile.write('Number of points in trend: NA\n')
            txtfile.write('Trend value: NA\n')
        txtfile.write('Number of points: '+str(np_last)+'\n')
        txtfile.write('Value: '+str(value_last)+'\n')

        txtfile.close()

    return fname


def write_last_state(datetime_last, fname):
    """
    writes SwissMetNet data in format datetime,avg_value, std_value

    Parameters
    ----------
    datetime_last : datetime object
        date and time of the last state
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    try:
        with open(fname, 'w', newline='') as txtfile:
            txtfile.write(datetime_last.strftime('%Y%m%d%H%M%S'))
            txtfile.close()

            return fname
    except EnvironmentError:
        warn('Unable to write on file '+fname)
        return None


def write_smn(datetime_vec, value_avg_vec, value_std_vec, fname):
    """
    writes SwissMetNet data in format datetime,avg_value, std_value

    Parameters
    ----------
    datetime_vec : datetime array
        array containing the measurement time
    value_avg_vec : float array
        array containing the average value
    value_std_vec : float array
        array containing the standard deviation
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    nvalues = len(value_avg_vec)
    with open(fname, 'w', newline='') as csvfile:
        fieldnames = ['datetime', 'avg', 'std']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i in range(nvalues):
            writer.writerow({
                'datetime': datetime_vec[i].strftime('%Y%m%d%H%M%S'),
                'avg': value_avg_vec[i],
                'std': value_std_vec[i]})

        csvfile.close()

    return fname


def write_rhi_profile(hvec, data, nvalid_vec, labels, fname, datatype=None,
                      timeinfo=None, sector=None):
    """
    writes the values of an RHI profile in a text file

    Parameters
    ----------
    hvec : float array
        array containing the alitude in m MSL
    data : list of float array
        the quantities at each altitude
    nvalid_vec : int array
        number of valid data points used to compute the quantiles
    labels : list of strings
        label specifying the quantitites in data
    fname : str
        file name where to store the data
    datatype : str
        the data type
    timeinfo : datetime object
        time of the rhi profile
    sector : dict
        dictionary specying the sector limits

    Returns
    -------
    fname : str
        the name of the file where data has been written

    """
    nvalues = len(hvec)
    data_aux = list()
    for i in range(len(labels)):
        data_aux.append(data[i].filled(fill_value=get_fillvalue()))
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Height Profile\n')
        csvfile.write('# ==============\n')
        csvfile.write('#\n')
        if datatype is None:
            csvfile.write('# Datatype (Unit) : Not specified\n')
        else:
            csvfile.write('# Datatype (Unit) : '+datatype+'\n')
        csvfile.write('# Fill value : '+str(get_fillvalue())+'\n')
        if timeinfo is None:
            csvfile.write('# Time     : Not specified\n')
        else:
            csvfile.write('# Time     : ' +
                          timeinfo.strftime('%Y-%m-%d %H:%M:%S')+'\n')
        if sector is None:
            csvfile.write('# Sector specification: None\n')
        else:
            csvfile.write('# Sector specification:\n')
            csvfile.write('#   Azimuth         : ' + str(sector['az']) +
                          ' deg\n')
            csvfile.write('#   Range start     : ' +
                          str(sector['rmin'])+' m\n')
            csvfile.write('#   Range stop      : ' +
                          str(sector['rmax'])+' m\n')

        fieldnames = ['Altitude [m MSL]']
        fieldnames.extend(labels)
        fieldnames.append('N valid')
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for j in range(nvalues):
            data_dict = {fieldnames[0]: hvec[j]}
            for i in range(len(labels)):
                data_dict.update({fieldnames[i+1]: data_aux[i][j]})
            data_dict.update({fieldnames[-1]: nvalid_vec[j]})
            writer.writerow(data_dict)

        csvfile.close()

    return fname


def write_field_coverage(quantiles, values, ele_start, ele_stop, azi_start,
                         azi_stop, threshold, nvalid_min, datatype, timeinfo,
                         fname):
    """
    writes the quantiles of the coverage on a particular sector

    Parameters
    ----------
    quantiles : datetime array
        array containing the quantiles computed
    values : float array
        quantile value
    ele_start, ele_stop, azi_start, azi_stop : float
        The limits of the sector
    threshold : float
        The minimum value to consider the data valid
    nvalid_min : int
        the minimum number of points to consider that there are values in a
        ray
    datatype : str
        data type and units
    timeinfo : datetime object
        the time stamp of the data
    fname : str
        name of the file where to write the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Quantiles of the field coverage\n')
        csvfile.write('# Datatype (Unit) : '+datatype+'\n')
        csvfile.write('# Time : '+timeinfo.strftime('%Y-%m-%d %H:%M:%S')+'\n')
        csvfile.write('# Sector specification:\n')
        csvfile.write('#   Azimuth start   : '+str(azi_start)+' deg\n')
        csvfile.write('#   Azimuth stop    : '+str(azi_stop)+' deg\n')
        csvfile.write('#   Elevation start : '+str(ele_start)+' deg\n')
        csvfile.write('#   Elevation stop  : '+str(ele_stop)+' deg\n')
        csvfile.write('# Threshold : '+str(threshold)+'\n')
        csvfile.write('# Minimum number of valid gates per ray : ' +
                      str(nvalid_min)+'\n')
        fieldnames = ['Quantile [%]', 'Rain extension [m]']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i, quant in enumerate(quantiles):
            writer.writerow({
                'Quantile [%]': quant,
                'Rain extension [m]': values[i]})

        csvfile.close()

    return fname


def write_cdf(quantiles, values, ntot, nnan, nclut, nblocked, nprec_filter,
              noutliers, ncdf, fname, use_nans=False, nan_value=0.,
              filterprec=[], vismin=None, sector=None, datatype=None,
              timeinfo=None):
    """
    writes a cumulative distribution function

    Parameters
    ----------
    quantiles : datetime array
        array containing the measurement time
    values : float array
        array containing the average value
    fname : float array
        array containing the standard deviation
    sector : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    hydrotype_list = ['NC', 'DS', 'CR', 'LR', 'GR', 'RN', 'VI', 'WS', 'MH',
                      'IH/HDG']
    nvalues = len(values)
    with open(fname, 'w', newline='') as txtfile:
        txtfile.write('Statistical analysis\n')
        txtfile.write('====================\n\n')
        if datatype is None:
            txtfile.write('Datatype (Unit) : Not specified\n')
        else:
            txtfile.write('Datatype (Unit) : '+datatype+'\n')
        if timeinfo is None:
            txtfile.write('Time     : Not specified\n')
        else:
            txtfile.write('Time     : ' +
                          timeinfo.strftime('%Y-%m-%d %H:%M:%S')+'\n')
        if sector is None:
            txtfile.write('Sector specification: None\n')
        else:
            txtfile.write('Sector specification:\n')
            if sector['rmin'] is None:
                txtfile.write('  Range start     : Not specified\n')
            else:
                txtfile.write('  Range start     : ' +
                              str(sector['rmin'])+' m\n')

            if sector['rmax'] is None:
                txtfile.write('  Range stop      : Not specified\n')
            else:
                txtfile.write('  Range stop      : ' +
                              str(sector['rmax'])+' m\n')

            if sector['azmin'] is None:
                txtfile.write('  Azimuth start   : Not specified\n')
            else:
                txtfile.write('  Azimuth start   : ' +
                              str(sector['azmin'])+' deg\n')

            if sector['azmax'] is None:
                txtfile.write('  Azimuth stop    : Not specified\n')
            else:
                txtfile.write('  Azimuth stop    : ' +
                              str(sector['azmax'])+' deg\n')

            if sector['elmin'] is None:
                txtfile.write('  Elevation start : Not specified\n')
            else:
                txtfile.write('  Elevation start : ' +
                              str(sector['elmin'])+' deg\n')

            if sector['elmax'] is None:
                txtfile.write('  Elevation stop  : Not specified\n')
            else:
                txtfile.write('  Elevation stop  : ' +
                              str(sector['elmax'])+' deg\n')

            if sector['hmin'] is None:
                txtfile.write('  Height start    : Not specified\n')
            else:
                txtfile.write('  Height start    : ' +
                              str(sector['hmin'])+' m\n')

            if sector['hmax'] is None:
                txtfile.write('  Height stop     : Not specified\n')
            else:
                txtfile.write('  Height stop     : ' +
                              str(sector['hmax'])+' m\n')
            txtfile.write('')
        txtfile.write('Total number of gates in sector      : ' +
                      str(ntot)+'\n')
        txtfile.write('Number of gates with no value (NaNs) : ' +
                      str(nnan)+'\n')
        if use_nans:
            txtfile.write('  NaNs are set to          : '+str(nan_value)+'\n')
        else:
            txtfile.write('  NaNs are ignored!\n')
        if nclut == -1:
            txtfile.write('Clutter contaminated gates           : ' +
                          'Not checked\n')
        else:
            txtfile.write('Clutter contaminated gates           : ' +
                          str(nclut)+'\n')
        if nblocked == -1:
            txtfile.write('Blocked gates                        : ' +
                          'Not checked\n')
        else:
            txtfile.write('Blocked gates (vismin = '+str(int(vismin)) +
                          ') : '+str(nblocked)+'\n')
        if nprec_filter == -1:
            txtfile.write('Filtered precipitation gates         : None\n')
        else:
            txtfile.write('Filtered precipitation gates         : ' +
                          str(nprec_filter)+'\n')
            txtfile.write('  precipitation types filtered: ')
            for ind_hydro in filterprec:
                txtfile.write(hydrotype_list[ind_hydro]+' ')
            txtfile.write('\n')
        txtfile.write('Number of outliers                   : ' +
                      str(noutliers)+'\n')
        txtfile.write('Number of gates used for histogram   : '+str(ncdf) +
                      ' ('+str(int(ncdf*100./ntot))+'%)\n\n')
        txtfile.write('Quantiles\n')
        txtfile.write('=========\n')
        for i in range(nvalues):
            txtfile.write('  Quantile_'+str(int(quantiles[i]))+' = ' +
                          '{:5.1f}'.format(values[i])+'\n')

        txtfile.close()

    return fname


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
    if not filelist:
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


def write_ts_cum(dataset, fname):
    """
    writes time series accumulation of data

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
    nvalues = len(dataset['time'])
    radar_value = dataset['radar_value'].filled(fill_value=get_fillvalue())
    sensor_value = dataset['sensor_value'].filled(fill_value=get_fillvalue())
    np_radar = dataset['np_radar']
    np_sensor = dataset['np_sensor']

    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Precipitation accumulation data file\n')
        csvfile.write('# Comment lines are preceded by "#"\n')
        csvfile.write('# Description: \n')
        csvfile.write('# Time series of a precipitation accumulation of ' +
                      'weather radar data and another sensor over a ' +
                      'fixed location.\n')
        csvfile.write(
            '# Location [lon, lat, alt]: ' +
            str(dataset['point_coordinates_WGS84_lon_lat_alt']) + '\n')
        csvfile.write(
            '# Nominal antenna coordinates used [az, el, r]: ' +
            str(dataset['antenna_coordinates_az_el_r'])+'\n')
        csvfile.write('# sensor type: ' + dataset['sensor']+'\n')
        csvfile.write('# sensor ID: ' + dataset['sensorid']+'\n')
        csvfile.write('# Data: Precipitation accumulation over ' +
                      str(dataset['cum_time'])+' s\n')
        csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
        csvfile.write(
            '# Start: ' +
            dataset['time'][0].strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
        csvfile.write('#\n')

        fieldnames = ['date', 'np_radar', 'radar_value', 'np_sensor',
                      'sensor_value']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i in range(nvalues):
            writer.writerow({
                'date': dataset['time'][i],
                'np_radar': np_radar[i],
                'radar_value': radar_value[i],
                'np_sensor': np_sensor[i],
                'sensor_value': sensor_value[i]})

        csvfile.close()

    return fname


def write_monitoring_ts(start_time, np_t, values, quantiles, datatype, fname,
                        rewrite=False):
    """
    writes time series of data

    Parameters
    ----------
    start_time : datetime object or array of date time objects
        the time of the monitoring
    np_t : int or array of ints
        the total number of points
    values: float array with 3 elements of array of arrays
        the values at certain quantiles
    quantiles: float array with 3 elements
        the quantiles computed
    fname : str
        file name where to store the data
    rewrite : bool
        if True a new file is created

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    nvalues = np.size(start_time)
    if nvalues == 1:
        start_time_aux = np.asarray([start_time])
        np_t_aux = np.asarray([np_t])
        values_aux = np.asarray([values.filled(fill_value=get_fillvalue())])
    else:
        start_time_aux = np.asarray(start_time)
        values_aux = values.filled(fill_value=get_fillvalue())
        np_t_aux = np_t

    if rewrite:
        file_exists = False
    else:
        filelist = glob.glob(fname)
        if not filelist:
            file_exists = False
        else:
            file_exists = True

    if not file_exists:
        with open(fname, 'w', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)

            csvfile.write(
                '# Weather radar monitoring timeseries data file\n' +
                '# Comment lines are preceded by "#"\n' +
                '# Description: \n' +
                '# Time series of a monitoring of weather radar data.\n' +
                '# Quantiles: '+str(quantiles[1])+', '+str(quantiles[0])+', ' +
                str(quantiles[2])+' percent.\n' +
                '# Data: '+generate_field_name_str(datatype)+'\n' +
                '# Fill Value: '+str(get_fillvalue())+'\n' +
                '# Start: '+start_time_aux[0].strftime(
                    '%Y-%m-%d %H:%M:%S UTC')+'\n' +
                '#\n')

            fieldnames = ['date', 'NP', 'central_quantile', 'low_quantile',
                          'high_quantile']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i in range(nvalues):
                writer.writerow({
                    'date': start_time_aux[i].strftime('%Y%m%d%H%M%S'),
                    'NP': np_t_aux[i],
                    'central_quantile': values_aux[i, 1],
                    'low_quantile': values_aux[i, 0],
                    'high_quantile': values_aux[i, 2]})

            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)

            fieldnames = ['date', 'NP', 'central_quantile', 'low_quantile',
                          'high_quantile']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i in range(nvalues):
                writer.writerow({
                    'date': start_time_aux[i].strftime('%Y%m%d%H%M%S'),
                    'NP': np_t_aux[i],
                    'central_quantile': values_aux[i, 1],
                    'low_quantile': values_aux[i, 0],
                    'high_quantile': values_aux[i, 2]})
            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()
    return fname


def write_excess_gates(excess_dict, fname):
    """
    Writes the position and values of gates that have a frequency of
    occurrence higher than a particular threshold

    Parameters
    ----------
    excess_dict : dict
        dictionary containing the gates parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    ngates = len(excess_dict['ray_ind'])
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Gates exceeding '+str(excess_dict['quant_min']) +
                      ' percentile data file\n')
        csvfile.write('# Comment lines are preceded by "#"\n')
        csvfile.write(
            '# Data collection start time: ' +
            excess_dict['starttime'].strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
        csvfile.write(
            '# Data collection end time: ' +
            excess_dict['endtime'].strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
        csvfile.write('# Number of gates in file: '+str(ngates)+'\n')
        csvfile.write('#\n')

        fieldnames = [
            'ray_ind', 'rng_ind', 'ele', 'azi', 'rng', 'nsamples',
            'occurrence', 'freq_occu']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i in range(ngates):
            writer.writerow({
                'ray_ind': excess_dict['ray_ind'][i],
                'rng_ind': excess_dict['rng_ind'][i],
                'ele': excess_dict['ele'][i],
                'azi': excess_dict['azi'][i],
                'rng': excess_dict['rng'][i],
                'nsamples': excess_dict['nsamples'][i],
                'occurrence': excess_dict['occurrence'][i],
                'freq_occu': excess_dict['freq_occu'][i]})

        csvfile.close()

    return fname


def write_intercomp_scores_ts(start_time, stats, field_name, fname,
                              rad1_name='RADAR001', rad2_name='RADAR002',
                              rewrite=False):
    """
    writes time series of radar intercomparison scores

    Parameters
    ----------
    start_time : datetime object or array of date time objects
        the time of the intercomparison
    stats : dict
        dictionary containing the statistics
    field_name : str
        The name of the field
    fname : str
        file name where to store the data
    rad1_name, rad2_name : str
        Name of the radars intercompared
    rewrite : bool
        if True a new file is created

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    nvalues = np.size(start_time)

    meanbias = stats['meanbias'].filled(fill_value=get_fillvalue())
    medianbias = stats['medianbias'].filled(fill_value=get_fillvalue())
    quant25bias = stats['quant25bias'].filled(fill_value=get_fillvalue())
    quant75bias = stats['quant75bias'].filled(fill_value=get_fillvalue())
    modebias = stats['modebias'].filled(fill_value=get_fillvalue())
    corr = stats['corr'].filled(fill_value=get_fillvalue())
    slope = stats['slope'].filled(fill_value=get_fillvalue())
    intercep = stats['intercep'].filled(fill_value=get_fillvalue())
    intercep_slope_1 = stats['intercep_slope_1'].filled(
        fill_value=get_fillvalue())

    if nvalues == 1:
        start_time_aux = np.asarray([start_time])
        meanbias = np.asarray([meanbias])
        medianbias = np.asarray([medianbias])
        quant25bias = np.asarray([quant25bias])
        quant75bias = np.asarray([quant75bias])
        modebias = np.asarray([modebias])
        corr = np.asarray([corr])
        slope = np.asarray([slope])
        intercep = np.asarray([intercep])
        intercep_slope_1 = np.asarray([intercep_slope_1])
        np_t = np.asarray([stats['npoints']])
    else:
        start_time_aux = np.asarray(start_time)
        np_t = stats['npoints']

    if rewrite:
        file_exists = False
    else:
        filelist = glob.glob(fname)
        if not filelist:
            file_exists = False
        else:
            file_exists = True

    if not file_exists:
        with open(fname, 'w', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)

            csvfile.write(
                '# Weather radar intercomparison scores timeseries file\n' +
                '# Comment lines are preceded by "#"\n' +
                '# Description: \n' +
                '# Time series of the intercomparison between two radars.\n' +
                '# Radar 1: '+rad1_name+'\n' +
                '# Radar 2: '+rad2_name+'\n' +
                '# Field name: '+field_name+'\n' +
                '# Fill Value: '+str(get_fillvalue())+'\n' +
                '# Start: '+start_time_aux[0].strftime(
                    '%Y-%m-%d %H:%M:%S UTC')+'\n' +
                '#\n')

            fieldnames = ['date', 'NP', 'mean_bias', 'median_bias',
                          'quant25_bias', 'quant75_bias', 'mode_bias', 'corr',
                          'slope_of_linear_regression',
                          'intercep_of_linear_regression',
                          'intercep_of_linear_regression_of_slope_1']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()

            for i in range(nvalues):
                writer.writerow({
                    'date': start_time_aux[i].strftime('%Y%m%d%H%M%S'),
                    'NP': np_t[i],
                    'mean_bias': meanbias[i],
                    'median_bias': medianbias[i],
                    'quant25_bias': quant25bias[i],
                    'quant75_bias': quant75bias[i],
                    'mode_bias': modebias[i],
                    'corr': corr[i],
                    'slope_of_linear_regression': slope[i],
                    'intercep_of_linear_regression': intercep[i],
                    'intercep_of_linear_regression_of_slope_1': (
                        intercep_slope_1[i])
                    })

            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)

            fieldnames = ['date', 'NP', 'mean_bias', 'median_bias',
                          'quant25_bias', 'quant75_bias', 'mode_bias', 'corr',
                          'slope_of_linear_regression',
                          'intercep_of_linear_regression',
                          'intercep_of_linear_regression_of_slope_1']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i in range(nvalues):
                writer.writerow({
                    'date': start_time_aux[i].strftime('%Y%m%d%H%M%S'),
                    'NP': np_t[i],
                    'mean_bias': meanbias[i],
                    'median_bias': medianbias[i],
                    'quant25_bias': quant25bias[i],
                    'quant75_bias': quant75bias[i],
                    'mode_bias': modebias[i],
                    'corr': corr[i],
                    'slope_of_linear_regression': slope[i],
                    'intercep_of_linear_regression': intercep[i],
                    'intercep_of_linear_regression_of_slope_1': (
                        intercep_slope_1[i])
                    })
            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()

    return fname


def write_colocated_gates(coloc_gates, fname):
    """
    Writes the position of gates colocated with two radars

    Parameters
    ----------
    coloc_gates : dict
        dictionary containing the colocated gates parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    ngates = len(coloc_gates['rad1_ele'])
    with open(fname, 'w', newline='') as csvfile:
        csvfile.write('# Colocated radar gates data file\n')
        csvfile.write('# Comment lines are preceded by "#"\n')
        csvfile.write('#\n')

        fieldnames = [
            'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele', 'rad1_azi', 'rad1_rng',
            'rad2_ray_ind', 'rad2_rng_ind', 'rad2_ele', 'rad2_azi', 'rad2_rng']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i in range(ngates):
            writer.writerow({
                'rad1_ray_ind': coloc_gates['rad1_ray_ind'][i],
                'rad1_rng_ind': coloc_gates['rad1_rng_ind'][i],
                'rad1_ele': coloc_gates['rad1_ele'][i],
                'rad1_azi': coloc_gates['rad1_azi'][i],
                'rad1_rng': coloc_gates['rad1_rng'][i],
                'rad2_ray_ind': coloc_gates['rad2_ray_ind'][i],
                'rad2_rng_ind': coloc_gates['rad2_rng_ind'][i],
                'rad2_ele': coloc_gates['rad2_ele'][i],
                'rad2_azi': coloc_gates['rad2_azi'][i],
                'rad2_rng': coloc_gates['rad2_rng'][i]})

        csvfile.close()

    return fname


def write_colocated_data(coloc_data, fname):
    """
    Writes the data of gates colocated with two radars

    Parameters
    ----------
    coloc_data : dict
        dictionary containing the colocated data parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    filelist = glob.glob(fname)
    ngates = len(coloc_data['rad1_ele'])
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Colocated radar gates data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('#\n')

            fieldnames = [
                'rad1_time', 'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele',
                'rad1_azi', 'rad1_rng', 'rad1_val', 'rad2_time',
                'rad2_ray_ind', 'rad2_rng_ind', 'rad2_ele', 'rad2_azi',
                'rad2_rng', 'rad2_val']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i in range(ngates):
                writer.writerow({
                    'rad1_time': (
                        coloc_data['rad1_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad1_ray_ind': coloc_data['rad1_ray_ind'][i],
                    'rad1_rng_ind': coloc_data['rad1_rng_ind'][i],
                    'rad1_ele': coloc_data['rad1_ele'][i],
                    'rad1_azi': coloc_data['rad1_azi'][i],
                    'rad1_rng': coloc_data['rad1_rng'][i],
                    'rad1_val': coloc_data['rad1_val'][i],
                    'rad2_time': (
                        coloc_data['rad2_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad2_ray_ind': coloc_data['rad2_ray_ind'][i],
                    'rad2_rng_ind': coloc_data['rad2_rng_ind'][i],
                    'rad2_ele': coloc_data['rad2_ele'][i],
                    'rad2_azi': coloc_data['rad2_azi'][i],
                    'rad2_rng': coloc_data['rad2_rng'][i],
                    'rad2_val': coloc_data['rad2_val'][i]})

            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'rad1_time', 'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele',
                'rad1_azi', 'rad1_rng', 'rad1_val', 'rad2_time',
                'rad2_ray_ind', 'rad2_rng_ind', 'rad2_ele', 'rad2_azi',
                'rad2_rng', 'rad2_val']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i in range(ngates):
                writer.writerow({
                    'rad1_time': (
                        coloc_data['rad1_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad1_ray_ind': coloc_data['rad1_ray_ind'][i],
                    'rad1_rng_ind': coloc_data['rad1_rng_ind'][i],
                    'rad1_ele': coloc_data['rad1_ele'][i],
                    'rad1_azi': coloc_data['rad1_azi'][i],
                    'rad1_rng': coloc_data['rad1_rng'][i],
                    'rad1_val': coloc_data['rad1_val'][i],
                    'rad2_time': (
                        coloc_data['rad2_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad2_ray_ind': coloc_data['rad2_ray_ind'][i],
                    'rad2_rng_ind': coloc_data['rad2_rng_ind'][i],
                    'rad2_ele': coloc_data['rad2_ele'][i],
                    'rad2_azi': coloc_data['rad2_azi'][i],
                    'rad2_rng': coloc_data['rad2_rng'][i],
                    'rad2_val': coloc_data['rad2_val'][i]})
            csvfile.close()

    return fname


def write_colocated_data_time_avg(coloc_data, fname):
    """
    Writes the time averaged data of gates colocated with two radars

    Parameters
    ----------
    coloc_data : dict
        dictionary containing the colocated data parameters
    fname : str
        file name where to store the data

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    filelist = glob.glob(fname)
    ngates = len(coloc_data['rad1_ele'])
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Colocated radar gates data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('#\n')

            fieldnames = [
                'rad1_time', 'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele',
                'rad1_azi', 'rad1_rng', 'rad1_dBZavg', 'rad1_PhiDPavg',
                'rad1_Flagavg', 'rad2_time', 'rad2_ray_ind', 'rad2_rng_ind',
                'rad2_ele', 'rad2_azi', 'rad2_rng', 'rad2_dBZavg',
                'rad2_PhiDPavg', 'rad2_Flagavg']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i in range(ngates):
                writer.writerow({
                    'rad1_time': (
                        coloc_data['rad1_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad1_ray_ind': coloc_data['rad1_ray_ind'][i],
                    'rad1_rng_ind': coloc_data['rad1_rng_ind'][i],
                    'rad1_ele': coloc_data['rad1_ele'][i],
                    'rad1_azi': coloc_data['rad1_azi'][i],
                    'rad1_rng': coloc_data['rad1_rng'][i],
                    'rad1_dBZavg': coloc_data['rad1_dBZavg'][i],
                    'rad1_PhiDPavg': coloc_data['rad1_PhiDPavg'][i],
                    'rad1_Flagavg': coloc_data['rad1_Flagavg'][i],
                    'rad2_time': (
                        coloc_data['rad2_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad2_ray_ind': coloc_data['rad2_ray_ind'][i],
                    'rad2_rng_ind': coloc_data['rad2_rng_ind'][i],
                    'rad2_ele': coloc_data['rad2_ele'][i],
                    'rad2_azi': coloc_data['rad2_azi'][i],
                    'rad2_rng': coloc_data['rad2_rng'][i],
                    'rad2_dBZavg': coloc_data['rad2_dBZavg'][i],
                    'rad2_PhiDPavg': coloc_data['rad2_PhiDPavg'][i],
                    'rad2_Flagavg': coloc_data['rad2_Flagavg'][i]})

            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'rad1_time', 'rad1_ray_ind', 'rad1_rng_ind', 'rad1_ele',
                'rad1_azi', 'rad1_rng', 'rad1_dBZavg', 'rad1_PhiDPavg',
                'rad1_Flagavg', 'rad2_time', 'rad2_ray_ind', 'rad2_rng_ind',
                'rad2_ele', 'rad2_azi', 'rad2_rng', 'rad2_dBZavg',
                'rad2_PhiDPavg', 'rad2_Flagavg']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i in range(ngates):
                writer.writerow({
                    'rad1_time': (
                        coloc_data['rad1_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad1_ray_ind': coloc_data['rad1_ray_ind'][i],
                    'rad1_rng_ind': coloc_data['rad1_rng_ind'][i],
                    'rad1_ele': coloc_data['rad1_ele'][i],
                    'rad1_azi': coloc_data['rad1_azi'][i],
                    'rad1_rng': coloc_data['rad1_rng'][i],
                    'rad1_dBZavg': coloc_data['rad1_dBZavg'][i],
                    'rad1_PhiDPavg': coloc_data['rad1_PhiDPavg'][i],
                    'rad1_Flagavg': coloc_data['rad1_Flagavg'][i],
                    'rad2_time': (
                        coloc_data['rad2_time'][i].strftime('%Y%m%d%H%M%S')),
                    'rad2_ray_ind': coloc_data['rad2_ray_ind'][i],
                    'rad2_rng_ind': coloc_data['rad2_rng_ind'][i],
                    'rad2_ele': coloc_data['rad2_ele'][i],
                    'rad2_azi': coloc_data['rad2_azi'][i],
                    'rad2_rng': coloc_data['rad2_rng'][i],
                    'rad2_dBZavg': coloc_data['rad2_dBZavg'][i],
                    'rad2_PhiDPavg': coloc_data['rad2_PhiDPavg'][i],
                    'rad2_Flagavg': coloc_data['rad2_Flagavg'][i]})
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
    if not filelist:
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
                writer.writerow({
                    'time': sun_hits['time'][i].strftime(
                        '%Y-%m-%d %H:%M:%S.%f'),
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
                writer.writerow({
                    'time': sun_hits['time'][i].strftime(
                        '%Y-%m-%d %H:%M:%S.%f'),
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
    sf_h = sun_retrieval['sf_h'].filled(fill_value=get_fillvalue())

    el_width_v = sun_retrieval['el_width_v'].filled(fill_value=get_fillvalue())
    az_width_v = sun_retrieval['az_width_v'].filled(fill_value=get_fillvalue())
    el_bias_v = sun_retrieval['el_bias_v'].filled(fill_value=get_fillvalue())
    az_bias_v = sun_retrieval['az_bias_v'].filled(fill_value=get_fillvalue())
    dBmv_sun_est = sun_retrieval['dBmv_sun_est'].filled(
        fill_value=get_fillvalue())
    std_dBmv_sun_est = sun_retrieval['std(dBmv_sun_est)'].filled(
        fill_value=get_fillvalue())
    sf_v = sun_retrieval['sf_v'].filled(fill_value=get_fillvalue())

    zdr_sun_est = sun_retrieval['ZDR_sun_est'].filled(
        fill_value=get_fillvalue())
    std_zdr_sun_est = sun_retrieval['std(ZDR_sun_est)'].filled(
        fill_value=get_fillvalue())
    sf_ref = sun_retrieval['sf_ref'].filled(
        fill_value=get_fillvalue())
    ref_time = 'None'
    if sun_retrieval['ref_time'] is not None:
        ref_time = sun_retrieval['ref_time'].strftime('%Y%m%d%H%M%S')

    filelist = glob.glob(fname)
    if not filelist:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar sun retrievals data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write('#\n')

            fieldnames = [
                'first_hit_time', 'last_hit_time',
                'nhits_h', 'el_width_h', 'az_width_h',
                'el_bias_h', 'az_bias_h', 'dBm_sun_est', 'std(dBm_sun_est)',
                'sf_h',
                'nhits_v', 'el_width_v', 'az_width_v',
                'el_bias_v', 'az_bias_v', 'dBmv_sun_est', 'std(dBmv_sun_est)',
                'sf_v',
                'nhits_zdr', 'ZDR_sun_est', 'std(ZDR_sun_est)', 'sf_ref',
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
                 'sf_h': sf_h,
                 'nhits_v': sun_retrieval['nhits_v'],
                 'el_width_v': el_width_v,
                 'az_width_v': az_width_v,
                 'el_bias_v': el_bias_v,
                 'az_bias_v': az_bias_v,
                 'dBmv_sun_est': dBmv_sun_est,
                 'std(dBmv_sun_est)': std_dBmv_sun_est,
                 'sf_v': sf_v,
                 'nhits_zdr': sun_retrieval['nhits_zdr'],
                 'ZDR_sun_est': zdr_sun_est,
                 'std(ZDR_sun_est)': std_zdr_sun_est,
                 'sf_ref': sf_ref,
                 'ref_time': ref_time})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'first_hit_time', 'last_hit_time',
                'nhits_h', 'el_width_h', 'az_width_h',
                'el_bias_h', 'az_bias_h', 'dBm_sun_est', 'std(dBm_sun_est)',
                'sf_h',
                'nhits_v', 'el_width_v', 'az_width_v',
                'el_bias_v', 'az_bias_v', 'dBmv_sun_est', 'std(dBmv_sun_est)',
                'sf_v',
                'nhits_zdr', 'ZDR_sun_est', 'std(ZDR_sun_est)', 'sf_ref',
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
                 'sf_h': sf_h,
                 'nhits_v': sun_retrieval['nhits_v'],
                 'el_width_v': el_width_v,
                 'az_width_v': az_width_v,
                 'el_bias_v': el_bias_v,
                 'az_bias_v': az_bias_v,
                 'dBmv_sun_est': dBmv_sun_est,
                 'std(dBmv_sun_est)': std_dBmv_sun_est,
                 'sf_v': sf_v,
                 'nhits_zdr': sun_retrieval['nhits_zdr'],
                 'ZDR_sun_est': zdr_sun_est,
                 'std(ZDR_sun_est)': std_zdr_sun_est,
                 'sf_ref': sf_ref,
                 'ref_time': ref_time})
            csvfile.close()

    return fname
