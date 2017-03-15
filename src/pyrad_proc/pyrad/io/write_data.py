"""
pyrad.io.write_data
====================

Functions for writing pyrad output data

.. autosummary::
    :toctree: generated/

    write_smn
    write_cdf
    write_ts_polar_data
    write_ts_cum
    write_monitoring_ts
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

from pyart.config import get_fillvalue, get_metadata

from .io_aux import generate_field_name_str


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
            writer.writerow(
                {'datetime': datetime_vec[i].strftime('%Y%m%d%H%M%S'),
                 'avg': value_avg_vec[i],
                 'std': value_std_vec[i]})
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
            for i in range(len(filterprec)):
                txtfile.write(hydrotype_list[filterprec[i]]+' ')
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
            writer.writerow(
                {'date': dataset['time'][i],
                 'np_radar': np_radar[i],
                 'radar_value': radar_value[i],
                 'np_sensor': np_sensor[i],
                 'sensor_value': sensor_value[i]})
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


def write_intercomp_scores_ts(start_time, stats, field_name, fname,
                              rad1_name='RADAR001', rad2_name='RADAR002'):
    """
    writes time series of radar intercomparison scores

    Parameters
    ----------
    start_time : datetime object
        the time of the intercomparison
    stats : dict
        dictionary containing the statistics
    field_name : str
        The name of the field
    fname : str
        file name where to store the data
    rad1_name, rad2_name : str
        Name of the radars intercompared

    Returns
    -------
    fname : str
        the name of the file where data has written

    """
    meanbias = stats['meanbias'].filled(fill_value=get_fillvalue())
    medianbias = stats['medianbias'].filled(fill_value=get_fillvalue())
    modebias = stats['modebias'].filled(fill_value=get_fillvalue())
    corr = stats['corr'].filled(fill_value=get_fillvalue())
    slope = stats['slope'].filled(fill_value=get_fillvalue())
    intercep = stats['intercep'].filled(fill_value=get_fillvalue())
    intercep_slope_1 = stats['intercep_slope_1'].filled(
        fill_value=get_fillvalue())

    filelist = glob.glob(fname)
    if len(filelist) == 0:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Weather radar intercomparison scores ' +
                          'timeseries file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('# Description: \n')
            csvfile.write(
                '# Time series of the intercomparison between two radars.\n')
            csvfile.write('# Radar 1: '+rad1_name+'\n')
            csvfile.write('# Radar 2: '+rad2_name+'\n')
            csvfile.write('# Fill Value: '+str(get_fillvalue())+'\n')
            csvfile.write(
                '# Start: ' +
                start_time.strftime('%Y-%m-%d %H:%M:%S UTC')+'\n')
            csvfile.write('#\n')

            fieldnames = ['date', 'NP', 'mean_bias', 'median_bias',
                          'mode_bias', 'corr', 'slope_of_linear_regression',
                          'intercep_of_linear_regression',
                          'intercep_of_linear_regression_of_slope_1']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()

            writer.writerow(
                {'date': start_time.strftime('%Y%m%d%H%M%S'),
                 'NP': stats['npoints'],
                 'mean_bias': meanbias,
                 'median_bias': medianbias,
                 'mode_bias': modebias,
                 'corr': corr,
                 'slope_of_linear_regression': slope,
                 'intercep_of_linear_regression': intercep,
                 'intercep_of_linear_regression_of_slope_1': intercep_slope_1
                 })
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = ['date', 'NP', 'mean_bias', 'median_bias',
                          'mode_bias', 'corr', 'slope_of_linear_regression',
                          'intercep_of_linear_regression',
                          'intercep_of_linear_regression_of_slope_1']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writerow(
                {'date': start_time.strftime('%Y%m%d%H%M%S'),
                 'NP': stats['npoints'],
                 'mean_bias': meanbias,
                 'median_bias': medianbias,
                 'mode_bias': modebias,
                 'corr': corr,
                 'slope_of_linear_regression': slope,
                 'intercep_of_linear_regression': intercep,
                 'intercep_of_linear_regression_of_slope_1': intercep_slope_1
                 })
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
            'rad1_ele', 'rad1_azi', 'rad1_rng',
            'rad2_ele', 'rad2_azi', 'rad2_rng']
        writer = csv.DictWriter(csvfile, fieldnames)
        writer.writeheader()
        for i in range(ngates):
            writer.writerow(
                {'rad1_ele': coloc_gates['rad1_ele'][i],
                 'rad1_azi': coloc_gates['rad1_azi'][i],
                 'rad1_rng': coloc_gates['rad1_rng'][i],
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
    if len(filelist) == 0:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Colocated radar gates data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('#\n')

            fieldnames = [
                'rad1_ele', 'rad1_azi', 'rad1_rng', 'rad1_val',
                'rad2_ele', 'rad2_azi', 'rad2_rng', 'rad2_val']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i in range(ngates):
                writer.writerow(
                    {'rad1_ele': coloc_data['rad1_ele'][i],
                     'rad1_azi': coloc_data['rad1_azi'][i],
                     'rad1_rng': coloc_data['rad1_rng'][i],
                     'rad1_val': coloc_data['rad1_val'][i],
                     'rad2_ele': coloc_data['rad2_ele'][i],
                     'rad2_azi': coloc_data['rad2_azi'][i],
                     'rad2_rng': coloc_data['rad2_rng'][i],
                     'rad2_val': coloc_data['rad2_val'][i]})
            csvfile.close()
    else:
        with open(fname, 'a', newline='') as csvfile:
            fieldnames = [
                'rad1_ele', 'rad1_azi', 'rad1_rng', 'rad1_val',
                'rad2_ele', 'rad2_azi', 'rad2_rng', 'rad2_val']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i in range(ngates):
                writer.writerow(
                    {'rad1_ele': coloc_data['rad1_ele'][i],
                     'rad1_azi': coloc_data['rad1_azi'][i],
                     'rad1_rng': coloc_data['rad1_rng'][i],
                     'rad1_val': coloc_data['rad1_val'][i],
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
    if len(filelist) == 0:
        with open(fname, 'w', newline='') as csvfile:
            csvfile.write('# Colocated radar gates data file\n')
            csvfile.write('# Comment lines are preceded by "#"\n')
            csvfile.write('#\n')

            fieldnames = [
                'rad1_ele', 'rad1_azi', 'rad1_rng',
                'rad1_dBZavg', 'rad1_PhiDPavg', 'rad1_Flagavg',
                'rad2_ele', 'rad2_azi', 'rad2_rng',
                'rad2_dBZavg', 'rad2_PhiDPavg', 'rad2_Flagavg']
            writer = csv.DictWriter(csvfile, fieldnames)
            writer.writeheader()
            for i in range(ngates):
                writer.writerow(
                    {'rad1_ele': coloc_data['rad1_ele'][i],
                     'rad1_azi': coloc_data['rad1_azi'][i],
                     'rad1_rng': coloc_data['rad1_rng'][i],
                     'rad1_dBZavg': coloc_data['rad1_dBZavg'][i],
                     'rad1_PhiDPavg': coloc_data['rad1_PhiDPavg'][i],
                     'rad1_dBZavg': coloc_data['rad1_dBZavg'][i],
                     'rad1_Flagavg': coloc_data['rad1_Flagavg'][i],
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
                'rad1_ele', 'rad1_azi', 'rad1_rng',
                'rad1_dBZavg', 'rad1_PhiDPavg', 'rad1_Flagavg',
                'rad2_ele', 'rad2_azi', 'rad2_rng',
                'rad2_dBZavg', 'rad2_PhiDPavg', 'rad2_Flagavg']
            writer = csv.DictWriter(csvfile, fieldnames)
            for i in range(ngates):
                writer.writerow(
                    {'rad1_ele': coloc_data['rad1_ele'][i],
                     'rad1_azi': coloc_data['rad1_azi'][i],
                     'rad1_rng': coloc_data['rad1_rng'][i],
                     'rad1_dBZavg': coloc_data['rad1_dBZavg'][i],
                     'rad1_PhiDPavg': coloc_data['rad1_PhiDPavg'][i],
                     'rad1_dBZavg': coloc_data['rad1_dBZavg'][i],
                     'rad1_Flagavg': coloc_data['rad1_Flagavg'][i],
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
