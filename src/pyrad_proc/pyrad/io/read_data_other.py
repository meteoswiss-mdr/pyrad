"""
pyrad.io.read_data_other
========================

Functions for reading auxiliary data

.. autosummary::
    :toctree: generated/

    read_last_state
    read_status
    read_rad4alp_cosmo
    read_rad4alp_vis
    read_excess_gates
    read_colocated_gates
    read_colocated_data
    read_colocated_data_time_avg
    read_timeseries
    read_ts_cum
    read_monitoring_ts
    read_monitoring_ts_old
    read_intercomp_scores_ts
    read_intercomp_scores_ts_old
    read_intercomp_scores_ts_old_v0
    read_selfconsistency
    read_antenna_pattern
"""

import os
import glob
import datetime
import csv
import xml.etree.ElementTree as et
from warnings import warn
import fcntl
import time

import numpy as np

from pyart.config import get_fillvalue, get_metadata

from .io_aux import get_fieldname_pyart


def read_last_state(fname):
    """
    Reads a file containing the date of acquisition of the last volume
    processed

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    last_state : datetime object
        the date

    """
    try:
        with open(fname, 'r', newline='') as txtfile:
            line = txtfile.readline()
            txtfile.close()
            if not line:
                warn('File '+fname+' is empty.')
                return None
            try:
                last_state = datetime.datetime.strptime(
                    line, '%Y%m%d%H%M%S')
                return last_state
            except ValueError:
                warn('File '+fname+' does not contain a valid date.')
                return None
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None


def read_status(voltime, cfg, ind_rad=0):
    """
    Reads rad4alp xml status file.

    Parameters
    ----------
    voltime : datetime object
        volume scan time

    cfg: dictionary of dictionaries
        configuration info to figure out where the data is
    ind_rad: int
        radar index

    Returns
    -------
    root : root element object
        The information contained in the status file

    """
    if cfg['RadarName'] is None:
        raise ValueError(
            'ERROR: Radar Name not specified in config file. \
            Unable to read status data')

    dayinfo = voltime.strftime('%y%j')
    timeinfo = voltime.strftime('%H%M')
    basename = 'ST'+cfg['RadarName'][ind_rad]+dayinfo
    if cfg['path_convention'] == 'RT':
        datapath = cfg['datapath'][ind_rad]+'ST'+cfg['RadarName'][ind_rad]+'/'
    else:
        datapath = cfg['datapath'][ind_rad]+dayinfo+'/'+basename+'/'
    filename = glob.glob(datapath+basename+timeinfo+'*.xml')
    if not filename:
        warn('rad4alp status file '+datapath+basename+timeinfo +
             '*.xml not found')
        return None
    try:
        root = et.parse(filename[0]).getroot()
        return root
    except et.ParseError as ee:
        warn(str(ee))
        warn('Unable to read file '+filename[0])

        return None


def read_rad4alp_cosmo(fname, datatype, ngates=0):
    """
    Reads rad4alp COSMO data binary file.

    Parameters
    ----------
    fname : str
        name of the file to read
    datatype : str
        name of the data type
    ngates : int
        maximum number of range gates per ray. If larger than 0 the radar
        field will be cut accordingly.

    Returns
    -------
    field : dictionary
        The data field

    """
    try:
        bindata = np.fromfile(fname, dtype='uint8', count=-1)
        nbins = bindata.size
        naz = 360
        nr = int(nbins/naz)
        bindata = np.reshape(bindata, (naz, nr))
        mask = bindata == 0

        if datatype == 'TEMP':
            field_name = get_fieldname_pyart(datatype)
            field_data = np.ma.masked_where(
                mask, (bindata-1).astype(float)*0.5-87.)

            field = get_metadata(field_name)
            field['data'] = field_data
            if ngates > 0:
                field['data'] = field['data'][:, :ngates]
            return field
        elif datatype == 'ISO0':
            field_name = get_fieldname_pyart(datatype)
            field_data = np.ma.masked_where(mask, (bindata-1).astype(float))

            field = get_metadata(field_name)
            field['data'] = field_data
            if ngates > 0:
                field['data'] = field['data'][:, :ngates]
            return field
        else:
            warn('Unknown COSMO data type '+datatype)
            return None
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None


def read_rad4alp_vis(fname, datatype):
    """
    Reads rad4alp visibility data binary file.

    Parameters
    ----------
    fname : str
        name of the file to read

    datatype : str
        name of the data type

    Returns
    -------
    field_list : list of dictionaries
        A data field. Each element of the list corresponds to one elevation

    """
    if datatype != 'VIS':
        warn('Unknown DEM data type '+datatype)
        return None

    header_size = 64
    nel = 20
    naz = 360
    nr = [492, 420, 492, 324, 366, 324, 292, 324, 280, 242,
          222, 200, 174, 150, 124, 100, 82, 68, 60, 54]

    try:
        with open(fname, 'rb') as visfile:
            field_list = list()
            # skip header
            visfile.seek(header_size, os.SEEK_SET)

            bindata = np.fromfile(visfile, dtype='uint8', count=-1)

            sweep_start_index = 0
            for i in range(nel):
                sweep_end_index = sweep_start_index+naz*nr[i]-1
                field_data = np.reshape(
                    bindata[sweep_start_index:sweep_end_index+1],
                    (naz, nr[i])).astype(float)
                sweep_start_index = sweep_end_index+1

                field_name = get_fieldname_pyart(datatype)
                field = get_metadata(field_name)
                field['data'] = field_data
                field_list.append(field)

            visfile.close()

            return field_list

    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None


def read_excess_gates(fname):
    """
    Reads a csv files containing the position of gates exceeding
    a given percentile of frequency of occurrence

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
    rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng : tupple
        A tupple with the data read. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)
            ray_ind = np.empty(nrows, dtype=int)
            rng_ind = np.empty(nrows, dtype=int)
            ele = np.empty(nrows, dtype=float)
            azi = np.empty(nrows, dtype=float)
            rng = np.empty(nrows, dtype=float)
            nsamples = np.empty(nrows, dtype=int)
            occurrence = np.empty(nrows, dtype=int)
            freq_occu = np.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            i = 0
            for row in reader:
                ray_ind[i] = int(row['ray_ind'])
                rng_ind[i] = int(row['rng_ind'])
                ele[i] = float(row['ele'])
                azi[i] = float(row['azi'])
                rng[i] = float(row['rng'])
                nsamples[i] = int(row['nsamples'])
                occurrence[i] = int(row['occurrence'])
                freq_occu[i] = float(row['freq_occu'])
                i += 1

            csvfile.close()

            return (ray_ind, rng_ind, ele, azi, rng, nsamples, occurrence,
                    freq_occu)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None, None, None, None


def read_colocated_gates(fname):
    """
    Reads a csv files containing the position of colocated gates

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
    rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng : tupple
        A tupple with the data read. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)
            rad1_ray_ind = np.empty(nrows, dtype=int)
            rad1_rng_ind = np.empty(nrows, dtype=int)
            rad1_ele = np.empty(nrows, dtype=float)
            rad1_azi = np.empty(nrows, dtype=float)
            rad1_rng = np.empty(nrows, dtype=float)
            rad2_ray_ind = np.empty(nrows, dtype=int)
            rad2_rng_ind = np.empty(nrows, dtype=int)
            rad2_ele = np.empty(nrows, dtype=float)
            rad2_azi = np.empty(nrows, dtype=float)
            rad2_rng = np.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            i = 0
            for row in reader:
                rad1_ray_ind[i] = int(row['rad1_ray_ind'])
                rad1_rng_ind[i] = int(row['rad1_rng_ind'])
                rad1_ele[i] = float(row['rad1_ele'])
                rad1_azi[i] = float(row['rad1_azi'])
                rad1_rng[i] = float(row['rad1_rng'])
                rad2_ray_ind[i] = int(row['rad2_ray_ind'])
                rad2_rng_ind[i] = int(row['rad2_rng_ind'])
                rad2_ele[i] = float(row['rad2_ele'])
                rad2_azi[i] = float(row['rad2_azi'])
                rad2_rng[i] = float(row['rad2_rng'])
                i += 1

            csvfile.close()

            return (rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
                    rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None, None, None, None, None, None


def read_colocated_data(fname):
    """
    Reads a csv files containing colocated data

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    rad1_time, rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
    rad1_val, rad2_time, rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi,
    rad2_rng, rad2_val : tupple
        A tupple with the data read. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)
            rad1_time = np.empty(nrows, dtype=datetime.datetime)
            rad1_ray_ind = np.empty(nrows, dtype=int)
            rad1_rng_ind = np.empty(nrows, dtype=int)
            rad1_ele = np.empty(nrows, dtype=float)
            rad1_azi = np.empty(nrows, dtype=float)
            rad1_rng = np.empty(nrows, dtype=float)
            rad1_val = np.empty(nrows, dtype=float)
            rad2_time = np.empty(nrows, dtype=datetime.datetime)
            rad2_ray_ind = np.empty(nrows, dtype=int)
            rad2_rng_ind = np.empty(nrows, dtype=int)
            rad2_ele = np.empty(nrows, dtype=float)
            rad2_azi = np.empty(nrows, dtype=float)
            rad2_rng = np.empty(nrows, dtype=float)
            rad2_val = np.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            i = 0
            for row in reader:
                rad1_time[i] = datetime.datetime.strptime(
                    row['rad1_time'], '%Y%m%d%H%M%S')
                rad1_ray_ind[i] = int(row['rad1_ray_ind'])
                rad1_rng_ind[i] = int(row['rad1_rng_ind'])
                rad1_ele[i] = float(row['rad1_ele'])
                rad1_azi[i] = float(row['rad1_azi'])
                rad1_rng[i] = float(row['rad1_rng'])
                rad1_val[i] = float(row['rad1_val'])
                rad2_time[i] = datetime.datetime.strptime(
                    row['rad2_time'], '%Y%m%d%H%M%S')
                rad2_ray_ind[i] = int(row['rad2_ray_ind'])
                rad2_rng_ind[i] = int(row['rad2_rng_ind'])
                rad2_ele[i] = float(row['rad2_ele'])
                rad2_azi[i] = float(row['rad2_azi'])
                rad2_rng[i] = float(row['rad2_rng'])
                rad2_val[i] = float(row['rad2_val'])
                i += 1

            csvfile.close()

            return (rad1_time, rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi,
                    rad1_rng, rad1_val, rad2_time, rad2_ray_ind, rad2_rng_ind,
                    rad2_ele, rad2_azi, rad2_rng, rad2_val)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None, None, None, None)


def read_colocated_data_time_avg(fname):
    """
    Reads a csv files containing time averaged colocated data

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    rad1_time, rad1_ray_ind, rad1_rng_ind, rad1_ele , rad1_azi, rad1_rng,
    rad1_val, rad2_time, rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi,
    rad2_rng, rad2_val : tupple
        A tupple with the data read. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)
            rad1_time = np.empty(nrows, dtype=datetime.datetime)
            rad1_ray_ind = np.empty(nrows, dtype=int)
            rad1_rng_ind = np.empty(nrows, dtype=int)
            rad1_ele = np.empty(nrows, dtype=float)
            rad1_azi = np.empty(nrows, dtype=float)
            rad1_rng = np.empty(nrows, dtype=float)
            rad1_dBZavg = np.empty(nrows, dtype=float)
            rad1_PhiDPavg = np.empty(nrows, dtype=float)
            rad1_Flagavg = np.empty(nrows, dtype=float)
            rad2_time = np.empty(nrows, dtype=datetime.datetime)
            rad2_ray_ind = np.empty(nrows, dtype=int)
            rad2_rng_ind = np.empty(nrows, dtype=int)
            rad2_ele = np.empty(nrows, dtype=float)
            rad2_azi = np.empty(nrows, dtype=float)
            rad2_rng = np.empty(nrows, dtype=float)
            rad2_dBZavg = np.empty(nrows, dtype=float)
            rad2_PhiDPavg = np.empty(nrows, dtype=float)
            rad2_Flagavg = np.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            i = 0
            for row in reader:
                rad1_time[i] = datetime.datetime.strptime(
                    row['rad1_time'], '%Y%m%d%H%M%S')
                rad1_ray_ind[i] = int(row['rad1_ray_ind'])
                rad1_rng_ind[i] = int(row['rad1_rng_ind'])
                rad1_ele[i] = float(row['rad1_ele'])
                rad1_azi[i] = float(row['rad1_azi'])
                rad1_rng[i] = float(row['rad1_rng'])
                rad1_dBZavg[i] = float(row['rad1_dBZavg'])
                rad1_PhiDPavg[i] = float(row['rad1_PhiDPavg'])
                rad1_Flagavg[i] = float(row['rad1_Flagavg'])
                rad2_time[i] = datetime.datetime.strptime(
                    row['rad2_time'], '%Y%m%d%H%M%S')
                rad2_ray_ind[i] = int(row['rad2_ray_ind'])
                rad2_rng_ind[i] = int(row['rad2_rng_ind'])
                rad2_ele[i] = float(row['rad2_ele'])
                rad2_azi[i] = float(row['rad2_azi'])
                rad2_rng[i] = float(row['rad2_rng'])
                rad2_dBZavg[i] = float(row['rad2_dBZavg'])
                rad2_PhiDPavg[i] = float(row['rad2_PhiDPavg'])
                rad2_Flagavg[i] = float(row['rad2_Flagavg'])
                i += 1

            csvfile.close()

            return (rad1_time, rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi,
                    rad1_rng, rad1_dBZavg, rad1_PhiDPavg, rad1_Flagavg,
                    rad2_time, rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi,
                    rad2_rng, rad2_dBZavg, rad2_PhiDPavg, rad2_Flagavg)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None, None, None, None, None, None, None, None)


def read_timeseries(fname):
    """
    Reads a time series contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date , value : tupple
        A datetime object array containing the time and a numpy masked array
        containing the value. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)
            value = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            i = 0
            date = list()
            for row in reader:
                date.append(datetime.datetime.strptime(
                    row['date'], '%Y-%m-%d %H:%M:%S.%f'))
                value[i] = float(row['value'])
                i += 1

            value = np.ma.masked_values(value, get_fillvalue())

            csvfile.close()

            return date, value
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None


def read_ts_cum(fname):
    """
    Reads a time series of precipitation accumulation contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date, np_radar, radar_value, np_sensor, sensor_value : tupple
        The data read

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)
            np_radar = np.zeros(nrows, dtype=int)
            radar_value = np.ma.empty(nrows, dtype=float)
            np_sensor = np.zeros(nrows, dtype=int)
            sensor_value = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            i = 0
            date = list()
            for row in reader:
                date.append(datetime.datetime.strptime(
                    row['date'], '%Y-%m-%d %H:%M:%S'))
                np_radar[i] = int(row['np_radar'])
                radar_value[i] = float(row['radar_value'])
                np_sensor[i] = int(row['np_sensor'])
                sensor_value[i] = float(row['sensor_value'])
                i += 1

            radar_value = np.ma.masked_values(radar_value, get_fillvalue())
            sensor_value = np.ma.masked_values(sensor_value, get_fillvalue())

            csvfile.close()

            return date, np_radar, radar_value, np_sensor, sensor_value
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None


def read_monitoring_ts(fname, sort_by_date=False):
    """
    Reads a monitoring time series contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file
    sort_by_date : bool
        if True, the read data is sorted by date prior to exit

    Returns
    -------
    date , np_t, central_quantile, low_quantile, high_quantile : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)

            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)
            np_t = np.zeros(nrows, dtype=int)
            central_quantile = np.ma.empty(nrows, dtype=float)
            low_quantile = np.ma.empty(nrows, dtype=float)
            high_quantile = np.ma.empty(nrows, dtype=float)
            date = np.empty(nrows, dtype=datetime.datetime)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#')
                )
            i = 0
            for row in reader:
                date[i] = datetime.datetime.strptime(
                    row['date'], '%Y%m%d%H%M%S')
                np_t[i] = int(row['NP'])
                central_quantile[i] = float(row['central_quantile'])
                low_quantile[i] = float(row['low_quantile'])
                high_quantile[i] = float(row['high_quantile'])
                i += 1

            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()

            central_quantile = np.ma.masked_values(
                central_quantile, get_fillvalue())
            low_quantile = np.ma.masked_values(
                low_quantile, get_fillvalue())
            high_quantile = np.ma.masked_values(
                high_quantile, get_fillvalue())

            if sort_by_date:
                ind = np.argsort(date)
                date = date[ind]
                np_t = np_t[ind]
                central_quantile = central_quantile[ind]
                low_quantile = low_quantile[ind]
                high_quantile = high_quantile[ind]

            return date, np_t, central_quantile, low_quantile, high_quantile
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None


def read_monitoring_ts_old(fname):
    """
    Reads an old format of the monitoring time series contained in a text file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date , np_t, central_quantile, low_quantile, high_quantile : tupple
        The read data in the current format. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith('#')),
                fieldnames=['Date [YYJJJ]', 'Npoints', 'Value'])
            nrows = sum(1 for row in reader)
            np_t = np.zeros(nrows, dtype=int)
            central_quantile = np.ma.empty(nrows, dtype=float)
            low_quantile = np.ma.empty(nrows, dtype=float)
            low_quantile[:] = np.ma.masked
            high_quantile = np.ma.empty(nrows, dtype=float)
            high_quantile[:] = np.ma.masked

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith('#')),
                fieldnames=['Date [YYJJJ]', 'Npoints', 'Value']
                )
            i = 0
            date = list()
            for row in reader:
                date.append(datetime.datetime.strptime(
                    row['Date [YYJJJ]'], '%y%j'))
                np_t[i] = int(row['Npoints'])
                central_quantile[i] = float(row['Value'])
                i += 1

            central_quantile = np.ma.masked_invalid(central_quantile)

            csvfile.close()

            return date, np_t, central_quantile, low_quantile, high_quantile
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None


def read_intercomp_scores_ts(fname, sort_by_date=False):
    """
    Reads a radar intercomparison scores csv file

    Parameters
    ----------
    fname : str
        path of time series file
    sort_by_date : bool
        if True, the read data is sorted by date prior to exit

    Returns
    -------
    date_vec, np_vec, meanbias_vec, medianbias_vec, quant25bias_vec,
    quant75bias_vec, modebias_vec, corr_vec, slope_vec, intercep_vec,
    intercep_slope1_vec : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            while True:
                try:
                    fcntl.flock(csvfile, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except OSError as e:
                    if e.errno != errno.EAGAIN:
                        raise
                    else:
                        time.sleep(0.1)
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)

            np_vec = np.zeros(nrows, dtype=int)
            meanbias_vec = np.ma.empty(nrows, dtype=float)
            medianbias_vec = np.ma.empty(nrows, dtype=float)
            quant25bias_vec = np.ma.empty(nrows, dtype=float)
            quant75bias_vec = np.ma.empty(nrows, dtype=float)
            modebias_vec = np.ma.empty(nrows, dtype=float)
            corr_vec = np.ma.empty(nrows, dtype=float)
            slope_vec = np.ma.empty(nrows, dtype=float)
            intercep_vec = np.ma.empty(nrows, dtype=float)
            intercep_slope1_vec = np.ma.empty(nrows, dtype=float)
            date_vec = np.empty(nrows, dtype=datetime.datetime)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#')
                )
            i = 0
            for row in reader:
                date_vec[i] = datetime.datetime.strptime(
                    row['date'], '%Y%m%d%H%M%S')
                np_vec[i] = int(row['NP'])
                meanbias_vec[i] = float(row['mean_bias'])
                medianbias_vec[i] = float(row['median_bias'])
                quant25bias_vec[i] = float(row['quant25_bias'])
                quant75bias_vec[i] = float(row['quant75_bias'])
                modebias_vec[i] = float(row['mode_bias'])
                corr_vec[i] = float(row['corr'])
                slope_vec[i] = float(row['slope_of_linear_regression'])
                intercep_vec[i] = float(row['intercep_of_linear_regression'])
                intercep_slope1_vec[i] = float(
                    row['intercep_of_linear_regression_of_slope_1'])
                i += 1

            fcntl.flock(csvfile, fcntl.LOCK_UN)
            csvfile.close()

            meanbias_vec = np.ma.masked_values(
                meanbias_vec, get_fillvalue())
            medianbias_vec = np.ma.masked_values(
                medianbias_vec, get_fillvalue())
            quant25bias_vec = np.ma.masked_values(
                quant25bias_vec, get_fillvalue())
            quant75bias_vec = np.ma.masked_values(
                quant75bias_vec, get_fillvalue())
            modebias_vec = np.ma.masked_values(
                modebias_vec, get_fillvalue())
            corr_vec = np.ma.masked_values(
                corr_vec, get_fillvalue())
            slope_vec = np.ma.masked_values(
                slope_vec, get_fillvalue())
            intercep_vec = np.ma.masked_values(
                intercep_vec, get_fillvalue())
            intercep_slope1_vec = np.ma.masked_values(
                intercep_slope1_vec, get_fillvalue())

            if sort_by_date:
                ind = np.argsort(date_vec)
                date_vec = date_vec[ind]
                np_vec = np_vec[ind]
                meanbias_vec = meanbias_vec[ind]
                medianbias_vec = medianbias_vec[ind]
                quant25bias_vec = quant25bias_vec[ind]
                quant75bias_vec = quant75bias_vec[ind]
                modebias_vec = modebias_vec[ind]
                corr_vec = corr_vec[ind]
                slope_vec = slope_vec[ind]
                intercep_vec = intercep_vec[ind]
                intercep_slope1_vec = intercep_slope1_vec[ind]

            return (date_vec, np_vec, meanbias_vec, medianbias_vec,
                    quant25bias_vec, quant75bias_vec, modebias_vec, corr_vec,
                    slope_vec, intercep_vec, intercep_slope1_vec)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None)


def read_intercomp_scores_ts_old(fname):
    """
    Reads a radar intercomparison scores csv file in old format

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date_vec, np_vec, meanbias_vec, medianbias_vec, quant25bias_vec,
    quant75bias_vec, modebias_vec, corr_vec, slope_vec, intercep_vec,
    intercep_slope1_vec : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith('#')),
                fieldnames=['date', 'NP', 'mode_bias', 'median_bias',
                            'mean_bias', 'corr', 'slope_of_linear_regression',
                            'intercep_of_linear_regression'],
                dialect='excel-tab')
            nrows = sum(1 for row in reader)

            np_vec = np.zeros(nrows, dtype=int)
            meanbias_vec = np.ma.empty(nrows, dtype=float)
            medianbias_vec = np.ma.empty(nrows, dtype=float)
            quant25bias_vec = np.ma.empty(nrows, dtype=float)
            quant25bias_vec[:] = np.ma.masked
            quant75bias_vec = np.ma.empty(nrows, dtype=float)
            quant75bias_vec[:] = np.ma.masked
            modebias_vec = np.ma.empty(nrows, dtype=float)
            corr_vec = np.ma.empty(nrows, dtype=float)
            slope_vec = np.ma.empty(nrows, dtype=float)
            intercep_vec = np.ma.empty(nrows, dtype=float)
            intercep_slope1_vec = np.ma.empty(nrows, dtype=float)
            intercep_slope1_vec[:] = np.ma.masked

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith('#')),
                fieldnames=['date', 'NP', 'mode_bias', 'median_bias',
                            'mean_bias', 'corr', 'slope_of_linear_regression',
                            'intercep_of_linear_regression'],
                dialect='excel-tab')
            i = 0
            date_vec = list()
            for row in reader:
                date_vec.append(datetime.datetime.strptime(
                    row['date'], '%y%j'))
                np_vec[i] = int(row['NP'])
                meanbias_vec[i] = float(row['mean_bias'])
                medianbias_vec[i] = float(row['median_bias'])
                modebias_vec[i] = float(row['mode_bias'])
                corr_vec[i] = float(row['corr'])
                slope_vec[i] = float(row['slope_of_linear_regression'])
                intercep_vec[i] = float(row['intercep_of_linear_regression'])
                i += 1

            meanbias_vec = np.ma.masked_invalid(meanbias_vec)
            medianbias_vec = np.ma.masked_invalid(medianbias_vec)
            modebias_vec = np.ma.masked_invalid(modebias_vec)
            corr_vec = np.ma.masked_invalid(corr_vec)
            slope_vec = np.ma.masked_invalid(slope_vec)
            intercep_vec = np.ma.masked_invalid(intercep_vec)

            csvfile.close()

            return (date_vec, np_vec, meanbias_vec, medianbias_vec,
                    quant25bias_vec, quant75bias_vec, modebias_vec, corr_vec,
                    slope_vec, intercep_vec, intercep_slope1_vec)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None)


def read_intercomp_scores_ts_old_v0(fname, corr_min=0.6, np_min=9):
    """
    Reads a radar intercomparison scores csv file in the oldest format

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date_vec, np_vec, meanbias_vec, medianbias_vec, quant25bias_vec,
    quant75bias_vec, modebias_vec, corr_vec, slope_vec, intercep_vec,
    intercep_slope1_vec : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith('#')),
                fieldnames=['date', 'ele', 'NP', 'mean_bias', 'corr'],
                dialect='excel-tab')
            nrows = sum(1 for row in reader)

            np_aux_vec = np.zeros(nrows, dtype=int)
            meanbias_aux_vec = np.ma.empty(nrows, dtype=float)
            corr_aux_vec = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith('#')),
                fieldnames=['date', 'ele', 'NP', 'mean_bias', 'corr'],
                dialect='excel-tab')
            i = 0
            date_aux_vec = list()
            for row in reader:
                date_aux_vec.append(datetime.datetime.strptime(
                    row['date'], '%y%j'))
                np_aux_vec[i] = int(row['NP'])
                meanbias_aux_vec[i] = float(row['mean_bias'])
                corr_aux_vec[i] = float(row['corr'])
                i += 1

            meanbias_aux_vec = np.ma.masked_invalid(meanbias_aux_vec)
            corr_aux_vec = np.ma.masked_invalid(corr_aux_vec)

            csvfile.close()

            date_aux_vec = np.asarray(date_aux_vec)
            date_vec, unique_ind = np.unique(date_aux_vec, return_index=True)
            nelements = len(date_vec)

            np_vec = np.zeros(nelements, dtype=int)
            meanbias_vec = np.ma.empty(nelements, dtype=float)
            corr_vec = np.ma.empty(nelements, dtype=float)
            for i in range(nelements-1):
                ind_aux = np.arange(unique_ind[i], unique_ind[i+1])
                np_aux = np_aux_vec[ind_aux]
                meanbias_aux = meanbias_aux_vec[ind_aux]
                corr_aux = corr_aux_vec[ind_aux]

                ind_valid = np.where(np.logical_and(
                    corr_aux > corr_min, np_aux > np_min))
                np_aux = np_aux[ind_valid]
                meanbias_aux = meanbias_aux[ind_valid]
                corr_aux = corr_aux[ind_valid]

                np_vec[i] = np.sum(np_aux, dtype=int)
                if np_vec[i] == 0:
                    continue
                meanbias_vec[i] = np.sum(np_aux*meanbias_aux)/np_vec[i]
                corr_vec[i] = np.sum(np_aux*corr_aux)/np_vec[i]

            # last date
            ind_aux = np.arange(unique_ind[-1], len(date_aux_vec))
            np_aux = np_aux_vec[ind_aux]
            meanbias_aux = meanbias_aux_vec[ind_aux]
            corr_aux = corr_aux_vec[ind_aux]

            ind_valid = np.where(np.logical_and(
                corr_aux > corr_min, np_aux > np_min))
            np_aux = np_aux[ind_valid]
            meanbias_aux = meanbias_aux[ind_valid]
            corr_aux = corr_aux[ind_valid]
            np_vec[-1] = np.sum(np_aux, dtype=int)
            if np_vec[-1] > 0:
                meanbias_vec[-1] = np.sum(np_aux*meanbias_aux)/np_vec[-1]
                corr_vec[-1] = np.sum(np_aux*corr_aux)/np_vec[-1]

            medianbias_vec = np.ma.empty(nelements, dtype=float)
            medianbias_vec[:] = np.ma.masked
            quant25bias_vec = np.ma.empty(nelements, dtype=float)
            quant25bias_vec[:] = np.ma.masked
            quant75bias_vec = np.ma.empty(nelements, dtype=float)
            quant75bias_vec[:] = np.ma.masked
            modebias_vec = np.ma.empty(nelements, dtype=float)
            modebias_vec[:] = np.ma.masked
            slope_vec = np.ma.empty(nelements, dtype=float)
            slope_vec[:] = np.ma.masked
            intercep_vec = np.ma.empty(nelements, dtype=float)
            intercep_vec[:] = np.ma.masked
            intercep_slope1_vec = np.ma.empty(nelements, dtype=float)
            intercep_slope1_vec[:] = np.ma.masked

            return (date_vec, np_vec, meanbias_vec, medianbias_vec,
                    quant25bias_vec, quant75bias_vec, modebias_vec, corr_vec,
                    slope_vec, intercep_vec, intercep_slope1_vec)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None)


def read_selfconsistency(fname):
    """
    Reads a self-consistency table with Zdr, Kdp/Zh columns

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    zdr, kdpzh : arrays
        The read values

    """
    try:
        with open(fname, 'r', newline='', encoding='utf-8', errors='ignore') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                csvfile, fieldnames=['zdr', 'kdpzh'], dialect='excel-tab')
            nrows = sum(1 for row in reader)

            zdr_kdpzh_table = np.empty((2, nrows), dtype='float32')

            # now read the data
            csvfile.seek(0)

            reader = csv.DictReader(
                csvfile, fieldnames=['zdr', 'kdpzh'], dialect='excel-tab')
            i = 0

            for row in reader:
                zdr_kdpzh_table[0, i] = float(row['zdr'])
                zdr_kdpzh_table[1, i] = float(row['kdpzh'])
                i += 1

            csvfile.close()

            return zdr_kdpzh_table
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None


def read_antenna_pattern(fname, linear=False, twoway=False):
    """
    Read antenna pattern from file

    Parameters
    ----------
    fname : str
        path of the antenna pattern file
    linear : boolean
        if true the antenna pattern is given in linear units
    twoway : boolean
        if true the attenuation is two-way

    Returns
    -------
    pattern : dict
        dictionary with the fields angle and attenuation

    """

    try:
        pfile = open(fname, "r")
    except EnvironmentError as ee:
        warn(str(ee))
        raise Exception("ERROR: Could not find|open antenna file '" +
                        fname+"'")

    # Get array length
    linenum = 0
    for line in pfile:
        if line.startswith("# CRC:"):
            break
        linenum += 1

    pattern = {
        'angle': np.empty(linenum),
        'attenuation': np.empty(linenum)
    }

    pfile.seek(0)
    cnt = 0
    for line in pfile:
        if line.startswith("# CRC:"):
            break
        line = line.strip()
        fields = line.split(',')
        pattern['angle'][cnt] = float(fields[0])
        pattern['attenuation'][cnt] = float(fields[1])
        cnt += 1

    pfile.close()

    if twoway:
        pattern['attenuation'] = 2. * pattern['attenuation']

    if linear:
        pattern['attenuation'] = 10.**(pattern['attenuation']/10.)

    return pattern
