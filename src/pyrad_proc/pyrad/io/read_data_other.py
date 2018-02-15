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
    read_colocated_gates
    read_colocated_data
    read_colocated_data_time_avg
    read_timeseries
    read_lightning
    read_ts_cum
    read_monitoring_ts
    read_intercomp_scores_ts
    read_sun_hits_multiple_days
    read_sun_hits
    read_sun_retrieval
    read_solar_flux
    get_sensor_data
    read_smn
    read_smn2
    read_disdro_scattering
    read_selfconsistency
    read_antenna_pattern
"""

import os
import glob
import datetime
import csv
import xml.etree.ElementTree as et
from warnings import warn
from copy import deepcopy

import numpy as np

from pyart.config import get_fillvalue, get_metadata

from .io_aux import get_save_dir, make_filename
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
            if len(line) == 0:
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
    if len(filename) == 0:
        warn('rad4alp status file '+datapath+basename+timeinfo +
             '*.xml not found')
        return None
    root = et.parse(filename[0]).getroot()

    return root


def read_rad4alp_cosmo(fname, datatype):
    """
    Reads rad4alp COSMO data binary file.

    Parameters
    ----------
    fname : str
        name of the file to read

    datatype : str
        name of the data type

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
            return field
        elif datatype == 'ISO0':
            field_name = get_fieldname_pyart(datatype)
            field_data = np.ma.masked_where(mask, (bindata-1).astype(float))

            field = get_metadata(field_name)
            field['data'] = field_data
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
    rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng, rad1_val,
    rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng, rad2_val :
        tupple
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
            rad1_val = np.empty(nrows, dtype=float)
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
                rad1_ray_ind[i] = int(row['rad1_ray_ind'])
                rad1_rng_ind[i] = int(row['rad1_rng_ind'])
                rad1_ele[i] = float(row['rad1_ele'])
                rad1_azi[i] = float(row['rad1_azi'])
                rad1_rng[i] = float(row['rad1_rng'])
                rad1_val[i] = float(row['rad1_val'])
                rad2_ray_ind[i] = int(row['rad2_ray_ind'])
                rad2_rng_ind[i] = int(row['rad2_rng_ind'])
                rad2_ele[i] = float(row['rad2_ele'])
                rad2_azi[i] = float(row['rad2_azi'])
                rad2_rng[i] = float(row['rad2_rng'])
                rad2_val[i] = float(row['rad2_val'])
                i += 1

            csvfile.close()

            return (rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
                    rad1_val, rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi,
                    rad2_rng, rad2_val)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None, None)


def read_colocated_data_time_avg(fname):
    """
    Reads a csv files containing time averaged colocated data

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    rad1_ray_ind, rad1_rng_ind, rad1_ele , rad1_azi, rad1_rng, rad1_val,
    rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng, rad2_val :
        tupple
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
            rad1_dBZavg = np.empty(nrows, dtype=float)
            rad1_PhiDPavg = np.empty(nrows, dtype=float)
            rad1_Flagavg = np.empty(nrows, dtype=float)
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
                rad1_ray_ind[i] = int(row['rad1_ray_ind'])
                rad1_rng_ind[i] = int(row['rad1_rng_ind'])
                rad1_ele[i] = float(row['rad1_ele'])
                rad1_azi[i] = float(row['rad1_azi'])
                rad1_rng[i] = float(row['rad1_rng'])
                rad1_dBZavg[i] = float(row['rad1_dBZavg'])
                rad1_PhiDPavg[i] = float(row['rad1_PhiDPavg'])
                rad1_Flagavg[i] = float(row['rad1_Flagavg'])
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

            return (rad1_ray_ind, rad1_rng_ind, rad1_ele, rad1_azi, rad1_rng,
                    rad1_dBZavg, rad1_PhiDPavg, rad1_Flagavg,
                    rad2_ray_ind, rad2_rng_ind, rad2_ele, rad2_azi, rad2_rng,
                    rad2_dBZavg, rad2_PhiDPavg, rad2_Flagavg)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None, None, None, None, None, None)


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


def read_lightning(fname, filter_data=True):
    """
    Reads lightning data contained in a text file. The file has the following
    fields:
        flashnr: (0 is for noise)
        UTC seconds of the day
        Time within flash (in seconds)
        Latitude (decimal degrees)
        Longitude (decimal degrees)
        Altitude (m MSL)
        Power (dBm)

    Parameters
    ----------
    fname : str
        path of time series file
    filter_data : Boolean
        if True filter noise (flashnr = 0)

    Returns
    -------
    flashnr, time, time_in_flash, lat, lon, alt, dBm : tupple
        A tupple containing the read values. None otherwise

    """
    try:
        # get date from file name
        bfile = os.path.basename(fname)
        datetimestr = bfile[0:6]
        fdatetime = datetime.datetime.strptime(datetimestr, '%y%m%d')

        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                csvfile, fieldnames=['flashnr', 'time', 'time_in_flash',
                                     'lat', 'lon', 'alt', 'dBm'],
                delimiter=' ')
            nrows = sum(1 for row in reader)

            flashnr = np.ma.empty(nrows, dtype=int)
            time_in_flash = np.ma.empty(nrows, dtype=float)
            lat = np.ma.empty(nrows, dtype=float)
            lon = np.ma.empty(nrows, dtype=float)
            alt = np.ma.empty(nrows, dtype=float)
            dBm = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                csvfile, fieldnames=['flashnr', 'time', 'time_in_flash',
                                     'lat', 'lon', 'alt', 'dBm'],
                delimiter=' ')
            i = 0
            time = list()
            for row in reader:
                flashnr[i] = int(row['flashnr'])
                time.append(fdatetime+datetime.timedelta(
                    seconds=float(row['time'])))
                time_in_flash[i] = float(row['time_in_flash'])
                lat[i] = float(row['lat'])
                lon[i] = float(row['lon'])
                alt[i] = float(row['alt'])
                dBm[i] = float(row['dBm'])

                i += 1

            time = np.array(time)

            if filter_data:
                flashnr_aux = deepcopy(flashnr)
                flashnr = flashnr[flashnr_aux > 0]
                time = time[flashnr_aux > 0]
                time_in_flash = time_in_flash[flashnr_aux > 0]
                lat = lat[flashnr_aux > 0]
                lon = lon[flashnr_aux > 0]
                alt = alt[flashnr_aux > 0]
                dBm = dBm[flashnr_aux > 0]

            csvfile.close()

            return flashnr, time, time_in_flash, lat, lon, alt, dBm
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None, None, None


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


def read_monitoring_ts(fname):
    """
    Reads a monitoring time series contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date , np_t, central_quantile, low_quantile, high_quantile : tupple
        The read data. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)
            np_t = np.zeros(nrows, dtype=int)
            central_quantile = np.ma.empty(nrows, dtype=float)
            low_quantile = np.ma.empty(nrows, dtype=float)
            high_quantile = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#')
                )
            i = 0
            date = list()
            for row in reader:
                date.append(datetime.datetime.strptime(
                    row['date'], '%Y%m%d%H%M%S'))
                np_t[i] = int(row['NP'])
                central_quantile[i] = float(row['central_quantile'])
                low_quantile[i] = float(row['low_quantile'])
                high_quantile[i] = float(row['high_quantile'])
                i += 1

            central_quantile = np.ma.masked_values(
                central_quantile, get_fillvalue())
            low_quantile = np.ma.masked_values(
                low_quantile, get_fillvalue())
            high_quantile = np.ma.masked_values(
                high_quantile, get_fillvalue())

            csvfile.close()

            return date, np_t, central_quantile, low_quantile, high_quantile
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None


def read_intercomp_scores_ts(fname):
    """
    Reads a radar intercomparison scores csv file

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

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#')
                )
            i = 0
            date_vec = list()
            for row in reader:
                date_vec.append(datetime.datetime.strptime(
                    row['date'], '%Y%m%d%H%M%S'))
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

            csvfile.close()

            return (date_vec, np_vec, meanbias_vec, medianbias_vec,
                    quant25bias_vec, quant75bias_vec, modebias_vec, corr_vec,
                    slope_vec, intercep_vec, intercep_slope1_vec)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None)


def read_sun_hits_multiple_days(cfg, time_ref, nfiles=1):
    """
    Reads sun hits data from multiple file sources

    Parameters
    ----------
    cfg : dict
        dictionary with configuration data to find out the right file
    time_ref : datetime object
        reference time
    nfiles : int
        number of files to read

    Returns
    -------
    date, ray, nrng, rad_el, rad_az, sun_el, sun_az, ph, ph_std, nph, nvalh,
    pv, pv_std, npv, nvalv, zdr, zdr_std, nzdr, nvalzdr : tupple
        Each parameter is an array containing a time series of information on
        a variable

    """
    timeinfo = time_ref - datetime.timedelta(days=nfiles-1)

    date = list()
    ray = list()
    nrng = list()
    rad_el = list()
    rad_az = list()
    sun_el = list()
    sun_az = list()
    ph = np.ma.array(list())
    ph_std = np.ma.array(list())
    nph = list()
    nvalh = list()
    pv = np.ma.array(list())
    pv_std = np.ma.array(list())
    npv = list()
    nvalv = list()
    zdr = np.ma.array(list())
    zdr_std = np.ma.array(list())
    nzdr = list()
    nvalzdr = list()

    for i in range(nfiles):
        savedir = get_save_dir(
            cfg['basepath'], cfg['procname'], cfg['dsname'],
            cfg['sun_hits_dir'], timeinfo=timeinfo, create_dir=False)

        fname = make_filename(
            'info', cfg['type'], 'detected', ['csv'],
            timeinfo=timeinfo, timeformat='%Y%m%d')

        (date_aux, ray_aux, nrng_aux, rad_el_aux, rad_az_aux, sun_el_aux,
         sun_az_aux, ph_aux, ph_std_aux, nph_aux, nvalh_aux, pv_aux,
         pv_std_aux, npv_aux, nvalv_aux, zdr_aux, zdr_std_aux, nzdr_aux,
         nvalzdr_aux) = read_sun_hits(savedir+fname[0])

        if date_aux is None:
            return (None, None, None, None, None, None, None, None, None,
                    None, None, None, None, None, None, None, None, None,
                    None)

        date += date_aux
        ray += list(ray_aux)
        nrng += list(nrng_aux)
        rad_el += list(rad_el_aux)
        rad_az += list(rad_az_aux)
        sun_el += list(sun_el_aux)
        sun_az += list(sun_az_aux)
        ph = np.ma.append(ph, ph_aux)
        ph_std = np.ma.append(ph_std, ph_std_aux)
        nph += list(nph_aux)
        nvalh += list(nvalh_aux)
        pv = np.ma.append(pv, pv_aux)
        pv_std = np.ma.append(pv_std, pv_std_aux)
        npv += list(npv_aux)
        nvalv += list(nvalv_aux)
        zdr = np.ma.append(zdr, zdr_aux)
        zdr_std = np.ma.append(zdr_std, zdr_std_aux)
        nzdr += list(nzdr_aux)
        nvalzdr += list(nvalzdr_aux)

        timeinfo += datetime.timedelta(days=1)

    ray = np.squeeze(np.array(ray))
    nrng = np.squeeze(np.array(nrng))
    rad_el = np.squeeze(np.array(rad_el))
    rad_az = np.squeeze(np.array(rad_az))
    sun_el = np.squeeze(np.array(sun_el))
    sun_az = np.squeeze(np.array(sun_az))
    nph = np.squeeze(np.array(nph))
    nvalh = np.squeeze(np.array(nvalh))
    npv = np.squeeze(np.array(npv))
    nvalv = np.squeeze(np.array(nvalv))
    nzdr = np.squeeze(np.array(nzdr))
    nvalzdr = np.squeeze(np.array(nvalzdr))

    return (date, ray, nrng, rad_el, rad_az, sun_el, sun_az, ph, ph_std, nph,
            nvalh, pv, pv_std, npv, nvalv, zdr, zdr_std, nzdr, nvalzdr)


def read_sun_hits(fname):
    """
    Reads sun hits data contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date, ray, nrng, rad_el, rad_az, sun_el, sun_az, ph, ph_std, nph, nvalh,
    pv, pv_std, npv, nvalv, zdr, zdr_std, nzdr, nvalzdr : tupple
        Each parameter is an array containing a time series of information on
        a variable

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)

            ray = np.empty(nrows, dtype=int)
            nrng = np.empty(nrows, dtype=int)
            rad_el = np.empty(nrows, dtype=float)
            rad_az = np.empty(nrows, dtype=float)
            sun_el = np.empty(nrows, dtype=float)
            sun_az = np.empty(nrows, dtype=float)
            ph = np.ma.empty(nrows, dtype=float)
            ph_std = np.ma.empty(nrows, dtype=float)
            nph = np.empty(nrows, dtype=int)
            nvalh = np.empty(nrows, dtype=int)
            pv = np.ma.empty(nrows, dtype=float)
            pv_std = np.ma.empty(nrows, dtype=float)
            npv = np.empty(nrows, dtype=int)
            nvalv = np.empty(nrows, dtype=int)
            zdr = np.ma.empty(nrows, dtype=float)
            zdr_std = np.ma.empty(nrows, dtype=float)
            nzdr = np.empty(nrows, dtype=int)
            nvalzdr = np.empty(nrows, dtype=int)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            i = 0
            date = list()
            for row in reader:
                date.append(datetime.datetime.strptime(
                    row['time'], '%Y-%m-%d %H:%M:%S.%f'))
                ray[i] = int(row['ray'])
                nrng[i] = int(row['NPrng'])
                rad_el[i] = float(row['rad_el'])
                rad_az[i] = float(row['rad_az'])
                sun_el[i] = float(row['sun_el'])
                sun_az[i] = float(row['sun_az'])
                ph[i] = float(row['dBm_sun_hit'])
                ph_std[i] = float(row['std(dBm_sun_hit)'])
                nph[i] = int(row['NPh'])
                nvalh[i] = int(row['NPhval'])
                pv[i] = float(row['dBmv_sun_hit'])
                pv_std[i] = float(row['std(dBmv_sun_hit)'])
                npv[i] = int(row['NPv'])
                nvalv[i] = int(row['NPvval'])
                zdr[i] = float(row['ZDR_sun_hit'])
                zdr_std[i] = float(row['std(ZDR_sun_hit)'])
                nzdr[i] = int(row['NPzdr'])
                nvalzdr[i] = int(row['NPzdrval'])

                i += 1

            ph = np.ma.masked_values(ph, get_fillvalue())
            ph_std = np.ma.masked_values(ph_std, get_fillvalue())
            pv = np.ma.masked_values(pv, get_fillvalue())
            pv_std = np.ma.masked_values(pv_std, get_fillvalue())
            zdr = np.ma.masked_values(zdr, get_fillvalue())
            zdr_std = np.ma.masked_values(zdr_std, get_fillvalue())

            csvfile.close()

            return (date, ray, nrng, rad_el, rad_az, sun_el, sun_az,
                    ph, ph_std, nph, nvalh, pv, pv_std, npv, nvalv,
                    zdr, zdr_std, nzdr, nvalzdr)

    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None, None, None, None, None, None, None, None, None)


def read_sun_retrieval(fname):
    """
    Reads sun retrieval data contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    first_hit_time, last_hit_time, nhits_h, el_width_h, az_width_h, el_bias_h,
    az_bias_h, dBm_sun_est, std_dBm_sun_est, sf_h,
    nhits_v, el_width_v, az_width_v, el_bias_v, az_bias_v, dBmv_sun_est,
    std_dBmv_sun_est, sf_v,
    nhits_zdr, zdr_sun_est, std_zdr_sun_est,
    sf_ref, ref_time : tupple
        Each parameter is an array containing a time series of information on
        a variable

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))
            nrows = sum(1 for row in reader)

            nhits_h = np.empty(nrows, dtype=int)
            el_width_h = np.ma.empty(nrows, dtype=float)
            az_width_h = np.ma.empty(nrows, dtype=float)
            el_bias_h = np.ma.empty(nrows, dtype=float)
            az_bias_h = np.ma.empty(nrows, dtype=float)
            dBm_sun_est = np.ma.empty(nrows, dtype=float)
            std_dBm_sun_est = np.ma.empty(nrows, dtype=float)
            sf_h = np.ma.empty(nrows, dtype=float)

            nhits_v = np.empty(nrows, dtype=int)
            el_width_v = np.ma.empty(nrows, dtype=float)
            az_width_v = np.ma.empty(nrows, dtype=float)
            el_bias_v = np.ma.empty(nrows, dtype=float)
            az_bias_v = np.ma.empty(nrows, dtype=float)
            dBmv_sun_est = np.ma.empty(nrows, dtype=float)
            std_dBmv_sun_est = np.ma.empty(nrows, dtype=float)
            sf_v = np.ma.empty(nrows, dtype=float)

            nhits_zdr = np.empty(nrows, dtype=int)
            zdr_sun_est = np.ma.empty(nrows, dtype=float)
            std_zdr_sun_est = np.ma.empty(nrows, dtype=float)

            sf_ref = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                row for row in csvfile if not row.startswith('#'))

            i = 0
            first_hit_time = list()
            last_hit_time = list()
            ref_time = list()
            for row in reader:
                first_hit_time.append(datetime.datetime.strptime(
                    row['first_hit_time'], '%Y%m%d%H%M%S'))
                last_hit_time.append(datetime.datetime.strptime(
                    row['last_hit_time'], '%Y%m%d%H%M%S'))

                nhits_h[i] = int(row['nhits_h'])
                el_width_h[i] = float(row['el_width_h'])
                az_width_h[i] = float(row['az_width_h'])
                el_bias_h[i] = float(row['el_bias_h'])
                az_bias_h[i] = float(row['az_bias_h'])
                dBm_sun_est[i] = float(row['dBm_sun_est'])
                std_dBm_sun_est[i] = float(row['std(dBm_sun_est)'])
                sf_h[i] = float(row['sf_h'])

                nhits_v[i] = int(row['nhits_v'])
                el_width_v[i] = float(row['el_width_v'])
                az_width_v[i] = float(row['az_width_v'])
                el_bias_v[i] = float(row['el_bias_v'])
                az_bias_v[i] = float(row['az_bias_v'])
                dBmv_sun_est[i] = float(row['dBmv_sun_est'])
                std_dBmv_sun_est[i] = float(row['std(dBmv_sun_est)'])
                sf_v[i] = float(row['sf_v'])

                nhits_zdr[i] = int(row['nhits_zdr'])
                zdr_sun_est[i] = float(row['ZDR_sun_est'])
                std_zdr_sun_est[i] = float(row['std(ZDR_sun_est)'])

                sf_ref[i] = float(row['sf_ref'])
                if row['ref_time'] == 'None':
                    ref_time.append(None)
                else:
                    ref_time.append(datetime.datetime.strptime(
                        row['ref_time'], '%Y%m%d%H%M%S'))

                i += 1

            el_width_h = np.ma.masked_values(el_width_h, get_fillvalue())
            az_width_h = np.ma.masked_values(az_width_h, get_fillvalue())
            el_bias_h = np.ma.masked_values(el_bias_h, get_fillvalue())
            az_bias_h = np.ma.masked_values(az_bias_h, get_fillvalue())
            dBm_sun_est = np.ma.masked_values(dBm_sun_est, get_fillvalue())
            std_dBm_sun_est = np.ma.masked_values(
                std_dBm_sun_est, get_fillvalue())
            sf_h = np.ma.masked_values(sf_h, get_fillvalue())

            el_width_v = np.ma.masked_values(el_width_v, get_fillvalue())
            az_width_v = np.ma.masked_values(az_width_v, get_fillvalue())
            el_bias_v = np.ma.masked_values(el_bias_v, get_fillvalue())
            az_bias_v = np.ma.masked_values(az_bias_v, get_fillvalue())
            dBmv_sun_est = np.ma.masked_values(dBmv_sun_est, get_fillvalue())
            std_dBmv_sun_est = np.ma.masked_values(
                std_dBmv_sun_est, get_fillvalue())
            sf_v = np.ma.masked_values(sf_v, get_fillvalue())

            zdr_sun_est = np.ma.masked_values(zdr_sun_est, get_fillvalue())
            std_zdr_sun_est = np.ma.masked_values(
                std_zdr_sun_est, get_fillvalue())

            sf_ref = np.ma.masked_values(sf_ref, get_fillvalue())

            csvfile.close()

            return (first_hit_time, last_hit_time, nhits_h,
                    el_width_h, az_width_h, el_bias_h, az_bias_h,
                    dBm_sun_est, std_dBm_sun_est, sf_h,
                    nhits_v, el_width_v, az_width_v, el_bias_v, az_bias_v,
                    dBmv_sun_est, std_dBmv_sun_est, sf_v,
                    nhits_zdr, zdr_sun_est, std_zdr_sun_est, sf_ref, ref_time)

    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None, None, None, None, None, None, None, None, None, None,
                None, None, None)


def read_solar_flux(fname):
    """
    Reads solar flux data from the DRAO observatory in Canada

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    flux_datetime : datetime array
        the date and time of the solar flux retrievals
    flux_value : array
        the observed solar flux

    """
    try:
        with open(fname, 'r', newline='') as txtfile:
            # skip the first two lines
            for i in range(2):
                next(txtfile)

            # first count the lines
            reader = csv.DictReader(
                txtfile, delimiter=' ', skipinitialspace=True, fieldnames=[
                    'fluxdate', 'fluxtime', 'fluxjulian', 'fluxcarrington',
                    'fluxobsflux', 'fluxadjflux', 'fluxursi'])
            nrows = sum(1 for row in reader)

            flux_value = np.empty(nrows, dtype=float)

            # now read the data
            txtfile.seek(0)

            # skip the first two lines
            for i in range(2):
                next(txtfile)

            reader = csv.DictReader(
                txtfile, delimiter=' ', skipinitialspace=True, fieldnames=[
                    'fluxdate', 'fluxtime', 'fluxjulian', 'fluxcarrington',
                    'fluxobsflux', 'fluxadjflux', 'fluxursi'])

            i = 0
            flux_datetime = list()
            for row in reader:
                flux_datetime.append(datetime.datetime.strptime(
                    row['fluxdate']+row['fluxtime'], '%Y%m%d%H%M%S'))
                flux_value[i] = float(row['fluxobsflux'])

                i += 1

            txtfile.close()

            return flux_datetime, flux_value

    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None


def get_sensor_data(date, datatype, cfg):
    """
    Gets data from a point measurement sensor (rain gauge or disdrometer)

    Parameters
    ----------
    date : datetime object
        measurement date

    datatype : str
        name of the data type to read

    cfg : dictionary
        dictionary containing sensor information

    Returns
    -------
    sensordate , sensorvalue, label, period : tupple
        date, value, type of sensor and measurement period

    """
    if cfg['sensor'] == 'rgage':
        datapath = cfg['smnpath']+date.strftime('%Y%m')+'/'
        datafile = date.strftime('%Y%m%d')+'_' + cfg['sensorid']+'.csv'
        (id, sensordate, pressure, temp,
         rh, sensorvalue, wspeed, wdir) = read_smn(datapath+datafile)
        if sensordate is None:
            return None, None, None, None
        label = 'RG'
        period = (sensordate[1]-sensordate[0]).total_seconds()
    elif cfg['sensor'] == 'disdro':
        datapath = cfg['disdropath']
        datafile = ('DSDfiltpolvar-'+cfg['sensorid']+'_' +
                    date.strftime('%Y%m%d')+'_Xband_temp' + cfg['temp'] +
                    '_elev'+cfg['elev']+'.txt')
        (sensordate, prectype, lwc, rr, zh, zv, zdr, ldr, ah, av,
         adiff, kdp, detaco, rhohv) = read_disdro_scattering(
            datapath+datafile)
        if sensordate is None:
            return None, None, None, None
        label = 'Disdro'
        period = (sensordate[1]-sensordate[0]).total_seconds()
        if datatype == 'RR':
            sensorvalue = rr
        elif (datatype == 'dBZ') or (datatype == 'dBZc'):
            sensorvalue = zh
    else:
        warn('Unknown sensor: '+cfg['sensor'])
        return None, None, None, None

    return sensordate, sensorvalue, label, period


def read_smn(fname):
    """
    Reads SwissMetNet data contained in a csv file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    id, date , pressure, temp, rh, precip, wspeed, wdir : tupple
        The read values

    """
    fill_value = 10000000.0
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(csvfile)
            nrows = sum(1 for row in reader)
            id = np.ma.empty(nrows, dtype='float32')
            pressure = np.ma.empty(nrows, dtype='float32')
            temp = np.ma.empty(nrows, dtype='float32')
            rh = np.ma.empty(nrows, dtype='float32')
            precip = np.ma.empty(nrows, dtype='float32')
            wspeed = np.ma.empty(nrows, dtype='float32')
            wdir = np.ma.empty(nrows, dtype='float32')

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(csvfile)
            i = 0
            date = list()
            for row in reader:
                id[i] = float(row['StationID'])
                date.append(datetime.datetime.strptime(
                    row['DateTime'], '%Y%m%d%H%M%S'))
                pressure[i] = float(row['AirPressure'])
                temp[i] = float(row['2mTemperature'])
                rh[i] = float(row['RH'])
                precip[i] = float(row['Precipitation'])
                wspeed[i] = float(row['Windspeed'])
                wdir[i] = float(row['Winddirection'])
                i += 1

            pressure = np.ma.masked_values(pressure, fill_value)
            temp = np.ma.masked_values(temp, fill_value)
            rh = np.ma.masked_values(rh, fill_value)
            precip = np.ma.masked_values(precip, fill_value)
            wspeed = np.ma.masked_values(wspeed, fill_value)
            wdir = np.ma.masked_values(wdir, fill_value)

            # convert precip from mm/10min to mm/h
            precip *= 6.

            csvfile.close()

            return id, date, pressure, temp, rh, precip, wspeed, wdir
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None, None, None, None


def read_smn2(fname):
    """
    Reads SwissMetNet data contained in a csv file with format
    station,time,value

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    id, date , value : tupple
        The read values

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # skip the first 2 lines
            next(csvfile)
            next(csvfile)

            # first count the lines
            reader = csv.DictReader(
                csvfile, fieldnames=['StationID', 'DateTime', 'Value'])
            nrows = sum(1 for row in reader)

            if nrows == 0:
                warn('Empty file '+fname)
                return None, None, None
            id = np.ma.empty(nrows, dtype=int)
            value = np.ma.empty(nrows, dtype=float)

            # now read the data
            csvfile.seek(0)

            # skip the first 2 lines
            next(csvfile)
            next(csvfile)

            reader = csv.DictReader(
                csvfile, fieldnames=['StationID', 'DateTime', 'Value'])
            i = 0
            date = list()
            for row in reader:
                id[i] = float(row['StationID'])
                date.append(datetime.datetime.strptime(
                    row['DateTime'], '%Y%m%d%H%M%S'))
                value[i] = float(row['Value'])
                i += 1

            csvfile.close()

            return id, date, value
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None


def read_disdro_scattering(fname):
    """
    Reads scattering parameters computed from disdrometer data contained in a
    text file

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    date, preciptype, lwc, rr, zh, zv, zdr, ldr, ah, av, adiff, kdp, deltaco,
    rhohv : tupple
        The read values

    """
    try:
        with open(fname, 'r', newline='', encoding='utf-8', errors='ignore') as csvfile:
            # skip the first line
            next(csvfile)

            # first count the lines
            reader = csv.DictReader(
                csvfile, fieldnames=['date', 'preciptype', 'lwc', 'rr', 'zh',
                                     'zv', 'zdr', 'ldr', 'ah', 'av', 'adiff',
                                     'kdp', 'deltaco', 'rhohv'],
                dialect='excel-tab')
            nrows = sum(1 for row in reader)

            preciptype = np.ma.empty(nrows, dtype='float32')
            lwc = np.ma.empty(nrows, dtype='float32')
            rr = np.ma.empty(nrows, dtype='float32')
            zh = np.ma.empty(nrows, dtype='float32')
            zv = np.ma.empty(nrows, dtype='float32')
            zdr = np.ma.empty(nrows, dtype='float32')
            ldr = np.ma.empty(nrows, dtype='float32')
            ah = np.ma.empty(nrows, dtype='float32')
            av = np.ma.empty(nrows, dtype='float32')
            adiff = np.ma.empty(nrows, dtype='float32')
            kdp = np.ma.empty(nrows, dtype='float32')
            deltaco = np.ma.empty(nrows, dtype='float32')
            rhohv = np.ma.empty(nrows, dtype='float32')

            # now read the data
            csvfile.seek(0)
            # skip the first line
            next(csvfile)

            reader = csv.DictReader(
                csvfile, fieldnames=['date', 'preciptype', 'lwc', 'rr', 'zh',
                                     'zv', 'zdr', 'ldr', 'ah', 'av', 'adiff',
                                     'kdp', 'deltaco', 'rhohv'],
                dialect='excel-tab')
            i = 0
            date = list()
            for row in reader:
                date.append(datetime.datetime.strptime(
                    row['date'], '%Y-%m-%d %H:%M:%S'))
                preciptype[i] = float(row['preciptype'])
                lwc[i] = float(row['lwc'])
                rr[i] = float(row['rr'])
                zh[i] = float(row['zh'])
                zv[i] = float(row['zv'])
                zdr[i] = float(row['zdr'])
                ldr[i] = float(row['ldr'])
                ah[i] = float(row['ah'])
                av[i] = float(row['av'])
                adiff[i] = float(row['adiff'])
                kdp[i] = float(row['kdp'])
                deltaco[i] = float(row['deltaco'])
                rhohv[i] = float(row['rhohv'])
                i += 1

            csvfile.close()

            return (date, preciptype, lwc, rr, zh, zv, zdr, ldr, ah, av,
                    adiff, kdp, deltaco, rhohv)
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None, None, None)


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
        if (line.startswith("# CRC:")):
            break
        linenum += 1

    pattern = {
        'angle': np.empty(linenum),
        'attenuation': np.empty(linenum)
    }

    pfile.seek(0)
    cnt = 0
    for line in pfile:
        if (line.startswith("# CRC:")):
            break
        line = line.strip()
        fields = line.split(',')
        pattern['angle'][cnt] = float(fields[0])
        pattern['attenuation'][cnt] = float(fields[1])
        cnt += 1

    pfile.close()

    if (twoway):
        pattern['attenuation'] = 2. * pattern['attenuation']

    if (linear):
        pattern['attenuation'] = 10.**(pattern['attenuation']/10.)

    return pattern
