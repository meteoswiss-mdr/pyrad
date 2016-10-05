"""
pyrad.io.read_data
====================

Functions for reading pyrad input data, i.e. radar files

.. autosummary::
    :toctree: generated/

    read_status
    read_rad4alp_cosmo
    read_timeseries
    read_sun_hits
    get_sensor_data
    read_smn
    read_disdro_scattering
    read_selfconsistency

"""

import glob
import datetime
import csv
import xml.etree.ElementTree as et
from warnings import warn

import numpy as np

from pyart.config import get_fillvalue, get_metadata


def read_status(voltime, cfg):
    """
    Reads rad4alp xml status file.

    Parameters
    ----------
    voltime : datetime object
        volume scan time

    cfg: dictionary of dictionaries
        configuration info to figure out where the data is

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
    basename = 'ST'+cfg['RadarName']+dayinfo
    datapath = cfg['datapath']+dayinfo+'/'+basename+'/'
    filename = glob.glob(datapath+basename+timeinfo+'*.xml')
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
            field_name = get_fieldname_rainbow(datatype)
            field_data = np.ma.masked_where(
                mask, (bindata-1).astype(float)*0.5-87.)

            field = get_metadata(field_name)
            field['data'] = field_data
            return field
        elif datatype == 'ISO0':
            field_name = get_fieldname_rainbow(datatype)
            field_data = np.ma.masked_where(mask, (bindata-1).astype(float))

            field = get_metadata(field_name)
            field['data'] = field_data
            return field
        else:
            warn('WARNING: Unknown COSMO data type '+datatype)
            return None

    except EnvironmentError:
        warn('WARNING: Unable to read file '+fname)
        return None


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

            return date, value
    except EnvironmentError:
        warn('WARNING: Unable to read file '+fname)
        return None, None


def read_sun_hits(fname):
    """
    Reads sun hits data contained in a csv file

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

            return (date, ray, nrng, rad_el, rad_az, sun_el, sun_az,
                    ph, ph_std, nph, nvalh, pv, pv_std, npv, nvalv,
                    zdr, zdr_std, nzdr, nvalzdr)

    except EnvironmentError:
        warn('WARNING: Unable to read file '+fname)
        return (None, None, None, None, None, None, None, None, None, None,
                None, None, None, None, None, None, None, None, None)


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
        warn('WARNING: Unknown sensor: '+cfg['sensor'])
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

            # convert precip from mm/10min to mm/h
            precip *= 6.

            return id, date, pressure, temp, rh, precip, wspeed, wdir
    except EnvironmentError:
        warn('WARNING: Unable to read file '+fname)
        return None, None, None, None, None, None, None, None


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
    id, date , pressure, temp, rh, precip, wspeed, wdir : arrays
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

            return (date, preciptype, lwc, rr, zh, zv, zdr, ldr, ah, av,
                    adiff, kdp, deltaco, rhohv)
    except EnvironmentError:
        warn('WARNING: Unable to read file '+fname)
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

            return zdr_kdpzh_table
    except EnvironmentError:
        warn('WARNING: Unable to read file '+fname)
        return None
