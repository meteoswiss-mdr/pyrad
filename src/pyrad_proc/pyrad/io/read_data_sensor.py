"""
pyrad.io.read_data_other
========================

Functions for reading data from other sensors

.. autosummary::
    :toctree: generated/

    read_lightning
    read_lightning_traj
    get_sensor_data
    read_smn
    read_smn2
    read_disdro_scattering

"""

import os
import datetime
import csv
from warnings import warn
from copy import deepcopy

import numpy as np


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
    flashnr, time_data, time_in_flash, lat, lon, alt, dBm : tupple
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
            time_data = list()
            for row in reader:
                flashnr[i] = int(row['flashnr'])
                time_data.append(fdatetime+datetime.timedelta(
                    seconds=float(row['time'])))
                time_in_flash[i] = float(row['time_in_flash'])
                lat[i] = float(row['lat'])
                lon[i] = float(row['lon'])
                alt[i] = float(row['alt'])
                dBm[i] = float(row['dBm'])

                i += 1

            time_data = np.array(time_data)

            if filter_data:
                flashnr_aux = deepcopy(flashnr)
                flashnr = flashnr[flashnr_aux > 0]
                time_data = time_data[flashnr_aux > 0]
                time_in_flash = time_in_flash[flashnr_aux > 0]
                lat = lat[flashnr_aux > 0]
                lon = lon[flashnr_aux > 0]
                alt = alt[flashnr_aux > 0]
                dBm = dBm[flashnr_aux > 0]

            csvfile.close()

            return flashnr, time_data, time_in_flash, lat, lon, alt, dBm
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None, None, None


def read_lightning_traj(fname):
    """
    Reads lightning trajectory data contained in a csv file. The file has the following
    fields:
        Date
        UTC [seconds since midnight]
        # Flash
        Flash Power (dBm)
        Value at flash
        Mean value in a 3x3x3 polar box
        Min value in a 3x3x3 polar box
        Max value in a 3x3x3 polar box
        # valid values in the polar box

    Parameters
    ----------
    fname : str
        path of time series file

    Returns
    -------
    time_flash, flashnr, dBm, val_at_flash, val_mean, val_min, val_max,
    nval : tupple
        A tupple containing the read values. None otherwise

    """
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith('#')),
                fieldnames=['Date', 'UTC', 'flashnr', 'dBm', 'at_flash',
                            'mean', 'min', 'max', 'nvalid'],
                delimiter=',')
            nrows = sum(1 for row in reader)

            time_flash = np.empty(nrows, dtype=datetime.datetime)
            flashnr = np.empty(nrows, dtype=int)
            dBm = np.empty(nrows, dtype=float)
            val_at_flash = np.ma.empty(nrows, dtype=float)
            val_mean = np.ma.empty(nrows, dtype=float)
            val_min = np.ma.empty(nrows, dtype=float)
            val_max = np.ma.empty(nrows, dtype=float)
            nval = np.empty(nrows, dtype=int)

            # now read the data
            csvfile.seek(0)
            reader = csv.DictReader(
                (row for row in csvfile if not row.startswith('#')),
                fieldnames=['Date', 'UTC', 'flashnr', 'dBm', 'at_flash',
                            'mean', 'min', 'max', 'nvalid'],
                delimiter=',')
            i = 0
            for row in reader:
                date_flash_aux = datetime.datetime.strptime(
                    row['Date'], '%d-%b-%Y')
                time_flash_aux = float(row['UTC'])
                time_flash[i] = date_flash_aux+datetime.timedelta(
                    seconds=time_flash_aux)

                flashnr[i] = int(float(row['flashnr']))
                dBm[i] = float(row['dBm'])
                val_at_flash[i] = float(row['at_flash'])
                val_mean[i] = float(row['mean'])
                val_min[i] = float(row['min'])
                val_max[i] = float(row['max'])
                nval[i] = int(float(row['nvalid']))

                i += 1

            csvfile.close()

            val_at_flash = np.ma.masked_invalid(val_at_flash)
            val_mean = np.ma.masked_invalid(val_mean)
            val_min = np.ma.masked_invalid(val_min)
            val_max = np.ma.masked_invalid(val_max)

            return (time_flash, flashnr, dBm, val_at_flash, val_mean, val_min,
                    val_max, nval)

    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None, None, None, None, None, None, None, None


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
        (sensor_id, sensordate, pressure, temp,
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
    smn_id, date , pressure, temp, rh, precip, wspeed, wdir : tupple
        The read values

    """
    fill_value = 10000000.0
    try:
        with open(fname, 'r', newline='') as csvfile:
            # first count the lines
            reader = csv.DictReader(csvfile)
            nrows = sum(1 for row in reader)
            smn_id = np.ma.empty(nrows, dtype='float32')
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
                smn_id[i] = float(row['StationID'])
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

            return smn_id, date, pressure, temp, rh, precip, wspeed, wdir
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
    smn_id, date , value : tupple
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
            smn_id = np.ma.empty(nrows, dtype=int)
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
                smn_id[i] = float(row['StationID'])
                date.append(datetime.datetime.strptime(
                    row['DateTime'], '%Y%m%d%H%M%S'))
                value[i] = float(row['Value'])
                i += 1

            csvfile.close()

            return smn_id, date, value
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
