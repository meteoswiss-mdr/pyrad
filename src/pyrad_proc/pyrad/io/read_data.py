"""
pyrad.io.read_data
====================

Functions for reading pyrad input data, i.e. radar files

.. autosummary::
    :toctree: generated/

    get_data
    read_status
    read_rad4alp_cosmo
    read_timeseries
    get_sensor_data
    read_smn
    read_disdro_scattering
    find_cosmo_file
    find_rad4alpcosmo_file
    get_datatypemetranet
    get_fieldname_rainbow
    get_file_list
    get_datatypefields
    get_datasetfields
    get_datetime

"""

import glob
import datetime
import os
import csv
import xml.etree.ElementTree as et
from warnings import warn

import numpy as np

import pyart
import wradlib as wrl


def get_data(voltime, datatypesdescr, cfg):
    """
    Reads pyrad input data.

    Parameters
    ----------
    voltime : datetime object
        volume scan time

    datatypesdescr : list
        list of radar field types to read.
        Format : [radar file type]:[datatype]

    cfg: dictionary of dictionaries
        configuration info to figure out where the data is

    Returns
    -------
    radar : Radar
        radar object

    """

    datatype_rainbow = list()
    datatype_rad4alp = list()
    datatype_cosmo = list()
    datatype_rad4alpcosmo = list()
    for datatypedescr in datatypesdescr:
        datagroup, datatype, dataset, product = get_datatypefields(
            datatypedescr)
        if datagroup == 'RAINBOW':
            datatype_rainbow.append(datatype)
        elif datagroup == 'RAD4ALP':
            datatype_rad4alp.append(datatype)
        elif datagroup == 'COSMO':
            datatype_cosmo.append(datatype)
        elif datagroup == 'RAD4ALPCOSMO':
            datatype_rad4alpcosmo.append(datatype)

    ndatatypes_rainbow = len(datatype_rainbow)
    ndatatypes_rad4alp = len(datatype_rad4alp)
    ndatatypes_cosmo = len(datatype_cosmo)
    ndatatypes_rad4alpcosmo = len(datatype_rad4alpcosmo)
    radar = None

    if ndatatypes_rainbow > 0:
        # read master data type
        # if radar object does not exist yet create it
        # otherwise add the new field
        daydir = voltime.strftime('%Y-%m-%d')
        datapath = cfg['datapath']+cfg['ScanList'][0]+daydir+'/'

        fdatetime = voltime.strftime('%Y%m%d%H%M%S')+'00'

        if (datatype_rainbow[0] != 'Nh') and (datatype_rainbow[0] != 'Nv'):
            filename = glob.glob(datapath+fdatetime+datatype_rainbow[0]+'.*')
        elif datatype_rainbow[0] == 'Nh':
            filename = glob.glob(datapath+fdatetime+'dBZ.*')
        else:
            filename = glob.glob(datapath+fdatetime+'dBZv.*')

        # create radar object
        radar = pyart.aux_io.read_rainbow_wrl(filename[0])
        if (datatype_rainbow[0] == 'Nh') or (datatype_rainbow[0] == 'Nv'):
            rbf = wrl.io.read_Rainbow(filename[0], loaddata=False)
            # check the number of slices
            nslices = int(rbf['volume']['scan']['pargroup']['numele'])
            if nslices > 1:
                single_slice = False
                common_slice_info = rbf['volume']['scan']['slice'][0]
            else:
                single_slice = True
                common_slice_info = rbf['volume']['scan']['slice']

            if datatype_rainbow[0] == 'Nh':
                noisedBZ1km_h = float(
                    common_slice_info['noise_power_dbz'])
                noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                    radar.nrays, noisedBZ1km_h, radar.range['data'], 1.,
                    noise_field='noisedBZ_hh')
                radar.fields = dict()
                radar.add_field('noisedBZ_hh', noisedBZ_h)
            else:
                noisedBZ1km_v = float(
                    common_slice_info['noise_power_dbz_dpv'])
                noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                    radar.nrays, noisedBZ1km_v, radar.range['data'], 1.,
                    noise_field='noisedBZ_vv')
                radar.fields = dict()
                radar.add_field('noisedBZ_vv', noisedBZ_v)

        # add other fields in the same scan
        for i in range(1, ndatatypes_rainbow):
            if (datatype_rainbow[i] != 'Nh') and (datatype_rainbow[i] != 'Nv'):
                filename = glob.glob(
                    datapath+fdatetime+datatype_rainbow[i]+'.*')
            elif datatype_rainbow[i] == 'Nh':
                filename = glob.glob(datapath+fdatetime+'dBZ.*')
            else:
                filename = glob.glob(datapath+fdatetime+'dBZv.*')

            radar_aux = pyart.aux_io.read_rainbow_wrl(filename[0])
            if (datatype_rainbow[i] == 'Nh') or (datatype_rainbow[i] == 'Nv'):
                rbf = wrl.io.read_Rainbow(filename[0], loaddata=False)
                # check the number of slices
                nslices = int(rbf['volume']['scan']['pargroup']['numele'])
                if nslices > 1:
                    single_slice = False
                    common_slice_info = rbf['volume']['scan']['slice'][0]
                else:
                    single_slice = True
                    common_slice_info = rbf['volume']['scan']['slice']

                if datatype_rainbow[i] == 'Nh':
                    noisedBZ1km_h = float(
                        common_slice_info['noise_power_dbz'])
                    noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                        radar_aux.nrays, noisedBZ1km_h,
                        radar_aux.range['data'], 1.,
                        noise_field='noisedBZ_hh')
                    radar_aux.fields = dict()
                    radar_aux.add_field('noisedBZ_hh', noisedBZ_h)
                else:
                    noisedBZ1km_v = float(
                        common_slice_info['noise_power_dbz_dpv'])
                    noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                        radar_aux.nrays, noisedBZ1km_v,
                        radar_aux.range['data'], 1.,
                        noise_field='noisedBZ_vv')
                    radar_aux.fields = dict()
                    radar_aux.add_field('noisedBZ_vv', noisedBZ_v)

            for field_name in radar_aux.fields.keys():
                break
            field_data = radar_aux.fields[field_name]['data']
            field_metadata = pyart.config.get_metadata(field_name)
            field_metadata['data'] = field_data
            radar.add_field(field_name, field_metadata)

        # merge scans into a single radar instance
        nscans = len(cfg['ScanList'])
        if nscans > 1:
            if (datatype_rainbow[0] == 'Nh') or (datatype_rainbow[0] == 'Nv'):
                datadescriptor = 'RAINBOW:dBZ'
            else:
                datadescriptor = 'RAINBOW:'+datatype_rainbow[0]
            endtime = voltime+datetime.timedelta(minutes=cfg['ScanPeriod'])
            for i in range(1, nscans):
                filelist = get_file_list(
                    cfg['ScanList'][i], datadescriptor, voltime, endtime, cfg)
                scantime = get_datetime(filelist[0], datadescriptor)
                daydir = scantime.strftime('%Y-%m-%d')
                datapath = cfg['datapath']+cfg['ScanList'][i]+daydir+'/'
                fdatetime = scantime.strftime('%Y%m%d%H%M%S')+'00'

                if ((datatype_rainbow[0] != 'Nh') and
                        (datatype_rainbow[0] != 'Nv')):
                    filename = glob.glob(
                        datapath+fdatetime+datatype_rainbow[0]+'.*')
                elif datatype_rainbow[0] == 'Nh':
                    filename = glob.glob(datapath+fdatetime+'dBZ.*')
                else:
                    filename = glob.glob(datapath+fdatetime+'dBZv.*')

                radar_aux = pyart.aux_io.read_rainbow_wrl(filename[0])
                if ((datatype_rainbow[0] == 'Nh') or
                        (datatype_rainbow[0] == 'Nv')):
                    rbf = wrl.io.read_Rainbow(filename[0], loaddata=False)
                    # check the number of slices
                    nslices = int(rbf['volume']['scan']['pargroup']['numele'])
                    if nslices > 1:
                        single_slice = False
                        common_slice_info = rbf['volume']['scan']['slice'][0]
                    else:
                        single_slice = True
                        common_slice_info = rbf['volume']['scan']['slice']

                    if datatype_rainbow[0] == 'Nh':
                        noisedBZ1km_h = float(
                            common_slice_info['noise_power_dbz'])
                        noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                            radar_aux.nrays, noisedBZ1km_h,
                            radar_aux.range['data'], 1.,
                            noise_field='noisedBZ_hh')
                        radar_aux.add_field('noisedBZ_hh', noisedBZ_h)
                    else:
                        noisedBZ1km_v = float(
                            common_slice_info['noise_power_dbz_dpv'])
                        noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                            radar_aux.nrays, noisedBZ1km_v,
                            radar_aux.range['data'], 1.,
                            noise_field='noisedBZ_vv')
                        radar_aux.add_field('noisedBZ_vv', noisedBZ_v)

                # add other fields in the same scan
                for j in range(1, ndatatypes_rainbow):
                    if ((datatype_rainbow[j] != 'Nh') and
                            (datatype_rainbow[j] != 'Nv')):
                        filename = glob.glob(
                            datapath+fdatetime+datatype_rainbow[j]+'.*')
                    elif datatype_rainbow[j] == 'Nh':
                        filename = glob.glob(datapath+fdatetime+'dBZ.*')
                    else:
                        filename = glob.glob(datapath+fdatetime+'dBZv.*')

                    radar_aux2 = pyart.aux_io.read_rainbow_wrl(filename[0])
                    if ((datatype_rainbow[j] == 'Nh') or
                            (datatype_rainbow[j] == 'Nv')):
                        rbf = wrl.io.read_Rainbow(filename[0], loaddata=False)
                        # check the number of slices
                        nslices = int(
                            rbf['volume']['scan']['pargroup']['numele'])
                        if nslices > 1:
                            single_slice = False
                            common_slice_info = (
                                rbf['volume']['scan']['slice'][0])
                        else:
                            single_slice = True
                            common_slice_info = rbf['volume']['scan']['slice']

                        if datatype_rainbow[j] == 'Nh':
                            noisedBZ1km_h = float(
                                common_slice_info['noise_power_dbz'])
                            noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                                radar_aux2.nrays, noisedBZ1km_h,
                                radar_aux2.range['data'], 1.,
                                noise_field='noisedBZ_hh')
                            radar_aux2.fields = dict()
                            radar_aux2.add_field('noisedBZ_hh', noisedBZ_h)
                        else:
                            noisedBZ1km_v = float(
                                common_slice_info['noise_power_dbz_dpv'])
                            noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                                radar_aux2.nrays, noisedBZ1km_v,
                                radar_aux2.range['data'], 1.,
                                noise_field='noisedBZ_vv')
                            radar_aux2.fields = dict()
                            radar_aux2.add_field('noisedBZ_vv', noisedBZ_v)

                    for field_name in radar_aux2.fields.keys():
                        break
                    field_data = radar_aux2.fields[field_name]['data']
                    field_metadata = pyart.config.get_metadata(field_name)
                    field_metadata['data'] = field_data
                    radar_aux.add_field(field_name, field_metadata)

                radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    elif ndatatypes_rad4alp > 0:
        if (cfg['RadarRes'] is None) or (cfg['RadarName'] is None):
            raise ValueError(
                'ERROR: Radar Name and Resolution \
                not specified in config file. Unable to load rad4alp data')

        metranet_field_names = dict()
        for datatype in datatype_rad4alp:
            if (datatype != 'Nh') and (datatype != 'Nv'):
                metranet_field_names.update(get_datatypemetranet(datatype))

        dayinfo = voltime.strftime('%y%j')
        timeinfo = voltime.strftime('%H%M')
        basename = 'P'+cfg['RadarRes']+cfg['RadarName']+dayinfo
        datapath = cfg['datapath']+dayinfo+'/'+basename+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+cfg['ScanList'][0])
        radar = pyart.aux_io.read_metranet(
            filename[0], field_names=metranet_field_names)

        # create secondary moments
        if ('Nh' in datatype_rad4alp) or ('Nv' in datatype_rad4alp):
            # read radar information in status file
            root = read_status(voltime, cfg)
            sweep_number = int(cfg['ScanList'][0])-1

            if 'Nh' in datatype_rad4alp:
                noise_h_vec = root.findall(
                    "./sweep/RADAR/STAT/CALIB/noisepower_frontend_h_inuse")
                rconst_h_vec = root.findall(
                    "./sweep/RADAR/STAT/CALIB/rconst_h")

                noisedBADU_h = 10.*np.log10(
                    float(noise_h_vec[sweep_number].attrib['value']))
                rconst_h = float(rconst_h_vec[sweep_number].attrib['value'])

                noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                    radar.nrays, noisedBADU_h+rconst_h, radar.range['data'],
                    100., noise_field='noisedBZ_hh')

                radar.add_field('noisedBZ_hh', noisedBZ_h)

            if 'Nv' in datatype_rad4alp:
                noise_v_vec = root.findall(
                    "./sweep/RADAR/STAT/CALIB/noisepower_frontend_v_inuse")
                rconst_v_vec = root.findall(
                    "./sweep/RADAR/STAT/CALIB/rconst_v")

                noisedBADU_v = 10.*np.log10(
                    float(noise_v_vec[sweep_number].attrib['value']))
                rconst_v = float(rconst_v_vec[sweep_number].attrib['value'])

                noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                    radar.nrays, noisedBADU_v+rconst_v, radar.range['data'],
                    100., noise_field='noisedBZ_vv')

                radar.add_field('noisedBZ_vv', noisedBZ_v)

        nelevs = len(cfg['ScanList'])
        # merge the elevations into a single radar instance
        for i in range(1, nelevs):
            filename = glob.glob(
                datapath+basename+timeinfo+'*.'+cfg['ScanList'][i])
            radar_aux = pyart.aux_io.read_metranet(
                filename[0], field_names=metranet_field_names)

            if ('Nh' in datatype_rad4alp) or ('Nv' in datatype_rad4alp):
                sweep_number = int(cfg['ScanList'][i])-1

                if 'Nh' in datatype_rad4alp:
                    noisedBADU_h = 10.*np.log10(
                        float(noise_h_vec[sweep_number].attrib['value']))
                    rconst_h = float(
                        rconst_h_vec[sweep_number].attrib['value'])

                    noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                        radar_aux.nrays, noisedBADU_h+rconst_h,
                        radar_aux.range['data'], 100.,
                        noise_field='noisedBZ_hh')

                    radar_aux.add_field('noisedBZ_hh', noisedBZ_h)

                if 'Nv' in datatype_rad4alp:
                    noisedBADU_v = 10.*np.log10(
                        float(noise_v_vec[sweep_number].attrib['value']))
                    rconst_v = float(
                        rconst_v_vec[sweep_number].attrib['value'])

                    noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                        radar_aux.nrays, noisedBADU_v+rconst_v,
                        radar_aux.range['data'], 100.,
                        noise_field='noisedBZ_vv')

                    radar_aux.add_field('noisedBZ_vv', noisedBZ_v)

            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    # add COSMO files to the radar field
    if ndatatypes_cosmo > 0:
        # look for COSMO data
        filename_list = list()
        for i in range(ndatatypes_cosmo):
            filename = find_cosmo_file(
                voltime, datatype_cosmo[i], cfg, cfg['ScanList'][0])
            if filename is not None:
                filename_list.append(filename)

        nfiles_valid = len(filename_list)
        if nfiles_valid > 0:
            if radar is None:
                radar = pyart.aux_io.read_rainbow_wrl(filename_list[0])
            else:
                radar_aux = pyart.aux_io.read_rainbow_wrl(filename_list[0])
                for field_name in radar_aux.fields.keys():
                    break
                field_data = radar_aux.fields[field_name]['data']
                field_metadata = pyart.config.get_metadata(field_name)
                field_metadata['data'] = field_data
                radar.add_field(field_name, field_metadata)

            # add other COSMO fields in the same scan
            for i in range(1, nfiles_valid):
                radar_aux = pyart.aux_io.read_rainbow_wrl(filename_list[i])
                for field_name in radar_aux.fields.keys():
                    break
                field_data = radar_aux.fields[field_name]['data']
                field_metadata = pyart.config.get_metadata(field_name)
                field_metadata['data'] = field_data
                radar.add_field(field_name, field_metadata)

        # merge scans into a single radar instance
        nscans = len(cfg['ScanList'])
        if nscans > 1:
            endtime = voltime+datetime.timedelta(minutes=cfg['ScanPeriod'])
            for i in range(1, nscans):
                radar_aux = None
                # look for COSMO data
                filename = find_cosmo_file(
                    voltime, datatype_cosmo[0], cfg, cfg['ScanList'][i])
                if filename is not None:
                    radar_aux = pyart.aux_io.read_rainbow_wrl(filename)

                # add other fields in the same scan
                for j in range(1, ndatatypes_cosmo):
                    # look for COSMO data
                    filename = find_cosmo_file(
                        voltime, datatype_cosmo[j], cfg, cfg['ScanList'][i])
                    if filename is not None:
                        radar_aux2 = pyart.aux_io.read_rainbow_wrl(filename)

                        if radar_aux is not None:
                            for field_name in radar_aux2.fields.keys():
                                break
                            field_data = radar_aux2.fields[field_name]['data']
                            field_metadata = (
                                pyart.config.get_metadata(field_name))
                            field_metadata['data'] = field_data
                            radar_aux.add_field(field_name, field_metadata)
                        else:
                            radar_aux = radar_aux2

                if radar_aux is not None:
                    radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    elif ndatatypes_rad4alpcosmo > 0:
        if (cfg['RadarRes'] is None) or (cfg['RadarName'] is None):
            raise ValueError(
                'ERROR: Radar Name and Resolution \
                not specified in config file. \
                Unable to load rad4alp COSMO data')

        for i in range(ndatatypes_rad4alpcosmo):
            # create the radar object where to store the data
            # taking as reference the metranet polar file
            metranet_field_names = dict()
            metranet_field_names.update(
                get_datatypemetranet('dBZ'))

            dayinfo = voltime.strftime('%y%j')
            timeinfo = voltime.strftime('%H%M')
            basename = (
                'P'+cfg['RadarRes']+cfg['RadarName']+dayinfo)
            datapath = cfg['datapath']+dayinfo+'/'+basename+'/'
            filename = glob.glob(
                datapath+basename+timeinfo+'*.'+cfg['ScanList'][0])
            radar_aux = pyart.aux_io.read_metranet(
                filename[0], field_names=metranet_field_names)
            radar_aux.fields = dict()

            # look for rad4alp COSMO data
            filename = find_rad4alpcosmo_file(
                voltime, datatype_rad4alpcosmo[i], cfg, cfg['ScanList'][0])
            if filename is not None:
                cosmo_field = read_rad4alp_cosmo(
                    filename, datatype_rad4alpcosmo[i])
                if cosmo_field is not None:
                    radar_aux.add_field(
                        get_fieldname_rainbow(datatype_rad4alpcosmo[i]),
                        cosmo_field)

                    # add other elevations
                    for j in range(1, len(cfg['ScanList'])):
                        # create the radar object where to store the data
                        # taking as reference the metranet polar file
                        metranet_field_names = dict()
                        metranet_field_names.update(
                            get_datatypemetranet('dBZ'))

                        dayinfo = voltime.strftime('%y%j')
                        timeinfo = voltime.strftime('%H%M')
                        basename = (
                            'P'+cfg['RadarRes']+cfg['RadarName']+dayinfo)
                        datapath = cfg['datapath']+dayinfo+'/'+basename+'/'
                        filename = glob.glob(
                            datapath+basename+timeinfo+'*.' +
                            cfg['ScanList'][j])
                        radar_aux2 = pyart.aux_io.read_metranet(
                            filename[0], field_names=metranet_field_names)
                        radar_aux2.fields = dict()

                        # look for rad4alp COSMO data
                        filename = find_rad4alpcosmo_file(
                            voltime, datatype_rad4alpcosmo[i], cfg,
                            cfg['ScanList'][j])
                        if filename is not None:
                            cosmo_field = read_rad4alp_cosmo(
                                filename, datatype_rad4alpcosmo[i])
                            if cosmo_field is not None:
                                radar_aux2.add_field(
                                    get_fieldname_rainbow(
                                        datatype_rad4alpcosmo[i]),
                                    cosmo_field)

                                radar_aux = pyart.util.radar_utils.join_radar(
                                    radar_aux, radar_aux2)
                            else:
                                radar_aux = None
                        else:
                            radar_aux = None
                    if radar_aux is not None:
                        if radar is None:
                            radar = radar_aux
                        else:
                            for field_name in radar_aux.fields.keys():
                                break

                            field_data = radar_aux.fields[field_name]['data']
                            field_metadata = pyart.config.get_metadata(
                                field_name)
                            field_metadata['data'] = field_data
                            radar.add_field(field_name, field_metadata)

    return radar


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
        fill_value = pyart.config.get_fillvalue()

        if datatype == 'TEMP':
            field_name = get_fieldname_rainbow(datatype)
            field_data = np.ma.masked_where(
                mask, (bindata-1).astype(float)*0.5-87.)
            field_data.set_fill_value(fill_value)
            field_data.data[mask] = fill_value

            field = pyart.config.get_metadata(field_name)
            field['data'] = field_data
            return field
        elif datatype == 'ISO0':
            field_name = get_fieldname_rainbow(datatype)
            field_data = np.ma.masked_where(mask, (bindata-1).astype(float))
            field_data.set_fill_value(fill_value)
            field_data.data[mask] = fill_value

            field = pyart.config.get_metadata(field_name)
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
            value = np.ma.empty(nrows, dtype='float32')

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

            fill_value = pyart.config.get_fillvalue()
            value = np.ma.masked_where(value == fill_value, value)
            value.set_fill_value(fill_value)

            return date, value
    except EnvironmentError:
        warn('WARNING: Unable to read file '+fname)
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
    id, date , pressure, temp, rh, precip, wspeed, wdir : tupple
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


def find_cosmo_file(voltime, datatype, cfg, scanid):
    """
    Search a COSMO file

    Parameters
    ----------
    voltime : datetime object
        volume scan time

    datatype : type of COSMO data to look for

    cfg: dictionary of dictionaries
        configuration info to figure out where the data is

    scanid: str
        name of the scan

    Returns
    -------
    fname : str
        Name of COSMO file if it exists. None otherwise

    """
    # hour rounded date-time
    fdatetime = voltime.strftime('%Y%m%d%H')+'000000'

    # initial run time to look for
    hvol = int(voltime.strftime('%H'))
    runhour0 = int(hvol/cfg['CosmoRunFreq'])*cfg['CosmoRunFreq']
    runtime0 = voltime.replace(hour=runhour0, minute=0, second=0)

    # look for cosmo file
    found = False
    nruns_to_check = int((cfg['CosmoForecasted']-1)/cfg['CosmoRunFreq'])
    for i in range(nruns_to_check):
        runtime = runtime0-datetime.timedelta(hours=i * cfg['CosmoRunFreq'])
        runtimestr = runtime.strftime('%Y%m%d%H')+'000000'

        daydir = runtime.strftime('%Y-%m-%d')
        datapath = cfg['cosmopath']+datatype+'/'+scanid+daydir+'/'

        search_name = (
            datapath+datatype+'_RUN'+runtimestr+'_DX50'+fdatetime+'.*')
        print('Looking for file: '+search_name)
        fname = glob.glob(search_name)
        if len(fname) > 0:
            found = True
            break

    if not found:
        warn('WARNING: Unable to get COSMO '+datatype+' information')
        return None
    else:
        return fname[0]


def find_rad4alpcosmo_file(voltime, datatype, cfg, scanid):
    """
    Search a COSMO file

    Parameters
    ----------
    voltime : datetime object
        volume scan time

    datatype : type of COSMO data to look for

    cfg: dictionary of dictionaries
        configuration info to figure out where the data is

    Returns
    -------
    fname : str
        Name of COSMO file if it exists. None otherwise

    scanid: str
        name of the scan

    """
    # hour rounded date-time
    fdatetime = voltime.strftime('%y%j%H')+'00'

    # initial run time to look for
    hvol = int(voltime.strftime('%H'))
    runhour0 = int(hvol/cfg['CosmoRunFreq'])*cfg['CosmoRunFreq']
    runtime0 = voltime.replace(hour=runhour0, minute=0, second=0)

    # look for cosmo file
    found = False
    nruns_to_check = int((cfg['CosmoForecasted']-1)/cfg['CosmoRunFreq'])
    id = 'P'+cfg['RadarRes']+cfg['RadarName']
    for i in range(nruns_to_check):
        runtime = runtime0-datetime.timedelta(hours=i * cfg['CosmoRunFreq'])
        runtimestr = runtime.strftime('%y%j%H')+'00'

        daydir = runtime.strftime('%y%j')
        datapath = cfg['cosmopath']+datatype+'/'+id+'/'+daydir+'/'

        search_name = (
            datapath+datatype+'_RUN'+runtimestr+'_'+id+fdatetime+'.'+scanid +
            '.bin')
        print('Looking for file: '+search_name)
        fname = glob.glob(search_name)
        if len(fname) > 0:
            found = True
            break

    if not found:
        warn('WARNING: Unable to get COSMO '+datatype+' information')
        return None
    else:
        return fname[0]


def get_datatypemetranet(datatype):
    """
    maps de config file radar data type name into the corresponding metranet
    data type  name and Py-ART field name

    Parameters
    ----------
    datatype : str
        config file radar data type name

    Returns
    -------
    metranet type : dict
        dictionary containing the metranet data type name and its
        corresponding Py-ART field name

    """
    if datatype == 'dBZ':
        datatype_metranet = 'ZH'
        field_name = 'reflectivity'
    elif datatype == 'dBZv':
        datatype_metranet = 'ZV'
        field_name = 'reflectivity_vv'
    elif datatype == 'ZDR':
        datatype_metranet = 'ZDR'
        field_name = 'differential_reflectivity'
    elif datatype == 'uRhoHV':
        datatype_metranet = 'RHO'
        field_name = 'uncorrected_cross_correlation_ratio'
    elif datatype == 'uPhiDP':
        datatype_metranet = 'PHI'
        field_name = 'uncorrected_differential_phase'
    elif datatype == 'V':
        datatype_metranet = 'VEL'
        field_name = 'velocity'
    elif datatype == 'W':
        datatype_metranet = 'WID'
        field_name = 'spectrum_width'
    else:
        raise ValueError(
            'ERROR: Metranet fields do not contain datatype '+datatype)

    return {datatype_metranet: field_name}


def get_fieldname_rainbow(datatype):
    """
    maps de config file radar data type name into the corresponding rainbow
    Py-ART field name

    Parameters
    ----------
    datatype : str
        config file radar data type name

    Returns
    -------
    field_name : str
        Py-ART field name

    """
    if datatype == 'dBZ':
        field_name = 'reflectivity'
    elif datatype == 'dBZv':
        field_name = 'reflectivity_vv'
    elif datatype == 'Nh':
        field_name = 'noisedBZ_hh'
    elif datatype == 'Nv':
        field_name = 'noisedBZ_vv'
    elif datatype == 'SNRh':
        field_name = 'signal_to_noise_ratio_hh'
    elif datatype == 'SNRv':
        field_name = 'signal_to_noise_ratio_vv'
    elif datatype == 'dBZc':
        field_name = 'corrected_reflectivity'
    elif datatype == 'ZDR':
        field_name = 'differential_reflectivity'
    elif datatype == 'ZDRc':
        field_name = 'corrected_differential_reflectivity'
    elif datatype == 'RhoHV':
        field_name = 'cross_correlation_ratio'
    elif datatype == 'uRhoHV':
        field_name = 'uncorrected_cross_correlation_ratio'
    elif datatype == 'RhoHVc':
        field_name = 'corrected_cross_correlation_ratio'
    elif datatype == 'L':
        field_name = 'logarithmic_cross_correlation_ratio'
    elif datatype == 'CDR':
        field_name = 'circular_depolarization_ratio'
    elif datatype == 'PhiDP':
        field_name = 'differential_phase'
    elif datatype == 'uPhiDP':
        field_name = 'uncorrected_differential_phase'
    elif datatype == 'PhiDPc':
        field_name = 'corrected_differential_phase'
    elif datatype == 'KDP':
        field_name = 'specific_differential_phase'
    elif datatype == 'KDPc':
        field_name = 'corrected_specific_differential_phase'
    elif datatype == 'V':
        field_name = 'velocity'
    elif datatype == 'W':
        field_name = 'spectrum_width'
    elif datatype == 'Ah':
        field_name = 'specific_attenuation'
    elif datatype == 'Ahc':
        field_name = 'corrected_specific_attenuation'
    elif datatype == 'Adp':
        field_name = 'specific_differential_attenuation'
    elif datatype == 'dBZc':
        field_name = 'corrected_reflectivity'
    elif datatype == 'ZDRc':
        field_name = 'corrected_differential_reflectivity'
    elif datatype == 'TEMP':
        field_name = 'temperature'
    elif datatype == 'ISO0':
        field_name = 'iso0'
    elif datatype == 'echoID':
        field_name = 'radar_echo_id'
    elif datatype == 'RR':
        field_name = 'radar_estimated_rain_rate'
    elif datatype == 'hydro':
        field_name = 'radar_echo_classification'
    else:
        raise ValueError('ERROR: Unknown data type '+datatype)

    return field_name


def get_file_list(scan, datadescriptor, starttime, endtime, cfg):
    """
    gets the list of files with a time period

    Parameters
    ----------
    scan : str
        scan name

    datadescriptor : str
        radar field type. Format : [radar file type]:[datatype]

    startime : datetime object
        start of time period

    endtime : datetime object
        end of time period

    cfg: dictionary of dictionaries
        configuration info to figure out where the data is

    Returns
    -------
    radar : Radar
        radar object

    """
    ndays = int(np.ceil(((endtime-starttime).total_seconds())/(3600.*24.)))
    datagroup, datatype, dataset, product = get_datatypefields(datadescriptor)

    if (datatype == 'Nh') or (datatype == 'Nv'):
        datatype = 'dBZ'

    t_filelist = []
    for i in range(ndays):
        if datagroup == 'RAINBOW':
            daydir = (
                starttime+datetime.timedelta(days=i)).strftime('%Y-%m-%d')
            dayinfo = (starttime+datetime.timedelta(days=i)).strftime('%Y%m%d')
            datapath = cfg['datapath']+scan+daydir+'/'
            dayfilelist = glob.glob(datapath+dayinfo+'*'+datatype+'.*')
            for filename in dayfilelist:
                t_filelist.append(filename)
        elif datagroup == 'RAD4ALP':
            dayinfo = (starttime+datetime.timedelta(days=i)).strftime('%y%j')
            basename = 'P'+cfg['RadarRes']+cfg['RadarName']+dayinfo
            datapath = cfg['datapath']+dayinfo+'/'+basename+'/'
            dayfilelist = glob.glob(datapath+basename+'*.'+scan)
            for filename in dayfilelist:
                t_filelist.append(filename)
        elif datagroup == 'SAVED':
            print('caca')

    filelist = []
    for filename in t_filelist:
        filenamestr = str(filename)
        fdatetime = get_datetime(filenamestr, datadescriptor)
        if (fdatetime >= starttime) and (fdatetime <= endtime):
            filelist.append(filenamestr)

    return sorted(filelist)


def get_datatypefields(datadescriptor):
    """
    splits the data type descriptor and provides each individual member

    Parameters
    ----------
    datadescriptor : str
        radar field type. Format : [radar file type]:[datatype]

    Returns
    -------
    datagroup : str
        data type group, i.e. RAINBOW, RAD4ALP, SAVED, COSMO, ...

    datatype : str
        data type, i.e. dBZ, ZDR, ISO0, ...

    dataset : str
        dataset type (for saved data only)
    product : str
        product type (for saved data only)

    """
    descrfields = datadescriptor.split(':')
    if len(descrfields) == 1:
        datagroup = 'RAINBOW'
        datatype = descrfields[0]
        dataset = None
        product = None
    else:
        datagroup = descrfields[0]
        if datagroup == 'SAVED':
            descrfields2 = descrfields[1].split(',')
            datatype = descrfields2[0]
            dataset = descrfields2[1]
            product = descrfields2[2]
        else:
            datatype = descrfields[1]
            dataset = None
            product = None

    return datagroup, datatype, dataset, product


def get_datasetfields(datasetdescr):
    """
    splits the dataset type descriptor and provides each individual member

    Parameters
    ----------
    datasetdescr : str
        dataset type. Format : [processing level]:[dataset type]

    Returns
    -------
    proclevel : str
        dataset processing level

    dataset : str
        dataset type, i.e. dBZ, ZDR, ISO0, ...

    """
    descrfields = datasetdescr.split(':')
    if len(descrfields) == 1:
        proclevel = 'l0'
        dataset = descrfields[0]
    else:
        proclevel = descrfields[0]
        dataset = descrfields[1]

    return proclevel, dataset


def get_datetime(fname, datadescriptor):
    """
    gets date and time from file name

    Parameters
    ----------
    fname : file name

    datadescriptor : str
        radar field type. Format : [radar file type]:[datatype]

    Returns
    -------
    fdatetime : datetime object
        date and time in file name

    """

    bfile = os.path.basename(fname)
    datagroup, datatype, dataset, product = get_datatypefields(datadescriptor)
    if datagroup == 'RAINBOW':
        datetimestr = bfile[0:14]
        fdatetime = datetime.datetime.strptime(datetimestr, '%Y%m%d%H%M%S')
    elif datagroup == 'RAD4ALP':
        datetimestr = bfile[3:12]
        fdatetime = datetime.datetime.strptime(datetimestr, '%y%j%H%M')
    elif datagroup == 'SAVED':
        print('caca')
    return fdatetime
