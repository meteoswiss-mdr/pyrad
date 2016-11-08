"""
pyrad.io.read_data_radar
========================

Functions for reading radar data files

.. autosummary::
    :toctree: generated/

    get_data
    merge_scans_rainbow
    merge_scans_dem
    merge_scans_rad4alp
    merge_scans_cosmo
    merge_scans_cosmo_rad4alp
    merge_scans_dem_rad4alp
    merge_fields_rainbow
    merge_fields_cosmo
    merge_fields_dem
    get_data_rainbow
    get_data_rad4alp
    interpol_field

"""

import glob
import datetime
import os
from warnings import warn

import numpy as np
from scipy.interpolate import RegularGridInterpolator

import pyart
import wradlib as wrl

from .read_data_other import read_status, read_rad4alp_cosmo, read_rad4alp_vis

from .io_aux import get_datatype_metranet, get_fieldname_pyart, get_file_list
from .io_aux import get_datatype_fields, get_datetime
from .io_aux import find_cosmo_file, find_rad4alpcosmo_file


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
    datatype_dem = list()
    datatype_rad4alpdem = list()
    for datatypedescr in datatypesdescr:
        datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datagroup == 'RAINBOW':
            datatype_rainbow.append(datatype)
        elif datagroup == 'RAD4ALP':
            datatype_rad4alp.append(datatype)
        elif datagroup == 'COSMO':
            datatype_cosmo.append(datatype)
        elif datagroup == 'RAD4ALPCOSMO':
            datatype_rad4alpcosmo.append(datatype)
        elif datagroup == 'DEM':
            datatype_dem.append(datatype)
        elif datagroup == 'RAD4ALPDEM':
            datatype_rad4alpdem.append(datatype)

    ndatatypes_rainbow = len(datatype_rainbow)
    ndatatypes_rad4alp = len(datatype_rad4alp)
    ndatatypes_cosmo = len(datatype_cosmo)
    ndatatypes_rad4alpcosmo = len(datatype_rad4alpcosmo)
    ndatatypes_dem = len(datatype_dem)
    ndatatypes_rad4alpdem = len(datatype_rad4alpdem)

    radar = None
    if ndatatypes_rainbow > 0:
        radar = merge_scans_rainbow(
            cfg['datapath'], cfg['ScanList'], voltime, cfg['ScanPeriod'],
            datatype_rainbow, cfg)

    elif ndatatypes_rad4alp > 0:
        radar = merge_scans_rad4alp(
            cfg['datapath'], cfg['ScanList'], cfg['RadarName'],
            cfg['RadarRes'], voltime, datatype_rad4alp, cfg)

    # add COSMO files to the radar field
    if ndatatypes_cosmo > 0:
        radar_aux = merge_scans_cosmo(voltime, datatype_cosmo, cfg)

        if radar is None:
            radar = radar_aux
        else:
            if radar_aux is not None:
                if ((np.allclose(
                        radar.azimuth['data'], radar_aux.azimuth['data'],
                        atol=0.5, equal_nan=True)) and
                        (np.allclose(
                            radar.elevation['data'],
                            radar_aux.elevation['data'], atol=0.5,
                            equal_nan=True))):
                    for field_name in radar_aux.fields.keys():
                        radar.add_field(
                            field_name, radar_aux.fields[field_name])
                else:
                    for field_name in radar_aux.fields.keys():
                        field_interp = interpol_field(
                            radar, radar_aux, field_name)
                        radar.add_field(field_name, field_interp)

    elif ndatatypes_rad4alpcosmo > 0:
        if (cfg['RadarRes'] is None) or (cfg['RadarName'] is None):
            raise ValueError(
                'ERROR: Radar Name and Resolution ' +
                'not specified in config file. ' +
                'Unable to load rad4alp COSMO data')

        for i in range(ndatatypes_rad4alpcosmo):
            radar_aux = merge_scans_cosmo_rad4alp(
                voltime, datatype_rad4alpcosmo[i], cfg)
            if radar is None:
                radar = radar_aux
            else:
                if radar_aux is not None:
                    for field_name in radar_aux.fields.keys():
                        radar.add_field(
                            field_name, radar_aux.fields[field_name])

    # add DEM files to the radar field
    if ndatatypes_dem > 0:
        radar_aux = merge_scans_dem(
            cfg['dempath'], cfg['ScanList'], datatype_dem, cfg)

        if radar is None:
            radar = radar_aux
        else:
            if radar_aux is not None:
                if ((np.allclose(
                        radar.azimuth['data'], radar_aux.azimuth['data'],
                        atol=0.5, equal_nan=True)) and
                        (np.allclose(
                            radar.elevation['data'],
                            radar_aux.elevation['data'], atol=0.5,
                            equal_nan=True))):
                    for field_name in radar_aux.fields.keys():
                        radar.add_field(
                            field_name, radar_aux.fields[field_name])
                else:
                    for field_name in radar_aux.fields.keys():
                        field_interp = interpol_field(
                            radar, radar_aux, field_name)
                        radar.add_field(field_name, field_interp)

    elif ndatatypes_rad4alpdem > 0:
        if (cfg['RadarRes'] is None) or (cfg['RadarName'] is None):
            raise ValueError(
                'ERROR: Radar Name and Resolution ' +
                'not specified in config file. ' +
                'Unable to load rad4alp DEM data')

        if cfg['RadarRes'] != 'L':
            raise ValueError(
                'ERROR: DEM files only available for rad4alp PL data')

        for i in range(ndatatypes_rad4alpdem):
            radar_aux = merge_scans_dem_rad4alp(
                voltime, datatype_rad4alpdem[i], cfg)
            if radar is None:
                radar = radar_aux
            else:
                if radar_aux is not None:
                    for field_name in radar_aux.fields.keys():
                        radar.add_field(
                            field_name, radar_aux.fields[field_name])

    return radar


def merge_scans_rainbow(basepath, scan_list, voltime, scan_period,
                        datatype_list, cfg):
    """
    merge rainbow scans

    Parameters
    ----------
    basepath : str
        base path of rad4alp radar data
    scan_list : list
        list of scans
    voltime: datetime object
        reference time of the scan
    scan_period : float
        time from reference time where to look for other scans data
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = merge_fields_rainbow(
        basepath, scan_list[0], voltime, datatype_list)

    # merge scans into a single radar instance
    nscans = len(scan_list)
    if nscans > 1:
        if (datatype_list[0] == 'Nh') or (datatype_list[0] == 'Nv'):
            datadescriptor = 'RAINBOW:dBZ'
        else:
            datadescriptor = 'RAINBOW:'+datatype_list[0]
        endtime = voltime+datetime.timedelta(minutes=scan_period)
        for i in range(1, nscans):
            filelist = get_file_list(
                scan_list[i], datadescriptor, voltime, endtime, cfg)
            scantime = get_datetime(filelist[0], datadescriptor)

            radar_aux = merge_fields_rainbow(
                basepath, scan_list[i], scantime, datatype_list)

            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_dem(basepath, scan_list, datatype_list, cfg):
    """
    merge rainbow scans

    Parameters
    ----------
    basepath : str
        base path of rad4alp radar data
    scan_list : list
        list of scans
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = merge_fields_dem(
        basepath, scan_list[0], datatype_list)

    # merge scans into a single radar instance
    nscans = len(scan_list)
    if nscans > 1:
        for i in range(1, nscans):
            radar_aux = merge_fields_dem(
                basepath, scan_list[i], datatype_list)

            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_rad4alp(basepath, scan_list, radar_name, radar_res, voltime,
                        datatype_list, cfg):
    """
    merge rad4alp data.

    Parameters
    ----------
    basepath : str
        base path of rad4alp radar data
    scan_list : list
        list of scans (001 to 020)
    radar_name : str
        radar_name (A, D, L, ...)
    radar_res : str
        radar resolution (H or L)
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary

    Returns
    -------
    radar : Radar
        radar object

    """
    if (radar_name is None) or (radar_res is None):
        raise ValueError(
            'ERROR: Radar Name and Resolution not specified in config file.' +
            ' Unable to load rad4alp data')

    dayinfo = voltime.strftime('%y%j')
    timeinfo = voltime.strftime('%H%M')
    basename = 'P'+radar_res+radar_name+dayinfo
    datapath = basepath+dayinfo+'/'+basename+'/'
    filename = glob.glob(datapath+basename+timeinfo+'*.'+scan_list[0])

    radar = get_data_rad4alp(filename[0], datatype_list, scan_list[0], cfg)

    nelevs = len(scan_list)
    # merge the elevations into a single radar instance
    for i in range(1, nelevs):
        filename = glob.glob(datapath+basename+timeinfo+'*.'+scan_list[i])
        radar_aux = get_data_rad4alp(
            filename[0], datatype_list, scan_list[i], cfg)

        radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_cosmo(voltime, datatype_list, cfg):
    """
    merge rainbow scans

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = None
    ndatatypes = len(datatype_list)
    # look for COSMO data
    filename_list = list()
    for i in range(ndatatypes):
        filename = find_cosmo_file(
            voltime, datatype_list[i], cfg, cfg['ScanList'][0])
        if filename is not None:
            filename_list.append(filename)

    nfiles_valid = len(filename_list)

    if nfiles_valid > 0:
        radar = merge_fields_cosmo(filename_list)

    # merge scans into a single radar instance
    nscans = len(cfg['ScanList'])
    if nscans > 1:
        endtime = voltime+datetime.timedelta(minutes=cfg['ScanPeriod'])
        for i in range(1, nscans):
            filename_list = list()
            for j in range(ndatatypes):
                filename = find_cosmo_file(
                    voltime, datatype_list[j], cfg, cfg['ScanList'][i])
                if filename is not None:
                    filename_list.append(filename)

            nfiles_valid = len(filename_list)
            if nfiles_valid > 0:
                radar_aux = merge_fields_cosmo(filename_list)

                if radar is None:
                    radar = radar_aux
                else:
                    radar = pyart.util.radar_utils.join_radar(
                        radar, radar_aux)

    return radar


def merge_scans_cosmo_rad4alp(voltime, datatype, cfg):
    """
    merge cosmo rad4alp scans. If data for all the scans cannot be retrieved
    returns None

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype : str
        name of the data type to read
    cfg : dict
        configuration dictionary

    Returns
    -------
    radar : Radar
        radar object

    """
    # look for rad4alp COSMO data. Data must be present in all scans
    # to consider the volume valid
    filename_list = list()
    for i in range(len(cfg['ScanList'])):
        filename = find_rad4alpcosmo_file(
            voltime, datatype, cfg, cfg['ScanList'][i])
        if filename is None:
            return None
        filename_list.append(filename)

    # create the radar object where to store the data
    # taking as reference the metranet polar file
    dayinfo = voltime.strftime('%y%j')
    timeinfo = voltime.strftime('%H%M')
    basename = 'P'+cfg['RadarRes']+cfg['RadarName']+dayinfo
    datapath = cfg['datapath']+dayinfo+'/'+basename+'/'
    filename = glob.glob(datapath+basename+timeinfo+'*.'+cfg['ScanList'][0])

    radar = get_data_rad4alp(filename[0], ['dBZ'], cfg['ScanList'][0], cfg)
    radar.fields = dict()

    cosmo_field = read_rad4alp_cosmo(filename_list[0], datatype)
    if cosmo_field is None:
        return None

    radar.add_field(get_fieldname_pyart(datatype), cosmo_field)

    # add the other scans
    for i in range(1, len(cfg['ScanList'])):
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+cfg['ScanList'][i])
        radar_aux = get_data_rad4alp(
            filename[0], ['dBZ'], cfg['ScanList'][i], cfg)
        radar_aux.fields = dict()

        cosmo_field = read_rad4alp_cosmo(filename_list[i], datatype)
        if cosmo_field is None:
            return None

        radar_aux.add_field(get_fieldname_pyart(datatype), cosmo_field)

        radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_dem_rad4alp(voltime, datatype, cfg):
    """
    merge cosmo rad4alp scans. If data for all the scans cannot be retrieved
    returns None

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype : str
        name of the data type to read
    cfg : dict
        configuration dictionary

    Returns
    -------
    radar : Radar
        radar object

    """
    # read visibility data file
    vis_list = read_rad4alp_vis(
        cfg['dempath']+cfg['RadarName']+'_visib_volume_40', datatype)
    if vis_list is None:
        return None

    # create the radar object where to store the data
    # taking as reference the metranet polar file
    dayinfo = voltime.strftime('%y%j')
    timeinfo = voltime.strftime('%H%M')
    basename = 'P'+cfg['RadarRes']+cfg['RadarName']+dayinfo
    datapath = cfg['datapath']+dayinfo+'/'+basename+'/'
    filename = glob.glob(datapath+basename+timeinfo+'*.'+cfg['ScanList'][0])

    radar = get_data_rad4alp(filename[0], ['dBZ'], cfg['ScanList'][0], cfg)
    radar.fields = dict()

    # add visibility data for first scan
    radar.add_field(
        get_fieldname_pyart(datatype), vis_list[int(cfg['ScanList'][0])-1])

    # add the other scans
    for i in range(1, len(cfg['ScanList'])):
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+cfg['ScanList'][i])
        radar_aux = get_data_rad4alp(
            filename[0], ['dBZ'], cfg['ScanList'][i], cfg)
        radar_aux.fields = dict()

        radar_aux.add_field(
            get_fieldname_pyart(datatype),
            vis_list[int(cfg['ScanList'][i])-1])

        radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_fields_rainbow(basepath, scan_name, voltime, datatype_list):
    """
    merge Rainbow fields into a single radar object.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    scan_name: str
        name of the scan
    voltime : datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get

    Returns
    -------
    radar : Radar
        radar object

    """
    datapath = basepath+scan_name+voltime.strftime('%Y-%m-%d')+'/'
    fdatetime = voltime.strftime('%Y%m%d%H%M%S')+'00'

    if (datatype_list[0] != 'Nh') and (datatype_list[0] != 'Nv'):
        filename = glob.glob(datapath+fdatetime+datatype_list[0]+'.*')
    elif datatype_list[0] == 'Nh':
        filename = glob.glob(datapath+fdatetime+'dBZ.*')
    else:
        filename = glob.glob(datapath+fdatetime+'dBZv.*')

    # create radar object
    radar = get_data_rainbow(filename[0], datatype_list[0])

    # add other fields in the same scan
    for i in range(1, len(datatype_list)):
        if (datatype_list[i] != 'Nh') and (datatype_list[i] != 'Nv'):
            filename = glob.glob(datapath+fdatetime+datatype_list[i]+'.*')
        elif datatype_list[i] == 'Nh':
            filename = glob.glob(datapath+fdatetime+'dBZ.*')
        else:
            filename = glob.glob(datapath+fdatetime+'dBZv.*')

        radar_aux = get_data_rainbow(filename[0], datatype_list[i])

        for field_name in radar_aux.fields.keys():
            break

        try:
            radar.add_field(field_name, radar_aux.fields[field_name])
        except (ValueError, KeyError):
            warn('Unable to add field '+field_name+' to radar object')

    return radar


def merge_fields_dem(basepath, scan_name, datatype_list):
    """
    merge DEM fields into a single radar object.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    scan_name: str
        name of the scan
    datatype_list : list
        lists of data types to get

    Returns
    -------
    radar : Radar
        radar object

    """
    datapath = basepath+datatype_list[0]+'/'+scan_name
    scan_name_aux = scan_name.partition('/')[0]
    filename = glob.glob(datapath+datatype_list[0]+'_'+scan_name_aux)

    # create radar object
    radar = get_data_rainbow(filename[0], datatype_list[0])

    # add other fields in the same scan
    for i in range(1, len(datatype_list)):
        datapath = basepath+datatype_list[i]+'/'+scan_name+'/'
        filename = glob.glob(datapath+datatype_list[i]+'_'+scan_name_aux)
        radar_aux = get_data_rainbow(filename[0], datatype_list[i])

        for field_name in radar_aux.fields.keys():
            break

        try:
            radar.add_field(field_name, radar_aux.fields[field_name])
        except (ValueError, KeyError):
            warn('Unable to add field '+field_name+' to radar object')

    return radar


def merge_fields_cosmo(filename_list):
    """
    merge COSMO fields in Rainbow file format

    Parameters
    ----------
    filename_list : str
        list of file paths where to find the data

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = pyart.aux_io.read_rainbow_wrl(filename_list[0])

    # add other COSMO fields in the same scan
    for i in range(1, len(filename_list)):
        radar_aux = pyart.aux_io.read_rainbow_wrl(filename_list[i])
        for field_name in radar_aux.fields.keys():
            break
        radar.add_field(field_name, radar_aux.fields[field_name])

    return radar


def get_data_rainbow(filename, datatype):
    """
    gets rainbow radar data

    Parameters
    ----------
    filename : str
        name of file containing rainbow data
    datatype : str
        field name

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = pyart.aux_io.read_rainbow_wrl(filename)
    if (datatype == 'Nh') or (datatype == 'Nv'):
        rbf = wrl.io.read_Rainbow(filename, loaddata=False)
        # check the number of slices
        nslices = int(rbf['volume']['scan']['pargroup']['numele'])
        if nslices > 1:
            single_slice = False
            common_slice_info = rbf['volume']['scan']['slice'][0]
        else:
            single_slice = True
            common_slice_info = rbf['volume']['scan']['slice']

        if datatype[0] == 'Nh':
            noisedBZ1km_h = float(common_slice_info['noise_power_dbz'])
            noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                radar.nrays, noisedBZ1km_h, radar.range['data'], 1.,
                noise_field='noisedBZ_hh')
            radar.fields = dict()
            radar.add_field('noisedBZ_hh', noisedBZ_h)
        else:
            noisedBZ1km_v = float(common_slice_info['noise_power_dbz_dpv'])
            noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                radar.nrays, noisedBZ1km_v, radar.range['data'], 1.,
                noise_field='noisedBZ_vv')
            radar.fields = dict()
            radar.add_field('noisedBZ_vv', noisedBZ_v)

    return radar


def get_data_rad4alp(filename, datatype_list, scan_name, cfg):
    """
    gets rad4alp radar data

    Parameters
    ----------
    filename : str
        name of file containing rainbow data
    datatype_list : list of strings
        list of data fields to get
    scan_name : str
        name of the elevation (001 to 020)
    cfg : dict
        configuration dictionary

    Returns
    -------
    radar : Radar
        radar object

    """
    metranet_field_names = dict()
    for datatype in datatype_list:
        if (datatype != 'Nh') and (datatype != 'Nv'):
            metranet_field_names.update(get_datatype_metranet(datatype))

    radar = pyart.aux_io.read_metranet(
        filename, field_names=metranet_field_names)

    # create secondary moments
    if ('Nh' in datatype_list) or ('Nv' in datatype_list):
        # read radar information in status file
        root = read_status(voltime, cfg)
        sweep_number = int(scan_name)-1

        if 'Nh' in datatype_list:
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

        if 'Nv' in datatype_list:
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

    return radar


def interpol_field(radar_dest, radar_orig, field_name):
    """
    interpolates field field_name contained in radar_orig to the grid in
    radar_dest

    Parameters
    ----------
    radar_dest : radar object
        the destination radar
    radar_orig : radar object
        the radar object containing the original field
    field_name: str
        name of the field to interpolate

    Returns
    -------
    field_dest : dict
        interpolated field and metadata

    """
    fill_value = pyart.config.get_fillvalue()
    field_orig_data = radar_orig.fields[field_name]['data'].filled(
        fill_value=fill_value)
    field_dest = radar_orig.fields[field_name]
    field_dest['data'] = np.ma.empty((radar_dest.nrays, radar_dest.ngates))
    field_dest['data'][:] = np.ma.masked

    for sweep in range(radar_dest.nsweeps):
        sweep_start_orig = radar_orig.sweep_start_ray_index['data'][sweep]
        sweep_end_orig = radar_orig.sweep_end_ray_index['data'][sweep]

        sweep_start_dest = radar_dest.sweep_start_ray_index['data'][sweep]
        sweep_end_dest = radar_dest.sweep_end_ray_index['data'][sweep]

        if radar_dest.scan_type == 'ppi':
            angle_old = radar_orig.azimuth['data'][
                sweep_start_orig:sweep_end_orig+1]
            angle_new = radar_dest.azimuth['data'][
                sweep_start_dest:sweep_end_dest+1]
        elif radar_dest.scan_type == 'rhi':
            angle_old = radar_orig.elevation['data'][
                sweep_start_orig:sweep_end_orig+1]
            angle_new = radar_dest.elevation['data'][
                sweep_start_dest:sweep_end_dest+1]

        interpol_func = RegularGridInterpolator(
            (angle_old, radar_orig.range['data']),
            field_orig_data[sweep_start_orig:sweep_end_orig+1, :],
            method='nearest', bounds_error=False, fill_value=fill_value)

        # interpolate data to radar_dest grid
        angv, rngv = np.meshgrid(
            angle_new, radar_dest.range['data'], indexing='ij')

        field_dest_sweep = interpol_func((angv, rngv))
        field_dest_sweep = np.ma.masked_where(
            field_dest_sweep == fill_value, field_dest_sweep)

        field_dest['data'][sweep_start_dest:sweep_end_dest+1, :] = (
            field_dest_sweep)

    return field_dest
