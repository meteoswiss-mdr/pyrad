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
    merge_scans_hydro_rad4alp
    merge_fields_rainbow
    merge_fields_cfradial
    merge_fields_dem
    merge_fields_cosmo
    get_data_rainbow
    get_data_rad4alp
    add_field
    interpol_field

"""

import sys
import glob
import datetime
import os
from warnings import warn
from copy import deepcopy

import numpy as np
from scipy.interpolate import RegularGridInterpolator

import pyart
try:
    import wradlib as wrl
    _WRADLIB_AVAILABLE = True
except Exception:
    _WRADLIB_AVAILABLE = False

from .read_data_other import read_status, read_rad4alp_cosmo, read_rad4alp_vis
from .read_data_mxpol import pyrad_MXPOL, pyrad_MCH

from .io_aux import get_datatype_metranet, get_fieldname_pyart, get_file_list
from .io_aux import get_datatype_fields, get_datetime, map_hydro
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
    datatype_cfradial = list()
    dataset_cfradial = list()
    product_cfradial = list()
    datatype_cosmo = list()
    datatype_rad4alpcosmo = list()
    datatype_dem = list()
    datatype_rad4alpdem = list()
    datatype_rad4alphydro = list()
    datatype_mxpol = list()
    for datatypedescr in datatypesdescr:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datagroup == 'RAINBOW':
            datatype_rainbow.append(datatype)
        elif datagroup == 'RAD4ALP':
            datatype_rad4alp.append(datatype)
        elif datagroup == 'CFRADIAL':
            datatype_cfradial.append(datatype)
            dataset_cfradial.append(dataset)
            product_cfradial.append(product)
        elif datagroup == 'COSMO':
            datatype_cosmo.append(datatype)
        elif datagroup == 'RAD4ALPCOSMO':
            datatype_rad4alpcosmo.append(datatype)
        elif datagroup == 'DEM':
            datatype_dem.append(datatype)
        elif datagroup == 'RAD4ALPDEM':
            datatype_rad4alpdem.append(datatype)
        elif datagroup == 'RAD4ALPHYDRO':
            datatype_rad4alphydro.append(datatype)
        elif datagroup == 'MXPOL':
            datatype_mxpol.append(datatype)

    ind_rad = int(radarnr[5:8])-1

    ndatatypes_rainbow = len(datatype_rainbow)
    ndatatypes_rad4alp = len(datatype_rad4alp)
    ndatatypes_cfradial = len(datatype_cfradial)
    ndatatypes_cosmo = len(datatype_cosmo)
    ndatatypes_rad4alpcosmo = len(datatype_rad4alpcosmo)
    ndatatypes_dem = len(datatype_dem)
    ndatatypes_rad4alpdem = len(datatype_rad4alpdem)
    ndatatypes_rad4alphydro = len(datatype_rad4alphydro)
    ndatatypes_mxpol = len(datatype_mxpol)

    radar = None
    if ndatatypes_rainbow > 0 and _WRADLIB_AVAILABLE:
        radar = merge_scans_rainbow(
            cfg['datapath'][ind_rad], cfg['ScanList'][ind_rad], voltime,
            cfg['ScanPeriod'], datatype_rainbow, cfg, radarnr=radarnr)

    elif ndatatypes_rad4alp > 0:
        radar = merge_scans_rad4alp(
            cfg['datapath'][ind_rad], cfg['ScanList'][ind_rad],
            cfg['RadarName'][ind_rad], cfg['RadarRes'][ind_rad], voltime,
            datatype_rad4alp, cfg, ind_rad=ind_rad)

    if ndatatypes_cfradial > 0:
        radar_aux = merge_fields_cfradial(
            cfg['loadbasepath'][ind_rad], cfg['loadname'][ind_rad], voltime,
            datatype_cfradial, dataset_cfradial, product_cfradial)
        radar = add_field(radar, radar_aux)

    if ndatatypes_mxpol > 0:
        radar = merge_scans_mxpol(
            cfg['datapath'][ind_rad], cfg['ScanList'][ind_rad], voltime,
            datatype_mxpol, cfg, ind_rad=ind_rad)

    # add COSMO files to the radar field
    if ndatatypes_cosmo > 0 and _WRADLIB_AVAILABLE:
        radar_aux = merge_scans_cosmo(
            voltime, datatype_cosmo, cfg, ind_rad=ind_rad)
        radar = add_field(radar, radar_aux)

    elif ndatatypes_rad4alpcosmo > 0:
        if ((cfg['RadarRes'][ind_rad] is None) or
                (cfg['RadarName'][ind_rad] is None)):
            raise ValueError(
                'ERROR: Radar Name and Resolution ' +
                'not specified in config file. ' +
                'Unable to load rad4alp COSMO data')

        for i in range(ndatatypes_rad4alpcosmo):
            radar_aux = merge_scans_cosmo_rad4alp(
                voltime, datatype_rad4alpcosmo[i], cfg, ind_rad=ind_rad)
            if radar is None:
                radar = radar_aux
            else:
                if radar_aux is not None:
                    for field_name in radar_aux.fields.keys():
                        radar.add_field(
                            field_name, radar_aux.fields[field_name])

    # add DEM files to the radar field
    if ndatatypes_dem > 0 and _WRADLIB_AVAILABLE:
        radar_aux = merge_scans_dem(
            cfg['dempath'][ind_rad], cfg['ScanList'][ind_rad], datatype_dem)
        radar = add_field(radar, radar_aux)

    elif ndatatypes_rad4alpdem > 0:
        if ((cfg['RadarRes'][ind_rad] is None) or
                (cfg['RadarName'][ind_rad] is None)):
            raise ValueError(
                'ERROR: Radar Name and Resolution ' +
                'not specified in config file. ' +
                'Unable to load rad4alp DEM data')

        if cfg['RadarRes'][ind_rad] != 'L':
            raise ValueError(
                'ERROR: DEM files only available for rad4alp PL data. ' +
                'Current radar '+cfg['RadarName'][ind_rad] +
                cfg['RadarRes'][ind_rad])

        for i in range(ndatatypes_rad4alpdem):
            radar_aux = merge_scans_dem_rad4alp(
                voltime, datatype_rad4alpdem[i], cfg, ind_rad=ind_rad)
            if radar is None:
                radar = radar_aux
            else:
                if radar_aux is not None:
                    for field_name in radar_aux.fields.keys():
                        radar.add_field(
                            field_name, radar_aux.fields[field_name])

    if ndatatypes_rad4alphydro > 0:
        if ((cfg['RadarRes'][ind_rad] is None) or
                (cfg['RadarName'][ind_rad] is None)):
            raise ValueError(
                'ERROR: Radar Name and Resolution ' +
                'not specified in config file. ' +
                'Unable to load rad4alp hydro data')

        for i in range(ndatatypes_rad4alphydro):
            radar_aux = merge_scans_hydro_rad4alp(
                voltime, datatype_rad4alphydro[i], cfg, ind_rad=ind_rad)
            if radar is None:
                radar = radar_aux
            else:
                if radar_aux is not None:
                    for field_name in radar_aux.fields.keys():
                        radar.add_field(
                            field_name, radar_aux.fields[field_name])

    # if it is specified, get the position from the config file
    if 'RadarPosition' in cfg:
        if 'latitude' in cfg['RadarPosition']:
            radar.latitude['data'][0] = (
                cfg['RadarPosition']['latitude'][ind_rad])
        if 'longitude' in cfg['RadarPosition']:
            radar.longitude['data'][0] = (
                cfg['RadarPosition']['longitude'][ind_rad])
        if 'altitude' in cfg['RadarPosition']:
            radar.altitude['data'][0] = (
                cfg['RadarPosition']['altitude'][ind_rad])
        radar.init_gate_longitude_latitude()
        radar.init_gate_altitude()

    return radar


def merge_scans_rainbow(basepath, scan_list, voltime, scan_period,
                        datatype_list, cfg, radarnr='RADAR001'):
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
    radarnr : str
        radar identifier number

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
            datadescriptor = radarnr+':RAINBOW:dBZ'
        else:
            datadescriptor = radarnr+':RAINBOW:'+datatype_list[0]
        endtime = voltime+datetime.timedelta(minutes=scan_period)
        for i in range(1, nscans):
            filelist = get_file_list(datadescriptor, voltime, endtime,
                                     cfg, scan=scan_list[i])

            if (len(filelist) == 0):
                print("ERROR: No data file found for scan '%s' "
                      "between %s and %s" % (scan_list[i], voltime, endtime),
                      file=sys.stderr)
                continue
            scantime = get_datetime(filelist[0], datadescriptor)

            radar_aux = merge_fields_rainbow(
                basepath, scan_list[i], scantime, datatype_list)

            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_dem(basepath, scan_list, datatype_list, radarnr='RADAR001'):
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
    radarnr : str
        radar identifier number

    Returns
    -------
    radar : Radar
        radar object

    """
    ind_rad = int(radarnr[5:8])-1
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
                        datatype_list, cfg, ind_rad=0):
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
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    if (radar_name is None) or (radar_res is None):
        raise ValueError(
            'ERROR: Radar Name and Resolution not specified in config file.' +
            ' Unable to load rad4alp data')

    radar = None
    dayinfo = voltime.strftime('%y%j')
    timeinfo = voltime.strftime('%H%M')
    basename = 'M'+radar_res+radar_name+dayinfo
    if cfg['path_convention'] == 'LTE':
        yy = dayinfo[0:2]
        dy = dayinfo[2:]
        subf = 'M'+radar_res+radar_name+yy+'hdf'+dy
        datapath = basepath+subf+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            subf = 'P'+radar_res+radar_name+yy+'hdf'+dy
            datapath = basepath+subf+'/'
    else:
        datapath = basepath+dayinfo+'/'+basename+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+dayinfo+'/'+basename+'/'
    filename = glob.glob(datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
    if not filename:
        warn('No file found in '+datapath+basename+timeinfo+'*.'+scan_list[0])
    else:
        radar = get_data_rad4alp(
            filename[0], datatype_list, scan_list[0], cfg, ind_rad=ind_rad)

    nelevs = len(scan_list)
    # merge the elevations into a single radar instance
    for i in range(1, nelevs):
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[i]+'*')
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.' +
                 scan_list[i])
        else:
            radar_aux = get_data_rad4alp(
                filename[0], datatype_list, scan_list[i], cfg, ind_rad=ind_rad)

            if radar is None:
                radar = radar_aux
            else:
                radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_mxpol(basepath, scan_list, voltime, datatype_list, cfg,
                      ind_rad=0):
    """
    merge rad4alp data.

    Parameters
    ----------
    basepath : str
        base path of mxpol radar data
    scan_list : list
        list of scans, in the case of mxpol, the elevation or azimuth denoted
        as 005 or 090 (for 5 or 90 degrees elevation) or 330 (for 330 degrees
        azimuth respectively)
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = None
    if cfg['path_convention'] == 'LTE':
        sub1 = str(voltime.year)
        sub2 = voltime.strftime('%m')
        sub3 = voltime.strftime('%d')
        dayinfo = voltime.strftime('%Y%m%d')
        timeinfo = voltime.strftime('%H%M')
        datapath = cfg['datapath'][ind_rad]+'/'+sub1+'/'+sub2+'/'+sub3+'/'
        scanname = 'MXPol-polar-'+dayinfo+'-'+timeinfo+'*-'
        filename = glob.glob(datapath+scanname+scan_list[0]+'*')
    else:
        daydir = voltime.strftime('%Y-%m-%d')
        dayinfo = voltime.strftime('%Y%m%d')
        timeinfo = voltime.strftime('%H%M')
        datapath = cfg['datapath'][ind_rad]+scan_list[0]+'/'+daydir+'/'
        if (not os.path.isdir(datapath)):
            warn("WARNING: Unknown datapath '%s'" % datapath)
            return None
        filename = glob.glob(
            datapath+'MXPol-polar-'+dayinfo+'-'+timeinfo+'*-' +
            scan_list[0]+'.nc')
    if not filename:
        warn('No file found matching '+datapath+scanname+scan_list[0]+'*')
    else:
        radar = get_data_mxpol(
            filename[0], datatype_list, scan_list[0], cfg, ind_rad=ind_rad)

    nelevs = len(scan_list)
    # merge the elevations into a single radar instance
    for i in range(1, nelevs):
        if cfg['path_convention'] == 'LTE':
            sub1 = str(voltime.year)
            sub2 = voltime.strftime('%m')
            sub3 = voltime.strftime('%d')
            dayinfo = voltime.strftime('%Y%m%d')
            timeinfo = voltime.strftime('%H%M')
            datapath = cfg['datapath'][ind_rad]+'/'+sub1+'/'+sub2+'/'+sub3+'/'
            scanname = 'MXPol-polar-'+dayinfo+'-'+timeinfo+'*-'
            filename = glob.glob(datapath+scanname+scan_list[i]+'*')
        else:
            daydir = voltime.strftime('%Y-%m-%d')
            dayinfo = voltime.strftime('%Y%m%d')
            timeinfo = voltime.strftime('%H%M')
            datapath = cfg['datapath'][ind_rad]+scan_list[i]+'/'+daydir+'/'
            if (not os.path.isdir(datapath)):
                warn("WARNING: Unknown datapath '%s'" % datapath)
                return None
            filename = glob.glob(
                datapath+'MXPol-polar-'+dayinfo+'-'+timeinfo+'*-' +
                scan_list[i]+'.nc')
        if not filename:
            warn('No file found in '+datapath+scanname+scan_list[i])
        else:
            radar_aux = get_data_mxpol(
                filename[0], datatype_list, scan_list[i], cfg, ind_rad=ind_rad)

            if radar is None:
                radar = radar_aux
            else:
                radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_cosmo(voltime, datatype_list, cfg, ind_rad=0):
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
    ind_rad : int
        radar index

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
            voltime, datatype_list[i], cfg, cfg['ScanList'][ind_rad][0],
            ind_rad=ind_rad)
        if filename is not None:
            filename_list.append(filename)

    nfiles_valid = len(filename_list)

    if nfiles_valid > 0:
        radar = merge_fields_cosmo(filename_list)

    # merge scans into a single radar instance
    nscans = len(cfg['ScanList'][ind_rad])
    if nscans > 1:
        endtime = voltime+datetime.timedelta(minutes=cfg['ScanPeriod'])
        for i in range(1, nscans):
            filename_list = list()
            for j in range(ndatatypes):
                filename = find_cosmo_file(
                    voltime, datatype_list[j], cfg,
                    cfg['ScanList'][ind_rad][i], ind_rad=ind_rad)
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


def merge_scans_cosmo_rad4alp(voltime, datatype, cfg, ind_rad=0):
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
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    # look for rad4alp COSMO data. Data must be present in all scans
    # to consider the volume valid
    filename_list = list()
    for i in range(len(cfg['ScanList'][ind_rad])):
        filename = find_rad4alpcosmo_file(
            voltime, datatype, cfg, cfg['ScanList'][ind_rad][i])
        if filename is None:
            return None
        filename_list.append(filename)

    # create the radar object where to store the data
    # taking as reference the metranet polar file
    dayinfo = voltime.strftime('%y%j')
    timeinfo = voltime.strftime('%H%M')
    basename = 'P'+cfg['RadarRes'][ind_rad]+cfg['RadarName'][ind_rad]+dayinfo
    datapath = cfg['datapath'][ind_rad]+dayinfo+'/'+basename+'/'
    filename = glob.glob(
        datapath+basename+timeinfo+'*.'+cfg['ScanList'][ind_rad][0])

    radar = None
    if not filename:
        warn('No file found in '+datapath+basename+timeinfo+'*.' +
             cfg['ScanList'][ind_rad][0])
    else:
        radar = get_data_rad4alp(
            filename[0], ['dBZ'], cfg['ScanList'][ind_rad][0], cfg,
            ind_rad=ind_rad)
        radar.fields = dict()

        cosmo_field = read_rad4alp_cosmo(filename_list[0], datatype)
        if cosmo_field is None:
            return None

        radar.add_field(get_fieldname_pyart(datatype), cosmo_field)

    # add the other scans
    for i in range(1, len(cfg['ScanList'][ind_rad])):
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+cfg['ScanList'][ind_rad][i])
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.' +
                 cfg['ScanList'][ind_rad][i])
        else:
            radar_aux = get_data_rad4alp(
                filename[0], ['dBZ'], cfg['ScanList'][ind_rad][i], cfg,
                ind_rad=ind_rad)
            radar_aux.fields = dict()

            cosmo_field = read_rad4alp_cosmo(filename_list[i], datatype)
            if cosmo_field is None:
                return None

            radar_aux.add_field(get_fieldname_pyart(datatype), cosmo_field)

            if radar is None:
                radar = radar_aux
            else:
                radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_dem_rad4alp(voltime, datatype, cfg, ind_rad=0):
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
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    # read visibility data file
    vis_list = read_rad4alp_vis(
        cfg['dempath'][ind_rad]+cfg['RadarName'][ind_rad]+'_visib_volume_40',
        datatype)
    if vis_list is None:
        return None

    # create the radar object where to store the data
    # taking as reference the metranet polar file
    radar = None
    dayinfo = voltime.strftime('%y%j')
    timeinfo = voltime.strftime('%H%M')
    radar_res = cfg['RadarRes'][ind_rad]
    radar_name = cfg['RadarName'][ind_rad]
    basepath = cfg['datapath'][ind_rad]
    scan_list = cfg['ScanList'][ind_rad]

    basename = 'M'+radar_res+radar_name+dayinfo
    if cfg['path_convention'] == 'LTE':
        yy = dayinfo[0:2]
        dy = dayinfo[2:]
        subf = 'M'+radar_res+radar_name+yy+'hdf'+dy
        datapath = basepath+subf+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            subf = 'P'+radar_res+radar_name+yy+'hdf'+dy
            datapath = basepath+subf+'/'
    else:
        datapath = basepath+dayinfo+'/'+basename+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+dayinfo+'/'+basename+'/'
    filename = glob.glob(datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
    if not filename:
        warn('No file found in '+datapath+basename+timeinfo+'*.'+scan_list[0])
    else:
        radar = get_data_rad4alp(
            filename[0], ['dBZ'], scan_list[0], cfg, ind_rad=ind_rad)
        radar.fields = dict()

        # add visibility data for first scan
        radar.add_field(
            get_fieldname_pyart(datatype), vis_list[int(scan_list[0])-1])

    # add the other scans
    for i in range(1, len(scan_list)):
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[i])
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.' +
                 scan_list[i])
            continue
        radar_aux = get_data_rad4alp(
            filename[0], ['dBZ'], scan_list[i], cfg, ind_rad=ind_rad)
        radar_aux.fields = dict()

        radar_aux.add_field(
            get_fieldname_pyart(datatype), vis_list[int(scan_list[i])-1])
        if radar is None:
            radar = radar_aux
        else:
            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_hydro_rad4alp(voltime, datatype, cfg, ind_rad=0):
    """
    merge rad4alp hydrometeor classification scans. If data for all the scans
    cannot be retrieved returns None

    Parameters
    ----------
    voltime: datetime object
        reference time of the scan
    datatype : str
        name of the data type to read
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    radar = None
    dayinfo = voltime.strftime('%y%j')
    timeinfo = voltime.strftime('%H%M')
    radar_res = cfg['RadarRes'][ind_rad]
    radar_name = cfg['RadarName'][ind_rad]
    basepath = cfg['datapath'][ind_rad]
    scan_list = cfg['ScanList'][ind_rad]

    # read hydrometeor classification data file for first scan
    basename_hydro = 'YM'+radar_name+dayinfo
    datapath_hydro = basepath+dayinfo+'/'+basename_hydro+'/'
    filename_hydro = glob.glob(datapath_hydro+basename_hydro+timeinfo+'*.' +
                               str(800+int(scan_list[0]))+'*')
    if not filename_hydro:
        warn('No file found in '+datapath_hydro+basename_hydro+timeinfo+'*.' +
             str(800+int(scan_list[0])))
        return None

    hydro_obj = pyart.aux_io.read_product(
        filename_hydro[0], physic_value=False, masked_array=True)

    if hydro_obj is None:
        warn('Unable to read file '+filename_hydro)
        return None

    hydro_field = get_fieldname_pyart(datatype)
    hydro_dict = pyart.config.get_metadata(hydro_field)
    hydro_dict['data'] = map_hydro(hydro_obj.data)

    # create the radar object where to store the data
    # taking as reference the metranet polar file
    basename = 'M'+radar_res+radar_name+dayinfo
    if cfg['path_convention'] == 'LTE':
        yy = dayinfo[0:2]
        dy = dayinfo[2:]
        subf = 'M'+radar_res+radar_name+yy+'hdf'+dy
        datapath = basepath+subf+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            subf = 'P'+radar_res+radar_name+yy+'hdf'+dy
            datapath = basepath+subf+'/'
    else:
        datapath = basepath+dayinfo+'/'+basename+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+dayinfo+'/'+basename+'/'
    filename = glob.glob(datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
    if not filename:
        warn('No file found in '+datapath+basename+timeinfo+'*.'+scan_list[0])
    else:
        radar = get_data_rad4alp(
            filename[0], ['dBZ'], scan_list[0], cfg, ind_rad=ind_rad)
        print(np.shape(radar.fields['reflectivity']['data']))
        radar.fields = dict()

        # add hydrometeor classification data for first scan
        print(np.shape(hydro_dict['data']))

        radar.add_field(hydro_field, hydro_dict)

    # add the other scans
    for i in range(1, len(scan_list)):
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[i]+'*')
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.' +
                 scan_list[i])
            continue

        radar_aux = get_data_rad4alp(
            filename[0], ['dBZ'], scan_list[i], cfg, ind_rad=ind_rad)
        radar_aux.fields = dict()

        # read hydrometeor classification data file for other scans
        filename_hydro = glob.glob(datapath_hydro+basename_hydro+timeinfo +
                                   '*.'+str(800+int(scan_list[i]))+'*')
        if not filename_hydro:
            warn('No file found in '+datapath_hydro+basename_hydro+timeinfo +
                 '*.'+str(800+int(scan_list[i])))
            continue

        hydro_obj = pyart.aux_io.read_product(
            filename_hydro[0], physic_value=False, masked_array=True)

        if hydro_obj is None:
            warn('Unable to read file '+filename_hydro)
            continue

        hydro_field = get_fieldname_pyart(datatype)
        hydro_dict = pyart.config.get_metadata(hydro_field)
        hydro_dict['data'] = map_hydro(hydro_obj.data)

        radar_aux.add_field(hydro_field, hydro_dict)
        if radar is None:
            radar = radar_aux
        else:
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
    radar = None
    if not filename:
        warn('No file found in '+datapath+fdatetime+datatype_list[0]+'.*')
    else:
        radar = get_data_rainbow(filename[0], datatype_list[0])

    # add other fields in the same scan
    for i in range(1, len(datatype_list)):
        if (datatype_list[i] != 'Nh') and (datatype_list[i] != 'Nv'):
            filename = glob.glob(datapath+fdatetime+datatype_list[i]+'.*')
        elif datatype_list[i] == 'Nh':
            filename = glob.glob(datapath+fdatetime+'dBZ.*')
        else:
            filename = glob.glob(datapath+fdatetime+'dBZv.*')
        if not filename:
            warn('No file found in '+datapath+fdatetime+datatype_list[i]+'.*')
        else:
            radar_aux = get_data_rainbow(filename[0], datatype_list[i])
            if radar is None:
                radar = radar_aux
            else:
                for field_name in radar_aux.fields.keys():
                    break
                try:
                    radar.add_field(field_name, radar_aux.fields[field_name])
                except (ValueError, KeyError) as ee:
                    warn("Unable to add field '"+field_name+"' to radar object"
                         ": (%s)" % str(ee))

    return radar


def merge_fields_cfradial(basepath, loadname, voltime, datatype_list,
                          dataset_list, product_list):
    """
    merge CF/Radial fields into a single radar object.

    Parameters
    ----------
    basepath : str
        name of the base path where to find the data
    loadname: str
        name of the saving directory
    voltime : datetime object
        reference time of the scan
    datatype_list : list
        list of data types to get
    dataset_list : list
        list of datasets that produced the data type to get.
        Used to get path.
    product_list : list
        list of products. Used to get path

    Returns
    -------
    radar : Radar
        radar object

    """
    datapath = (basepath+loadname+'/'+voltime.strftime('%Y-%m-%d')+'/' +
                dataset_list[0]+'/'+product_list[0]+'/')
    fdatetime = voltime.strftime('%Y%m%d%H%M%S')
    filename = glob.glob(datapath+fdatetime+'*'+datatype_list[0]+'.nc')

    # create radar object
    radar = None
    if not filename:
        warn('No file found in '+datapath+fdatetime+'*'+datatype_list[0] +
             '.nc')
    else:
        radar = pyart.io.read_cfradial(filename[0])

    # add other fields in the same scan
    for i in range(1, len(datatype_list)):
        datapath = (basepath+loadname+'/'+voltime.strftime('%Y-%m-%d')+'/' +
                    dataset_list[i]+'/'+product_list[i]+'/')
        filename = glob.glob(datapath+fdatetime+'*'+datatype_list[i]+'.nc')
        if not filename:
            warn('No file found in '+datapath+fdatetime+'*'+datatype_list[0] +
                 '.nc')
        else:
            radar_aux = pyart.io.read_cfradial(filename[0])
            if radar is None:
                radar = radar_aux
            else:
                add_field(radar, radar_aux)

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
    radar = None
    if not filename:
        warn('No file found in '+datapath+datatype_list[0]+'_'+scan_name_aux)
    else:
        # create radar object
        radar = get_data_rainbow(filename[0], datatype_list[0])

    # add other fields in the same scan
    for i in range(1, len(datatype_list)):
        datapath = basepath+datatype_list[i]+'/'+scan_name+'/'
        filename = glob.glob(datapath+datatype_list[i]+'_'+scan_name_aux)
        if not filename:
            warn('No file found in '+datapath+datatype_list[i]+'_' +
                 scan_name_aux)
        else:
            radar_aux = get_data_rainbow(filename[0], datatype_list[i])
            if radar is None:
                radar = radar_aux
            else:
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

        if datatype == 'Nh':
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


def get_data_rad4alp(filename, datatype_list, scan_name, cfg, ind_rad=0):
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
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    metranet_field_names = dict()
    for datatype in datatype_list:
        if (datatype != 'Nh') and (datatype != 'Nv'):
            metranet_field_names.update(get_datatype_metranet(datatype))

    if cfg['path_convention'] == 'LTE':
        radar = pyrad_MCH(filename, field_names=metranet_field_names)
    else:
        radar = pyart.aux_io.read_metranet(
            filename, field_names=metranet_field_names)

    # create secondary moments
    if ('Nh' in datatype_list) or ('Nv' in datatype_list):
        # read radar information in status file
        voltime = get_datetime(filename, 'RAD4ALP:dBZ')
        root = read_status(voltime, cfg, ind_rad=ind_rad)

        if root is not None:
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


def get_data_mxpol(filename, datatype_list, scan_name, cfg, ind_rad=0):
    """
    gets MXPol radar data

    Parameters
    ----------
    filename : str
        name of file containing MXPol data
    datatype_list : list of strings
        list of data fields to get
    scan_name : list
        list of scans, in the case of mxpol, the elevation or azimuth denoted
        as 005 or 090 (for 5 or 90 degrees elevation) or 330 (for 330 degrees
        azimuth respectively)
    cfg : dict
        configuration dictionary
    ind_rad : int
        radar index

    Returns
    -------
    radar : Radar
        radar object

    """
    field_names = dict()
    for datatype in datatype_list:
        if (datatype != 'Nh') and (datatype != 'Nv'):
            field_names.update(get_datatype_metranet(datatype))
    radar = pyrad_MXPOL(filename, field_names=field_names)

    # create secondary moments (TODO)
    if ('Nh' in datatype_list) or ('Nv' in datatype_list):
        pass

    return radar


def add_field(radar_dest, radar_orig):
    """
    adds the fields from orig radar into dest radar. If they are not in the
    same grid, interpolates them to dest grid

    Parameters
    ----------
    radar_dest : radar object
        the destination radar
    radar_orig : radar object
        the radar object containing the original field

    Returns
    -------
    field_dest : dict
        interpolated field and metadata

    """
    if radar_dest is None:
        radar_dest = radar_orig
    else:
        if radar_orig is not None:
            if radar_dest.nrays == radar_orig.nrays:
                if ((np.allclose(
                        radar_dest.azimuth['data'],
                        radar_orig.azimuth['data'],
                        atol=0.5, equal_nan=True)) and
                        (np.allclose(
                            radar_dest.elevation['data'],
                            radar_orig.elevation['data'],
                            atol=0.5, equal_nan=True))):
                    for field_name in radar_orig.fields.keys():
                        radar_dest.add_field(
                            field_name, radar_orig.fields[field_name])
                else:
                    for field_name in radar_orig.fields.keys():
                        field_interp = interpol_field(
                            radar_dest, radar_orig, field_name)
                        radar_dest.add_field(field_name, field_interp)
            else:
                for field_name in radar_orig.fields.keys():
                    field_interp = interpol_field(
                        radar_dest, radar_orig, field_name)
                    radar_dest.add_field(field_name, field_interp)

    return radar_dest


def interpol_field(radar_dest, radar_orig, field_name, fill_value=None):
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
    if fill_value is None:
        fill_value = pyart.config.get_fillvalue()

    field_orig_data = radar_orig.fields[field_name]['data'].filled(
        fill_value=fill_value)
    field_dest = deepcopy(radar_orig.fields[field_name])
    field_dest['data'] = np.ma.empty((radar_dest.nrays, radar_dest.ngates))
    field_dest['data'][:] = np.ma.masked

    for sweep in range(radar_dest.nsweeps):
        sweep_start_orig = radar_orig.sweep_start_ray_index['data'][sweep]
        sweep_end_orig = radar_orig.sweep_end_ray_index['data'][sweep]

        sweep_start_dest = radar_dest.sweep_start_ray_index['data'][sweep]
        sweep_end_dest = radar_dest.sweep_end_ray_index['data'][sweep]

        if radar_dest.scan_type == 'ppi':
            angle_old = np.sort(radar_orig.azimuth['data'][
                sweep_start_orig:sweep_end_orig+1])
            ind_ang = np.argsort(radar_orig.azimuth['data'][
                sweep_start_orig:sweep_end_orig+1])
            angle_new = radar_dest.azimuth['data'][
                sweep_start_dest:sweep_end_dest+1]
        elif radar_dest.scan_type == 'rhi':
            angle_old = np.sort(radar_orig.elevation['data'][
                sweep_start_orig:sweep_end_orig+1])
            ind_ang = np.argsort(radar_orig.azimuth['data'][
                sweep_start_orig:sweep_end_orig+1])
            angle_new = radar_dest.elevation['data'][
                sweep_start_dest:sweep_end_dest+1]

        field_orig_sweep_data = field_orig_data[
            sweep_start_orig:sweep_end_orig+1, :]
        interpol_func = RegularGridInterpolator(
            (angle_old, radar_orig.range['data']),
            field_orig_sweep_data[ind_ang], method='nearest',
            bounds_error=False, fill_value=fill_value)

        # interpolate data to radar_dest grid
        angv, rngv = np.meshgrid(
            angle_new, radar_dest.range['data'], indexing='ij')

        field_dest_sweep = interpol_func((angv, rngv))
        field_dest_sweep = np.ma.masked_where(
            field_dest_sweep == fill_value, field_dest_sweep)

        field_dest['data'][sweep_start_dest:sweep_end_dest+1, :] = (
            field_dest_sweep)

    return field_dest
