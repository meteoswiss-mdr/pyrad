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
	merge_scans_odim
    merge_scans_cosmo
    merge_scans_cosmo_rad4alp
    merge_scans_dem_rad4alp
    merge_scans_hydro_rad4alp
    merge_fields_rainbow
    merge_fields_pyrad
    merge_fields_dem
    merge_fields_cosmo
    get_data_rainbow
    get_data_rad4alp
    get_data_odim
    add_field
    interpol_field

"""

import glob
import datetime
import os
from warnings import warn
from copy import deepcopy

import numpy as np
from scipy.interpolate import RegularGridInterpolator

try:
    import wradlib as wrl
    _WRADLIB_AVAILABLE = True
except Exception:
    _WRADLIB_AVAILABLE = False

import pyart

from .read_data_other import read_status, read_rad4alp_cosmo, read_rad4alp_vis
from .read_data_mxpol import pyrad_MXPOL, pyrad_MCH

from .io_aux import get_datatype_metranet, get_fieldname_pyart, get_file_list
from .io_aux import get_datatype_odim, find_date_in_file_name
from .io_aux import get_datatype_fields, get_datetime, map_hydro, map_Doppler
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
        Format : [radarnr]:[datagroup]:[datatype],[dataset],[product]
        'dataset' is only specified for data groups 'ODIM',
        'CFRADIAL' and 'ODIMPYRAD'. 'product' is only specified for data
        groups 'CFRADIAL' and 'ODIMPYRAD'
        The data group specifies the type file from which data is extracted.
        It can be:
            'RAINBOW': Propietary Leonardo format
            'COSMO': COSMO model data saved in Rainbow file format
            'DEM': Visibility data saved in Rainbow file format

            'RAD4ALP': METRANET format used for the operational MeteoSwiss
                data. To find out which datatype to use to match a particular
                METRANET field name check the function 'get_datatype_metranet'
                in pyrad/io/io_aux.py
            'RAD4ALPCOSMO': COSMO model data saved in a binary file format.
                Used by operational MeteoSwiss radars
            'RAD4ALPDEM': Visibility data saved in a binary format used by
                operational MeteoSwiss radars
            'RAD4ALPHYDRO': Used to read the MeteoSwiss operational
                hydrometeor classification
            'RAD4ALPDOPPLER': Used to read the MeteoSwiss operational
                dealiased Doppler velocity

            'ODIM': Generic ODIM file format. For such types 'dataset'
                specifies the directory and file name date convention.
                Example: ODIM:dBZ,D{%Y-%m-%d}-F{%Y%m%d%H%M%S}. To find out
                which datatype to use to match a particular ODIM field name
                check the function 'get_datatype_odim' in pyrad/io/io_aux.py

            'MXPOL': MXPOL (EPFL) data written in a netcdf file

            'CFRADIAL': CFRadial format with the naming convention and
                directory structure in which Pyrad saves the data. For such
                datatypes 'dataset' specifies the directory where the dataset
                is stored and 'product' specifies the directroy where the
                product is stored.
                Example: CFRADIAL:dBZc,Att_ZPhi,SAVEVOL_dBZc
            'ODIMPYRAD': ODIM file format with the naming convention and
                directory structure in which Pyrad saves the data.  For such
                datatypes 'dataset' specifies the directory where the dataset
                is stored and 'product' specifies the directroy where the
                product is stored.
                Example: ODIMPYRAD:dBZc,Att_ZPhi,SAVEVOL_dBZc
        'RAINBOW', 'RAD4ALP', 'ODIM' and 'MXPOL' are primary data file sources
        and they cannot be mixed for the same radar. It is also the case for
        their complementary data files, i.e. 'COSMO' and 'RAD4ALPCOSMO', etc.
        'CFRADIAL' and 'ODIMPYRAD' are secondary data file sources and they
        can be combined with any other datagroup type.
        For a list of accepted datatypes and how they map to the Py-ART name
        convention check function 'get_field_name_pyart' in pyrad/io/io_aux.py
    cfg: dictionary of dictionaries
        configuration info to figure out where the data is

    Returns
    -------
    radar : Radar
        radar object

    """
    datatype_rainbow = list()
    datatype_rad4alp = list()
    datatype_odim = list()
    dataset_odim = list()
    datatype_cfradial = list()
    dataset_cfradial = list()
    product_cfradial = list()
    datatype_odimpyrad = list()
    dataset_odimpyrad = list()
    product_odimpyrad = list()
    datatype_cosmo = list()
    datatype_rad4alpcosmo = list()
    datatype_dem = list()
    datatype_rad4alpdem = list()
    datatype_rad4alphydro = list()
    datatype_rad4alpDoppler = list()
    datatype_mxpol = list()
    for datatypedescr in datatypesdescr:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if datagroup == 'RAINBOW':
            datatype_rainbow.append(datatype)
        elif datagroup == 'RAD4ALP':
            datatype_rad4alp.append(datatype)
        elif datagroup == 'ODIM':
            datatype_odim.append(datatype)
            dataset_odim.append(dataset)
        elif datagroup == 'CFRADIAL':
            datatype_cfradial.append(datatype)
            dataset_cfradial.append(dataset)
            product_cfradial.append(product)
        elif datagroup == 'ODIMPYRAD':
            datatype_odimpyrad.append(datatype)
            dataset_odimpyrad.append(dataset)
            product_odimpyrad.append(product)
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
        elif datagroup == 'RAD4ALPDOPPLER':
            datatype_rad4alpDoppler.append(datatype)
        elif datagroup == 'MXPOL':
            datatype_mxpol.append(datatype)

    ind_rad = int(radarnr[5:8])-1

    ndatatypes_rainbow = len(datatype_rainbow)
    ndatatypes_rad4alp = len(datatype_rad4alp)
    ndatatypes_odim = len(datatype_odim)
    ndatatypes_cfradial = len(datatype_cfradial)
    ndatatypes_odimpyrad = len(datatype_odimpyrad)
    ndatatypes_cosmo = len(datatype_cosmo)
    ndatatypes_rad4alpcosmo = len(datatype_rad4alpcosmo)
    ndatatypes_dem = len(datatype_dem)
    ndatatypes_rad4alpdem = len(datatype_rad4alpdem)
    ndatatypes_rad4alphydro = len(datatype_rad4alphydro)
    ndatatypes_rad4alpDoppler = len(datatype_rad4alpDoppler)
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

    elif ndatatypes_odim > 0:
        try:
            radar_name = cfg['RadarName'][ind_rad]
            radar_res = cfg['RadarRes'][ind_rad]
        except TypeError:
            radar_name = None
            radar_res = None
        radar = merge_scans_odim(
            cfg['datapath'][ind_rad], cfg['ScanList'][ind_rad], radar_name, radar_res,
            voltime, datatype_odim, dataset_odim, cfg, ind_rad=ind_rad)

    elif ndatatypes_mxpol > 0:
        radar = merge_scans_mxpol(
            cfg['datapath'][ind_rad], cfg['ScanList'][ind_rad], voltime,
            datatype_mxpol, cfg)

    if ndatatypes_cfradial > 0:
        radar_aux = merge_fields_pyrad(
            cfg['loadbasepath'][ind_rad], cfg['loadname'][ind_rad], voltime,
            datatype_cfradial, dataset_cfradial, product_cfradial,
            rmax=cfg['rmax'])
        radar = add_field(radar, radar_aux)
    if ndatatypes_odimpyrad > 0:
        radar_aux = merge_fields_pyrad(
            cfg['loadbasepath'][ind_rad], cfg['loadname'][ind_rad], voltime,
            datatype_odimpyrad, dataset_odimpyrad, product_odimpyrad,
            rmax=cfg['rmax'], termination='.h5')
        radar = add_field(radar, radar_aux)

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

        for dt_rad4alpcosmo in datatype_rad4alpcosmo:
            radar_aux = merge_scans_cosmo_rad4alp(
                voltime, dt_rad4alpcosmo, cfg, ind_rad=ind_rad)
            if radar is None:
                radar = radar_aux
                continue
            if radar_aux is None:
                continue
            for field_name in radar_aux.fields.keys():
                try:
                    radar.add_field(
                        field_name, radar_aux.fields[field_name])
                except (ValueError, KeyError) as ee:
                    warn("Unable to add field '"+field_name +
                         "' to radar object"
                         ": (%s)" % str(ee))

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

        for dt_rad4alpdem in datatype_rad4alpdem:
            radar_aux = merge_scans_dem_rad4alp(
                voltime, dt_rad4alpdem, cfg, ind_rad=ind_rad)
            if radar is None:
                radar = radar_aux
                continue
            if radar_aux is None:
                continue
            for field_name in radar_aux.fields.keys():
                try:
                    radar.add_field(
                        field_name, radar_aux.fields[field_name])
                except (ValueError, KeyError) as ee:
                    warn("Unable to add field '"+field_name +
                         "' to radar object"
                         ": (%s)" % str(ee))

    if ndatatypes_rad4alphydro > 0:
        if ((cfg['RadarRes'][ind_rad] is None) or
                (cfg['RadarName'][ind_rad] is None)):
            raise ValueError(
                'ERROR: Radar Name and Resolution ' +
                'not specified in config file. ' +
                'Unable to load rad4alp hydro data')

        for dt_rad4alphydro in datatype_rad4alphydro:
            radar_aux = merge_scans_hydro_rad4alp(
                voltime, dt_rad4alphydro, cfg, ind_rad=ind_rad)
            if radar is None:
                radar = radar_aux
                continue
            if radar_aux is None:
                continue
            if radar_aux is not None:
                for field_name in radar_aux.fields.keys():
                    try:
                        radar.add_field(
                            field_name, radar_aux.fields[field_name])
                    except (ValueError, KeyError) as ee:
                        warn("Unable to add field '"+field_name +
                             "' to radar object"
                             ": (%s)" % str(ee))

    if ndatatypes_rad4alpDoppler > 0:
        if ((cfg['RadarRes'][ind_rad] is None) or
                (cfg['RadarName'][ind_rad] is None)):
            raise ValueError(
                'ERROR: Radar Name and Resolution ' +
                'not specified in config file. ' +
                'Unable to load rad4alp dealiased Doppler data')

        for dt_rad4alpDoppler in datatype_rad4alpDoppler:
            radar_aux = merge_scans_Doppler_rad4alp(
                voltime, dt_rad4alpDoppler, cfg, ind_rad=ind_rad)
            if radar is None:
                radar = radar_aux
                continue
            if radar_aux is None:
                continue
            if radar_aux is not None:
                for field_name in radar_aux.fields.keys():
                    try:
                        radar.add_field(
                            field_name, radar_aux.fields[field_name])
                    except (ValueError, KeyError) as ee:
                        warn("Unable to add field '"+field_name +
                             "' to radar object"
                             ": (%s)" % str(ee))

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
        for scan in scan_list[1:]:
            filelist = get_file_list(datadescriptor, voltime, endtime,
                                     cfg, scan=scan)

            if not filelist:
                warn("ERROR: No data file found for scan '%s' "
                     "between %s and %s" % (scan, voltime, endtime))
                continue
            scantime = get_datetime(filelist[0], datadescriptor)

            radar_aux = merge_fields_rainbow(
                basepath, scan, scantime, datatype_list)

            if radar_aux is None:
                continue

            if radar is None:
                radar = radar_aux
            else:
                radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    if radar is None:
        return radar

    # keep only PPIs within elevation limits
    if cfg['elmin'] != -600. or cfg['elmax'] != 600.:
        if radar.scan_type == 'ppi':
            ind_sweeps = np.where(np.logical_and(
                radar.fixed_angle['data'] > cfg['elmin'],
                radar.fixed_angle['data'] < cfg['elmax']))[0]
            if ind_sweeps.size == 0:
                warn('Elevation angles outside of specified angle range. ' +
                     'Min angle: '+str(cfg['elmin']) +
                     ' Max angle: '+str(cfg['elmax']))
                return None
            radar = radar.extract_sweeps(ind_sweeps)

    return radar


def merge_scans_dem(basepath, scan_list, datatype_list):
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
    radar = merge_fields_dem(
        basepath, scan_list[0], datatype_list)

    # merge scans into a single radar instance
    nscans = len(scan_list)
    if nscans > 1:
        for scan in scan_list[1:]:
            radar_aux = merge_fields_dem(basepath, scan, datatype_list)
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
    elif cfg['path_convention'] == 'MCH':
        datapath = basepath+dayinfo+'/'+basename+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+dayinfo+'/'+basename+'/'
    else:
        datapath = basepath+'M'+radar_res+radar_name+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+'P'+radar_res+radar_name+'/'

    filename = glob.glob(datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
    if not filename:
        warn('No file found in '+datapath+basename+timeinfo+'*.'+scan_list[0])
    else:
        radar = get_data_rad4alp(
            filename[0], datatype_list, scan_list[0], cfg, ind_rad=ind_rad)

    if len(scan_list) == 1:
        return radar

    # merge the elevations into a single radar instance
    for scan in scan_list[1:]:
        filename = glob.glob(datapath+basename+timeinfo+'*.'+scan+'*')
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.'+scan)
        else:
            radar_aux = get_data_rad4alp(
                filename[0], datatype_list, scan, cfg, ind_rad=ind_rad)
            if radar_aux is None:
                continue

            if radar is None:
                radar = radar_aux
            else:
                radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar

def merge_scans_odim(basepath, scan_list, radar_name, radar_res, voltime,
                     datatype_list, dataset_list, cfg, ind_rad=0):
    """
    merge odim data.

    Parameters
    ----------
    basepath : str
        base path of odim radar data
    scan_list : list
        list of scans (h5)
    voltime: datetime object
        reference time of the scan
    datatype_list : list
        lists of data types to get
    dataset_list : list
        list of datasets. Used to get path
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
    if radar_name is not None and radar_res is not None:
        basename = 'M'+radar_res+radar_name+dayinfo
    if cfg['path_convention'] == 'LTE':
        yy = dayinfo[0:2]
        dy = dayinfo[2:]
        subf = 'M'+radar_res+radar_name+yy+'hdf'+dy
        datapath = basepath+subf+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            subf = 'P'+radar_res+radar_name+yy+'hdf'+dy
            datapath = basepath+subf+'/'
    elif cfg['path_convention'] == 'MCH':
        datapath = basepath+dayinfo+'/'+basename+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+dayinfo+'/'+basename+'/'
    elif cfg['path_convention'] == 'ODIM':
        fpath_strf = dataset_list[0][dataset_list[0].find("D")+2:dataset_list[0].find("F")-2]
        fdate_strf = dataset_list[0][dataset_list[0].find("F")+2:-1]
        datapath = (basepath+voltime.strftime(fpath_strf)+'/')
        filenames = glob.glob(datapath+'*'+scan_list[0]+'*')
        filename = []
        for filename_aux in filenames:
            fdatetime = find_date_in_file_name(
                filename_aux, date_format=fdate_strf)
            if fdatetime == voltime:
                filename = [filename_aux]
    else:
        datapath = basepath+'M'+radar_res+radar_name+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+'P'+radar_res+radar_name+'/'
            filename = glob.glob(datapath+basename+timeinfo+'*'+scan_list[0] + '*')
    if not filename:
        warn('No file found in '+datapath[0]+basename+timeinfo+'*.h5')
    else:
        radar = get_data_odim(
            filename[0], datatype_list, scan_list[0], cfg, ind_rad=ind_rad)

    if len(scan_list) == 1:
        return radar

    # merge the elevations into a single radar instance
    for scan in scan_list[1:]:
        if cfg['path_convention'] == 'ODIM':
            filenames = glob.glob(datapath+'*'+scan+'*')
            filename = []
            for filename_aux in filenames:
                fdatetime = find_date_in_file_name(
                    filename_aux, date_format=fdate_strf)
                if fdatetime == voltime:
                    filename = [filename_aux]
                    break
        else:
            filename = glob.glob(datapath+basename+timeinfo+'*'+scan+'*')
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.'+scan)
        else:
            radar_aux = get_data_odim(
                filename[0], datatype_list, scan, cfg, ind_rad=ind_rad)
            if radar_aux is None:
                continue

            if radar is None:
                radar = radar_aux
            else:
                radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_mxpol(basepath, scan_list, voltime, datatype_list, cfg):
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
        datapath = basepath+'/'+sub1+'/'+sub2+'/'+sub3+'/'
        scanname = 'MXPol-polar-'+dayinfo+'-'+timeinfo+'*-'
        filename = glob.glob(datapath+scanname+scan_list[0]+'*')
    else:
        daydir = voltime.strftime('%Y-%m-%d')
        dayinfo = voltime.strftime('%Y%m%d')
        timeinfo = voltime.strftime('%H%M')
        datapath = basepath+scan_list[0]+'/'+daydir+'/'
        if not os.path.isdir(datapath):
            warn("WARNING: Unknown datapath '%s'" % datapath)
            return None
        filename = glob.glob(
            datapath+'MXPol-polar-'+dayinfo+'-'+timeinfo+'*-' +
            scan_list[0]+'.nc')
    if not filename:
        warn('No file found matching '+datapath+scanname+scan_list[0]+'*')
    else:
        radar = get_data_mxpol(filename[0], datatype_list)

    if len(scan_list) == 1:
        return radar

    # merge the elevations into a single radar instance
    for scan in scan_list[1:]:
        if cfg['path_convention'] == 'LTE':
            sub1 = str(voltime.year)
            sub2 = voltime.strftime('%m')
            sub3 = voltime.strftime('%d')
            dayinfo = voltime.strftime('%Y%m%d')
            timeinfo = voltime.strftime('%H%M')
            datapath = basepath+'/'+sub1+'/'+sub2+'/'+sub3+'/'
            scanname = 'MXPol-polar-'+dayinfo+'-'+timeinfo+'*-'
            filename = glob.glob(datapath+scanname+scan+'*')
        else:
            daydir = voltime.strftime('%Y-%m-%d')
            dayinfo = voltime.strftime('%Y%m%d')
            timeinfo = voltime.strftime('%H%M')
            datapath = basepath+scan+'/'+daydir+'/'
            if not os.path.isdir(datapath):
                warn("WARNING: Unknown datapath '%s'" % datapath)
                return None
            filename = glob.glob(
                datapath+'MXPol-polar-'+dayinfo+'-'+timeinfo+'*-'+scan+'.nc')
        if not filename:
            warn('No file found in '+datapath+scanname+scan)
        else:
            radar_aux = get_data_mxpol(filename[0], datatype_list)

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
    # look for COSMO data
    filename_list = list()
    for datatype in datatype_list:
        filename = find_cosmo_file(
            voltime, datatype, cfg, cfg['ScanList'][ind_rad][0],
            ind_rad=ind_rad)
        if filename is not None:
            filename_list.append(filename)

    nfiles_valid = len(filename_list)

    if nfiles_valid > 0:
        radar = merge_fields_cosmo(filename_list)

    # merge scans into a single radar instance
    nscans = len(cfg['ScanList'][ind_rad])
    if nscans > 1:
        for scan in cfg['ScanList'][ind_rad][1:]:
            filename_list = list()
            for datatype in datatype_list:
                filename = find_cosmo_file(
                    voltime, datatype, cfg, scan, ind_rad=ind_rad)
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
    for scan in cfg['ScanList'][ind_rad]:
        filename = find_rad4alpcosmo_file(voltime, datatype, cfg, scan)
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
        if radar is not None:
            radar.fields = dict()

            ngates = 0
            if cfg['rmax'] > 0:
                ngates = radar.ngates
            cosmo_field = read_rad4alp_cosmo(
                filename_list[0], datatype, ngates=ngates)
            if cosmo_field is None:
                return None

            radar.add_field(get_fieldname_pyart(datatype), cosmo_field)

    if len(cfg['ScanList'][ind_rad]) == 1:
        return radar

    # add the other scans
    for i, scan in enumerate(cfg['ScanList'][ind_rad][1:], start=1):
        filename = glob.glob(datapath+basename+timeinfo+'*.'+scan)
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.'+scan)
        else:
            radar_aux = get_data_rad4alp(
                filename[0], ['dBZ'], scan, cfg, ind_rad=ind_rad)
            if radar_aux is None:
                return None
            radar_aux.fields = dict()

            ngates = 0
            if cfg['rmax'] > 0:
                ngates = radar_aux.ngates
            cosmo_field = read_rad4alp_cosmo(
                filename_list[i], datatype, ngates=ngates)
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
    merge DEM rad4alp scans. If data for all the scans cannot be retrieved
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
    elif cfg['path_convention'] == 'MCH':
        datapath = basepath+dayinfo+'/'+basename+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+dayinfo+'/'+basename+'/'
    else:
        datapath = basepath+'M'+radar_res+radar_name+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+'P'+radar_res+radar_name+'/'

    filename = glob.glob(datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
    if not filename:
        warn('No file found in '+datapath+basename+timeinfo+'*.'+scan_list[0])
    else:
        radar = get_data_rad4alp(
            filename[0], ['dBZ'], scan_list[0], cfg, ind_rad=ind_rad)

        if radar is not None:
            radar.fields = dict()

            # add visibility data for first scan
            if cfg['rmax'] > 0:
                ngates = radar.ngates
                radar.add_field(
                    get_fieldname_pyart(datatype),
                    vis_list[int(scan_list[0])-1][:, :ngates])
            else:
                radar.add_field(
                    get_fieldname_pyart(datatype),
                    vis_list[int(scan_list[0])-1])

    if len(scan_list) == 1:
        return radar

    # add the other scans
    for scan in scan_list[1:]:
        filename = glob.glob(datapath+basename+timeinfo+'*.'+scan)
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.'+scan)
            continue
        radar_aux = get_data_rad4alp(
            filename[0], ['dBZ'], scan, cfg, ind_rad=ind_rad)
        if radar_aux is None:
            continue

        radar_aux.fields = dict()
        if cfg['rmax'] > 0:
            ngates = radar_aux.ngates
            radar_aux.add_field(
                get_fieldname_pyart(datatype),
                vis_list[int(scan_list[0])-1][:, :ngates])
        else:
            radar_aux.add_field(
                get_fieldname_pyart(datatype), vis_list[int(scan)-1])
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
    if cfg['path_convention'] == 'LTE':
        yy = dayinfo[0:2]
        dy = dayinfo[2:]
        subf = 'M'+radar_res+radar_name+yy+'hdf'+dy
        datapath_hydro = basepath+subf+'/'
    elif cfg['path_convention'] == 'MCH':
        datapath_hydro = basepath+dayinfo+'/'+basename_hydro+'/'
    else:
        datapath_hydro = basepath+'YM'+radar_name+'/'

    filename_hydro = glob.glob(datapath_hydro+basename_hydro+timeinfo+'*.' +
                               str(800+int(scan_list[0]))+'*')
    if not filename_hydro:
        warn('No file found in '+datapath_hydro+basename_hydro+timeinfo+'*.' +
             str(800+int(scan_list[0])))
        return None

    filename_hydro = filename_hydro[0]

    hydro_obj = pyart.aux_io.read_product(
        filename_hydro, physic_value=False, masked_array=True)

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
    elif cfg['path_convention'] == 'MCH':
        datapath = basepath+dayinfo+'/'+basename+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+dayinfo+'/'+basename+'/'
    else:
        datapath = basepath+'M'+radar_res+radar_name+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+'P'+radar_res+radar_name+'/'

    filename = glob.glob(datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
    if not filename:
        warn('No file found in '+datapath+basename+timeinfo+'*.'+scan_list[0])
    else:
        radar = get_data_rad4alp(
            filename[0], ['dBZ'], scan_list[0], cfg, ind_rad=ind_rad)
        if radar is not None:
            radar.fields = dict()

            # add hydrometeor classification data for first scan
            if cfg['rmax'] > 0.:
                ngates = radar.ngates
                hydro_dict['data'] = hydro_dict['data'][:, :ngates]
            radar.add_field(hydro_field, hydro_dict)

    if len(scan_list) == 1:
        return radar

    # add the other scans
    for scan in scan_list[1:]:
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan+'*')
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.' +
                 scan)
            continue

        radar_aux = get_data_rad4alp(
            filename[0], ['dBZ'], scan, cfg, ind_rad=ind_rad)
        if radar_aux is None:
            continue
        radar_aux.fields = dict()

        # read hydrometeor classification data file for other scans
        filename_hydro = glob.glob(datapath_hydro+basename_hydro+timeinfo +
                                   '*.'+str(800+int(scan))+'*')
        if not filename_hydro:
            warn('No file found in '+datapath_hydro+basename_hydro+timeinfo +
                 '*.'+str(800+int(scan)))
            continue

        filename_hydro = filename_hydro[0]

        hydro_obj = pyart.aux_io.read_product(
            filename_hydro, physic_value=False, masked_array=True)

        if hydro_obj is None:
            warn('Unable to read file '+filename_hydro)
            continue

        hydro_field = get_fieldname_pyart(datatype)
        hydro_dict = pyart.config.get_metadata(hydro_field)
        hydro_dict['data'] = map_hydro(hydro_obj.data)

        if cfg['rmax'] > 0.:
            ngates = radar_aux.ngates
            hydro_dict['data'] = hydro_dict['data'][:, :ngates]
        radar_aux.add_field(hydro_field, hydro_dict)
        if radar is None:
            radar = radar_aux
        else:
            radar = pyart.util.radar_utils.join_radar(radar, radar_aux)

    return radar


def merge_scans_Doppler_rad4alp(voltime, datatype, cfg, ind_rad=0):
    """
    merge rad4alp dealised Doppler velocity scans. If data for all the scans
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

    # read Dopplermeteor classification data file for first scan
    basename_Doppler = 'DV'+radar_name+dayinfo
    if cfg['path_convention'] == 'LTE':
        yy = dayinfo[0:2]
        dy = dayinfo[2:]
        subf = 'M'+radar_res+radar_name+yy+'hdf'+dy
        datapath_Doppler = basepath+subf+'/'
    elif cfg['path_convention'] == 'MCH':
        datapath_Doppler = basepath+dayinfo+'/'+basename_Doppler+'/'
    else:
        datapath_Doppler = basepath+'YM'+radar_name+'/'

    filename_Doppler = glob.glob(datapath_Doppler+basename_Doppler+timeinfo+'*.' +
                                 str(800+int(scan_list[0]))+'*')
    if not filename_Doppler:
        warn('No file found in '+datapath_Doppler+basename_Doppler+timeinfo+'*.' +
             str(800+int(scan_list[0])))
        return None

    filename_Doppler = filename_Doppler[0]

    Doppler_obj = pyart.aux_io.read_product(
        filename_Doppler, physic_value=False, masked_array=True)

    if Doppler_obj is None:
        warn('Unable to read file '+filename_Doppler)
        return None

    Doppler_field = get_fieldname_pyart(datatype)
    Doppler_dict = pyart.config.get_metadata(Doppler_field)
    Doppler_dict['data'] = map_Doppler(
        Doppler_obj.data, float(Doppler_obj.header['nyquist']))

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
    elif cfg['path_convention'] == 'MCH':
        datapath = basepath+dayinfo+'/'+basename+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+dayinfo+'/'+basename+'/'
    else:
        datapath = basepath+'M'+radar_res+radar_name+'/'
        filename = glob.glob(
            datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
        if not filename:
            basename = 'P'+radar_res+radar_name+dayinfo
            datapath = basepath+'P'+radar_res+radar_name+'/'

    filename = glob.glob(datapath+basename+timeinfo+'*.'+scan_list[0] + '*')
    if not filename:
        warn('No file found in '+datapath+basename+timeinfo+'*.'+scan_list[0])
    else:
        radar = get_data_rad4alp(
            filename[0], ['dBZ'], scan_list[0], cfg, ind_rad=ind_rad)
        if radar is not None:
            radar.fields = dict()

            # add Doppler data for first scan
            if cfg['rmax'] > 0.:
                ngates = radar.ngates
                Doppler_dict['data'] = Doppler_dict['data'][:, :ngates]
            radar.add_field(Doppler_field, Doppler_dict)

    if len(scan_list) == 1:
        return radar

    # add the other scans
    for scan in scan_list[1:]:
        filename = glob.glob(datapath+basename+timeinfo+'*.'+scan+'*')
        if not filename:
            warn('No file found in '+datapath+basename+timeinfo+'*.'+scan)
            continue

        radar_aux = get_data_rad4alp(
            filename[0], ['dBZ'], scan, cfg, ind_rad=ind_rad)
        if radar_aux is None:
            continue
        radar_aux.fields = dict()

        # read Dopplermeteor classification data file for other scans
        filename_Doppler = glob.glob(datapath_Doppler+basename_Doppler+timeinfo +
                                     '*.'+str(800+int(scan))+'*')
        if not filename_Doppler:
            warn('No file found in '+datapath_Doppler+basename_Doppler+timeinfo +
                 '*.'+str(800+int(scan)))
            continue

        filename_Doppler = filename_Doppler[0]

        Doppler_obj = pyart.aux_io.read_product(
            filename_Doppler, physic_value=False, masked_array=True)

        if Doppler_obj is None:
            warn('Unable to read file '+filename_Doppler)
            continue

        Doppler_field = get_fieldname_pyart(datatype)
        Doppler_dict = pyart.config.get_metadata(Doppler_field)
        Doppler_dict['data'] = map_Doppler(
            Doppler_obj.data, float(Doppler_obj.header['nyquist']))

        if cfg['rmax'] > 0.:
            ngates = radar_aux.ngates
            Doppler_dict['data'] = Doppler_dict['data'][:, :ngates]
        radar_aux.add_field(Doppler_field, Doppler_dict)
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

    if len(datatype_list) == 1:
        return radar

    # add other fields in the same scan
    for datatype in datatype_list[1:]:
        if datatype not in ('Nh', 'Nv'):
            filename = glob.glob(datapath+fdatetime+datatype+'.*')
        elif datatype == 'Nh':
            filename = glob.glob(datapath+fdatetime+'dBZ.*')
        else:
            filename = glob.glob(datapath+fdatetime+'dBZv.*')
        if not filename:
            warn('No file found in '+datapath+fdatetime+datatype+'.*')
        else:
            radar_aux = get_data_rainbow(filename[0], datatype)
            if radar_aux is None:
                continue

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


def merge_fields_pyrad(basepath, loadname, voltime, datatype_list,
                       dataset_list, product_list, rmax=0.,
                       termination='.nc'):
    """
    merge fields from Pyrad-generated files into a single radar object.
    Accepted file types are CFRadial and ODIM.

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
    rmax : float
        maximum range that will be kept.
    termination : str
        file termination type. Can be '.nc' or '.h5'

    Returns
    -------
    radar : Radar
        radar object

    """
    datapath = (basepath+loadname+'/'+voltime.strftime('%Y-%m-%d')+'/' +
                dataset_list[0]+'/'+product_list[0]+'/')
    fdatetime = voltime.strftime('%Y%m%d%H%M%S')
    filename = glob.glob(datapath+fdatetime+'*'+datatype_list[0]+termination)

    # create radar object
    radar = None
    if not filename:
        warn('No file found in '+datapath+fdatetime+'*'+datatype_list[0] +
             '.nc')
    else:
        if termination == '.nc':
            radar = pyart.io.read_cfradial(filename[0])
        else:
            radar = pyart.aux_io.read_odim_h5(filename[0])
        if rmax > 0.:
            radar.range['data'] = radar.range['data'][
                radar.range['data'] < rmax]
            radar.ngates = len(radar.range['data'])
            for field in radar.fields:
                radar.fields[field]['data'] = (
                    radar.fields[field]['data'][:, :radar.ngates])
            radar.gate_x['data'] = radar.gate_x['data'][:, :radar.ngates]
            radar.init_gate_x_y_z()
            radar.init_gate_longitude_latitude()
            radar.init_gate_altitude()

    if len(datatype_list) > 1:
        # add other fields in the same scan
        for i, dataset in enumerate(dataset_list[1:], start=1):
            datapath = (
                basepath+loadname+'/'+voltime.strftime('%Y-%m-%d')+'/' +
                dataset+'/'+product_list[i]+'/')
            filename = glob.glob(
                datapath+fdatetime+'*'+datatype_list[i]+termination)
            if not filename:
                warn('No file found in '+datapath+fdatetime+'*' +
                     datatype_list[i]+'.nc')
                continue

            if termination == '.nc':
                radar_aux = pyart.io.read_cfradial(filename[0])
            else:
                radar_aux = pyart.aux_io.read_odim_h5(filename[0])
            if rmax > 0.:
                radar_aux.range['data'] = radar_aux.range['data'][
                    radar_aux.range['data'] < rmax]
                radar_aux.ngates = len(radar_aux.range['data'])
                for field in radar_aux.fields:
                    radar_aux.fields[field]['data'] = (
                        radar_aux.fields[field]['data'][
                            :, :radar_aux.ngates])
                radar_aux.gate_x['data'] = (
                    radar_aux.gate_x['data'][:, :radar_aux.ngates])
                radar_aux.init_gate_x_y_z()
                radar_aux.init_gate_longitude_latitude()
                radar_aux.init_gate_altitude()
            if radar is None:
                radar = radar_aux
            else:
                add_field(radar, radar_aux)

    for field in radar.fields:
        radar.fields[field]['data'] = np.ma.asarray(
            radar.fields[field]['data'])

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

    if len(datatype_list) == 1:
        return radar

    # add other fields in the same scan
    for datatype in datatype_list:
        datapath = basepath+datatype+'/'+scan_name+'/'
        filename = glob.glob(datapath+datatype+'_'+scan_name_aux)
        if not filename:
            warn('No file found in '+datapath+datatype+'_' +
                 scan_name_aux)
        else:
            radar_aux = get_data_rainbow(filename[0], datatype)
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
    try:
        radar = pyart.aux_io.read_rainbow_wrl(filename_list[0])
    except OSError as ee:
        warn(str(ee))
        warn('Unable to read file '+filename_list[0])
        return None

    if radar is None:
        return None

    if len(filename_list) == 1:
        return radar

    # add other COSMO fields in the same scan
    for filename in filename_list:
        try:
            radar_aux = pyart.aux_io.read_rainbow_wrl(filename)
        except OSError as ee:
            warn(str(ee))
            warn('Unable to read file '+filename)
            continue
        if radar_aux is None:
            continue
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
    radar : Radar or None
        radar object if the reading of the data has been successful.
        None otherwise

    """
    try:
        radar = pyart.aux_io.read_rainbow_wrl(filename)
    except OSError as ee:
        warn(str(ee))
        warn('Unable to read file '+filename)
        return None
    if radar is None:
        return None

    if datatype in ('Nh', 'Nv'):
        try:
            with open(filename, 'rb') as fid:
                rbf = wrl.io.read_rainbow(fid, loaddata=True)
                fid.close()
        except OSError as ee:
            warn(str(ee))
            warn('Unable to read file '+filename)
            return None

        # check the number of slices
        nslices = int(rbf['volume']['scan']['pargroup']['numele'])
        if nslices > 1:
            common_slice_info = rbf['volume']['scan']['slice'][0]
        else:
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
        radar object. None if the reading has not been successful

    """
    metranet_field_names = dict()
    for datatype in datatype_list:
        if datatype not in ('Nh', 'Nv'):
            metranet_field_names.update(get_datatype_metranet(datatype))

    if cfg['path_convention'] == 'LTE':
        radar = pyrad_MCH(filename, field_names=metranet_field_names)
    else:
        try:
            radar = pyart.aux_io.read_metranet(
                filename, field_names=metranet_field_names, rmax=cfg['rmax'])
        except ValueError as ee:
            warn("Unable to read file '"+filename+": (%s)" % str(ee))
            return None

    if ('Nh' not in datatype_list) and ('Nv' not in datatype_list):
        return radar

    # create noise moments
    # read radar information in status file
    voltime = get_datetime(filename, 'RAD4ALP:dBZ')
    root = read_status(voltime, cfg, ind_rad=ind_rad)
    if root is None:
        return radar

    sweep_number = int(scan_name)-1
    if 'Nh' in datatype_list:
        found = False
        for sweep in root.findall('sweep'):
            sweep_number_file = (
                int(sweep.attrib['name'].split('.')[1])-1)
            if sweep_number_file == sweep_number:
                noise_h = sweep.find(
                    "./RADAR/STAT/CALIB/noisepower_frontend_h_inuse")
                rconst_h = sweep.find("./RADAR/STAT/CALIB/rconst_h")
                if noise_h is None or rconst_h is None:
                    warn('Horizontal channel noise power not ' +
                         'available for sweep '+scan_name)
                    break

                noisedBADU_h = 10.*np.log10(
                    float(noise_h.attrib['value']))
                rconst_h = float(rconst_h.attrib['value'])

                noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                    radar.nrays, noisedBADU_h+rconst_h,
                    radar.range['data'], 100.,
                    noise_field='noisedBZ_hh')

                radar.add_field('noisedBZ_hh', noisedBZ_h)

                found = True
        if not found:
            warn('Horizontal channel noise power not ' +
                 'available for sweep '+scan_name)

    if 'Nv' in datatype_list:
        found = False
        for sweep in root.findall('sweep'):
            sweep_number_file = (
                int(sweep.attrib['name'].split('.')[1])-1)
            if sweep_number_file == sweep_number:
                noise_v = sweep.find(
                    "./RADAR/STAT/CALIB/noisepower_frontend_v_inuse")
                rconst_v = sweep.find("./RADAR/STAT/CALIB/rconst_v")
                if noise_v is None or rconst_v is None:
                    warn('Vertical channel noise power not ' +
                         'available for sweep '+scan_name)
                    break

                noisedBADU_v = 10.*np.log10(
                    float(noise_v.attrib['value']))
                rconst_v = float(rconst_v.attrib['value'])

                noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                    radar.nrays, noisedBADU_v+rconst_v,
                    radar.range['data'], 100.,
                    noise_field='noisedBZ_vv')

                radar.add_field('noisedBZ_vv', noisedBZ_v)

                found = True
        if not found:
            warn('Horizontal channel noise power not ' +
                 'available for sweep '+scan_name)

    return radar


def get_data_odim(filename, datatype_list, scan_name, cfg, ind_rad=0):
    """
    gets ODIM radar data

    Parameters
    ----------
    filename : str
        name of file containing odim data
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
        radar object. None if the reading has not been successful

    """
    odim_field_names = dict()
    for datatype in datatype_list:
        if datatype not in ('Nh', 'Nv'):
            odim_field_names.update(get_datatype_odim(datatype))
        try:
            radar = pyart.aux_io.read_odim_h5(
                filename, field_names=odim_field_names)
        except ValueError as ee:
            warn("Unable to read file '"+filename+": (%s)" % str(ee))
            return None

    if ('Nh' not in datatype_list) and ('Nv' not in datatype_list):
        return radar

    # create noise moments
    # read radar information in status file
    voltime = get_datetime(filename, 'ODIM:dBZ')
    root = read_status(voltime, cfg, ind_rad=ind_rad)
    if root is None:
        return radar

    sweep_number = int(scan_name)-1
    if 'Nh' in datatype_list:
        found = False
        for sweep in root.findall('sweep'):
            sweep_number_file = (
                int(sweep.attrib['name'].split('.')[1])-1)
            if sweep_number_file == sweep_number:
                noise_h = sweep.find(
                    "./RADAR/STAT/CALIB/noisepower_frontend_h_inuse")
                rconst_h = sweep.find("./RADAR/STAT/CALIB/rconst_h")
                if noise_h is None or rconst_h is None:
                    warn('Horizontal channel noise power not ' +
                         'available for sweep '+scan_name)
                    break

                noisedBADU_h = 10.*np.log10(
                    float(noise_h.attrib['value']))
                rconst_h = float(rconst_h.attrib['value'])

                noisedBZ_h = pyart.retrieve.compute_noisedBZ(
                    radar.nrays, noisedBADU_h+rconst_h,
                    radar.range['data'], 100.,
                    noise_field='noisedBZ_hh')

                radar.add_field('noisedBZ_hh', noisedBZ_h)

                found = True
        if not found:
            warn('Horizontal channel noise power not ' +
                 'available for sweep '+scan_name)

    if 'Nv' in datatype_list:
        found = False
        for sweep in root.findall('sweep'):
            sweep_number_file = (
                int(sweep.attrib['name'].split('.')[1])-1)
            if sweep_number_file == sweep_number:
                noise_v = sweep.find(
                    "./RADAR/STAT/CALIB/noisepower_frontend_v_inuse")
                rconst_v = sweep.find("./RADAR/STAT/CALIB/rconst_v")
                if noise_v is None or rconst_v is None:
                    warn('Vertical channel noise power not ' +
                         'available for sweep '+scan_name)
                    break

                noisedBADU_v = 10.*np.log10(
                    float(noise_v.attrib['value']))
                rconst_v = float(rconst_v.attrib['value'])

                noisedBZ_v = pyart.retrieve.compute_noisedBZ(
                    radar.nrays, noisedBADU_v+rconst_v,
                    radar.range['data'], 100.,
                    noise_field='noisedBZ_vv')

                radar.add_field('noisedBZ_vv', noisedBZ_v)

                found = True
        if not found:
            warn('Horizontal channel noise power not ' +
                 'available for sweep '+scan_name)

    return radar


def get_data_mxpol(filename, datatype_list):
    """
    gets MXPol radar data

    Parameters
    ----------
    filename : str
        name of file containing MXPol data
    datatype_list : list of strings
        list of data fields to get

    Returns
    -------
    radar : Radar
        radar object

    """
    field_names = dict()
    for datatype in datatype_list:
        if datatype not in ('Nh', 'Nv'):
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


def interpol_field(radar_dest, radar_orig, field_name, fill_value=None,
                   ang_tol=0.5):
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
    fill_value: float
        The fill value
    ang_tol : float
        angle tolerance to determine whether the radar origin sweep is the
        radar destination sweep

    Returns
    -------
    field_dest : dict
        interpolated field and metadata

    """
    if radar_dest.nsweeps != radar_orig.nsweeps:
        warn('Number of sweeps in destination radar object different from ' +
             'origin radar object. Orig: '+str(radar_orig.nsweeps) +
             ' Dest : '+str(radar_dest.nsweeps))

    if fill_value is None:
        fill_value = pyart.config.get_fillvalue()

    field_orig_data = radar_orig.fields[field_name]['data'].filled(
        fill_value=fill_value)
    field_dest = deepcopy(radar_orig.fields[field_name])
    field_dest['data'] = np.ma.masked_all(
        (radar_dest.nrays, radar_dest.ngates))

    for sweep in range(radar_dest.nsweeps):
        sweep_start_dest = radar_dest.sweep_start_ray_index['data'][sweep]
        sweep_end_dest = radar_dest.sweep_end_ray_index['data'][sweep]
        fixed_angle = radar_dest.fixed_angle['data'][sweep]
        nrays_sweep = radar_dest.rays_per_sweep['data'][sweep]

        # look for nearest angle
        delta_ang = np.absolute(radar_orig.fixed_angle['data']-fixed_angle)
        ind_sweep_orig = np.argmin(delta_ang)

        if delta_ang[ind_sweep_orig] > ang_tol:
            warn('No fixed angle of origin radar object matches the fixed ' +
                 'angle of destination radar object for sweep nr ' +
                 str(sweep)+' with fixed angle '+str(fixed_angle)+'+/-' +
                 str(ang_tol))
            field_dest_sweep = np.ma.masked_all(
                (nrays_sweep, radar_dest.ngates), dtype=float)
        else:
            sweep_start_orig = radar_orig.sweep_start_ray_index['data'][
                ind_sweep_orig]
            sweep_end_orig = radar_orig.sweep_end_ray_index['data'][
                ind_sweep_orig]

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
                ind_ang = np.argsort(radar_orig.elevation['data'][
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
