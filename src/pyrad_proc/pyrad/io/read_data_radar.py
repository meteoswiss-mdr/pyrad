"""
pyrad.io.read_data
====================

Functions for reading pyrad input data, i.e. radar files

.. autosummary::
    :toctree: generated/

    get_data
    merge_scans_rainbow
    merge_scans_rad4alp
    merge_scans_cosmo
    merge_scans_cosmo_rad4alp
    merge_fields_rainbow
    merge_fields_cosmo
    get_data_rainbow
    get_data_rad4alp
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
from warnings import warn

import numpy as np

import pyart
import wradlib as wrl

from .read_data_aux import read_status, read_rad4alp_cosmo


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
                for field_name in radar_aux.fields.keys():
                    radar.add_field(field_name, radar_aux.fields[field_name])

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
    datatype_list : list
        lists of data types to get
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
    radar = merge_scans_rad4alp(
        cfg['datapath'], cfg['ScanList'][0], cfg['RadarName'],
        cfg['RadarRes'], voltime, get_datatypemetranet('dBZ'), cfg)
    radar.fields = dict()

    cosmo_field = read_rad4alp_cosmo(filename_list[0], datatype)
    if cosmo_field is None:
        return None

    radar.add_field(get_fieldname_rainbow(datatype), cosmo_field)

    # add the other scans
    for i in range(1, len(cfg['ScanList'])):
        radar_aux = merge_scans_rad4alp(
            cfg['datapath'], cfg['ScanList'][i], cfg['RadarName'],
            cfg['RadarRes'], voltime, get_datatypemetranet('dBZ'), cfg)
        radar_aux.fields = dict()

        cosmo_field = read_rad4alp_cosmo(filename_list[i], datatype)
        if cosmo_field is None:
            return None

        radar_aux.add_field(get_fieldname_rainbow(datatype), cosmo_field)

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
            metranet_field_names.update(get_datatypemetranet(datatype))

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
    elif datatype == 'dBuZ':
        field_name = 'unfiltered_reflectivity'
    elif datatype == 'dBZv':
        field_name = 'reflectivity_vv'
    elif datatype == 'dBuZv':
        field_name = 'unfiltered_reflectivity_vv'
    elif datatype == 'dBm':
        field_name = 'signal_power_hh'
    elif datatype == 'dBmv':
        field_name = 'signal_power_vv'
    elif datatype == 'dBm_sun_hit':
        field_name = 'sun_hit_power_h'
    elif datatype == 'dBmv_sun_hit':
        field_name = 'sun_hit_power_v'
    elif datatype == 'ZDR_sun_hit':
        field_name = 'sun_hit_differential_reflectivity'
    elif datatype == 'dBm_sun_est':
        field_name = 'sun_est_power_h'
    elif datatype == 'dBmv_sun_est':
        field_name = 'sun_est_power_v'
    elif datatype == 'ZDR_sun_est':
        field_name = 'sun_est_differential_reflectivity'
    elif datatype == 'sun_pos_h':
        field_name = 'sun_hit_h'
    elif datatype == 'sun_pos_v':
        field_name = 'sun_hit_v'
    elif datatype == 'sun_pos_zdr':
        field_name = 'sun_hit_zdr'
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
    elif datatype == 'ZDRu':
        field_name = 'unfiltered_differential_reflectivity'
    elif datatype == 'ZDRc':
        field_name = 'corrected_differential_reflectivity'
    elif datatype == 'RhoHV':
        field_name = 'cross_correlation_ratio'
    elif datatype == 'uRhoHV':
        field_name = 'uncorrected_cross_correlation_ratio'
    elif datatype == 'RhoHVc':
        field_name = 'corrected_cross_correlation_ratio'
    elif datatype == 'RhoHV_rain':
        field_name = 'cross_correlation_ratio_in_rain'
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
    elif datatype == 'PhiDP0':
        field_name = 'system_differential_phase'
    elif datatype == 'PhiDP0_bin':
        field_name = 'first_gate_differential_phase'
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
    elif datatype == 'dBZ_bias':
        field_name = 'reflectivity_bias'
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
