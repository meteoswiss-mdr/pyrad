"""
pyrad.io.io_aux
===============

Auxiliary functions for reading/writing files

.. autosummary::
    :toctree: generated/

    get_save_dir
    make_filename
    generate_field_name_str
    get_datatype_metranet
    get_fieldname_pyart
    get_file_list
    get_datatype_fields
    get_dataset_fields
    get_datetime
    find_cosmo_file
    find_rad4alpcosmo_file

"""

import os
import glob
import datetime
import csv
import xml.etree.ElementTree as et
from warnings import warn

import numpy as np

from pyart.config import get_metadata


def get_save_dir(basepath, procname, dsname, prdname, timeinfo=None,
                 timeformat='%Y-%m-%d', create_dir=True):
    """
    obtains the path to a product directory and eventually creates it

    Parameters
    ----------
    basepath : str
        product base path
    procname : str
        name of processing space
    dsname : str
        data set name
    prdname : str
        product name
    timeinfo : datetime
        time info to generate the date directory. If None there is no time
        format in the path
    timeformat : str
        Optional. The time format.
    create_dir : boolean
        If True creates the directory

    Returns
    -------
    savedir : str
        path to product

    """
    if timeinfo is None:
        savedir = basepath+procname+'/'+dsname+'/'+prdname+'/'
    else:
        daydir = timeinfo.strftime('%Y-%m-%d')
        savedir = basepath+procname+'/'+daydir+'/'+dsname+'/'+prdname+'/'

    if create_dir is False:
        return savedir

    if os.path.isdir(savedir):
        pass
    else:
        os.makedirs(savedir)

    return savedir


def make_filename(prdtype, dstype, dsname, ext, prdcfginfo=None,
                  timeinfo=None, timeformat='%Y%m%d%H%M%S'):
    """
    creates a product file name

    Parameters
    ----------
    timeinfo : datetime
        time info to generate the date directory

    prdtype : str
        product type, i.e. 'ppi', etc.

    dstype : str
        data set type, i.e. 'raw', etc.

    dsname : str
        data set name

    ext : str
        file name extension, i.e. 'png'

    prdcfginfo : str
        Optional. string to add product configuration information, i.e. 'el0.4'

    timeformat : str
        Optional. The time format

    Returns
    -------
    fname : str
        file name

    """
    if timeinfo is None:
        timeinfostr = ''
    else:
        timeinfostr = timeinfo.strftime(timeformat)+'_'

    if prdcfginfo is None:
        cfgstr = ''
    else:
        cfgstr = '_'+prdcfginfo

    fname = timeinfostr+prdtype+'_'+dstype+'_'+dsname+cfgstr+'.'+ext

    return fname


def generate_field_name_str(datatype):
    """
    Generates a field name in a nice to read format.

    Parameters
    ----------
    datatype : str
        The data type

    Returns
    -------
    field_str : str
        The field name

    """
    field_name = get_fieldname_pyart(datatype)
    field_dic = get_metadata(field_name)
    field_str = field_dic['standard_name'].replace('_', ' ')
    field_str = field_str[0].upper() + field_str[1:]
    field_str += ' ('+field_dic['units']+')'

    return field_str


def get_datatype_metranet(datatype):
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


def get_fieldname_pyart(datatype):
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
    elif datatype == 'dBZc':
        field_name = 'corrected_reflectivity'
    elif datatype == 'dBuZc':
        field_name = 'corrected_unfiltered_reflectivity'
    elif datatype == 'dBZv':
        field_name = 'reflectivity_vv'
    elif datatype == 'dBuZv':
        field_name = 'unfiltered_reflectivity_vv'
    elif datatype == 'dBuZvc':
        field_name = 'corrected_unfiltered_reflectivity_vv'
    elif datatype == 'dBZ_bias':
        field_name = 'reflectivity_bias'

    elif datatype == 'ZDR':
        field_name = 'differential_reflectivity'
    elif datatype == 'ZDRu':
        field_name = 'unfiltered_differential_reflectivity'
    elif datatype == 'ZDRc':
        field_name = 'corrected_differential_reflectivity'
    elif datatype == 'ZDRuc':
        field_name = 'corrected_unfiltered_differential_reflectivity'
    elif datatype == 'ZDR_rain':
        field_name = 'differential_reflectivity_in_rain'

    elif datatype == 'dBm':
        field_name = 'signal_power_hh'
    elif datatype == 'dBmv':
        field_name = 'signal_power_vv'
    elif datatype == 'Nh':
        field_name = 'noisedBZ_hh'
    elif datatype == 'Nv':
        field_name = 'noisedBZ_vv'
    elif datatype == 'SNRh':
        field_name = 'signal_to_noise_ratio_hh'
    elif datatype == 'SNRv':
        field_name = 'signal_to_noise_ratio_vv'

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

    elif datatype == 'TEMP':
        field_name = 'temperature'
    elif datatype == 'ISO0':
        field_name = 'iso0'

    elif datatype == 'VIS':
        field_name = 'visibility'

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
    datagroup, datatype, dataset, product = get_datatype_fields(datadescriptor)

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


def get_datatype_fields(datadescriptor):
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


def get_dataset_fields(datasetdescr):
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
    datagroup, datatype, dataset, product = get_datatype_fields(datadescriptor)
    if datagroup == 'RAINBOW':
        datetimestr = bfile[0:14]
        fdatetime = datetime.datetime.strptime(datetimestr, '%Y%m%d%H%M%S')
    elif datagroup == 'RAD4ALP':
        datetimestr = bfile[3:12]
        fdatetime = datetime.datetime.strptime(datetimestr, '%y%j%H%M')
    elif datagroup == 'SAVED':
        print('caca')
    return fdatetime


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
