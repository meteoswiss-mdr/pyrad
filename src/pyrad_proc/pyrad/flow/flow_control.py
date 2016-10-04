"""
pyrad.flow.flow_control
=======================

functions to control the Pyrad data processing flow

.. autosummary::
    :toctree: generated/

    main
    _create_cfg_dict
    _create_datacfg_dict
    _create_dscfg_dict
    _create_prdcfg_dict
    _get_datatype_list
    _get_datasets_list
    _get_masterfile_list
    _add_dataset
    _process_dataset

"""
from warnings import warn
import os

from ..io.read_data import get_datetime, get_data, get_file_list
from ..io.read_data import get_datasetfields, get_datatypefields
from ..io.config import read_config

from ..proc.process_dataset import get_process_type
from ..proc.process_product import get_product_type

from pyrad import proc


def main(cfgfile, starttime, endtime):
    """
    main flow control. Processes data over a given period of time

    Parameters
    ----------
    cfgfile : str
        path of the main config file
    starttime, endtime : datetime object
        start and end time of the data to be processed

    Returns
    -------
    None

    """
    cfg = _create_cfg_dict(cfgfile)
    datacfg = _create_datacfg_dict(cfg)

    datatypesdescr = _get_datatype_list(cfg)
    dataset_levels = _get_datasets_list(cfg)
    masterfilelist, masterdatatypedescr = _get_masterfile_list(
        cfg['ScanList'][0], datatypesdescr, starttime, endtime, datacfg)

    nvolumes = len(masterfilelist)
    if nvolumes == 0:
        raise Exception(
            'ERROR: Could not find any volume within the specified times ' +
            ' for master scan '+masterscan+' and master data type ' +
            masterdatatypedescr)
    print('Number of volumes to process: '+str(nvolumes)+'\n\n')

    # initial processing of the datasets
    print('\n\nInitializing datasets:')
    for level in sorted(dataset_levels):
        print('\nProcess level: '+level)
        for dataset in dataset_levels[level]:
            print('Processing dataset: '+dataset)
            _process_dataset(cfg, dataset, proc_status=0, radar=None,
                             voltime=None)

    # process all data files in file list
    for masterfile in masterfilelist:
        print('\n\nmaster file: '+os.path.basename(masterfile))
        voltime = get_datetime(masterfile, masterdatatypedescr)

        # get all raw data
        radar = get_data(voltime, datatypesdescr, datacfg)

        # process all data sets
        for level in sorted(dataset_levels):
            print('\nProcess level: '+level)
            for dataset in dataset_levels[level]:
                print('Processing dataset: '+dataset)
                _process_dataset(cfg, dataset, proc_status=1, radar=radar,
                                 voltime=voltime)

    # post-processing of the datasets
    print('\n\nPost-processing datasets')
    for level in sorted(dataset_levels):
        print('\nProcess level: '+level)
        for dataset in dataset_levels[level]:
            print('Processing dataset: '+dataset)
            _process_dataset(cfg, dataset, proc_status=2, radar=None,
                             voltime=voltime)

    print('\n\n\nThis is the end my friend! See you soon!')


def _create_cfg_dict(cfgfile):
    """
    creates a configuration dictionary

    Parameters
    ----------
    cfgfile : str
        path of the main config file

    Returns
    -------
    cfg : dict
        dictionary containing the configuration data

    """
    cfg = dict({'configFile': cfgfile})
    cfg = read_config(cfg['configFile'], cfg=cfg)
    cfg = read_config(cfg['locationConfigFile'], cfg=cfg)
    cfg = read_config(cfg['productConfigFile'], cfg=cfg)

    # fill in defaults
    if 'cosmopath' not in cfg:
        cfg.update({'cosmopath': None})
    if 'dempath' not in cfg:
        cfg.update({'dempath': None})
    if 'smnpath' not in cfg:
        cfg.update({'smnpath': None})
    if 'disdropath' not in cfg:
        cfg.update({'disdropath': None})
    if 'loadbasepath' not in cfg:
        cfg.update({'loadbasepath': None})
    if 'loadname' not in cfg:
        cfg.update({'loadname': None})
    if 'RadarName' not in cfg:
        cfg.update({'RadarName': None})
    if 'RadarRes' not in cfg:
        cfg.update({'RadarRes': None})
    if 'mflossh' not in cfg:
        cfg.update({'mflossh': None})
    if 'mflossv' not in cfg:
        cfg.update({'mflossv': None})
    if 'radconsth' not in cfg:
        cfg.update({'radconsth': None})
    if 'radconstv' not in cfg:
        cfg.update({'radconstv': None})
    if 'attg' not in cfg:
        cfg.update({'attg': None})
    if 'ScanPeriod' not in cfg:
        warn(
            'WARNING: Scan period not specified.' +
            'Assumed default value 5 min')
        cfg.update({'ScanPeriod': 5})
    if 'CosmoRunFreq' not in cfg:
        warn(
            'WARNING: COSMO run frequency not specified.' +
            'Assumed default value 3h')
        cfg.update({'CosmoRunFreq': 3})
    if 'CosmoForecasted' not in cfg:
        warn(
            'WARNING: Hours forecasted by COSMO not specified.' +
            'Assumed default value 7h (including analysis)')
        cfg.update({'CosmoForecasted': 7})

    return cfg


def _create_datacfg_dict(cfg):
    """
    creates a data configuration dictionary from a config dictionary

    Parameters
    ----------
    cfg : dict
        config dictionary

    Returns
    -------
    datacfg : dict
        data config dictionary

    """
    datacfg = dict({'datapath': cfg['datapath']})
    datacfg.update({'ScanList': cfg['ScanList']})
    datacfg.update({'cosmopath': cfg['cosmopath']})
    datacfg.update({'dempath': cfg['dempath']})
    datacfg.update({'loadbasepath': cfg['loadbasepath']})
    datacfg.update({'loadname': cfg['loadname']})
    datacfg.update({'RadarName': cfg['RadarName']})
    datacfg.update({'RadarRes': cfg['RadarRes']})
    datacfg.update({'ScanPeriod': int(cfg['ScanPeriod'])})
    datacfg.update({'CosmoRunFreq': int(cfg['CosmoRunFreq'])})
    datacfg.update({'CosmoForecasted': int(cfg['CosmoForecasted'])})

    return datacfg


def _create_dscfg_dict(cfg, dataset, voltime=None):
    """
    creates a dataset configuration dictionary

    Parameters
    ----------
    cfg : dict
        config dictionary
    dataset : str
        name of the dataset
    voltime : datetime object
        time of the dataset

    Returns
    -------
    dscfg : dict
        dataset config dictionary

    """
    dscfg = cfg[dataset]
    dscfg.update({'configpath': cfg['configpath']})
    dscfg.update({'mflossh': cfg['mflossh']})
    dscfg.update({'mflossv': cfg['mflossv']})
    dscfg.update({'radconsth': cfg['radconsth']})
    dscfg.update({'radconstv': cfg['radconstv']})
    dscfg.update({'attg': cfg['attg']})

    dscfg.update({'basepath': cfg['saveimgbasepath']})
    dscfg.update({'procname': cfg['name']})
    dscfg.update({'dsname': dataset})
    if 'MAKE_GLOBAL' not in dscfg:
            dscfg.update({'MAKE_GLOBAL': 0})
    if voltime is not None:
        dscfg.update({'timeinfo': voltime})

    return dscfg


def _create_prdcfg_dict(cfg, dataset, product, voltime=None):
    """
    creates a product configuration dictionary

    Parameters
    ----------
    cfg : dict
        config dictionary
    dataset : str
        name of the dataset used to create the product
    product : str
        name of the product
    voltime : datetime object
        time of the dataset

    Returns
    -------
    prdcfg : dict
        product config dictionary

    """
    prdcfg = cfg[dataset]['products'][product]
    prdcfg.update({'procname': cfg['name']})
    prdcfg.update({'basepath': cfg['saveimgbasepath']})
    prdcfg.update({'smnpath': cfg['smnpath']})
    prdcfg.update({'disdropath': cfg['disdropath']})
    prdcfg.update({'ScanPeriod': cfg['ScanPeriod']})
    prdcfg.update({'imgformat': cfg['imgformat']})
    prdcfg.update({'convertformat': cfg['convertformat']})
    prdcfg.update(
        {'ppiImageConfig': cfg['ppiImageConfig']})
    prdcfg.update(
        {'rhiImageConfig': cfg['rhiImageConfig']})
    prdcfg.update(
        {'sunhitsImageConfig': cfg['sunhitsImageConfig']})
    prdcfg.update({'dsname': dataset})
    prdcfg.update({'dstype': cfg[dataset]['type']})
    prdcfg.update({'prdname': product})

    if voltime is not None:
        prdcfg.update({'timeinfo': voltime})

    return prdcfg


def _get_datatype_list(cfg):
    """
    get list of unique input data types

    Parameters
    ----------
    cfg : dict
        config dictionary

    Returns
    -------
    datatypesdescr : list
        list of data type descriptors

    """
    datatypesdescr = set()

    for datasetdescr in cfg['dataSetList']:
        proclevel, dataset = get_datasetfields(datasetdescr)
        if 'datatype' in cfg[dataset]:
            if isinstance(cfg[dataset]['datatype'], str):
                datagroup, datatype, dataset_save, product_save = (
                    get_datatypefields(cfg[dataset]['datatype']))
                if datagroup != 'PROC':
                    datatypesdescr.add(cfg[dataset]['datatype'])
            else:
                for datatype in cfg[dataset]['datatype']:
                    datagroup, datatype_aux, dataset_save, product_save = (
                        get_datatypefields(datatype))
                    if datagroup != 'PROC':
                        datatypesdescr.add(datatype)

    datatypesdescr = list(datatypesdescr)

    return datatypesdescr


def _get_datasets_list(cfg):
    """
    get list of dataset at each processing level

    Parameters
    ----------
    cfg : dict
        config dictionary

    Returns
    -------
    dataset_levels : dict
        a dictionary containing the list of datasets at each processing level

    """
    dataset_levels = dict({'l0': list()})
    for datasetdescr in cfg['dataSetList']:
        proclevel, dataset = get_datasetfields(datasetdescr)
        if proclevel in dataset_levels:
            dataset_levels[proclevel].append(dataset)
        else:
            dataset_levels.update({proclevel: [dataset]})

    return dataset_levels


def _get_masterfile_list(masterscan, datatypesdescr, starttime, endtime,
                         datacfg):
    """
    get master file list

    Parameters
    ----------
    masterscan : str
        name of the master scan
    datatypesdescr : list
        list of unique data type descriptors
    starttime, endtime : datetime object
        start and end of processing period
    datacfg : dict
        data configuration dictionary

    Returns
    -------
    masterfilelist : list
        the list of master files
    masterdatatypedescr : str
        the master data type descriptor

    """
    masterdatatypedescr = None
    for datatypedescr in datatypesdescr:
        datagroup, datatype, dataset, product = get_datatypefields(
            datatypedescr)
        if (datagroup != 'COSMO') and (datagroup != 'RAD4ALPCOSMO'):
            masterdatatypedescr = datatypedescr
            break

    # if only data type is COSMO use dBZ as reference
    if masterdatatypedescr is None:
        for datatypedescr in datatypesdescr:
            datagroup, datatype, dataset, product = (
                get_datatypefields(datatypedescr))
            if datagroup == 'COSMO':
                masterdatatypedescr = 'RAINBOW:dBZ'
                break
            elif datagroup == 'RAD4ALPCOSMO':
                masterdatatypedescr = 'RAD4ALP:dBZ'
                break

    masterfilelist = get_file_list(
        masterscan, masterdatatypedescr, starttime, endtime, datacfg)

    return masterfilelist, masterdatatypedescr


def _add_dataset(new_dataset, radar, make_global=True):
    """
    adds a new field to an existing radar object

    Parameters
    ----------
    new_dataset : radar object
        the radar object containing the new fields
    radar : radar object
        the radar object containing the global data
    make_global : boolean
        if true a new field is added to the global data

    Returns
    -------
    None

    """
    if radar is None:
        return

    if not make_global:
        return

    for field in new_dataset.fields:
        print('Adding field: '+field)
        radar.add_field(
            field, new_dataset.fields[field],
            replace_existing=True)


def _process_dataset(cfg, dataset, proc_status=0, radar=None, voltime=None):
    """
    processes a dataset

    Parameters
    ----------
    cfg : dict
        configuration dictionary
    dataset : str
        name of the dataset to be processed
    proc_status : int
        status of the processing 0: Initialization 1: process of radar volume
        2: Final processing
    radar : radar object
        radar object containing the data to be processed
    voltime : datetime object
        reference time of the radar

    Returns
    -------
    None

    """
    dscfg = _create_dscfg_dict(cfg, dataset, voltime=voltime)
    proc_func_name, dsformat = get_process_type(dscfg['type'])
    proc_func = getattr(proc, proc_func_name)
    new_dataset = proc_func(proc_status, dscfg, radar=radar)
    if new_dataset is not None:
        _add_dataset(new_dataset, radar, make_global=dscfg['MAKE_GLOBAL'])

        # create the data set products
        if 'products' in cfg[dataset]:
            for product in cfg[dataset]['products']:
                prdcfg = _create_prdcfg_dict(
                    cfg, dataset, product, voltime=voltime)
                prod_func_name = get_product_type(dsformat)
                prod_func = getattr(proc, prod_func_name)
                result = prod_func(new_dataset, prdcfg)
