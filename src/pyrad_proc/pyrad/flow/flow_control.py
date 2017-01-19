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
    _warning_format

"""
from __future__ import print_function
import sys
import warnings
from warnings import warn
import traceback
import os
from datetime import datetime
from datetime import timedelta
import atexit
import inspect

from ..io.config import read_config
from ..io.read_data_radar import get_data
from ..io.io_aux import get_datetime, get_file_list, get_scan_list
from ..io.io_aux import get_dataset_fields, get_datatype_fields
from ..io.trajectory import Trajectory

from ..proc.process_aux import get_process_func
from ..prod.product_aux import get_prodgen_func

from pyrad import proc, prod
from pyrad import version as pyrad_version
from pyart import version as pyart_version


def main(cfgfile, starttime, endtime, infostr="", trajfile=""):
    """
    main flow control. Processes data over a given period of time

    Parameters
    ----------
    cfgfile : str
        path of the main config file
    starttime, endtime : datetime object
        start and end time of the data to be processed
    trajfile : str
        path to file describing the trajectory
    infostr : Information string about the actual data processing
              (e.g. 'RUN57'). This sting is added to product files.

    """

    print("- PYRAD version: %s (compiled %s by %s)" %
          (pyrad_version.version, pyrad_version.compile_date_time,
           pyrad_version.username))
    print("- PYART version: " + pyart_version.version)

    # Definie behaviour of warnings
    warnings.simplefilter('always')  # always print matching warnings
    # warnings.simplefilter('error')  # turn matching warnings into exceptions
    warnings.formatwarning = _warning_format  # define format

    cfg = _create_cfg_dict(cfgfile)
    datacfg = _create_datacfg_dict(cfg)

    # Get the plane trajectory
    if (len(trajfile) > 0):
        print("- Trajectory file: " + trajfile)
        try:
            traj = Trajectory(trajfile, starttime=starttime, endtime=endtime)
        except Exception as ee:
            warn(str(ee))
            sys.exit(1)

        # Derive start and end time (if not specified by arguments)
        if (starttime is None):
            scan_min = cfg['ScanPeriod'] * 1.1  # [min]
            starttime = traj.get_start_time() - timedelta(minutes=scan_min)
        if (endtime is None):
            scan_min = cfg['ScanPeriod'] * 1.1  # [min]
            endtime = traj.get_end_time() + timedelta(minutes=scan_min)
    else:
        traj = None

    print("- Start time: " + str(starttime))
    print("- End time: " + str(endtime))

    if (len(infostr) > 0):
        print('- Info string : ' + infostr)

    datatypesdescr_list = list()
    for i in range(1, cfg['NumRadars']+1):
        datatypesdescr_list.append(
            _get_datatype_list(cfg, radarnr='RADAR'+'{:03d}'.format(i)))

    dataset_levels = _get_datasets_list(cfg)
    masterfilelist, masterdatatypedescr, masterscan = _get_masterfile_list(
        datatypesdescr_list[0], starttime, endtime, datacfg,
        scan_list=cfg['ScanList'])

    nvolumes = len(masterfilelist)
    if nvolumes == 0:
        raise ValueError("ERROR: Could not find any valid volumes between " +
                         starttime.strftime('%Y-%m-%d %H:%M:%S') + " and " +
                         endtime.strftime('%Y-%m-%d %H:%M:%S') + " for " +
                         "master scan '" + str(masterscan) +
                         "' and master data type '" + masterdatatypedescr +
                         "'")
    print('- Number of volumes to process: ' + str(nvolumes))

    # initial processing of the datasets
    print('- Initializing datasets:')

    dscfg = dict()
    for level in sorted(dataset_levels):
        print('-- Process level: '+level)
        for dataset in dataset_levels[level]:
            print('--- Processing dataset: '+dataset)
            dscfg.update({dataset: _create_dscfg_dict(cfg, dataset)})
            try:
                result = _process_dataset(cfg, dscfg[dataset], proc_status=0,
                                          radar_list=None, voltime=None,
                                          trajectory=traj, runinfo=infostr)
            except Exception as ee:
                traceback.print_exc()
                continue

    # process all data files in file list
    for masterfile in masterfilelist:
        print('- master file: ' + os.path.basename(masterfile))
        master_voltime = get_datetime(masterfile, masterdatatypedescr)

        # get data of master radar
        radar_list = list()
        radar_list.append(
            get_data(master_voltime, datatypesdescr_list[0], datacfg))

        # get data of rest of radars
        for i in range(1, cfg['NumRadars']):
            filelist_ref, datatypedescr_ref, scan_ref = _get_masterfile_list(
                datatypesdescr_list[i],
                master_voltime-timedelta(seconds=cfg['TimeTol']),
                master_voltime+timedelta(seconds=cfg['TimeTol']),
                datacfg, scan_list=cfg['ScanList'])

            if len(filelist_ref) == 0:
                warn("ERROR: Could not find any valid volume for reference " +
                     "time " + master_voltime.strftime('%Y-%m-%d %H:%M:%S') +
                     ' and radar RADAR'+'{:03d}'.format(i+1))
                radar_list.append(None)
            else:
                voltime_ref = get_datetime(filelist_ref[0], datatypedescr_ref)
                radar_list.append(
                    get_data(voltime_ref, datatypesdescr_list[i], datacfg))

        # process all data sets
        for level in sorted(dataset_levels):
            print('-- Process level: '+level)
            for dataset in dataset_levels[level]:
                print('--- Processing dataset: '+dataset)
                try:
                    result = _process_dataset(cfg, dscfg[dataset],
                                              proc_status=1,
                                              radar_list=radar_list,
                                              voltime=master_voltime,
                                              trajectory=traj, runinfo=infostr)
                except Exception as ee:
                    traceback.print_exc()
                    continue

    # post-processing of the datasets
    print('- Post-processing datasets:')
    for level in sorted(dataset_levels):
        print('-- Process level: '+level)
        for dataset in dataset_levels[level]:
            print('--- Processing dataset: '+dataset)
            try:
                result = _process_dataset(cfg, dscfg[dataset], proc_status=2,
                                          radar_list=None, voltime=None,
                                          trajectory=traj, runinfo=infostr)
            except Exception as ee:
                traceback.print_exc()
                continue

    print('- This is the end my friend! See you soon!')


def _process_dataset(cfg, dscfg, proc_status=0, radar_list=None, voltime=None,
                     trajectory=None, runinfo=None):
    """
    processes a dataset

    Parameters
    ----------
    cfg : dict
        configuration dictionary
    dscfg : dict
        dataset specific configuration dictionary
    proc_status : int
        status of the processing 0: Initialization 1: process of radar volume
        2: Final processing
    radar : radar object
        radar object containing the data to be processed
    voltime : datetime object
        reference time of the radar
    trajectory : Trajectory object
        containing trajectory samples

    Returns
    -------
    0 if a new dataset has been created. None otherwise

    """

    dscfg['timeinfo'] = voltime
    try:
        proc_ds_func, dsformat = get_process_func(dscfg['type'],
                                                  dscfg['dsname'])
    except Exception as ee:
        raise

    if (type(proc_ds_func) is str):
        proc_ds_func = getattr(proc, proc_ds_func)

    # Create dataset
    if ('trajectory' in inspect.getfullargspec(proc_ds_func).args):
        new_dataset, ind_rad = proc_ds_func(proc_status, dscfg,
                                            radar_list=radar_list,
                                            trajectory=trajectory)
    else:
        new_dataset, ind_rad = proc_ds_func(proc_status, dscfg,
                                            radar_list=radar_list)

    if new_dataset is None:
        return None

    result = _add_dataset(new_dataset, radar_list, ind_rad,
                          make_global=dscfg['MAKE_GLOBAL'])

    try:
        prod_func = get_prodgen_func(dsformat, dscfg['dsname'],
                                     dscfg['type'])
    except Exception as ee:
        raise

    # create the data set products
    if 'products' in dscfg:
        for product in dscfg['products']:
            print('---- Processing product: ' + product)
            prdcfg = _create_prdcfg_dict(cfg, dscfg['dsname'], product,
                                         voltime, runinfo=runinfo)
            try:
                result = prod_func(new_dataset, prdcfg)
            except Exception as ee:
                traceback.print_exc()
                continue

    return 0


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
    try:
        print("- Main config file : %s" % cfgfile)
        cfg = read_config(cfg['configFile'], cfg=cfg)
        print("- Location config file : %s" % cfg['locationConfigFile'])
        cfg = read_config(cfg['locationConfigFile'], cfg=cfg)
        print("- Product config file : %s" % cfg['productConfigFile'])
        cfg = read_config(cfg['productConfigFile'], cfg=cfg)
    except Exception as ee:
        warn(str(ee))
        sys.exit(1)

    # check for mandatory config parameters
    param_must = ['name', 'configpath', 'saveimgbasepath', 'dataSetList']
    for param in param_must:
        if param not in cfg:
            raise Exception("ERROR config: Parameter '%s' undefined!" % param)

    # fill in defaults
    if 'NumRadars' not in cfg:
        cfg.update({'NumRadars': 1})
    if 'TimeTol' not in cfg:
        cfg.update({'TimeTol': 3600.})
    if 'ScanList' not in cfg:
        cfg.update({'ScanList': None})
    else:
        cfg.update({'ScanList': get_scan_list(cfg['ScanList'])})
    if 'datapath' not in cfg:
        cfg.update({'datapath': None})
    if 'cosmopath' not in cfg:
        cfg.update({'cosmopath': None})
    if 'psrpath' not in cfg:
        cfg.update({'psrpath': None})
    if 'colocgatespath' not in cfg:
        cfg.update({'colocgatespath': None})
    if 'dempath' not in cfg:
        cfg.update({'dempath': None})
    if 'smnpath' not in cfg:
        cfg.update({'smnpath': None})
    if 'disdropath' not in cfg:
        cfg.update({'disdropath': None})
    if 'solarfluxpath' not in cfg:
        cfg.update({'solarfluxpath': None})
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
    if 'lrxh' not in cfg:
        cfg.update({'lrxh': None})
    if 'lrxv' not in cfg:
        cfg.update({'lrxv': None})
    if 'lradomeh' not in cfg:
        cfg.update({'lradomeh': None})
    if 'lradomev' not in cfg:
        cfg.update({'lradomev': None})
    if 'AntennaGain' not in cfg:
        cfg.update({'AntennaGain': None})
    if 'attg' not in cfg:
        cfg.update({'attg': None})
    if 'ScanPeriod' not in cfg:
        warn('WARNING: Scan period not specified. ' +
             'Assumed default value 5 min')
        cfg.update({'ScanPeriod': 5})
    if 'CosmoRunFreq' not in cfg:
        warn('WARNING: COSMO run frequency not specified. ' +
             'Assumed default value 3h')
        cfg.update({'CosmoRunFreq': 3})
    if 'CosmoForecasted' not in cfg:
        warn('WARNING: Hours forecasted by COSMO not specified. ' +
             'Assumed default value 7h (including analysis)')
        cfg.update({'CosmoForecasted': 7})

    # Convert the following strings to string arrays
    strarr_list = ['datapath', 'cosmopath', 'dempath', 'loadbasepath',
                   'loadname', 'RadarName', 'RadarRes', 'ScanList',
                   'imgformat']
    for param in strarr_list:
        if (type(cfg[param]) is str):
            cfg[param] = [cfg[param]]

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
    datacfg.update({'ScanPeriod': cfg['ScanPeriod']})
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
    dscfg.update({'solarfluxpath': cfg['solarfluxpath']})
    dscfg.update({'colocgatespath': cfg['colocgatespath']})
    dscfg.update({'mflossh': cfg['mflossh']})
    dscfg.update({'mflossv': cfg['mflossv']})
    dscfg.update({'radconsth': cfg['radconsth']})
    dscfg.update({'radconstv': cfg['radconstv']})
    dscfg.update({'lrxh': cfg['lrxh']})
    dscfg.update({'lrxv': cfg['lrxv']})
    dscfg.update({'lradomeh': cfg['lradomeh']})
    dscfg.update({'lradomev': cfg['lradomev']})
    dscfg.update({'AntennaGain': cfg['AntennaGain']})
    dscfg.update({'attg': cfg['attg']})
    dscfg.update({'basepath': cfg['saveimgbasepath']})
    dscfg.update({'procname': cfg['name']})
    dscfg.update({'dsname': dataset})
    dscfg.update({'timeinfo': None})
    if ('par_azimuth_antenna' in cfg):
        dscfg.update({'par_azimuth_antenna': cfg['par_azimuth_antenna']})
    if ('par_elevation_antenna' in cfg):
        dscfg.update({'par_elevation_antenna': cfg['par_elevation_antenna']})

    # indicates the dataset has been initialized and aux data is available
    dscfg.update({'initialized': False})
    dscfg.update({'global_data': None})

    if 'MAKE_GLOBAL' not in dscfg:
        dscfg.update({'MAKE_GLOBAL': 0})

    # Convert the following strings to string arrays
    strarr_list = ['datatype']
    for param in strarr_list:
        if param in dscfg:
            if (type(dscfg[param]) is str):
                dscfg[param] = [dscfg[param]]

    return dscfg


def _create_prdcfg_dict(cfg, dataset, product, voltime, runinfo=None):
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

    # Ugly copying of dataset config parameters to product
    # config dict. Better: Make dataset config dict available to
    # the product generation.
    prdcfg = cfg[dataset]['products'][product]
    prdcfg.update({'procname': cfg['name']})
    prdcfg.update({'basepath': cfg['saveimgbasepath']})
    prdcfg.update({'smnpath': cfg['smnpath']})
    prdcfg.update({'disdropath': cfg['disdropath']})
    prdcfg.update({'ScanPeriod': cfg['ScanPeriod']})
    prdcfg.update({'imgformat': cfg['imgformat']})
    if 'ppiImageConfig' in cfg:
        prdcfg.update({'ppiImageConfig': cfg['ppiImageConfig']})
    if 'rhiImageConfig' in cfg:
        prdcfg.update({'rhiImageConfig': cfg['rhiImageConfig']})
    if 'sunhitsImageConfig' in cfg:
        prdcfg.update({'sunhitsImageConfig': cfg['sunhitsImageConfig']})
    prdcfg.update({'dsname': dataset})
    prdcfg.update({'dstype': cfg[dataset]['type']})
    prdcfg.update({'prdname': product})
    prdcfg.update({'timeinfo': voltime})
    prdcfg.update({'runinfo': runinfo})
    if 'dssavename' in cfg[dataset]:
        prdcfg.update({'dssavename': cfg[dataset]['dssavename']})

    return prdcfg


def _get_datatype_list(cfg, radarnr='RADAR001'):
    """
    get list of unique input data types

    Parameters
    ----------
    cfg : dict
        config dictionary
    radarnr : str
        radar number identifier

    Returns
    -------
    datatypesdescr : list
        list of data type descriptors

    """
    datatypesdescr = set()

    for datasetdescr in cfg['dataSetList']:
        proclevel, dataset = get_dataset_fields(datasetdescr)
        if 'datatype' in cfg[dataset]:
            if isinstance(cfg[dataset]['datatype'], str):
                (radarnr_descr, datagroup, datatype_aux, dataset_save,
                 product_save) = (
                    get_datatype_fields(cfg[dataset]['datatype']))
                if datagroup != 'PROC' and radarnr_descr == radarnr:
                    if ((dataset_save is None) and (product_save is None)):
                        datatypesdescr.add(datagroup+":"+datatype_aux)
                    else:
                        datatypesdescr.add(datagroup + ":" + datatype_aux +
                                           "," + dataset_save + "," +
                                           product_save)
            else:
                for datatype in cfg[dataset]['datatype']:
                    (radarnr_descr, datagroup, datatype_aux, dataset_save,
                     product_save) = (
                        get_datatype_fields(datatype))
                    if datagroup != 'PROC' and radarnr_descr == radarnr:
                        if ((dataset_save is None) and (product_save is None)):
                            datatypesdescr.add(datagroup+":"+datatype_aux)
                        else:
                            datatypesdescr.add(datagroup + ":" + datatype_aux +
                                               "," + dataset_save + "," +
                                               product_save)

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
        proclevel, dataset = get_dataset_fields(datasetdescr)
        if proclevel in dataset_levels:
            dataset_levels[proclevel].append(dataset)
        else:
            dataset_levels.update({proclevel: [dataset]})

    return dataset_levels


def _get_masterfile_list(datatypesdescr, starttime, endtime, datacfg,
                         scan_list=None):
    """
    get master file list

    Parameters
    ----------
    datatypesdescr : list
        list of unique data type descriptors
    starttime, endtime : datetime object
        start and end of processing period
    datacfg : dict
        data configuration dictionary
    scan_list : list
        list of scans

    Returns
    -------
    masterfilelist : list
        the list of master files
    masterdatatypedescr : str
        the master data type descriptor

    """
    masterdatatypedescr = None
    masterscan = None
    for datatypedescr in datatypesdescr:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        if ((datagroup != 'COSMO') and (datagroup != 'RAD4ALPCOSMO') and
                (datagroup != 'DEM') and (datagroup != 'RAD4ALPDEM')):
            masterdatatypedescr = datatypedescr
            if scan_list is not None:
                masterscan = scan_list[int(radarnr[5:8])-1][0]
            break

    # if data type is not radar use dBZ as reference
    if masterdatatypedescr is None:
        for datatypedescr in datatypesdescr:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            if datagroup == 'COSMO':
                masterdatatypedescr = radarnr+':RAINBOW:dBZ'
                if scan_list is not None:
                    masterscan = scan_list[int(radarnr[5:8])-1][0]
                break
            elif datagroup == 'RAD4ALPCOSMO':
                masterdatatypedescr = radarnr+':RAD4ALP:dBZ'
                if scan_list is not None:
                    masterscan = scan_list[int(radarnr[5:8])-1][0]
                break
            elif datagroup == 'DEM':
                masterdatatypedescr = radarnr+':RAINBOW:dBZ'
                if scan_list is not None:
                    masterscan = scan_list[int(radarnr[5:8])-1][0]
                break
            elif datagroup == 'RAD4ALPDEM':
                masterdatatypedescr = radarnr+':RAD4ALP:dBZ'
                if scan_list is not None:
                    masterscan = scan_list[int(radarnr[5:8])-1][0]
                break

    masterfilelist = get_file_list(
        masterdatatypedescr, starttime, endtime, datacfg,
        scan=masterscan)

    return masterfilelist, masterdatatypedescr, masterscan


def _add_dataset(new_dataset, radar_list, ind_rad, make_global=True):
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
    0 if successful. None otherwise

    """
    if radar_list is None:
        return None

    if not make_global:
        return None

    for field in new_dataset.fields:
        print('Adding field: '+field)
        radar_list[ind_rad].add_field(
            field, new_dataset.fields[field],
            replace_existing=True)
    return 0


def _warning_format(message, category, filename, lineno, file=None, line=None):
    return '%s (%s:%s)\n' % (message, filename, lineno)
