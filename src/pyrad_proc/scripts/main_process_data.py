#!/usr/bin/env python

"""
================================================
Pyrad: The MeteoSwiss Radar Processing framework
================================================

Welcome to Pyrad!

To run the processing framework type:
    python main_process_data.py \
[config_file] [process_start_time] [process_end_time]
Example:
    python main_process_data.py \
'/home/lom/users/fvj/pyrad/config/processing/paradiso_fvj_vol.txt' \
'20140523000000' '20140523001000'

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import argparse
import os

import pyart
import pyrad

print(__doc__)


if __name__ == '__main__':

    # parse the arguments
    parser = argparse.ArgumentParser(
        description='Entry to Pyrad processing framework')

    # positional arguments
    parser.add_argument(
        'proccfgfile', type=str, help='path to main configuration file')
    parser.add_argument(
        'starttime', type=str,
        help='starting time of the data to be processed')
    parser.add_argument(
        'endtime', type=str, help='end time of the data to be processed ')

    args = parser.parse_args()

    print('config file: '+args.proccfgfile)
    print('start time: '+args.starttime)
    print('end time: '+args.endtime)

    proc_starttime = datetime.datetime.strptime(
        args.starttime, '%Y%m%d%H%M%S')
    proc_endtime = datetime.datetime.strptime(args.endtime, '%Y%m%d%H%M%S')

    # create config dictionary
    cfg = dict({'configFile': args.proccfgfile})
    cfg = pyrad.io.read_config(cfg['configFile'], cfg=cfg)
    cfg = pyrad.io.read_config(cfg['locationConfigFile'], cfg=cfg)
    cfg = pyrad.io.read_config(cfg['productConfigFile'], cfg=cfg)

    # configuration dictionary to figure out where the data is
    datacfg = dict({'datapath': cfg['datapath']})
    datacfg.update({'ScanList': cfg['ScanList']})
    if 'cosmopath' in cfg:
        datacfg.update({'cosmopath': cfg['cosmopath']})
    else:
        datacfg.update({'cosmopath': None})
    if 'dempath' in cfg:
        datacfg.update({'dempath': cfg['dempath']})
    else:
        datacfg.update({'dempath': None})
    if 'loadbasepath' in cfg:
        datacfg.update({'loadbasepath': cfg['loadbasepath']})
    else:
        datacfg.update({'loadbasepath': None})
    if 'loadname' in cfg:
        datacfg.update({'loadname': cfg['loadname']})
    else:
        datacfg.update({'loadname': None})
    if 'RadarName' in cfg:
        datacfg.update({'RadarName': cfg['RadarName']})
    else:
        datacfg.update({'RadarName': None})
    if 'RadarRes' in cfg:
        datacfg.update({'RadarRes': cfg['RadarRes']})
    else:
        datacfg.update({'RadarRes': None})
    if 'ScanPeriod' in cfg:
        datacfg.update({'ScanPeriod': int(cfg['ScanPeriod'])})
    else:
        print(
            'WARNING: Scan period not specified. \
            Assumed default value 5 min')
        datacfg.update({'ScanPeriod': 5})
    if 'CosmoRunFreq' in cfg:
        datacfg.update({'CosmoRunFreq': int(cfg['CosmoRunFreq'])})
    else:
        print(
            'WARNING: COSMO run frequency not specified. \
            Assumed default value 3h')
        datacfg.update({'CosmoRunFreq': 3})
    if 'CosmoForecasted' in cfg:
        datacfg.update({'CosmoForecasted': int(cfg['CosmoForecasted'])})
    else:
        print(
            'WARNING: Hours forecasted by COSMO not specified. \
             Assumed default value 7h (including analysis)')
        datacfg.update({'CosmoForecasted': 7})

    # get unique initial data types list
    datatypesdescr = set()

    for datasetdescr in cfg['dataSetList']:
        proclevel, dataset = pyrad.io.get_datasetfields(datasetdescr)
        if isinstance(cfg[dataset]['datatype'], str):
            datagroup, datatype, dataset_save, product_save = (
                pyrad.io.get_datatypefields(cfg[dataset]['datatype']))
            if datagroup != 'PROC':
                datatypesdescr.add(cfg[dataset]['datatype'])
        else:
            for datatype in cfg[dataset]['datatype']:
                datagroup, datatype_aux, dataset_save, product_save = (
                    pyrad.io.get_datatypefields(datatype))
                if datagroup != 'PROC':
                    datatypesdescr.add(datatype)

    datatypesdescr = list(datatypesdescr)

    # get number of processing levels and datasets corresponding to
    # each processing level
    dataset_levels = dict({'l0': list()})
    for datasetdescr in cfg['dataSetList']:
        proclevel, dataset = pyrad.io.get_datasetfields(datasetdescr)
        if proclevel in dataset_levels:
            dataset_levels[proclevel].append(dataset)
        else:
            dataset_levels.update({proclevel: [dataset]})

    # get lists of files to process using as reference a master scan
    for datatypedescr in datatypesdescr:
        datagroup, datatype, dataset, product = pyrad.io.get_datatypefields(
            datatypedescr)
        if (datagroup != 'COSMO') and (datagroup != 'RAD4ALPCOSMO'):
            masterdatatypedescr = datatypedescr
            masterscan = cfg['ScanList'][0]
            break

    masterfilelist = pyrad.io.get_file_list(
        masterscan, masterdatatypedescr, proc_starttime, proc_endtime,
        datacfg)
    nvolumes = len(masterfilelist)
    if nvolumes == 0:
        raise Exception(
            "ERROR: Could not find any volume within the specified times")
    print('Number of volumes to process: '+str(nvolumes)+'\n\n')

    # initial processing of the datasets
    print('\n\nInitializing datasets:')
    for level in sorted(dataset_levels):
        print('\nProcess level: '+level)
        for dataset in dataset_levels[level]:
            print('Processing dataset: '+dataset)

            dscfg = cfg[dataset]
            if 'MAKE_GLOBAL' not in dscfg:
                    dscfg.update({'MAKE_GLOBAL': 0})

            proc_func_name, dsformat = pyrad.proc.get_process_type(
                dscfg['type'])
            proc_func = getattr(pyrad.proc, proc_func_name)
            new_dataset = proc_func(0, dscfg, radar=None)

            if new_dataset is not None:
                if dscfg['MAKE_GLOBAL'] == 1:
                        for field in new_dataset.fields:
                            print('Adding field: '+field)
                            radar.add_field(
                                field, new_dataset.fields[field],
                                replace_existing=True)

                # create the data set products
                if 'products' in cfg[dataset]:
                    for product in cfg[dataset]['products']:
                        prdcfg = cfg[dataset]['products'][product]
                        prdcfg.update({'procname': cfg['name']})
                        prdcfg.update({'basepath': cfg['saveimgbasepath']})
                        prdcfg.update({'imgformat': cfg['imgformat']})
                        prdcfg.update({'convertformat': cfg['convertformat']})
                        prdcfg.update(
                            {'ppiImageConfig': cfg['ppiImageConfig']})
                        prdcfg.update(
                            {'rhiImageConfig': cfg['rhiImageConfig']})
                        prdcfg.update({'dsname': dataset})
                        prdcfg.update({'dstype': cfg[dataset]['type']})
                        prdcfg.update({'prdname': product})
                        prdcfg.update({'timeinfo': voltime})

                        prod_func_name = pyrad.proc.get_product_type(dsformat)
                        prod_func = getattr(pyrad.proc, prod_func_name)
                        result = prod_func(new_dataset, prdcfg)

    # process all data files in file list
    for masterfile in masterfilelist:
        print('\n\nmaster file: '+os.path.basename(masterfile))
        voltime = pyrad.io.get_datetime(masterfile, masterdatatypedescr)

        # get all raw data
        radar = pyrad.io.get_data(voltime, datatypesdescr, datacfg)

        # process all data sets
        for level in sorted(dataset_levels):
            print('\nProcess level: '+level)
            for dataset in dataset_levels[level]:
                print('Processing dataset: '+dataset)

                dscfg = cfg[dataset]
                if 'MAKE_GLOBAL' not in dscfg:
                    dscfg.update({'MAKE_GLOBAL': 0})

                proc_func_name, dsformat = pyrad.proc.get_process_type(
                    dscfg['type'])
                proc_func = getattr(pyrad.proc, proc_func_name)
                new_dataset = proc_func(1, dscfg, radar=radar)

                if new_dataset is not None:
                    if dscfg['MAKE_GLOBAL'] == 1:
                        for field in new_dataset.fields:
                            print('Adding field: '+field)
                            radar.add_field(
                                field, new_dataset.fields[field],
                                replace_existing=True)

                    # create the data set products
                    if 'products' in cfg[dataset]:
                        for product in cfg[dataset]['products']:
                            prdcfg = cfg[dataset]['products'][product]
                            prdcfg.update({'procname': cfg['name']})
                            prdcfg.update(
                                {'basepath': cfg['saveimgbasepath']})
                            prdcfg.update({'imgformat': cfg['imgformat']})
                            prdcfg.update(
                                {'convertformat': cfg['convertformat']})
                            prdcfg.update(
                                {'ppiImageConfig': cfg['ppiImageConfig']})
                            prdcfg.update(
                                {'rhiImageConfig': cfg['rhiImageConfig']})
                            prdcfg.update({'dsname': dataset})
                            prdcfg.update({'dstype': cfg[dataset]['type']})
                            prdcfg.update({'prdname': product})
                            prdcfg.update({'timeinfo': voltime})

                            prod_func_name = pyrad.proc.get_product_type(
                                dsformat)
                            prod_func = getattr(pyrad.proc, prod_func_name)
                            result = prod_func(new_dataset, prdcfg)

    # post-processing of the datasets
    print('\n\nPost-processing datasets')
    for level in sorted(dataset_levels):
        print('\nProcess level: '+level)
        for dataset in dataset_levels[level]:
            print('Processing dataset: '+dataset)

            dscfg = cfg[dataset]
            if 'MAKE_GLOBAL' not in dscfg:
                    dscfg.update({'MAKE_GLOBAL': 0})

            proc_func_name, dsformat = pyrad.proc.get_process_type(
                dscfg['type'])
            proc_func = getattr(pyrad.proc, proc_func_name)
            new_dataset = proc_func(2, dscfg, radar=None)

            if new_dataset is not None:
                if dscfg['MAKE_GLOBAL'] == 1:
                    for field in new_dataset.fields:
                        print('Adding field: '+field)
                        radar.add_field(
                            field, new_dataset.fields[field],
                            replace_existing=True)

                # create the data set products
                if 'products' in cfg[dataset]:
                    for product in cfg[dataset]['products']:
                        prdcfg = cfg[dataset]['products'][product]
                        prdcfg.update({'procname': cfg['name']})
                        prdcfg.update({'basepath': cfg['saveimgbasepath']})
                        prdcfg.update({'imgformat': cfg['imgformat']})
                        prdcfg.update({'convertformat': cfg['convertformat']})
                        prdcfg.update(
                            {'ppiImageConfig': cfg['ppiImageConfig']})
                        prdcfg.update(
                            {'rhiImageConfig': cfg['rhiImageConfig']})
                        prdcfg.update({'dsname': dataset})
                        prdcfg.update({'dstype': cfg[dataset]['type']})
                        prdcfg.update({'prdname': product})
                        prdcfg.update({'timeinfo': voltime})

                        prod_func_name = pyrad.proc.get_product_type(dsformat)
                        prod_func = getattr(pyrad.proc, prod_func_name)
                        result = prod_func(new_dataset, prdcfg)

    print('\n\n\nThis is the end my friend! See you soon!')
