"""
pyrad.proc.process_cosmo
===========================

Functions to manage COSMO data

.. autosummary::
    :toctree: generated/

    process_cosmo_temp
    process_cosmo_temp_lookup_table
    process_cosmo_coord

"""

from copy import deepcopy
from warnings import warn
import time

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, find_raw_cosmo_file
from ..io.read_data_cosmo import read_cosmo_temp, read_cosmo_coord
from ..io.read_data_cosmo import cosmo2radar_data, cosmo2radar_coord
from ..io.read_data_cosmo import get_cosmo_field
from ..io.read_data_radar import interpol_field

from netCDF4 import num2date


def process_cosmo_temp(procstatus, dscfg, radar_list=None):
    """
    Gets COSMO temperature data and put it in radar coordinates

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type
        keep_in_memory : int. Dataset keyword
            if set keeps the COSMO data dict, the COSMO coordinates dict and
            the COSMO field in radar coordinates in memory
        regular_grid : int. Dataset keyword
            if set it is assume that the radar has a grid constant in time and
            there is no need to compute a new COSMO field if the COSMO
            temperature has not changed
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # debugging
    # start_time = time.time()

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    keep_in_memory = 0
    if 'keep_in_memory' in dscfg:
        keep_in_memory = dscfg['keep_in_memory']

    regular_grid = 0
    if 'regular_grid' in dscfg:
        regular_grid = dscfg['regular_grid']

    fname = find_raw_cosmo_file(dscfg['timeinfo'], 'TEMP', dscfg,
                                ind_rad=ind_rad)

    if fname is None:
        return None, None

    if keep_in_memory:
        if dscfg['initialized'] == 0:
            cosmo_coord = read_cosmo_coord(
                dscfg['cosmopath'][ind_rad] +
                'rad2cosmo1/cosmo-1_MDR_3D_const.nc')
            dscfg['global_data'] = {
                'cosmo_fname': None,
                'cosmo_temp': None,
                'cosmo_field': None,
                'time_index': None,
                'cosmo_coord': cosmo_coord}
            dscfg['initialized'] = 1

        cosmo_coord = dscfg['global_data']['cosmo_coord']
        if fname != dscfg['global_data']['cosmo_fname']:
            # debugging
            # start_time2 = time.time()
            cosmo_temp = read_cosmo_temp(fname, celsius=True)
            # print(" reading COSMO takes %s seconds " %
            #      (time.time() - start_time2))

            dscfg['global_data']['cosmo_temp'] = cosmo_temp
            dscfg['global_data']['cosmo_fname'] = fname
        else:
            print('raw COSMO data already in memory')
            cosmo_temp = dscfg['global_data']['cosmo_temp']
    else:
        cosmo_coord = read_cosmo_coord(
            dscfg['cosmopath'][ind_rad]+'rad2cosmo1/cosmo-1_MDR_3D_const.nc')

        # debugging
        # start_time2 = time.time()
        cosmo_temp = read_cosmo_temp(fname, celsius=True)
        # print(" reading COSMO takes %s seconds " %
        #      (time.time() - start_time2))

    dtcosmo = num2date(
        cosmo_temp['time']['data'][:], cosmo_temp['time']['units'])
    time_index = np.argmin(abs(dtcosmo-dscfg['timeinfo']))

    if keep_in_memory and regular_grid:
        if time_index != dscfg['global_data']['time_index']:
            cosmo_field = cosmo2radar_data(
                radar, cosmo_coord, cosmo_temp, time_index=time_index,
                field_name='temperature')
            dscfg['global_data']['time_index'] = time_index
            dscfg['global_data']['cosmo_field'] = cosmo_field
        else:
            print('COSMO field already in memory')
            cosmo_field = dscfg['global_data']['cosmo_field']
    else:
        cosmo_field = cosmo2radar_data(
            radar, cosmo_coord, cosmo_temp, time_index=time_index,
            field_name='temperature')

    # prepare for exit
    new_dataset = deepcopy(radar)
    new_dataset.fields = dict()
    new_dataset.add_field('temperature', cosmo_field)

    # debugging
    # print(" putting COSMO data in radar coordinates takes %s seconds " %
    #      (time.time() - start_time))

    return new_dataset, ind_rad


def process_cosmo_temp_lookup_table(procstatus, dscfg, radar_list=None):
    """
    Gets COSMO temperature data and put it in radar coordinates
    using look up tables computed or loaded when initializing

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type
        lookup_table : int. Dataset keyword
            if set a pre-computed look up table for the COSMO coordinates is
            loaded. Otherwise the look up table is computed taking the first
            radar object as reference
        regular_grid : int. Dataset keyword
            if set it is assume that the radar has a grid constant in time and
            therefore there is no need to interpolate the COSMO field in
            memory to the current radar grid
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # debugging
    # start_time = time.time()

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    regular_grid = 0
    if 'regular_grid' in dscfg:
        regular_grid = dscfg['regular_grid']

    lookup_table = 0
    if 'lookup_table' in dscfg:
        lookup_table = dscfg['lookup_table']

    fname = find_raw_cosmo_file(dscfg['timeinfo'], 'TEMP', dscfg,
                                ind_rad=ind_rad)

    if fname is None:
        return None, None

    if dscfg['initialized'] == 0:
        if lookup_table:
            savedir = dscfg['cosmopath'][ind_rad]+'rad2cosmo1/'
            fname_ind = 'rad2cosmo_cosmo_index_'+dscfg['procname']+'.nc'
            cosmo_radar = pyart.io.read_cfradial(savedir+fname_ind)
        else:
            cosmo_coord = read_cosmo_coord(
                dscfg['cosmopath'][ind_rad] +
                'rad2cosmo1/cosmo-1_MDR_3D_const.nc')
            cosmo_ind_field = cosmo2radar_coord(radar, cosmo_coord)
            cosmo_radar = deepcopy(radar)
            cosmo_radar.fields = dict()
            cosmo_radar.add_field('cosmo_index', cosmo_ind_field)

        dscfg['global_data'] = {
            'cosmo_fname': None,
            'cosmo_temp': None,
            'cosmo_field': None,
            'cosmo_radar': cosmo_radar,
            'time_index': None}
        dscfg['initialized'] = 1

    if fname != dscfg['global_data']['cosmo_fname']:
        # debugging
        # start_time2 = time.time()
        cosmo_temp = read_cosmo_temp(fname, celsius=True)
        # print(" reading COSMO takes %s seconds " %
        #      (time.time() - start_time2))
        dscfg['global_data']['cosmo_temp'] = cosmo_temp
        dscfg['global_data']['cosmo_fname'] = fname
    else:
        print('raw COSMO data already in memory')
        cosmo_temp = dscfg['global_data']['cosmo_temp']

    dtcosmo = num2date(
        cosmo_temp['time']['data'][:], cosmo_temp['time']['units'])
    time_index = np.argmin(abs(dtcosmo-dscfg['timeinfo']))

    if time_index != dscfg['global_data']['time_index']:
        # debugging
        # start_time3 = time.time()
        cosmo_field = get_cosmo_field(
            cosmo_temp,
            dscfg['global_data']['cosmo_radar'].fields['cosmo_index'],
            time_index=time_index,
            field_name='temperature')
        # print(" getting COSMO data takes %s seconds "
        #      % (time.time() - start_time3))

        dscfg['global_data']['time_index'] = time_index
        dscfg['global_data']['cosmo_field'] = cosmo_field
    else:
        print('COSMO field already in memory')
        cosmo_field = dscfg['global_data']['cosmo_field']

    # prepare for exit
    new_dataset = deepcopy(dscfg['global_data']['cosmo_radar'])
    new_dataset.add_field('temperature', cosmo_field)

    # interpolate to current radar grid
    if not regular_grid:
        cosmo_field_interp = interpol_field(radar, new_dataset, 'temperature')
        new_dataset = deepcopy(radar)
        new_dataset.fields = dict()
        new_dataset.add_field('temperature', cosmo_field_interp)

    # debugging
    # print(" putting COSMO data in radar coordinates takes %s seconds "
    #      % (time.time() - start_time))

    return new_dataset, ind_rad


def process_cosmo_coord(procstatus, dscfg, radar_list=None):
    """
    Gets the COSMO indices corresponding to each cosmo coordinates

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type
        cosmopath : string. General keyword
            path where to store the look up table
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : Radar
        radar object
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    start_time = time.time()

    if dscfg['initialized'] == 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    cosmo_coord = read_cosmo_coord(
        dscfg['cosmopath'][ind_rad]+'rad2cosmo1/cosmo-1_MDR_3D_const.nc')
    cosmo_ind_field = cosmo2radar_coord(
        radar, cosmo_coord, slice_xy=True, slice_z=False)

    # prepare for exit
    radar_obj = deepcopy(radar)
    radar_obj.fields = dict()
    radar_obj.add_field('cosmo_index', cosmo_ind_field)

    new_dataset = {
        'ind_rad': ind_rad,
        'radar_obj': radar_obj}

    dscfg['initialized'] = 1

    return new_dataset, ind_rad
