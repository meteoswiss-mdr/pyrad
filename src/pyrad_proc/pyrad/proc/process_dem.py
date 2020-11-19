"""
pyrad.proc.process_dem
======================

Functions to manage DEM data

.. autosummary::
    :toctree: generated/

    process_dem
    process_visibility

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_dem import read_idrisi_data, dem2radar_data

# from memory_profiler import profile


def process_dem(procstatus, dscfg, radar_list=None):
    """
    Gets DEM data and put it in radar coordinates

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
            the COSMO field in radar coordinates in memory. Default False
        regular_grid : int. Dataset keyword
            if set it is assume that the radar has a grid constant in time and
            there is no need to compute a new COSMO field if the COSMO
            data has not changed. Default False
        dem_field : str. Dataset keyword
            name of the DEM field to process
        demfile : str. Dataset keyword
            Name of the file containing the DEM data
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    # debugging
    # start_time = time.time()

    for datatypedescr in dscfg['datatype']:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    keep_in_memory = dscfg.get('keep_in_memory', 0)
    regular_grid = dscfg.get('regular_grid', 0)
    field_name = get_fieldname_pyart(dscfg['dem_field'])

    fname = dscfg['dempath'][ind_rad]+dscfg['demfile']

    if keep_in_memory:
        if dscfg['initialized'] == 0:
            dem_data = read_idrisi_data(fname, field_name)

            dscfg['global_data'] = {
                'dem_data': dem_data,
                'dem_field': None}

            if regular_grid:
                dscfg['global_data']['dem_field'] = dem2radar_data(
                    radar, dem_data, field_name)

            dscfg['initialized'] = 1

        dem_data = dscfg['global_data']['dem_data']
    else:
        dem_data = read_idrisi_data(fname, field_name)
        if dem_data is None:
            warn('DEM data not found')
            return None, None

    if regular_grid:
        print('DEM field already in memory')
        dem_field = dscfg['global_data']['dem_fields']
    else:
        dem_field = dem2radar_data(radar, dem_data, field_name=field_name)
        if dem_field is None:
            warn('Unable to obtain DEM fields')
            return None, None

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field(field_name, dem_field)

    return new_dataset, ind_rad


def process_visibility(procstatus, dscfg, radar_list=None):
    """
    Gets the visibility in percentage from the minimum visible elevation.
    Anything with elevation lower than the minimum visible elevation plus and
    offset is set to 0 while above is set to 100.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            arbitrary data type
        offset : float. Dataset keyword
            The offset above the minimum visibility that must be filtered
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
        if datatype == 'minvisel':
            minvisel_field = get_fieldname_pyart(datatype)
            break

    ind_rad = int(radarnr[5:8])-1
    if radar_list[ind_rad] is None:
        warn('No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    offset = dscfg.get('offset', 0.)

    minvisel_data = radar.fields[minvisel_field]['data']+offset
    ele_data = np.broadcast_to(
        radar.elevation['data'].reshape(radar.nrays, 1),
        (radar.nrays, radar.ngates))

    vis_dict = pyart.config.get_metadata('visibility')
    vis_dict['data'] = 100.*np.ma.greater_equal(
        ele_data, minvisel_data, dtype=float)

    # if a gate has visibility 0 all the subsequent gates in the ray
    # are set to 0
    for ray in range(radar.nrays):
        ind = np.where(vis_dict['data'][ray, :] == 0.)[0]
        if ind.size > 0:
            vis_dict['data'][ray, ind[0]:] = 0.

    # prepare for exit
    new_dataset = {'radar_out': deepcopy(radar)}
    new_dataset['radar_out'].fields = dict()
    new_dataset['radar_out'].add_field('visibility', vis_dict)

    return new_dataset, ind_rad
