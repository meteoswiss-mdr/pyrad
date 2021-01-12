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
from ..io.read_data_dem import read_dem, dem2radar_data

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
            dem_data = read_dem(fname, field_name)

            dscfg['global_data'] = {
                'dem_data': dem_data,
                'dem_field': None}

            if regular_grid:
                dscfg['global_data']['dem_field'] = dem2radar_data(
                    radar, dem_data, field_name)

            dscfg['initialized'] = 1

        dem_data = dscfg['global_data']['dem_data']
    else:
        dem_data = read_dem(fname, field_name)
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





def process_gecsx(procstatus, dscfg, radar_list=None):
    """
    Computes ground clutter RCS, radar visibility and many others using the
    GECSX algorithmn translated from IDL into python

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of string. Dataset keyword
            The input data types

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
        
    ind_rad = int(radarnr[5:8])-1
   
    fname = dscfg['dempath'][ind_rad] + dscfg['demfile']
    dem_data = read_dem(fname)

    # If no radar data is provided we create empty radar object from user
    # specification
    if not len(radar_list):
        ranges = np.arange(dscfg['range_resolution'] / 2, dscfg['rmax'],
                           dscfg['range_resolution'])
        azimuths = np.arange(dscfg['azmin'], dscfg['azmax'], 
                             dscfg['anglestep'])
        elevations = dscfg['antenna_elevations'] 
        
        radar = pyart.testing.make_empty_ppi_radar(len(ranges), len(azimuths), 
                                                   len(elevations))
        radar.latitude['data'] = np.array([dscfg['RadarPosition']
                                           ['latitude']])
        radar.latitude['data'] = np.array([dscfg['RadarPosition']
                                           ['longitude']])
        radar.altitude['data'] = np.array([dscfg['RadarPosition']
                                           ['altitude']])
        radar.azimuth['data'] = azimuths
        radar.range['data'] =  ranges
        radar.fixed_angle['data'] =  np.array(elevations)
        radar.elevation['data'] = np.array([len(azimuths) * [e] 
                                            for e in elevations]).ravel()
    
    # Create dict with radar specifications
    radar_specs = {}
    radar_specs['frequency'] = dscfg['frequency']
    radar_specs['loss'] = dscfg['lrxh'] + dscfg['mflossh']
    radar_specs['power'] = dscfg['txpwrh']
    radar_specs['tau'] = dscfg['pulse_width']
    radar_specs['beamwidth'] = dscfg['radar_beam_width_h']
    radar_specs['gain'] = dscfg['AntennaGainH']
    
    az_conv = dscfg.get('AzimTol', 0)
    ke = dscfg.get('refcorr', 4/3.)
    atm_att = dscfg.get('attg', 0.012)
    mosotti_kw = dscfg.get('mosotti_factor', 0.9644)
    sigma0_method = dscfg.get('sigma0_method', 'Gabella')
    raster_oversampling = dscfg.get('raster_oversampling', 1)
    verbose = dscfg.get('verbose', 1)
    clip = dscfg.get('clip', 1)
    daz = dscfg.get('az_discretization', 0.2)
    dr = dscfg.get('range_discretization', 100)
    
    (bent_terrain_altitude_dic, terrain_slope_dic, 
    terrain_aspect_dic, elevation_dic, min_vis_elevation_dic,
    min_vis_altitude_dic, visibility_dic, incident_angle_dic,
    effective_area_dic, sigma_0_dic, rcs_clutter_dic, dBm_clutter_dic,
    dBZ_clutter_dic, visibility_polar_dic) = pyart.retrieve.gecsx(radar, 
                                     radar_specs,
                                     dem_data,
                                     fill_value = None,      
                                     az_conv = az_conv ,
                                     dr = dr,
                                     daz = daz,
                                     ke = ke,
                                     atm_att = atm_att,
                                     mosotti_kw = mosotti_kw,
                                     raster_oversampling = raster_oversampling,
                                     sigma0_method = sigma0_method,
                                     clip = clip,
                                     return_pyart_objects = False,
                                     verbose = verbose)
    import pdb
    pdb.set_trace()
    # return new_dataset, ind_rad