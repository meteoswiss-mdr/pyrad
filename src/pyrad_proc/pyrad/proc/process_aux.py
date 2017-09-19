"""
pyrad.proc.process_aux
=============================

Auxiliary functions. Functions to determine the process type, pass raw data to
the product generation functions, save radar data and extract data at
determined points or regions of interest.

.. autosummary::
    :toctree: generated/

    get_process_func
    process_raw
    process_save_radar
    process_point_measurement
    process_grid
"""

from copy import deepcopy
from warnings import warn
import numpy as np
from netCDF4 import num2date

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from .process_traj import process_trajectory, process_traj_atplane, \
    process_traj_antenna_pattern

from netCDF4 import num2date, date2num


def get_process_func(dataset_type, dsname):
    """
    maps the dataset type into its processing function and data set format

    Parameters
    ----------
    dataset_type : str
        data set type, i.e. 'RAW', 'SAN', etc.
    dsname : str
        Name of dataset

    Returns
    -------
    func_name : str or function
        pyrad function used to process the data set type
    dsformat : str
        data set format, i.e.: 'VOL', etc.

    """

    dsformat = 'VOL'
    if dataset_type == 'RAW':
        func_name = process_raw
    elif dataset_type == 'GRID':
        func_name = process_grid
        dsformat = 'GRID'
    elif dataset_type == 'CDF':
        func_name = 'process_cdf'
    elif dataset_type == 'NCVOL':
        func_name = process_save_radar
    elif dataset_type == 'PWR':
        func_name = 'process_signal_power'
    elif dataset_type == 'SNR':
        func_name = 'process_snr'
    elif dataset_type == 'RHOHV_CORRECTION':
        func_name = 'process_correct_noise_rhohv'
    elif dataset_type == 'BIAS_CORRECTION':
        func_name = 'process_correct_bias'
    elif dataset_type == 'L':
        func_name = 'process_l'
    elif dataset_type == 'CDR':
        func_name = 'process_cdr'
    elif dataset_type == 'SAN':
        func_name = 'process_echo_id'
    elif dataset_type == 'ECHO_FILTER':
        func_name = 'process_echo_filter'
    elif dataset_type == 'SNR_FILTER':
        func_name = 'process_filter_snr'
    elif dataset_type == 'VIS_FILTER':
        func_name = 'process_filter_visibility'
    elif dataset_type == 'OUTLIER_FILTER':
        func_name = 'process_outlier_filter'
    elif dataset_type == 'PHIDP0_CORRECTION':
        func_name = 'process_correct_phidp0'
    elif dataset_type == 'PHIDP_SMOOTH_1W':
        func_name = 'process_smooth_phidp_single_window'
    elif dataset_type == 'PHIDP_SMOOTH_2W':
        func_name = 'process_smooth_phidp_double_window'
    elif dataset_type == 'PHIDP_KDP_VULPIANI':
        func_name = 'process_phidp_kdp_Vulpiani'
    elif dataset_type == 'PHIDP_KDP_KALMAN':
        func_name = 'process_phidp_kdp_Kalman'
    elif dataset_type == 'PHIDP_KDP_MAESAKA':
        func_name = 'process_phidp_kdp_Maesaka'
    elif dataset_type == 'PHIDP_KDP_LP':
        func_name = 'process_phidp_kdp_lp'
    elif dataset_type == 'KDP_LEASTSQUARE_1W':
        func_name = 'process_kdp_leastsquare_single_window'
    elif dataset_type == 'KDP_LEASTSQUARE_2W':
        func_name = 'process_kdp_leastsquare_double_window'
    elif dataset_type == 'ATTENUATION':
        func_name = 'process_attenuation'
    elif dataset_type == 'RAINRATE':
        func_name = 'process_rainrate'
    elif dataset_type == 'WIND_VEL':
        func_name = 'process_wind_vel'
    elif dataset_type == 'WINDSHEAR':
        func_name = 'process_windshear'
    elif dataset_type == 'HYDROCLASS':
        func_name = 'process_hydroclass'
    elif dataset_type == 'PHIDP0_ESTIMATE':
        func_name = 'process_estimate_phidp0'
    elif dataset_type == 'RHOHV_RAIN':
        func_name = 'process_rhohv_rain'
    elif dataset_type == 'ZDR_PREC':
        func_name = 'process_zdr_precip'
    elif dataset_type == 'SELFCONSISTENCY_KDP_PHIDP':
        func_name = 'process_selfconsistency_kdp_phidp'
    elif dataset_type == 'SELFCONSISTENCY_BIAS':
        func_name = 'process_selfconsistency_bias'
    elif dataset_type == 'COSMO':
        func_name = 'process_cosmo'
    elif dataset_type == 'COSMO_LOOKUP':
        func_name = 'process_cosmo_lookup_table'
    elif dataset_type == 'COSMO_COORD':
        func_name = 'process_cosmo_coord'
        dsformat = 'COSMO_COORD'
    elif dataset_type == 'TIME_AVG':
        func_name = 'process_time_avg'
        dsformat = 'TIMEAVG'
    elif dataset_type == 'WEIGHTED_TIME_AVG':
        func_name = 'process_weighted_time_avg'
        dsformat = 'TIMEAVG'
    elif dataset_type == 'FLAG_TIME_AVG':
        func_name = 'process_time_avg_flag'
        dsformat = 'TIMEAVG'
    elif dataset_type == 'COLOCATED_GATES':
        func_name = 'process_colocated_gates'
        dsformat = 'COLOCATED_GATES'
    elif dataset_type == 'INTERCOMP':
        func_name = 'process_intercomp'
        dsformat = 'INTERCOMP'
    elif dataset_type == 'INTERCOMP_TIME_AVG':
        func_name = 'process_intercomp_time_avg'
        dsformat = 'INTERCOMP'
    elif dataset_type == 'MONITORING':
        func_name = 'process_monitoring'
        dsformat = 'MONITORING'
    elif dataset_type == 'SUN_HITS':
        func_name = 'process_sun_hits'
        dsformat = 'SUN_HITS'
    elif dataset_type == 'POINT_MEASUREMENT':
        func_name = process_point_measurement
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ':
        func_name = process_trajectory
        dsformat = 'TRAJ_ONLY'
    elif dataset_type == 'TRAJ_ATPLANE':
        func_name = process_traj_atplane
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_ANTENNA_PATTERN':
        func_name = process_traj_antenna_pattern
        dsformat = 'TIMESERIES'
    else:
        raise ValueError("ERROR: Unknown dataset type '%s' of dataset '%s'"
                         % (dataset_type, dsname))

    return func_name, dsformat


def process_raw(procstatus, dscfg, radar_list=None):
    """
    dummy function that returns the initial input data set

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
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

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        break
    ind_rad = int(radarnr[5:8])-1
    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar')
        return None, None
    new_dataset = deepcopy(radar_list[ind_rad])

    return new_dataset, ind_rad


def process_save_radar(procstatus, dscfg, radar_list=None):
    """
    dummy function that allows to save the entire radar object

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
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

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        break
    ind_rad = int(radarnr[5:8])-1
    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar')
        return None, None
    new_dataset = deepcopy(radar_list[ind_rad])

    return new_dataset, ind_rad


def process_point_measurement(procstatus, dscfg, radar_list=None):
    """
    Obtains the radar data at a point measurement

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        latlon : boolean. Dataset keyword
            if True position is obtained from latitude, longitude information,
            otherwise position is obtained from antenna coordinates
            (range, azimuth, elevation).
        truealt : boolean. Dataset keyword
            if True the user input altitude is used to determine the point of
            interest.
            if False use the altitude at a given radar elevation ele over the
            point of interest.
        lon : float. Dataset keyword
            the longitude [deg]. Use when latlon is True.
        lat : float. Dataset keyword
            the latitude [deg]. Use when latlon is True.
        alt : float. Dataset keyword
            altitude [m MSL]. Use when latlon is True.
        ele : float. Dataset keyword
            radar elevation [deg]. Use when latlon is False or when latlon is
            True and truealt is False
        azi : float. Dataset keyword
            radar azimuth [deg]. Use when latlon is False
        rng : float. Dataset keyword
            range from radar [m]. Use when latlon is False
        AziTol : float. Dataset keyword
            azimuthal tolerance to determine which radar azimuth to use [deg]
        EleTol : float. Dataset keyword
            elevation tolerance to determine which radar elevation to use [deg]
        RngTol : float. Dataset keyword
            range tolerance to determine which radar bin to use [m]
    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and metadata of the point of interest
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        break
    field_name = get_fieldname_pyart(datatype)
    ind_rad = int(radarnr[5:8])-1

    if procstatus == 2:
        if dscfg['initialized'] == 0:
            return None, None

        # prepare for exit
        new_dataset = {
            'time': dscfg['global_data']['time'],
            'datatype': 'RR',
            'point_coordinates_WGS84_lon_lat_alt': (
                dscfg['global_data']['point_coordinates_WGS84_lon_lat_alt']),
            'antenna_coordinates_az_el_r': (
                dscfg['global_data']['antenna_coordinates_az_el_r']),
            'final': True}

        return new_dataset, ind_rad

    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if field_name not in radar.fields:
        warn('Unable to extract point measurement information. ' +
             'Field not available')
        return None, None

    projparams = dict()
    projparams.update({'proj': 'pyart_aeqd'})
    projparams.update({'lon_0': radar.longitude['data']})
    projparams.update({'lat_0': radar.latitude['data']})

    if dscfg['latlon']:
        lon = dscfg['lon']
        lat = dscfg['lat']
        alt = dscfg['alt']
        x, y = pyart.core.geographic_to_cartesian(lon, lat, projparams)

        if not dscfg['truealt']:
            ke = 4./3.  # constant for effective radius
            a = 6378100.  # earth radius
            re = a * ke  # effective radius

            elrad = dscfg['ele'] * np.pi / 180.
            r_ground = np.sqrt(x ** 2. + y ** 2.)
            r = r_ground / np.cos(elrad)
            alt_radar = radar.altitude['data']+np.sqrt(
                r ** 2. + re ** 2. + 2. * r * re * np.sin(elrad)) - re
            alt_radar = alt_radar[0]
        else:
            alt_radar = dscfg['alt']

        r, az, el = pyart.core.cartesian_to_antenna(
            x, y, alt_radar-radar.altitude['data'])
        r = r[0]
        az = az[0]
        el = el[0]
    else:
        r = dscfg['rng']
        az = dscfg['azi']
        el = dscfg['ele']

        x, y, alt = pyart.core.antenna_to_cartesian(r, az, el)
        lon, lat = pyart.core.cartesian_to_geographic(x, y, projparams)

    d_az = np.min(np.abs(radar.azimuth['data'] - az))
    if d_az > dscfg['AziTol']:
        warn(' No radar bin found for point (az, el, r):(' +
             str(az)+', '+str(el)+', '+str(r) +
             '). Minimum distance to radar azimuth '+str(d_az) +
             ' larger than tolerance')
        return None, None

    d_el = np.min(np.abs(radar.elevation['data'] - el))
    if d_el > dscfg['EleTol']:
        warn(' No radar bin found for point (az, el, r):(' +
             str(az)+', '+str(el)+', '+str(r) +
             '). Minimum distance to radar elevation '+str(d_el) +
             ' larger than tolerance')
        return None, None

    d_r = np.min(np.abs(radar.range['data'] - r))
    if d_r > dscfg['RngTol']:
        warn(' No radar bin found for point (az, el, r):(' +
             str(az)+', '+str(el)+', '+str(r) +
             '). Minimum distance to radar range bin '+str(d_r) +
             ' larger than tolerance')
        return None, None

    ind_ray = np.argmin(np.abs(radar.azimuth['data'] - az) +
                        np.abs(radar.elevation['data'] - el))
    ind_r = np.argmin(np.abs(radar.range['data'] - r))

    val = radar.fields[field_name]['data'].data[ind_ray, ind_r]
    time = num2date(radar.time['data'][ind_ray], radar.time['units'],
                    radar.time['calendar'])

    # initialize dataset
    if dscfg['initialized'] == 0:
        poi = {
            'point_coordinates_WGS84_lon_lat_alt': [lon, lat, alt],
            'antenna_coordinates_az_el_r': [az, el, r],
            'time': time}
        dscfg['global_data'] = poi
        dscfg['initialized'] = 1

    # prepare for exit
    new_dataset = dict()
    new_dataset.update({'value': val})
    new_dataset.update({'datatype': datatype})
    new_dataset.update({'time': time})
    new_dataset.update(
        {'point_coordinates_WGS84_lon_lat_alt': [lon, lat, alt]})
    new_dataset.update({'antenna_coordinates_az_el_r': [az, el, r]})
    new_dataset.update(
        {'used_antenna_coordinates_az_el_r': [radar.azimuth['data'][ind_ray],
         radar.elevation['data'][ind_ray],
         radar.range['data'][ind_r]]})
    new_dataset.update({'final': False})

    return new_dataset, ind_rad


def process_grid(procstatus, dscfg, radar_list=None):
    """
    Puts the radar data in a regular grid

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        gridconfig : dictionary. Dataset keyword
            Dictionary containing some or all of this keywords:
            xmin, xmax, ymin, ymax, zmin, zmax : floats
                minimum and maximum horizontal distance from grid origin [km]
                and minimum and maximum vertical distance from grid origin [m]
                Defaults -40, 40, -40, 40, 0., 10000.
            hres, vres : floats
                horizontal and vertical grid resolution [m]
                Defaults 1000., 500.
            latorig, lonorig, altorig : floats
                latitude and longitude of grid origin [deg] and altitude of
                grid origin [m MSL]
                Defaults the latitude, longitude and altitude of the radar
        wfunc : str
            the weighting function used to combine the radar gates close to a
            grid point. Possible values BARNES, CRESSMAN, NEAREST_NEIGHBOUR
            Default NEAREST_NEIGHBOUR
        roif_func : str
            the function used to compute the region of interest.
            Possible values: dist_beam, constant
        roi : float
             the (minimum) radius of the region of interest in m. Default half
             the largest resolution

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the gridded data
    ind_rad : int
        radar index

    """
    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
            datatypedescr)
        break
    field_name = get_fieldname_pyart(datatype)
    ind_rad = int(radarnr[5:8])-1

    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar')
        return None, None

    radar = radar_list[ind_rad]

    # default parameters
    xmin = -40.
    xmax = 40.
    ymin = -40.
    ymax = 40.
    zmin = 0.
    zmax = 10000.
    hres = 1000.
    vres = 500.
    lat = float(radar.latitude['data'])
    lon = float(radar.longitude['data'])
    alt = float(radar.altitude['data'])

    if 'gridConfig' in dscfg:
        if 'xmin' in dscfg['gridConfig']:
            xmin = dscfg['gridConfig']['xmin']
        if 'xmax' in dscfg['gridConfig']:
            xmax = dscfg['gridConfig']['xmax']
        if 'ymin' in dscfg['gridConfig']:
            ymin = dscfg['gridConfig']['ymin']
        if 'ymax' in dscfg['gridConfig']:
            ymax = dscfg['gridConfig']['ymax']
        if 'zmin' in dscfg['gridConfig']:
            zmin = dscfg['gridConfig']['zmin']
        if 'zmax' in dscfg['gridConfig']:
            zmax = dscfg['gridConfig']['zmax']
        if 'hres' in dscfg['gridConfig']:
            hres = dscfg['gridConfig']['hres']
        if 'vres' in dscfg['gridConfig']:
            vres = dscfg['gridConfig']['vres']
        if 'latorig' in dscfg['gridConfig']:
            lat = dscfg['gridConfig']['latorig']
        if 'lonorig' in dscfg['gridConfig']:
            lon = dscfg['gridConfig']['lonorig']
        if 'altorig' in dscfg['gridConfig']:
            alt = dscfg['gridConfig']['altorig']

    wfunc = 'NEAREST_NEIGHBOUR'
    if 'wfunc' in dscfg:
        wfunc = dscfg['wfunc']

    roi_func = 'dist_beam'
    if 'roi_func' in dscfg:
        roi_func = dscfg['roi_func']

    # number of grid points in cappi
    nz = int((zmax-zmin)/vres)+1
    ny = int((ymax-ymin)*1000./hres)+1
    nx = int((xmax-xmin)*1000./hres)+1

    min_radius = np.max([vres, hres])/2.
    if 'roi' in dscfg:
        min_radius = dscfg['roi']

    # parameters to determine the gates to use for each grid point
    beamwidth = 1.
    beam_spacing = 1.
    if 'radar_beam_width_h' in radar.instrument_parameters:
        beamwidth = radar.instrument_parameters[
            'radar_beam_width_h']['data'][0]

    if radar.ray_angle_res is not None:
        beam_spacing = radar.ray_angle_res['data'][0]

    # cartesian mapping
    grid = pyart.map.grid_from_radars(
        (radar,), gridding_algo='map_to_grid',
        weighting_function=wfunc,
        roi_func=roi_func, h_factor=1.0, nb=beamwidth, bsp=beam_spacing,
        min_radius=min_radius, constant_roi=min_radius,
        grid_shape=(nz, ny, nx),
        grid_limits=((zmin, zmax), (ymin*1000., ymax*1000.),
                     (xmin*1000., xmax*1000.)),
        grid_origin=(lat, lon), grid_origin_alt=alt,
        fields=[field_name])

    return grid, ind_rad
