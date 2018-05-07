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
    process_roi
    process_grid
    process_qvp
    process_time_height
"""

from copy import deepcopy
from warnings import warn
import numpy as np
from netCDF4 import num2date

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_sensor import read_trt_traj_data
from .process_traj import process_trajectory, process_traj_atplane
from .process_traj import process_traj_antenna_pattern, process_traj_lightning
from .process_traj import process_traj_trt
from ..util.radar_utils import find_rng_index, belongs_roi_indices


def get_process_func(dataset_type, dsname):
    """
    Maps the dataset type into its processing function and data set format
    associated.

    Parameters
    ----------
    dataset_type : str
        data set type, i.e. 'RAW', 'SAN', etc.
    dsname : str
        Name of dataset

    Returns
    -------
    func_name : str or processing function
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
    elif dataset_type == 'QVP':
        func_name = process_qvp
        dsformat = 'QVP'
    elif dataset_type == 'TIME_HEIGHT':
        func_name = process_time_height
        dsformat = 'QVP'
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
    elif dataset_type == 'CLT_TO_SAN':
        func_name = 'process_clt_to_echo_id'
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
    elif dataset_type == 'ML_DETECTION':
        func_name = 'process_melting_layer'
    elif dataset_type == 'PHIDP0_ESTIMATE':
        func_name = 'process_estimate_phidp0'
    elif dataset_type == 'RHOHV_RAIN':
        func_name = 'process_rhohv_rain'
    elif dataset_type == 'ZDR_PREC':
        func_name = 'process_zdr_precip'
    elif dataset_type == 'ZDR_SNOW':
        func_name = 'process_zdr_snow'
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
    elif dataset_type == 'HZT_COORD':
        func_name = 'process_hzt_coord'
        dsformat = 'COSMO_COORD'
    elif dataset_type == 'HZT':
        func_name = 'process_hzt'
    elif dataset_type == 'HZT_LOOKUP':
        func_name = 'process_hzt_lookup_table'
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
    elif dataset_type == 'GC_MONITORING':
        func_name = 'process_gc_monitoring'
        dsformat = 'MONITORING'
    elif dataset_type == 'OCCURRENCE':
        func_name = 'process_occurrence'
        dsformat = 'OCCURRENCE'
    elif dataset_type == 'OCCURRENCE_PERIOD':
        func_name = 'process_occurrence_period'
        dsformat = 'OCCURRENCE'
    elif dataset_type == 'SUN_HITS':
        func_name = 'process_sun_hits'
        dsformat = 'SUN_HITS'
    elif dataset_type == 'POINT_MEASUREMENT':
        func_name = process_point_measurement
        dsformat = 'TIMESERIES'
    elif dataset_type == 'ROI':
        func_name = process_roi
    elif dataset_type == 'TRAJ':
        func_name = process_trajectory
        dsformat = 'TRAJ_ONLY'
    elif dataset_type == 'TRAJ_ATPLANE':
        func_name = process_traj_atplane
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_ANTENNA_PATTERN':
        func_name = process_traj_antenna_pattern
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_LIGHTNING':
        func_name = process_traj_lightning
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_TRT':
        func_name = process_traj_trt
    else:
        raise ValueError("ERROR: Unknown dataset type '%s' of dataset '%s'"
                         % (dataset_type, dsname))

    return func_name, dsformat


def process_raw(procstatus, dscfg, radar_list=None):
    """
    Dummy function that returns the initial input data set

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
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    new_dataset = deepcopy(radar_list[ind_rad])

    return new_dataset, ind_rad


def process_save_radar(procstatus, dscfg, radar_list=None):
    """
    Dummy function that allows to save the entire radar object

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
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    new_dataset = deepcopy(radar_list[ind_rad])

    return new_dataset, ind_rad


def process_point_measurement(procstatus, dscfg, radar_list=None):
    """
    Obtains the radar data at a point location.

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
        dictionary containing the data and metadata at the point of interest
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

    if (radar_list is None) or (radar_list[ind_rad] is None):
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
        {'used_antenna_coordinates_az_el_r': [
            radar.azimuth['data'][ind_ray],
            radar.elevation['data'][ind_ray],
            radar.range['data'][ind_r]]})
    new_dataset.update({'final': False})

    return new_dataset, ind_rad


def process_roi(procstatus, dscfg, radar_list=None):
    """
    Obtains the radar data at a region of interest.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement

    radar_list : list of Radar objects
          Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the data and metadata at the point of interest
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

    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    if field_name not in radar.fields:
        warn('Unable to extract ROI information. ' +
             'Field not available')
        return None, None

    if 'trtfile' in dscfg:
        (traj_ID, yyyymmddHHMM, lon, lat, ell_L, ell_S, ell_or, area,
         vel_x, vel_y, det, RANKr, CG_n, CG_p, CG, CG_percent_p, ET45,
         ET45m, ET15, ET15m, VIL, maxH, maxHm, POH, RANK, Dvel_x,
         Dvel_y, cell_contour) = read_trt_traj_data(dscfg['trtfile'])

        time_tol = dscfg.get('TimeTol', 100.)
        dt = np.empty(yyyymmddHHMM.size, dtype=float)
        for i, time_traj in enumerate(yyyymmddHHMM):
            dt[i] = np.abs((dscfg['timeinfo'] - time_traj).total_seconds())
        if dt.min() > time_tol:
            warn('No TRT data for radar volume time')
            return None, None

        ind = np.argmin(dt)
        lon_roi = cell_contour[ind]['lon']
        lat_roi = cell_contour[ind]['lat']
    else:
        lon_roi = dscfg.get('lon_roi', None)
        lat_roi = dscfg.get('lat_roi', None)

        if lon_roi is None or lat_roi is None:
            warn('Undefined ROI')
            return None, None

    alt_min = dscfg.get('alt_min', None)
    alt_max = dscfg.get('alt_max', None)

    roi_dict = {
        'lon': lon_roi,
        'lat': lat_roi,
        'alt_min': alt_min,
        'alt_max': alt_max}

    # extract the data within the ROI boundaries
    inds_ray, inds_rng = np.indices(np.shape(radar.gate_longitude['data']))

    mask = np.logical_and(
        np.logical_and(
            radar.gate_latitude['data'] >= roi_dict['lat'].min(),
            radar.gate_latitude['data'] <= roi_dict['lat'].max()),
        np.logical_and(
            radar.gate_longitude['data'] >= roi_dict['lon'].min(),
            radar.gate_longitude['data'] <= roi_dict['lon'].max()))

    if alt_min is not None:
        mask[radar.gate_altitude['data'] < alt_min] = 0
    if alt_max is not None:
        mask[radar.gate_altitude['data'] > alt_max] = 0

    if np.all(mask == 0):
        warn('No values within ROI')
        return None, None

    inds_ray = inds_ray[mask]
    inds_rng = inds_rng[mask]

    # extract the data inside the ROI
    lat = radar.gate_latitude['data'][mask]
    lon = radar.gate_longitude['data'][mask]
    inds, is_roi = belongs_roi_indices(lat, lon, roi_dict)

    if is_roi == 'None':
        warn('No values within ROI')
        return None, None

    inds_ray = inds_ray[inds]
    inds_rng = inds_rng[inds]

    lat = lat[inds].T
    lon = lon[inds].T
    alt = radar.gate_altitude['data'][inds_ray, inds_rng].T

    # prepare new radar object output
    radar_roi = deepcopy(radar)

    radar_roi.range['data'] = radar.range['data'][inds_rng]
    radar_roi.ngates = inds_rng.size
    radar_roi.time['data'] = np.asarray([radar_roi.time['data'][0]])
    radar_roi.scan_type = 'roi'
    radar_roi.sweep_mode['data'] = np.array(['roi'])
    radar_roi.sweep_start_ray_index['data'] = np.array([0], dtype='int32')
    radar_roi.fixed_angle['data'] = np.array([], dtype='float64')
    radar_roi.sweep_number['data'] = np.array([0], dtype='int32')
    radar_roi.nsweeps = 1

    if radar.rays_are_indexed is not None:
        radar_roi.rays_are_indexed['data'] = np.array(
            [radar.rays_are_indexed['data'][0]])
    if radar.ray_angle_res is not None:
        radar_roi.ray_angle_res['data'] = np.array(
            [radar.ray_angle_res['data'][0]])

    radar_roi.sweep_end_ray_index['data'] = np.array([1], dtype='int32')
    radar_roi.rays_per_sweep = np.array([1], dtype='int32')
    radar_roi.azimuth['data'] = np.array([], dtype='float64')
    radar_roi.elevation['data'] = np.array([], dtype='float64')
    radar_roi.nrays = 1

    radar_roi.gate_longitude['data'] = lon
    radar_roi.gate_latitude['data'] = lat
    radar_roi.gate_altitude['data'] = alt

    radar_roi.gate_x['data'] = radar.gate_x['data'][inds_ray, inds_rng].T
    radar_roi.gate_y['data'] = radar.gate_y['data'][inds_ray, inds_rng].T
    radar_roi.gate_z['data'] = radar.gate_z['data'][inds_ray, inds_rng].T

    radar_roi.fields = dict()
    field_dict = deepcopy(radar.fields[field_name])
    field_dict['data'] = radar.fields[field_name]['data'][inds_ray, inds_rng].T
    radar_roi.add_field(field_name, field_dict)

    return radar_roi, ind_rad


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

    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None

    radar = radar_list[ind_rad]

    if field_name not in radar.fields:
        warn('Field name '+field_name+' not available in radar object')
        return None, None

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

    wfunc = dscfg.get('wfunc', 'NEAREST_NEIGHBOUR')
    roi_func = dscfg.get('roi_func', 'dist_beam')

    # number of grid points in cappi
    nz = int((zmax-zmin)/vres)+1
    ny = int((ymax-ymin)*1000./hres)+1
    nx = int((xmax-xmin)*1000./hres)+1

    min_radius = dscfg.get('roi', np.max([vres, hres])/2.)
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


def process_qvp(procstatus, dscfg, radar_list=None):
    """
    Computes quasi vertical profiles, by averaging over height levels
    PPI or RHI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement

        anglenr : int
            The sweep number to use. It assumes the radar volume consists on
            PPI scans

        hmax : float
            The maximum height to plot [m]. Default 10000.

        hres : float
            The height resolution [m]. Default 50

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            break
        field_name = get_fieldname_pyart(datatype)
        ind_rad = int(radarnr[5:8])-1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn('ERROR: No valid radar')
            return None, None

        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn('Field name '+field_name+' not available in radar object')
            return None, None

        # default parameters
        anglenr = dscfg.get('anglenr', 0)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 50.)

        radar_aux = deepcopy(radar)
        radar_aux = radar_aux.extract_sweeps([anglenr])

        # initialize dataset
        if dscfg['initialized'] == 0:
            qvp_aux = deepcopy(radar_aux)
            # prepare space for field
            qvp_aux.fields = dict()
            qvp_aux.add_field(
                field_name, deepcopy(radar_aux.fields[field_name]))
            qvp_aux.fields[field_name]['data'] = np.array([], dtype='float64')

            # fixed radar objects parameters
            qvp_aux.range['data'] = np.arange(hmax/hres)*hres+hres/2.
            qvp_aux.ngates = len(qvp_aux.range['data'])

            qvp_aux.time['units'] = pyart.io.make_time_unit_str(
                dscfg['timeinfo'])
            qvp_aux.time['data'] = np.array([], dtype='float64')
            qvp_aux.scan_type = 'qvp'
            qvp_aux.sweep_mode['data'] = np.array(['qvp'])
            qvp_aux.sweep_start_ray_index['data'] = np.array(
                [0], dtype='int32')

            # ray dependent radar objects parameters
            qvp_aux.sweep_end_ray_index['data'] = np.array([-1], dtype='int32')
            qvp_aux.rays_per_sweep = np.array([0], dtype='int32')
            qvp_aux.azimuth['data'] = np.array([], dtype='float64')
            qvp_aux.elevation['data'] = np.array([], dtype='float64')
            qvp_aux.nrays = 0

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_obj': qvp_aux})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1

        # modify metadata
        qvp = dscfg['global_data']['radar_obj']

        start_time = num2date(0, qvp.time['units'], qvp.time['calendar'])
        qvp.time['data'] = np.append(
            qvp.time['data'], (dscfg['timeinfo'] - start_time).total_seconds())
        qvp.sweep_end_ray_index['data'][0] += 1
        qvp.rays_per_sweep[0] += 1
        qvp.nrays += 1

        qvp.azimuth['data'] = np.ones((qvp.nrays, ), dtype='float64')*0.
        qvp.elevation['data'] = (
            np.ones((qvp.nrays, ), dtype='float64') *
            qvp.fixed_angle['data'][0])

        qvp.gate_longitude['data'] = (
            np.ones((qvp.nrays, qvp.ngates), dtype='float64') *
            qvp.longitude['data'][0])
        qvp.gate_latitude['data'] = (
            np.ones((qvp.nrays, qvp.ngates), dtype='float64') *
            qvp.latitude['data'][0])
        qvp.gate_altitude['data'] = np.broadcast_to(
            qvp.range['data'], (qvp.nrays, qvp.ngates))

        # compute QVP data
        values = np.ma.mean(radar_aux.fields[field_name]['data'], axis=0)
        # altitude corresponding to qvp grid:
        qvp_data = np.ma.zeros(qvp.ngates)
        qvp_data[:] = np.ma.masked
        for ind_r, h in enumerate(qvp.range['data']):
            ind_h = find_rng_index(
                radar_aux.gate_altitude['data'][0, :], h, rng_tol=hres/2.)
            if ind_h is None:
                continue
            qvp_data[ind_r] = values[ind_h]

        if np.size(qvp.fields[field_name]['data']) == 0:
            qvp.fields[field_name]['data'] = qvp_data.reshape(1, qvp.ngates)
        else:
            qvp.fields[field_name]['data'] = np.ma.concatenate(
                (qvp.fields[field_name]['data'],
                 qvp_data.reshape(1, qvp.ngates)))

        dscfg['global_data']['radar_obj'] = qvp

        new_dataset = dict()
        new_dataset.update({'radar_obj': qvp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            break

        ind_rad = int(radarnr[5:8])-1

        qvp = dscfg['global_data']['radar_obj']

        new_dataset = dict()
        new_dataset.update({'radar_obj': qvp})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad


def process_time_height(procstatus, dscfg, radar_list=None):
    """
    Produces time height radar objects at a point of interest defined by
    latitude and longitude. A time-height contains the evolution
    of the vertical structure of radar measurements above the location
    of interest.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        lat, lon : float
            latitude and longitude of the point of interest [deg]
        latlon_tol : float
            tolerance in latitude and longitude in deg. Default 0.0005
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 50

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            break
        field_name = get_fieldname_pyart(datatype)
        ind_rad = int(radarnr[5:8])-1

        if (radar_list is None) or (radar_list[ind_rad] is None):
            warn('ERROR: No valid radar')
            return None, None

        radar = radar_list[ind_rad]

        if field_name not in radar.fields:
            warn('Field name '+field_name+' not available in radar object')
            return None, None

        # default parameters
        lon = dscfg['lon']
        lat = dscfg['lat']
        latlon_tol = dscfg.get('latlon_tol', 0.0005)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 50.)

        radar_aux = deepcopy(radar)

        # initialize dataset
        if dscfg['initialized'] == 0:
            th_aux = deepcopy(radar_aux)
            # prepare space for field
            th_aux.fields = dict()
            th_aux.add_field(
                field_name, deepcopy(radar_aux.fields[field_name]))
            th_aux.fields[field_name]['data'] = np.array([], dtype='float64')

            # fixed radar objects parameters
            th_aux.range['data'] = np.arange(hmax/hres)*hres+hres/2.
            th_aux.ngates = len(th_aux.range['data'])

            th_aux.time['units'] = pyart.io.make_time_unit_str(
                dscfg['timeinfo'])
            th_aux.time['data'] = np.array([], dtype='float64')
            th_aux.scan_type = 'time_height'
            th_aux.sweep_mode['data'] = np.array(['time_height'])
            th_aux.sweep_start_ray_index['data'] = np.array(
                [0], dtype='int32')
            th_aux.fixed_angle['data'] = np.array([90.], dtype='float64')
            th_aux.sweep_number['data'] = np.array([0], dtype='int32')
            th_aux.nsweeps = 1

            if radar_aux.rays_are_indexed is not None:
                th_aux.rays_are_indexed['data'] = np.array(
                    [radar_aux.rays_are_indexed['data'][0]])

            if radar_aux.ray_angle_res is not None:
                th_aux.ray_angle_res['data'] = np.array(
                    [radar_aux.ray_angle_res['data'][0]])

            # ray dependent radar objects parameters
            th_aux.sweep_end_ray_index['data'] = np.array([-1], dtype='int32')
            th_aux.rays_per_sweep = np.array([0], dtype='int32')
            th_aux.azimuth['data'] = np.array([], dtype='float64')
            th_aux.elevation['data'] = np.array([], dtype='float64')
            th_aux.nrays = 0

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_obj': th_aux})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1

        # modify metadata
        th = dscfg['global_data']['radar_obj']

        start_time = num2date(0, th.time['units'], th.time['calendar'])
        th.time['data'] = np.append(
            th.time['data'], (dscfg['timeinfo'] - start_time).total_seconds())
        th.sweep_end_ray_index['data'][0] += 1
        th.rays_per_sweep[0] += 1
        th.nrays += 1

        th.azimuth['data'] = np.ones((th.nrays, ), dtype='float64')*0.
        th.elevation['data'] = np.ones((th.nrays, ), dtype='float64')*90.

        th.gate_longitude['data'] = (
            np.ones((th.nrays, th.ngates), dtype='float64')*lon)
        th.gate_latitude['data'] = (
            np.ones((th.nrays, th.ngates), dtype='float64')*lat)
        th.gate_altitude['data'] = np.broadcast_to(
            th.range['data'], (th.nrays, th.ngates))

        # find data
        th_data = np.ma.zeros(th.ngates)
        th_data[:] = np.ma.masked

        # find gates close to lat lon point
        inds = np.where(np.logical_and(
            np.logical_and(
                radar_aux.gate_latitude['data'][:, :] < lat+latlon_tol,
                radar_aux.gate_latitude['data'][:, :] > lat-latlon_tol),
            np.logical_and(
                radar_aux.gate_longitude['data'][:, :] < lon+latlon_tol,
                radar_aux.gate_longitude['data'][:, :] > lon-latlon_tol)))

        # find closest altitude
        if inds[0].size > 0:
            values = radar_aux.fields[field_name]['data'][inds].flatten()
            altitudes = radar_aux.gate_altitude['data'][inds].flatten()
            for ind_r, h in enumerate(th.range['data']):
                ind_h = find_rng_index(altitudes, h, rng_tol=hres/2.)
                if ind_h is None:
                    continue
                th_data[ind_r] = values[ind_h]
        else:
            warn('No data found at point lat '+str(lat)+' +- ' +
                 str(latlon_tol)+' lon '+str(lon)+' +- ' +
                 str(latlon_tol)+' deg')

        if np.size(th.fields[field_name]['data']) == 0:
            th.fields[field_name]['data'] = th_data.reshape(1, th.ngates)
        else:
            th.fields[field_name]['data'] = np.ma.concatenate(
                (th.fields[field_name]['data'],
                 th_data.reshape(1, th.ngates)))

        dscfg['global_data']['radar_obj'] = th

        new_dataset = dict()
        new_dataset.update({'radar_obj': th})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        for datatypedescr in dscfg['datatype']:
            radarnr, datagroup, datatype, dataset, product = (
                get_datatype_fields(datatypedescr))
            break

        ind_rad = int(radarnr[5:8])-1

        th = dscfg['global_data']['radar_obj']

        new_dataset = dict()
        new_dataset.update({'radar_obj': th})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad
