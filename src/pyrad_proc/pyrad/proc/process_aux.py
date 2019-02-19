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
    process_fixed_rng
    process_roi
    process_grid
    process_azimuthal_average

"""

from copy import deepcopy
from warnings import warn
import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.read_data_sensor import read_trt_traj_data
from ..util.radar_utils import belongs_roi_indices, get_target_elevations
from ..util.radar_utils import find_neighbour_gates, compute_directional_stats
from ..util.radar_utils import get_fixed_rng_data


def get_process_func(dataset_type, dsname):
    """
    Maps the dataset type into its processing function and data set format
    associated.

    Parameters
    ----------
    dataset_type : str
        The following is a list of data set types ordered by type of output
        dataset with the function they call. For details of what they do check
        the function documentation:
            'VOL' format output:
                'ATTENUATION': process_attenuation
                'AZI_AVG': process_azimuthal_average
                'BIAS_CORRECTION': process_correct_bias
                'BIRDS_ID': process_birds_id
                'BIRD_DENSITY': process_bird_density
                'CDF': process_cdf
                'CDR': process_cdr
                'CLT_TO_SAN': process_clt_to_echo_id
                'COSMO': process_cosmo
                'COSMO_LOOKUP': process_cosmo_lookup_table
                'DEALIAS_FOURDD': process_dealias_fourdd
                'DEALIAS_REGION': process_dealias_region_based
                'DEALIAS_UNWRAP': process_dealias_unwrap_phase
                'ECHO_FILTER': process_echo_filter
                'FIXED_RNG': process_fixed_rng
                'HYDROCLASS': process_hydroclass
                'HZT': process_hzt
                'HZT_LOOKUP': process_hzt_lookup_table
                'KDP_LEASTSQUARE_1W': process_kdp_leastsquare_single_window
                'KDP_LEASTSQUARE_2W': process_kdp_leastsquare_double_window
                'L': process_l
                'NCVOL': process_save_radar
                'OUTLIER_FILTER': process_outlier_filter
                'PHIDP0_CORRECTION': process_correct_phidp0
                'PHIDP0_ESTIMATE': process_estimate_phidp0
                'PHIDP_KDP_KALMAN': process_phidp_kdp_Kalman
                'PHIDP_KDP_LP': process_phidp_kdp_lp
                'PHIDP_KDP_VULPIANI': process_phidp_kdp_Vulpiani
                'PHIDP_SMOOTH_1W': process_smooth_phidp_single_window
                'PHIDP_SMOOTH_2W': process_smooth_phidp_double_window
                'PWR': process_signal_power
                'RAINRATE': process_rainrate
                'RAW': process_raw
                'RCS': process_rcs
                'RCS_PR': process_rcs_pr
                'RHOHV_CORRECTION': process_correct_noise_rhohv
                'RHOHV_RAIN': process_rhohv_rain
                'ROI': process_roi
                'SAN': process_echo_id
                'SELFCONSISTENCY_BIAS': process_selfconsistency_bias
                'SELFCONSISTENCY_KDP_PHIDP': process_selfconsistency_kdp_phidp
                'SNR': process_snr
                'SNR_FILTER': process_filter_snr
                'TRAJ_TRT' : process_traj_trt
                'VAD': process_vad
                'VEL_FILTER': process_filter_vel_diff
                'VIS_FILTER': process_filter_visibility
                'VOL_REFL': process_vol_refl
                'WIND_VEL': process_wind_vel
                'WINDSHEAR': process_windshear
                'ZDR_PREC': process_zdr_precip
                'ZDR_SNOW': process_zdr_snow
            'COLOCATED_GATES' format output:
                'COLOCATED_GATES': process_colocated_gates
            'COSMO_COORD' format output:
                'COSMO_COORD': process_cosmo_coord
                'HZT_COORD': process_hzt_coord
            'GRID' format output:
                'GRID': process_grid
            'INTERCOMP' format output:
                'INTERCOMP': process_intercomp
                'INTERCOMP_TIME_AVG': process_intercomp_time_avg
            'ML' format output:
                'ML_DETECTION': process_melting_layer
            'MONITORING' format output:
                'GC_MONITORING': process_gc_monitoring
                'MONITORING': process_monitoring
            'OCCURRENCE' format output:
                'OCCURRENCE': process_occurrence
                'OCCURRENCE_PERIOD': process_occurrence_period
                'TIMEAVG_STD': process_time_avg_std
            'QVP' format output:
                'EVP': process_evp
                'QVP': process_qvp
                'rQVP': process_rqvp
                'SVP': process_svp
                'TIME_HEIGHT': process_time_height
            'SPARSE_GRID' format output:
                'ZDR_COLUMN': process_zdr_column
            'SUN_HITS' format output:
                'SUN_HITS': process_sun_hits
            'TIMEAVG' format output:
                'FLAG_TIME_AVG': process_time_avg_flag
                'TIME_AVG': process_time_avg
                'WEIGHTED_TIME_AVG': process_weighted_time_avg
            'TIMESERIES' format output:
                'POINT_MEASUREMENT': 'process_point_measurement'
                'TRAJ_ANTENNA_PATTERN': process_traj_antenna_pattern
                'TRAJ_ATPLANE': process_traj_atplane
                'TRAJ_LIGHTNING': process_traj_lightning
            'TRAJ_ONLY' format output:
                'TRAJ': process_trajectory
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
    elif dataset_type == 'AZI_AVG':
        func_name = process_azimuthal_average
    elif dataset_type == 'GRID':
        func_name = process_grid
        dsformat = 'GRID'
    elif dataset_type == 'QVP':
        func_name = 'process_qvp'
        dsformat = 'QVP'
    elif dataset_type == 'rQVP':
        func_name = 'process_rqvp'
        dsformat = 'QVP'
    elif dataset_type == 'SVP':
        func_name = 'process_svp'
        dsformat = 'QVP'
    elif dataset_type == 'EVP':
        func_name = 'process_evp'
        dsformat = 'QVP'
    elif dataset_type == 'TIME_HEIGHT':
        func_name = 'process_time_height'
        dsformat = 'QVP'
    elif dataset_type == 'CDF':
        func_name = 'process_cdf'
    elif dataset_type == 'NCVOL':
        func_name = process_save_radar
    elif dataset_type == 'PWR':
        func_name = 'process_signal_power'
    elif dataset_type == 'RCS_PR':
        func_name = 'process_rcs_pr'
    elif dataset_type == 'RCS':
        func_name = 'process_rcs'
    elif dataset_type == 'SNR':
        func_name = 'process_snr'
    elif dataset_type == 'VOL_REFL':
        func_name = 'process_vol_refl'
    elif dataset_type == 'BIRD_DENSITY':
        func_name = 'process_bird_density'
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
    elif dataset_type == 'BIRDS_ID':
        func_name = 'process_birds_id'
    elif dataset_type == 'CLT_TO_SAN':
        func_name = 'process_clt_to_echo_id'
    elif dataset_type == 'ECHO_FILTER':
        func_name = 'process_echo_filter'
    elif dataset_type == 'ZDR_COLUMN':
        func_name = 'process_zdr_column'
        dsformat = 'SPARSE_GRID'
    elif dataset_type == 'SNR_FILTER':
        func_name = 'process_filter_snr'
    elif dataset_type == 'VEL_FILTER':
        func_name = 'process_filter_vel_diff'
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
    elif dataset_type == 'DEALIAS_FOURDD':
        func_name = 'process_dealias_fourdd'
    elif dataset_type == 'DEALIAS_REGION':
        func_name = 'process_dealias_region_based'
    elif dataset_type == 'DEALIAS_UNWRAP':
        func_name = 'process_dealias_unwrap_phase'
    elif dataset_type == 'WIND_VEL':
        func_name = 'process_wind_vel'
    elif dataset_type == 'VAD':
        func_name = 'process_vad'
    elif dataset_type == 'WINDSHEAR':
        func_name = 'process_windshear'
    elif dataset_type == 'HYDROCLASS':
        func_name = 'process_hydroclass'
    elif dataset_type == 'ML_DETECTION':
        func_name = 'process_melting_layer'
        dsformat = 'ML'
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
    elif dataset_type == 'TIMEAVG_STD':
        func_name = 'process_time_avg_std'
        dsformat = 'OCCURRENCE'
    elif dataset_type == 'OCCURRENCE_PERIOD':
        func_name = 'process_occurrence_period'
        dsformat = 'OCCURRENCE'
    elif dataset_type == 'SUN_HITS':
        func_name = 'process_sun_hits'
        dsformat = 'SUN_HITS'
    elif dataset_type == 'POINT_MEASUREMENT':
        func_name = 'process_point_measurement'
        dsformat = 'TIMESERIES'
    elif dataset_type == 'ROI':
        func_name = process_roi
    elif dataset_type == 'TRAJ':
        func_name = 'process_trajectory'
        dsformat = 'TRAJ_ONLY'
    elif dataset_type == 'TRAJ_ATPLANE':
        func_name = 'process_traj_atplane'
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_ANTENNA_PATTERN':
        func_name = 'process_traj_antenna_pattern'
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_LIGHTNING':
        func_name = 'process_traj_lightning'
        dsformat = 'TIMESERIES'
    elif dataset_type == 'TRAJ_TRT':
        func_name = 'process_traj_trt'
    elif dataset_type == 'FIXED_RNG':
        func_name = process_fixed_rng
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
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break
    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    new_dataset = {'radar_out': deepcopy(radar_list[ind_rad])}

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
    new_dataset : dict
        dictionary containing the output
    ind_rad : int
        radar index

    """

    if procstatus != 1:
        return None, None

    for datatypedescr in dscfg['datatype']:
        radarnr, _, _, _, _ = get_datatype_fields(datatypedescr)
        break
    ind_rad = int(radarnr[5:8])-1
    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    new_dataset = {'radar_out': deepcopy(radar_list[ind_rad])}

    return new_dataset, ind_rad


def process_fixed_rng(procstatus, dscfg, radar_list=None):
    """
    Obtains radar data at a fixed range

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : list of strings. Dataset keyword
            The fields we want to extract
        rng : float. Dataset keyword
            The fixed range [m]
        RngTol : float. Dataset keyword
            The tolerance between the nominal range and the radar range
        ele_min, ele_max, azi_min, azi_max : floats. Dataset keyword
            The azimuth and elevation limits of the data [deg]

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

    field_names = []
    for datatypedescr in dscfg['datatype']:
        radarnr, _, datatype, _, _ = get_datatype_fields(
            datatypedescr)
        field_names.append(get_fieldname_pyart(datatype))
    ind_rad = int(radarnr[5:8])-1

    if (radar_list is None) or (radar_list[ind_rad] is None):
        warn('ERROR: No valid radar')
        return None, None
    radar = radar_list[ind_rad]

    # user defined parameters
    rng_tol = dscfg.get('RngTol', 50.)
    ele_min = dscfg.get('ele_min', None)
    ele_max = dscfg.get('ele_max', None)
    azi_min = dscfg.get('azi_min', None)
    azi_max = dscfg.get('azi_max', None)

    radar_aux = get_fixed_rng_data(
        radar, field_names, dscfg['rng'], rng_tol=rng_tol, ele_min=ele_min,
        ele_max=ele_max, azi_min=azi_min, azi_max=azi_max)

    if radar_aux is None:
        return None

    new_dataset = {'radar_out': radar_aux}

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
        radarnr, _, datatype, _, _ = get_datatype_fields(
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
        (_, yyyymmddHHMM, lon, lat, _, _, _, _, _, _, _, _, _, _, _, _, _, _,
         _, _, _, _, _, _, _, _, _, cell_contour) = read_trt_traj_data(
             dscfg['trtfile'])

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
    new_dataset = {'radar_out': deepcopy(radar)}

    new_dataset['radar_out'].range['data'] = radar.range['data'][inds_rng]
    new_dataset['radar_out'].ngates = inds_rng.size
    new_dataset['radar_out'].time['data'] = np.asarray(
        [new_dataset['radar_out'].time['data'][0]])
    new_dataset['radar_out'].scan_type = 'roi'
    new_dataset['radar_out'].sweep_mode['data'] = np.array(['roi'])
    new_dataset['radar_out'].sweep_start_ray_index['data'] = np.array(
        [0], dtype='int32')
    new_dataset['radar_out'].fixed_angle['data'] = np.array(
        [], dtype='float64')
    new_dataset['radar_out'].sweep_number['data'] = np.array(
        [0], dtype='int32')
    new_dataset['radar_out'].nsweeps = 1

    if radar.rays_are_indexed is not None:
        new_dataset['radar_out'].rays_are_indexed['data'] = np.array(
            [radar.rays_are_indexed['data'][0]])
    if radar.ray_angle_res is not None:
        new_dataset['radar_out'].ray_angle_res['data'] = np.array(
            [radar.ray_angle_res['data'][0]])

    new_dataset['radar_out'].sweep_end_ray_index['data'] = np.array(
        [1], dtype='int32')
    new_dataset['radar_out'].rays_per_sweep = np.array([1], dtype='int32')
    new_dataset['radar_out'].azimuth['data'] = np.array([], dtype='float64')
    new_dataset['radar_out'].elevation['data'] = np.array([], dtype='float64')
    new_dataset['radar_out'].nrays = 1

    new_dataset['radar_out'].gate_longitude['data'] = lon
    new_dataset['radar_out'].gate_latitude['data'] = lat
    new_dataset['radar_out'].gate_altitude['data'] = alt

    new_dataset['radar_out'].gate_x['data'] = (
        radar.gate_x['data'][inds_ray, inds_rng].T)
    new_dataset['radar_out'].gate_y['data'] = (
        radar.gate_y['data'][inds_ray, inds_rng].T)
    new_dataset['radar_out'].gate_z['data'] = (
        radar.gate_z['data'][inds_ray, inds_rng].T)

    new_dataset['radar_out'].fields = dict()
    field_dict = deepcopy(radar.fields[field_name])
    field_dict['data'] = radar.fields[field_name]['data'][inds_ray, inds_rng].T
    new_dataset['radar_out'].add_field(field_name, field_dict)

    return new_dataset['radar_out'], ind_rad


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
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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


def process_azimuthal_average(procstatus, dscfg, radar_list=None):
    """
    Averages radar data in azimuth obtaining and RHI as a result

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
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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
    angle = dscfg.get('angle', None)
    delta_azi = dscfg.get('delta_azi', None)
    avg_type = dscfg.get('avg_type', 'mean')
    nvalid_min = dscfg.get('nvalid_min', 1)
    if avg_type not in ('mean', 'median'):
        warn('Unsuported statistics '+avg_type)
        return None, None

    if delta_azi == -1:
        delta_azi = None
    if angle == -1:
        angle = None

    radar_aux = deepcopy(radar)
    # transform radar into ppi over the required elevation
    if radar_aux.scan_type == 'rhi':
        target_elevations, el_tol = get_target_elevations(radar_aux)
        radar_ppi = pyart.util.cross_section_rhi(
            radar_aux, target_elevations, el_tol=el_tol)
    elif radar_aux.scan_type == 'ppi':
        radar_ppi = radar_aux
    else:
        warn('Error: unsupported scan type.')
        return None, None

    # range, metadata, radar position are the same as the original
    # time
    radar_rhi = deepcopy(radar)
    radar_rhi.fields = dict()
    radar_rhi.scan_type = 'rhi'
    radar_rhi.sweep_number['data'] = np.array([0])
    radar_rhi.sweep_mode['data'] = np.array(['rhi'])
    radar_rhi.fixed_angle['data'] = np.array([0])
    radar_rhi.sweep_start_ray_index['data'] = np.array([0])
    radar_rhi.sweep_end_ray_index['data'] = np.array([radar_ppi.nsweeps-1])
    radar_rhi.rays_per_sweep['data'] = np.array([radar_ppi.nsweeps])
    radar_rhi.azimuth['data'] = np.ones(radar_ppi.nsweeps)
    radar_rhi.elevation['data'] = radar_ppi.fixed_angle['data']
    radar_rhi.nrays = radar_ppi.fixed_angle['data'].size
    radar_rhi.nsweeps = 1
    radar_rhi.rays_are_indexed = None
    radar_rhi.ray_angle_res = None

    # average radar data
    field_dict = pyart.config.get_metadata(field_name)
    field_dict['data'] = np.ma.masked_all((radar_ppi.nsweeps, radar_ppi.ngates))
    if angle is None:
        fixed_angle = np.zeros(radar_ppi.nsweeps)
    for sweep in range(radar_ppi.nsweeps):
        radar_aux = deepcopy(radar_ppi)
        radar_aux = radar_aux.extract_sweeps([sweep])

        # find neighbouring gates to be selected
        inds_ray, inds_rng = find_neighbour_gates(
            radar_aux, angle, None, delta_azi=delta_azi, delta_rng=None)

        # keep only data we are interested in
        field_aux = radar_aux.fields[field_name]['data'][:, inds_rng]
        field_aux = field_aux[inds_ray, :]

        vals, _ = compute_directional_stats(
            field_aux, avg_type=avg_type, nvalid_min=nvalid_min, axis=0)

        field_dict['data'][sweep, :] = vals
        if angle is None:
            fixed_angle[sweep] = np.median(radar_aux.azimuth['data'][inds_ray])
    if angle is None:
        radar_rhi.fixed_angle['data'] = np.array([np.mean(fixed_angle)])
    else:
        radar_rhi.fixed_angle['data'] = np.array([angle])
    radar_rhi.azimuth['data'] *= radar_rhi.fixed_angle['data'][0]
    radar_rhi.add_field(field_name, field_dict)

    # prepare for exit
    new_dataset = {'radar_out': radar_rhi}

    return new_dataset, ind_rad
