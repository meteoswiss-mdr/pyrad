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
    process_trajectory
    process_traj_atplane

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart

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
    elif dataset_type == 'PHIDP0_CORRECTION':
        func_name = 'process_correct_phidp0'
    elif dataset_type == 'PHIDP_SMOOTH_1W':
        func_name = 'process_smooth_phidp_single_window'
    elif dataset_type == 'PHIDP_SMOOTH_2W':
        func_name = 'process_smooth_phidp_double_window'
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
    elif dataset_type == 'HYDROCLASS':
        func_name = 'process_hydroclass'
    elif dataset_type == 'PHIDP0_ESTIMATE':
        func_name = 'process_estimate_phidp0'
    elif dataset_type == 'RHOHV_RAIN':
        func_name = 'process_rhohv_rain'
    elif dataset_type == 'ZDR_RAIN':
        func_name = 'process_zdr_rain'
    elif dataset_type == 'SELFCONSISTENCY_KDP_PHIDP':
        func_name = 'process_selfconsistency_kdp_phidp'
    elif dataset_type == 'SELFCONSISTENCY_BIAS':
        func_name = 'process_selfconsistency_bias'
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

    return new_dataset, ind_rad


def process_trajectory(procstatus, dscfg, radar_list=None, trajectory=None):
    """
    Return trajectory

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
    radar_list : list of Radar objects
        Optional. list of radar objects
    trajectory : Trajectory object
        containing trajectory samples

    Returns
    -------
    new_dataset : Trajectory object
        radar object
    ind_rad : int
        None

    """

    if procstatus != 1:
        return None, None

    if (not dscfg['initialized']):
        if (trajectory is None):
            raise Exception("ERROR: Undefined trajectory for dataset '%s'"
                            % dscfg['dsname'])

        if (radar_list is not None):
            for radar in radar_list:
                rad = trajectory.add_radar(radar)
                trajectory.calculate_velocities(rad)
        else:
            warn('ERROR: No valid radar found')
            return None, None

        dscfg['initialized'] = True
        return trajectory, None

    return None, None


def process_traj_atplane(procstatus, dscfg, radar_list=None, trajectory=None):
    """
    Return time series according to trajectory

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration
    radar_list : list of Radar objects
        Optional. list of radar objects
    trajectory : Trajectory object
        containing trajectory samples

    Returns
    -------
    new_dataset : Trajectory object
        radar object
    ind_rad : int
        None

    """

    if procstatus == 0:
        # first call: nothing to do
        return None, None

    if procstatus == 2:
        # last call: do the products
        # XXX
        return None, None

    # Process
    radarnr, datagroup, datatype, dataset, product = get_datatype_fields(
        dscfg['datatype'][0])
    field_name = get_fieldname_pyart(datatype)
    ind_rad = int(radarnr[5:8])-1
    if ((radar_list is None) or (radar_list[ind_rad] is None)):
        warn('ERROR: No valid radar found')
        return None, None
    radar = radar_list[ind_rad]

    if field_name not in radar.fields:
        warn("Datatype '%s' not available in radar data" % field_name)
        return None, None

    ttask_start = radar.time['data'].min()
    dt_task_start = num2date(ttask_start, radar.time['units'],
                             radar.time['calendar'])
    if (not dscfg['initialized']):
        # init
        if (trajectory is None):
            raise Exception("ERROR: Undefined trajectory for dataset '%s'"
                            % dscfg['dsname'])

        rad_traj = trajectory.add_radar(radar)

        trajdict = dict({
                'radar_traj': rad_traj,
                'radar_old': None,
                'radar_old2': None,
                'last_task_start_dt': None})
        traj_ind = trajectory.get_samples_in_period(end=dt_task_start)

        dscfg['traj_atplane_dict'] = trajdict
        dscfg['initialized'] = True
    else:
        trajdict = dscfg['traj_atplane_dict']
        rad_traj = trajdict['radar_traj']
        traj_ind = trajectory.get_samples_in_period(
            start=trajdict['last_task_start_dt'], end=dt_task_start)

    traj_time_vec = date2num(trajectory.time_vector, radar.time['units'],
                    radar.time['calendar'])
        
    for tind in np.nditer(traj_ind):
        az = rad_traj.azimuth_vec[tind]
        el = rad_traj.elevation_vec[tind]
        rr = rad_traj.range_vec[tind]
        tt = traj_time_vec[tind]

        # Find closest azimuth and elevation ray
        daz = np.abs(radar.azimuth['data'] - az)
        dele = np.abs(radar.elevation['data'] - el)
        dangle = np.sqrt(daz**2 + dele**2)
        ray_ind = np.argmin(dangle)
        dangle_min = dangle[ray_ind]
        dt = tt - radar.time['data'][ray_ind]

        # XXX
        if (dangle_min > 3.0):
            # print("INFO: Trajectory sample out of angle")
            continue

        # XXX
        print("=== ", end="")
        print(tind, az, el, rr, tt)
        #print(ray_ind, radar.azimuth['data'][ray_ind],
        #      radar.elevation['data'][ray_ind], dt)

        if (trajdict['radar_old'] is not None):
            rad_old = trajdict['radar_old']

            daz_old = np.abs(rad_old.azimuth['data'] - az)
            dele_old = np.abs(rad_old.elevation['data'] - el)
            dangle_old = np.sqrt(daz_old**2 + dele_old**2)
            ray_ind_old = np.argmin(dangle_old)
            dangle_min_old = dangle[ray_ind_old]
            dt_old = tt - rad_old.time['data'][ray_ind_old]

            # XXX
            #print(ray_ind_old, rad_old.azimuth['data'][ray_ind_old],
            #      rad_old.elevation['data'][ray_ind_old], dt_old)

            # Find closest time

            if ((dt_old < 0.) and (trajdict['radar_old2'] is not None)):
                rad_old2 = trajdict['radar_old2']
                daz_old2 = np.abs(rad_old2.azimuth['data'] - az)
                dele_old2 = np.abs(rad_old2.elevation['data'] - el)
                dangle_old2 = np.sqrt(daz_old2**2 + dele_old2**2)
                ray_ind_old2 = np.argmin(dangle_old2)
                dangle_min_old2 = dangle[ray_ind_old2]
                dt_old2 = tt - rad_old2.time['data'][ray_ind_old2]

                # XXX
                #print(ray_ind_old2, rad_old2.azimuth['data'][ray_ind_old2],
                #      rad_old2.elevation['data'][ray_ind_old2], dt_old2)

                if ((np.abs(dt_old) < np.abs(dt_old2))):
                    radar_sel = rad_old
                    ray_sel = ray_ind_old
                else:
                    radar_sel = rad_old2
                    ray_sel = ray_ind_old2
            else:
                if ((np.abs(dt) < np.abs(dt_old))):
                    radar_sel = radar
                    ray_sel = ray_ind
                else:
                    radar_sel = rad_old
                    ray_sel = ray_ind_old
        else:
            radar_sel = radar
            ray_sel = ray_ind

        # Find closest range bin
        rr_ind = np.argmin(np.abs(radar_sel.range['data'] - rr))
        rr_min = radar_sel.range['data'][rr_ind]

        print(rr_ind, rr_min) # XXX

        # Get sample at bin
        val = radar_sel.fields[field_name]['data'].data[ray_sel, rr_ind]

        #print(val) # XXX

        # Get samples around this cell (3x3 box)

        #print(radar_sel.elevation['data']) #XXX
        el_vec = np.sort(np.unique(radar_sel.elevation['data'].round(decimals=1)))
        print(el_vec) #XXX

        #print(radar_sel.azimuth['data']) #XXX
        az_vec = np.sort(np.unique(radar_sel.azimuth['data'].round(decimals=1)))
        print(az_vec) #XXX


        # end loop over traj samples within period

    trajdict['last_task_start_dt'] = dt_task_start
    trajdict['radar_old2'] = trajdict['radar_old']
    trajdict['radar_old'] = radar

    return None, None
