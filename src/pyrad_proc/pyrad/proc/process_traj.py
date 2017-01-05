"""
pyrad.proc.process_traj
=============================

Trajectory functions. Functions to pass trajectory dataset data to
the product generation functions.

.. autosummary::
    :toctree: generated/

    process_trajectory
    process_traj_atplane
    process_traj_antenna_pattern
"""

from warnings import warn
import numpy as np
from netCDF4 import num2date, date2num

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.io_aux import get_field_unit, get_field_name
from ..io.timeseries import TimeSeries
from ..io.read_data_other import read_antenna_pattern


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
    trajectory : Trajectory object
        Object holding time series
    ind_rad : int
        radar index

    """

    ang_tol = 1.2

    if procstatus == 0:
        # first call: nothing to do
        return None, None

    if procstatus == 2:
        # last call: do the products
        if (not dscfg['initialized']):
            warn('ERROR: No trajectory dataset available!')
            return None, None
        trajdict = dscfg['traj_atplane_dict']
        return trajdict['ts'], trajdict['ind_rad']

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

        description = [
            "Description:",
            "Time series of a weather radar data type at the location",
            "of the plane.",
            "The time samples where the plane was out of the weather radar",
            "sector are NOT included in this file.",
            "NaN (Not a number): No rain detected at the plane location."
        ]

        # Tmp: define timeformat="%Y-%m-%d %H:%M:%S.%f"
        ts = TimeSeries(description, maxlength=trajectory.time_vector.size,
                        datatype=datatype)

        unit = get_field_unit(datatype)
        name = get_field_name(datatype)
        # Append empty series: Note: sequence matters!
        ts.add_dataseries("at plane", name, unit, color='b')
        ts.add_dataseries("Mean", name, unit, color='r')
        ts.add_dataseries("Min", name, unit, color='k', linestyle=':')
        ts.add_dataseries("Max", name, unit, color='k', linestyle=':')
        ts.add_dataseries("#Valid", "", "", plot=False)

        trajdict = dict({
                'radar_traj': rad_traj,
                'radar_old': None,
                'radar_old2': None,
                'last_task_start_dt': None,
                'ind_rad': ind_rad,
                'ts': ts})
        traj_ind = trajectory.get_samples_in_period(end=dt_task_start)

        dscfg['traj_atplane_dict'] = trajdict
        dscfg['initialized'] = True
    else:
        # init already done
        trajdict = dscfg['traj_atplane_dict']
        rad_traj = trajdict['radar_traj']
        traj_ind = trajectory.get_samples_in_period(
            start=trajdict['last_task_start_dt'], end=dt_task_start)
        ts = trajdict['ts']

    traj_time_vec = date2num(trajectory.time_vector, radar.time['units'],
                             radar.time['calendar'])

    for tind in np.nditer(traj_ind):
        az = rad_traj.azimuth_vec[tind]
        el = rad_traj.elevation_vec[tind]
        rr = rad_traj.range_vec[tind]
        tt = traj_time_vec[tind]

        # Find closest azimuth and elevation ray
        (radar_sel, ray_sel, rr_ind, el_vec_rnd, az_vec_rnd) = \
            _get_closests_bin(az, el, rr, tt, radar, trajdict)

        # Check if traj sample is within scan sector
        if (_sample_out_of_sector(az, el, rr, radar_sel, ray_sel,
                                  rr_ind, el_vec_rnd, az_vec_rnd)):
            continue

        # ========================================================================
        # Get sample at bin
        rdata = radar_sel.fields[field_name]['data'].data
        val = rdata[ray_sel, rr_ind]

        # ========================================================================
        # Get samples around this cell (3x3 box)
        el_res = np.median(np.diff(el_vec_rnd))
        az_res = np.median(np.diff(az_vec_rnd))

        d_el = np.abs(radar_sel.elevation['data'] - radar_sel.elevation['data'][ray_sel])
        d_az = np.abs(radar_sel.azimuth['data'] - radar_sel.azimuth['data'][ray_sel])
        cell_ind = np.where((d_el < el_res*ang_tol) & (d_az < az_res*ang_tol))

        rr_min = rr_ind-1
        if (rr_min < 0):
            rr_min = 0
        cell_vals = rdata[cell_ind, rr_min:rr_ind+2]

        nvals_tot = cell_vals.size

        # ========================================================================
        # Exclude undefined values

        fill_val = radar_sel.fields[field_name]['_FillValue']
        cell_vals = np.ma.masked_values(cell_vals, fill_val)
        nvals_valid = cell_vals.count()
        val_mean = np.mean(cell_vals)
        val_min = np.min(cell_vals)
        val_max = np.max(cell_vals)

        # ========================================================================
        # Add to time series

        ts.add_timesample(trajectory.time_vector[tind],
                          (val, val_mean, val_min, val_max, nvals_valid))

        # end loop over traj samples within period

    trajdict['last_task_start_dt'] = dt_task_start
    trajdict['radar_old2'] = trajdict['radar_old']
    trajdict['radar_old'] = radar

    return None, None


def process_traj_antenna_pattern(procstatus, dscfg, radar_list=None, trajectory=None):
    """
    Process a new array of data volumes considering a plane
    trajectory. As result a timeseries with the values
    transposed for a given antenna pattern is created.
    The result is created when the LAST flag is set.

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
    trajectory : Trajectory object
        Object holding time series
    ind_rad : int
        radar index
    """

    if procstatus == 0:
        # first call: nothing to do
        return None, None

    if procstatus == 2:
        # last call: do the products
        if (not dscfg['initialized']):
            warn('ERROR: No trajectory dataset available!')
            return None, None
        # XXX
        tadict = dscfg['traj_antenna_dict']
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

        # Check config
        if ('antennaType' not in dscfg):
            raise Exception("ERROR: Undefined 'antennaType' for dataset '%s'"
                            % dscfg['dsname'])
        if ('configpath' not in dscfg):
            raise Exception("ERROR: Undefined 'configpath' for dataset '%s'"
                            % dscfg['dsname'])

        if (dscfg['antennaType'] == 'AZIMUTH'):
            is_azimuth_antenna = True
            info = 'parAzAnt'
            description = [
                "Antenna: PAR Azimuth antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the PAR",
                "antenna.",
                "The time samples where the plane was out of the weather radar",
                "sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]
            if ('par_azimuth_antenna' not in dscfg):
                raise Exception("ERROR: Undefined 'par_azimuth_antenna' for dataset '%s'"
                                % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['par_azimuth_antenna']['elPatternFile']
            fixed_angle = dscfg['par_azimuth_antenna']['fixed_angle']
        elif (dscfg['antennaType'] == 'ELEVATION'):
            is_azimuth_antenna = False
            info = 'parElAnt'
            description = [
                "Antenna: PAR Elevation antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the PAR",
                "antenna.",
                "The time samples where the plane was out of the weather radar",
                "sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]

            if ('par_elevation_antenna' not in dscfg):
                raise Exception("ERROR: Undefined 'par_elevation_antenna' for dataset '%s'"
                                % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['par_elevation_antenna']['azPatternFile']
            fixed_angle = dscfg['par_elevation_antenna']['fixed_angle']
        else:
            raise Exception("ERROR: Unexpected antenna type '%s' for dataset '%s'"
                            % (dscfg['antennaType'], dscfg['dsname']))

        # Get antenna pattern and make weight vector
        try:
            antpattern = read_antenna_pattern(patternfile, linear=True, twoway=True)
        except:
            raise
        pattern_angles = antpattern['angle'] + fixed_angle
        if (not is_azimuth_antenna):
            pattern_angles[pattern_angles < 0] += 360.
        pattern_angles[pattern_angles >= 360.] -= 360.

        if (is_azimuth_antenna):
            scan_angles = np.sort(np.unique(radar.elevation['data'].round(decimals=1)))
        else:
            scan_angles = np.sort(np.unique(radar.azimuth['data'].round(decimals=1)))

        weightvec = np.empty(scan_angles.size, dtype=float)
        for kk in range(scan_angles.size):
            ind = np.argmin(np.abs(pattern_angles - scan_angles[kk]))
            weightvec[kk] = antpattern['attenuation'][ind]

        ts = TimeSeries(description, maxlength=trajectory.time_vector.size,
                        datatype=datatype)

        unit = get_field_unit(datatype)
        name = get_field_name(datatype)

        # XXX

        tadict = dict({
                'radar_traj': rad_traj,
                'radar_old': None,
                'radar_old2': None,
                'last_task_start_dt': None,
                'ind_rad': ind_rad,
                'is_azimuth_antenna': is_azimuth_antenna,
                'info': info,
                'weightvec': weightvec,
                'ts': ts})
        traj_ind = trajectory.get_samples_in_period(end=dt_task_start)

        dscfg['traj_antenna_dict'] = tadict
        dscfg['initialized'] = True
        # end init
    else:
        # init already done
        tadict = dscfg['traj_antenna_dict']
        rad_traj = tadict['radar_traj']
        is_azimuth_antenna = tadict['is_azimuth_antenna']
        weightvec = tadict['weightvec']
        ts = tadict['ts']

        traj_ind = trajectory.get_samples_in_period(
            start=tadict['last_task_start_dt'], end=dt_task_start)

    traj_time_vec = date2num(trajectory.time_vector, radar.time['units'],
                             radar.time['calendar'])

    for tind in np.nditer(traj_ind):
        az = rad_traj.azimuth_vec[tind]
        el = rad_traj.elevation_vec[tind]
        rr = rad_traj.range_vec[tind]
        tt = traj_time_vec[tind]

        # Find closest azimuth and elevation ray
        (radar_sel, ray_sel, rr_ind, el_vec_rnd, az_vec_rnd) = \
             _get_closests_bin(az, el, rr, tt, radar, tadict)

        # Check if traj sample is within scan sector
        if (_sample_out_of_sector(az, el, rr, radar_sel, ray_sel,
                                  rr_ind, el_vec_rnd, az_vec_rnd)):
            continue

        # XXX

        # end loop over traj samples within period

    tadict['last_task_start_dt'] = dt_task_start
    tadict['radar_old2'] = tadict['radar_old']
    tadict['radar_old'] = radar

    return None, None


def _get_closests_bin(az, el, rr, tt, radar, tdict):
    """
    """

    daz = np.abs(radar.azimuth['data'] - az)
    dele = np.abs(radar.elevation['data'] - el)
    dangle = np.sqrt(daz**2 + dele**2)
    ray_ind = np.argmin(dangle)
    dt = tt - radar.time['data'][ray_ind]

    if (tdict['radar_old'] is not None):
        rad_old = tdict['radar_old']

        daz_old = np.abs(rad_old.azimuth['data'] - az)
        dele_old = np.abs(rad_old.elevation['data'] - el)
        dangle_old = np.sqrt(daz_old**2 + dele_old**2)
        ray_ind_old = np.argmin(dangle_old)
        dt_old = tt - rad_old.time['data'][ray_ind_old]

        # Find closest time
        if ((dt_old < 0.) and (tdict['radar_old2'] is not None)):
            rad_old2 = tdict['radar_old2']
            daz_old2 = np.abs(rad_old2.azimuth['data'] - az)
            dele_old2 = np.abs(rad_old2.elevation['data'] - el)
            dangle_old2 = np.sqrt(daz_old2**2 + dele_old2**2)
            ray_ind_old2 = np.argmin(dangle_old2)
            dt_old2 = tt - rad_old2.time['data'][ray_ind_old2]

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

    el_vec_rnd = np.sort(np.unique(radar_sel.elevation['data'].round(decimals=1)))
    az_vec_rnd = np.sort(np.unique(radar_sel.azimuth['data'].round(decimals=1)))

    return (radar_sel, ray_sel, rr_ind, el_vec_rnd, az_vec_rnd)


def _sample_out_of_sector(az, el, rr, radar_sel, ray_sel, rr_ind,
                          el_vec_rnd, az_vec_rnd):
    """
    """

    # Check if sample is within sector
    rr_min = radar_sel.range['data'][rr_ind]
    range_res = radar_sel.range['data'][1] - radar_sel.range['data'][0]
    if (np.abs(rr_min - rr) > (2*range_res)):
        # print("INFO: Trajectory sample out of range")
        return True

    if (((az_vec_rnd[0] - az) > 3.0) or ((az - az_vec_rnd[-1]) > 3.0)):
        # print("INFO: Trajectory sample out of azimuth angle")
        return True

    if (((el_vec_rnd[0] - el) > 3.0) or ((el - el_vec_rnd[-1]) > 3.0)):
        # print("INFO: Trajectory sample out of elevation angle")
        return True

    return False
