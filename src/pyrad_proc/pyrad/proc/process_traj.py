"""
pyrad.proc.process_traj
=============================

Trajectory functions. Functions to pass trajectory dataset data to
the product generation functions.

.. autosummary::
    :toctree: generated/

    process_trajectory
    process_traj_lightning
    process_traj_atplane
    process_traj_antenna_pattern
    _get_closests_bin
    _sample_out_of_sector
    Target_Radar
"""

from warnings import warn
import sys
import gc

import numpy as np
from netCDF4 import num2date, date2num

from pyart.config import get_metadata
from pyart.core import Radar
from pyart.util import colocated_gates, intersection

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..io.io_aux import get_field_unit, get_field_name
from ..io.timeseries import TimeSeries
from ..io.read_data_other import read_antenna_pattern

from ..util.stat_utils import quantiles_weighted


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

    if not dscfg['initialized']:
        if trajectory is None:
            raise Exception("ERROR: Undefined trajectory for dataset '%s'"
                            % dscfg['dsname'])

        if radar_list is not None:
            for radar in radar_list:
                rad = trajectory.add_radar(radar)
                trajectory.calculate_velocities(rad)
        else:
            warn('ERROR: No valid radar found')
            return None, None

        dscfg['initialized'] = True
        return trajectory, None

    return None, None


def process_traj_lightning(procstatus, dscfg, radar_list=None,
                           trajectory=None):
    """
    Return time series according to lightning trajectory

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
        if not dscfg['initialized']:
            warn('ERROR: No trajectory dataset available!')
            return None, None
        trajdict = dscfg['traj_atplane_dict']

        dataset = {
            'ts': trajdict['ts'],
            'final': 1
        }
        return dataset, trajdict['ind_rad']

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

    if not dscfg['initialized']:
        # init
        if trajectory is None:
            raise Exception("ERROR: Undefined trajectory for dataset '%s'"
                            % dscfg['dsname'])

        rad_traj = trajectory.add_radar(radar)

        description = [
            "Description:",
            "Time series of a weather radar data type at the location",
            "of the lightning flash.",
            "The time samples where the flash was out of the weather radar",
            "sector are NOT included in this file.",
            "NaN (Not a number): No rain detected at the plane location."
        ]

        # Tmp: define timeformat="%Y-%m-%d %H:%M:%S.%f"
        ts = TimeSeries(description, maxlength=trajectory.time_vector.size,
                        datatype=datatype)

        unit = get_field_unit(datatype)
        name = get_field_name(datatype)
        # Append empty series: Note: sequence matters!
        ts.add_dataseries("#Flash", "", "", plot=False)
        ts.add_dataseries("Power", "", "dBm", plot=False)
        ts.add_dataseries("at_flash", name, unit, color='b')
        ts.add_dataseries("Mean", name, unit, color='r')
        ts.add_dataseries("Min", name, unit, color='k', linestyle=':')
        ts.add_dataseries("Max", name, unit, color='k', linestyle=':')
        ts.add_dataseries("#Valid", "", "", plot=False)

        data_is_log = False
        if 'data_is_log' in dscfg:
            data_is_log = (dscfg['data_is_log'] != 0)

        trajdict = dict({
            'radar_traj': rad_traj,
            'radar_old': None,
            'radar_old2': None,
            'last_task_start_dt': None,
            'ind_rad': ind_rad,
            'ts': ts,
            'data_is_log': data_is_log})
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
        data_is_log = trajdict['data_is_log']

    traj_time_vec = date2num(trajectory.time_vector, radar.time['units'],
                             radar.time['calendar'])

    if np.size(traj_ind) == 0:
        warn('No trajectory samples within current period')

        trajdict['radar_old2'] = trajdict['radar_old']
        trajdict['radar_old'] = radar
        return None, None

    az_list = []
    el_list = []
    rr_list = []
    tt_list = []
    for tind in np.nditer(traj_ind):
        az = rad_traj.azimuth_vec[tind]
        el = rad_traj.elevation_vec[tind]
        rr = rad_traj.range_vec[tind]
        tt = traj_time_vec[tind]
        dBm = trajectory.dBm[tind]
        flashnr = trajectory.flashnr_vec[tind]

        # Find closest azimuth and elevation ray
        (radar_sel, ray_sel, rr_ind, el_vec_rnd, az_vec_rnd) = \
            _get_closests_bin(az, el, rr, tt, radar, trajdict)

        # Check if traj sample is within scan sector
        if (_sample_out_of_sector(az, el, rr, radar_sel, ray_sel,
                                  rr_ind, el_vec_rnd, az_vec_rnd)):
            continue

        # =====================================================================
        # Get sample at bin
        rdata = radar_sel.fields[field_name]['data']
        val = rdata[ray_sel, rr_ind]

        # =====================================================================
        # Get samples around this cell (3x3 box)
        el_res = np.median(np.diff(el_vec_rnd))
        az_res = np.median(np.diff(az_vec_rnd))

        d_el = np.abs(radar_sel.elevation['data'] -
                      radar_sel.elevation['data'][ray_sel])
        d_az = np.abs(radar_sel.azimuth['data'] -
                      radar_sel.azimuth['data'][ray_sel])
        cell_ind = np.where((d_el < el_res*ang_tol) & (d_az < az_res*ang_tol))

        rr_min = rr_ind-1
        if rr_min < 0:
            rr_min = 0
        cell_vals = rdata[cell_ind, rr_min:rr_ind+2]

        # =====================================================================
        # Compute statistics and get number of valid data
        if data_is_log:
            vals_lin = np.ma.power(10., cell_vals/10.)
            val_mean = np.ma.mean(vals_lin)
            val_mean = 10. * np.ma.log10(val_mean)
        else:
            val_mean = np.ma.mean(cell_vals)
        val_min = np.ma.min(cell_vals)
        val_max = np.ma.max(cell_vals)

        nvals_valid = np.count_nonzero(
            np.logical_not(np.ma.getmaskarray(cell_vals)))
        # =====================================================================
        # Add to time series

        ts.add_timesample(
            trajectory.time_vector[tind],
            (flashnr, dBm, val, val_mean, val_min, val_max, nvals_valid))

        az_list.append(az)
        el_list.append(el)
        rr_list.append(rr)
        tt_list.append(trajectory.time_vector[tind])
        # end loop over traj samples within period

    # output radar volume and flash coordinates respect to radar
    radar_sel = radar
    if trajdict['radar_old'] is not None:
        radar_sel = trajdict['radar_old']

    dataset = {
        'azi_traj': np.asarray(az_list),
        'ele_traj': np.asarray(el_list),
        'rng_traj': np.asarray(rr_list),
        'time_traj': np.asarray(tt_list),
        'radar': radar_sel,
        'final': 0
    }

    # update trajectory dictionary
    trajdict['last_task_start_dt'] = dt_task_start
    trajdict['radar_old2'] = trajdict['radar_old']
    trajdict['radar_old'] = radar

    return dataset, ind_rad


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
        if not dscfg['initialized']:
            warn('ERROR: No trajectory dataset available!')
            return None, None
        trajdict = dscfg['traj_atplane_dict']

        dataset = {
            'ts': trajdict['ts'],
            'final': 1
        }
        return dataset, trajdict['ind_rad']

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
    if not dscfg['initialized']:
        # init
        if trajectory is None:
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

        data_is_log = False
        if 'data_is_log' in dscfg:
            data_is_log = (dscfg['data_is_log'] != 0)

        trajdict = dict({
            'radar_traj': rad_traj,
            'radar_old': None,
            'radar_old2': None,
            'last_task_start_dt': None,
            'ind_rad': ind_rad,
            'ts': ts,
            'data_is_log': data_is_log})
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
        data_is_log = trajdict['data_is_log']

    traj_time_vec = date2num(trajectory.time_vector, radar.time['units'],
                             radar.time['calendar'])

    if np.size(traj_ind) == 0:
        warn('No trajectory samples within current period')

        trajdict['radar_old2'] = trajdict['radar_old']
        trajdict['radar_old'] = radar
        return None, None

    az_list = []
    el_list = []
    rr_list = []
    tt_list = []
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

        # =====================================================================
        # Get sample at bin
        rdata = radar_sel.fields[field_name]['data']
        val = rdata[ray_sel, rr_ind]

        # =====================================================================
        # Get samples around this cell (3x3 box)
        el_res = np.median(np.diff(el_vec_rnd))
        az_res = np.median(np.diff(az_vec_rnd))

        d_el = np.abs(radar_sel.elevation['data'] -
                      radar_sel.elevation['data'][ray_sel])
        d_az = np.abs(radar_sel.azimuth['data'] -
                      radar_sel.azimuth['data'][ray_sel])
        cell_ind = np.where((d_el < el_res*ang_tol) & (d_az < az_res*ang_tol))

        rr_min = rr_ind-1
        if rr_min < 0:
            rr_min = 0
        cell_vals = rdata[cell_ind, rr_min:rr_ind+2]

        # =====================================================================
        # Compute statistics and get number of valid data
        if data_is_log:
            vals_lin = np.ma.power(10., cell_vals/10.)
            val_mean = np.ma.mean(vals_lin)
            val_mean = 10. * np.ma.log10(val_mean)
        else:
            val_mean = np.ma.mean(cell_vals)
        val_min = np.ma.min(cell_vals)
        val_max = np.ma.max(cell_vals)

        nvals_valid = np.count_nonzero(
            np.logical_not(np.ma.getmaskarray(cell_vals)))
        # =====================================================================
        # Add to time series

        ts.add_timesample(trajectory.time_vector[tind],
                          (val, val_mean, val_min, val_max, nvals_valid))

        az_list.append(az)
        el_list.append(el)
        rr_list.append(rr)
        tt_list.append(trajectory.time_vector[tind])
        # end loop over traj samples within period

    # output radar volume and antenna coordinates respect to radar
    radar_sel = radar
    if trajdict['radar_old'] is not None:
        radar_sel = trajdict['radar_old']

    dataset = {
        'azi_traj': np.asarray(az_list),
        'ele_traj': np.asarray(el_list),
        'rng_traj': np.asarray(rr_list),
        'time_traj': np.asarray(tt_list),
        'radar': radar_sel,
        'final': 0
    }

    # update trajectory dictionary
    trajdict['last_task_start_dt'] = dt_task_start
    trajdict['radar_old2'] = trajdict['radar_old']
    trajdict['radar_old'] = radar

    return dataset, ind_rad


def process_traj_antenna_pattern(procstatus, dscfg, radar_list=None,
                                 trajectory=None):
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
        if not dscfg['initialized']:
            warn('ERROR: No trajectory dataset available!')
            return None, None
        tadict = dscfg['traj_antenna_dict']
        dataset = {
            'ts': tadict['ts'],
            'final': 1
        }
        return dataset, tadict['ind_rad']

    # Tolerance in azimuth to find traj for the reference radar, also in range
    az_traj_tol = 10
    rg_traj_tol = 10000 #m

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
    if not dscfg['initialized']:
        # === init ============================================================
        if trajectory is None:
            raise Exception("ERROR: Undefined trajectory for dataset '%s'"
                            % dscfg['dsname'])

        # Check config
        if 'antennaType' not in dscfg:
            raise Exception("ERROR: Undefined 'antennaType' for dataset '%s'"
                            % dscfg['dsname'])
        if 'configpath' not in dscfg:
            raise Exception("ERROR: Undefined 'configpath' for dataset '%s'"
                            % dscfg['dsname'])

        if 'target_radar_pos' not in dscfg:
            radar_antenna_atsameplace = True
            target_radar = None

            warn('No target radar position specified. ' +
                 'The radars are assumed colocated')
            rad_traj = trajectory.add_radar(radar)
        else:
            radar_antenna_atsameplace = False

            # create dummy radar object with target radar specs
            latitude = get_metadata('latitude')
            longitude = get_metadata('longitude')
            altitude = get_metadata('altitude')

            latitude['data'] = np.array(
                [dscfg['target_radar_pos']['latitude']], dtype='float64')
            longitude['data'] = np.array(
                [dscfg['target_radar_pos']['longitude']], dtype='float64')
            altitude['data'] = np.array(
                [dscfg['target_radar_pos']['altitude']], dtype='float64')

            target_radar = Target_Radar(latitude, longitude, altitude)
            rad_traj = trajectory.add_radar(target_radar)

        if dscfg['antennaType'] == 'AZIMUTH':
            is_azimuth_antenna = True
            info = 'parAzAnt'
            description = [
                "Antenna: PAR Azimuth antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the PAR",
                "antenna.",
                "The time samples where the plane was out of the weather",
                "radar sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]
            if 'par_azimuth_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'par_azimuth_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['par_azimuth_antenna']['elPatternFile']
            fixed_angle = dscfg['par_azimuth_antenna']['fixed_angle']

        elif dscfg['antennaType'] == 'ELEVATION':
            is_azimuth_antenna = False
            info = 'parElAnt'
            description = [
                "Antenna: PAR Elevation antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the PAR",
                "antenna.",
                "The time samples where the plane was out of the weather ",
                "radar sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]

            if 'par_elevation_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'par_elevation_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['par_elevation_antenna']['azPatternFile']
            fixed_angle = dscfg['par_elevation_antenna']['fixed_angle']

        elif dscfg['antennaType'] == 'LOWBEAM':
            is_azimuth_antenna = True
            info = 'asrLowBeamAnt'
            description = [
                "Antenna: ASR low beam antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the ASR",
                "antenna.",
                "The time samples where the plane was out of the weather",
                "radar sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]
            if 'asr_lowbeam_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'asr_lowbeam_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['asr_lowbeam_antenna']['elPatternFile']
            fixed_angle = dscfg['asr_lowbeam_antenna']['fixed_angle']

        elif dscfg['antennaType'] == 'HIGHBEAM':
            is_azimuth_antenna = True
            info = 'asrHighBeamAnt'
            description = [
                "Antenna: ASR high beam antenna",
                "Description:",
                "Time series of a weather radar data type at the location",
                "of the plane weighted by the antenna pattern of the ASR",
                "antenna.",
                "The time samples where the plane was out of the weather",
                "radar sector are NOT included in this file.",
                "NaN (Not a number): No rain detected at the plane location."
            ]
            if 'asr_highbeam_antenna' not in dscfg:
                raise Exception("ERROR: Undefined 'asr_highbeam_antenna' for"
                                " dataset '%s'" % dscfg['dsname'])

            patternfile = dscfg['configpath'] + 'antenna/' \
                + dscfg['asr_highbeam_antenna']['elPatternFile']
            patternfile_low = dscfg['configpath'] + 'antenna/' \
                + dscfg['asr_lowbeam_antenna']['elPatternFile']
            fixed_angle = dscfg['asr_highbeam_antenna']['fixed_angle']
        else:
            raise Exception("ERROR: Unexpected antenna type '%s' for dataset"
                            " '%s'" % (dscfg['antennaType'], dscfg['dsname']))

        # Read dataset config parameters:
        nan_value = dscfg.get('nan_value', 0.)
        use_nans = False
        if 'use_nans' in dscfg:
            if dscfg['use_nans'] != 0:
                use_nans = True

        weight_threshold = dscfg.get('pattern_thres', 0.)

        do_all_ranges = False
        if 'range_all' in dscfg:
            if dscfg['range_all'] != 0:
                do_all_ranges = True

        # Config parameters for processing when the weather radar and the
        # antenna are not at the same place:
        max_altitude = dscfg.get('max_altitude', 12000.) # [m]
        rhi_resolution = dscfg.get('rhi_resolution', 0.5)  # [deg]
        latlon_tol = dscfg.get('latlon_tol', 0.04)  # [deg]
        alt_tol = dscfg.get('alt_tol', 1000.)  # [m]

        data_is_log = False
        if 'data_is_log' in dscfg:
            data_is_log = (dscfg['data_is_log'] != 0)

        # Get antenna pattern and make weight vector
        try:
            if info == 'asrHighBeamAnt':
                antpattern = read_antenna_pattern(
                    patternfile, linear=True, twoway=False)
                antpattern_low = read_antenna_pattern(
                    patternfile_low, linear=True, twoway=False)
                antpattern['attenuation'] *= antpattern_low['attenuation']
            else:
                antpattern = read_antenna_pattern(patternfile, linear=True,
                                                  twoway=True)
        except:
            raise

        pattern_angles = antpattern['angle'] + fixed_angle
        if not is_azimuth_antenna:
            pattern_angles[pattern_angles < 0] += 360.
        pattern_angles[pattern_angles >= 360.] -= 360.

        if radar_antenna_atsameplace:
            if is_azimuth_antenna:
                scan_angles = np.sort(np.unique(
                    radar.elevation['data'].round(decimals=1)))
            else:
                scan_angles = np.sort(np.unique(
                    radar.azimuth['data'].round(decimals=1)))
        else:
            scan_angles = np.arange(0, 90, rhi_resolution, dtype=float)

        weightvec = np.empty(scan_angles.size, dtype=float)
        for kk in range(scan_angles.size):
            ind = np.argmin(np.abs(pattern_angles - scan_angles[kk]))
            weightvec[kk] = antpattern['attenuation'][ind]

        ts = TimeSeries(description, maxlength=trajectory.time_vector.size,
                        datatype=datatype)

        unit = get_field_unit(datatype)
        name = get_field_name(datatype)
        # Quantiles of interest
        quantiles = [{"val": 0.1, "plot": False, "color": None, "ltype": None},
                     {"val": 0.2, "plot": True, "color": 'k', "ltype": ':'},
                     {"val": 0.3, "plot": False, "color": None, "ltype": None},
                     {"val": 0.4, "plot": False, "color": None, "ltype": None},
                     {"val": 0.5, "plot": True, "color": 'r', "ltype": None},
                     {"val": 0.6, "plot": False, "color": None, "ltype": None},
                     {"val": 0.7, "plot": False, "color": None, "ltype": None},
                     {"val": 0.8, "plot": True, "color": 'k', "ltype": ':'},
                     {"val": 0.9, "plot": False, "color": None, "ltype": None},
                     {"val": 0.95, "plot": False, "color": None,
                      "ltype": None}]

        ts.add_dataseries("Weighted average", name, unit, color='b')
        for qq in quantiles:
            label = "Quantile_%4.2f" % qq["val"]
            ts.add_dataseries(label, name, unit, plot=qq["plot"],
                              color=qq["color"], linestyle=qq["ltype"])
        ts.add_dataseries("#Valid", "", "", plot=False)

        quants = np.array([ee['val'] for ee in quantiles])

        # Persistent data structure
        tadict = dict({
            'radar_traj': rad_traj,
            'radar_old': None,
            'radar_old2': None,
            'target_radar': target_radar,
            'last_task_start_dt': None,
            'ind_rad': ind_rad,
            'is_azimuth_antenna': is_azimuth_antenna,
            'info': info,
            'scan_angles': scan_angles,
            'radar_antenna_atsameplace': radar_antenna_atsameplace,
            'weightvec': weightvec,
            'quantiles': quants,
            'use_nans': use_nans,
            'nan_value': nan_value,
            'weight_threshold': weight_threshold,
            'do_all_ranges': do_all_ranges,
            'max_altitude': max_altitude,
            'latlon_tol': latlon_tol,
            'alt_tol': alt_tol,
            'data_is_log': data_is_log,
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
        scan_angles = tadict['scan_angles']
        radar_antenna_atsameplace = tadict['radar_antenna_atsameplace']
        weightvec = tadict['weightvec']
        nan_value = tadict['nan_value']
        use_nans = tadict['use_nans']
        weight_threshold = tadict['weight_threshold']
        do_all_ranges = tadict['do_all_ranges']
        ts = tadict['ts']
        target_radar = tadict['target_radar']
        max_altitude = tadict['max_altitude']
        latlon_tol = tadict['latlon_tol']
        alt_tol = tadict['alt_tol']
        data_is_log = tadict['data_is_log']

        traj_ind = trajectory.get_samples_in_period(
            start=tadict['last_task_start_dt'], end=dt_task_start)

    traj_time_vec = date2num(trajectory.time_vector, radar.time['units'],
                             radar.time['calendar'])

    if np.size(traj_ind) == 0:
        warn('No trajectory samples within current period')
        return None, None

    for tind in np.nditer(traj_ind):
        az = rad_traj.azimuth_vec[tind]
        el = rad_traj.elevation_vec[tind]
        rr = rad_traj.range_vec[tind]
        tt = traj_time_vec[tind]

        # Select radar object, find closest azimuth and elevation ray
        (radar_sel, ray_sel, rr_ind, el_vec_rnd, az_vec_rnd) = \
            _get_closests_bin(az, el, rr, tt, radar, tadict)

        # Check if traj sample is within scan sector
        if (_sample_out_of_sector(az, el, rr, radar_sel, ray_sel,
                                  rr_ind, el_vec_rnd, az_vec_rnd)):
            continue

        if radar_antenna_atsameplace:
            # =====================================================================
            # Radar and scanning antenna are at the SAME place
            # ===================================================================

            # =====================================================================
            # Get sample at bin

            if is_azimuth_antenna:
                angles = radar_sel.azimuth['data']
                angles_scan = radar_sel.elevation['data']
                ray_angle = radar_sel.azimuth['data'][ray_sel]
            else:
                angles = radar_sel.elevation['data']
                angles_scan = radar_sel.azimuth['data']
                ray_angle = radar_sel.elevation['data'][ray_sel]

            d_angle = np.abs(angles - ray_angle)
            ray_inds = np.where(d_angle < 0.09)[0]
            angles_sortind = np.argsort(angles_scan[ray_inds])

            ray_inds = ray_inds[angles_sortind]
            angles_sorted = angles_scan[ray_inds]

            # Set default values
            avg = None
            qvals = np.array([None] * tadict['quantiles'].size)
            nvals_valid = None

            if ((scan_angles.size != angles_sorted.size) or
                    (np.max(np.abs(scan_angles - angles_sorted)) > 0.1)):
                print("WARNING: Scan angle mismatch!", file=sys.stderr)
                ts.add_timesample(trajectory.time_vector[tind],
                                  (np.concatenate([[avg], qvals,
                                                   [nvals_valid]])))
                continue

            if do_all_ranges:
                rr_ind_min = 0
            else:
                rr_ind_min = rr_ind

            w_vec = tadict['weightvec']
            rdata = radar_sel.fields[field_name]['data']
            values = rdata[ray_inds, rr_ind_min:rr_ind+1]

        else:
            # ===================================================================
            # Radar and scanning antenna are NOT at the same place
            # ===================================================================

            # Create dummy Radar object with fix az and r, and all elevations

            n_rays = scan_angles.size

            r_time = {'data': [tt * n_rays]}
            r_range = {'data': [rr]}
            r_azimuth = get_metadata('azimuth')
            r_azimuth['data'] = np.ones(n_rays) * az
            r_elevation = {'data': scan_angles}
            r_sweep_number = {'data': [0]}
            r_fields = {'colocated_gates': get_metadata('colocated_gates')}
            r_fields['colocated_gates']['data'] = np.ma.ones((n_rays, 1),
                                                             dtype=int)

            r_radar = Radar(r_time, r_range, r_fields, None, 'rhi',
                            target_radar.latitude, target_radar.longitude,
                            target_radar.altitude,
                            r_sweep_number, None, None, None, None,
                            r_azimuth, r_elevation)

            # flag regions with colocated usable data in r_radar
            ##r_ind_invalid = r_radar.gate_altitude['data'] > max_altitude
            ##r_radar.fields['colocated_gates']['data'][r_ind_invalid] = 0

            # Find minimum and maximum azimuth in trajectory in this
            # time step.
            azmin = az-az_traj_tol
            azmax = az+az_traj_tol
            rmin = rr-rg_traj_tol
            rmax = rr+rg_traj_tol
            if azmin < 0:
                azmin = azmin+360
            if azmax > 360:
                azmax = azmax-360

            # flag regions with colocated usable data in radar_sel
            gate_coloc_radar_sel = intersection(
                radar_sel, r_radar, h_tol=alt_tol, latlon_tol=latlon_tol,
                vol_d_tol=None, vismin=None, hmin=None, hmax=max_altitude,
                rmin=rmin, rmax=rmax, elmin=None, elmax=None, azmin=azmin,
                azmax=azmax, visib_field=None,
                intersec_field='colocated_gates')
            radar_sel.add_field('colocated_gates', gate_coloc_radar_sel,
                                replace_existing=True)

            (colgates, r_radar_colg) = colocated_gates(r_radar, radar_sel,
                                                       h_tol=alt_tol,
                                                       latlon_tol=latlon_tol)

            w_ind = np.where(r_radar_colg['data'] != 0)[0]
            w_vec = tadict['weightvec'][w_ind]

            rdata = radar_sel.fields[field_name]['data']
            values = rdata[colgates['rad2_ray_ind'], colgates['rad2_rng_ind']]

        if use_nans:
            values_ma = np.ma.getmaskarray(values)
            values[values_ma] = nan_value

        try:
            (avg, qvals, nvals_valid) = quantiles_weighted(
                values,
                weight_vector=w_vec,
                quantiles=tadict['quantiles'],
                weight_threshold=weight_threshold,
                data_is_log=data_is_log)
        except Exception as ee:
            warn(str(ee))
            continue

        ts.add_timesample(trajectory.time_vector[tind],
                          (np.concatenate([[avg], qvals, [nvals_valid]])))

        # end loop over traj samples within period

    tadict['last_task_start_dt'] = dt_task_start
    tadict['radar_old2'] = tadict['radar_old']
    tadict['radar_old'] = radar

    # Collect garbage
    gc.collect()


    return None, None


def _get_closests_bin(az, el, rr, tt, radar, tdict):
    """
    """

    daz = np.abs(radar.azimuth['data'] - az)
    dele = np.abs(radar.elevation['data'] - el)
    dangle = np.sqrt(daz**2 + dele**2)
    ray_ind = np.argmin(dangle)
    dt = tt - radar.time['data'][ray_ind]

    if tdict['radar_old'] is not None:
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

            if np.abs(dt_old) < np.abs(dt_old2):
                radar_sel = rad_old
                ray_sel = ray_ind_old
            else:
                radar_sel = rad_old2
                ray_sel = ray_ind_old2
        else:
            if np.abs(dt) < np.abs(dt_old):
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

    el_vec_rnd = np.sort(np.unique(
        radar_sel.elevation['data'].round(decimals=1)))
    az_vec_rnd = np.sort(np.unique(
        radar_sel.azimuth['data'].round(decimals=1)))

    return (radar_sel, ray_sel, rr_ind, el_vec_rnd, az_vec_rnd)


def _sample_out_of_sector(az, el, rr, radar_sel, ray_sel, rr_ind,
                          el_vec_rnd, az_vec_rnd):
    """
    """

    # Check if sample is within sector
    rr_min = radar_sel.range['data'][rr_ind]
    range_res = radar_sel.range['data'][1] - radar_sel.range['data'][0]
    if np.abs(rr_min - rr) > (2*range_res):
        # print("INFO: Trajectory sample out of range")
        return True

    if (((az_vec_rnd[0] - az) > 3.0) or ((az - az_vec_rnd[-1]) > 3.0)):
        # print("INFO: Trajectory sample out of azimuth angle")
        return True

    if (((el_vec_rnd[0] - el) > 3.0) or ((el - el_vec_rnd[-1]) > 3.0)):
        # print("INFO: Trajectory sample out of elevation angle")
        return True

    return False


class Target_Radar:
    """

    """
    def __init__(self, latitude, longitude, altitude):
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
