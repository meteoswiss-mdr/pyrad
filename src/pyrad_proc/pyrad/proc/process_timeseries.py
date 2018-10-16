"""
pyrad.proc.process_timeseries
=============================

Functions to obtain time series of radar data

.. autosummary::
    :toctree: generated/

    process_point_measurement
    process_qvp
    process_rqvp
    process_evp
    process_svp
    process_time_height
"""

from copy import deepcopy
from warnings import warn
import numpy as np
from scipy.interpolate import interp1d
from netCDF4 import num2date

import pyart

from ..io.io_aux import get_datatype_fields, get_fieldname_pyart
from ..util.radar_utils import find_nearest_gate, find_neighbour_gates
from ..util.radar_utils import project_to_vertical, compute_directional_stats
from ..util.radar_utils import get_target_elevations


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
        radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
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


def process_qvp(procstatus, dscfg, radar_list=None):
    """
    Computes quasi vertical profiles, by averaging over height levels
    PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        angle : int or float
            If the radar object contains a PPI volume, the sweep number to
            use, if it contains an RHI volume the elevation angle.
            Default 0.
        ang_tol : float
            If the radar object contains an RHI volume, the tolerance in the
            elevation angle for the conversion into PPI
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 50
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to accept average. Default 30.
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Ryzhkov A., Zhang P., Reeves H., Kumjian M., Tschallener T., Trömel S.,
    Simmer C. 2016: Quasi-Vertical Profiles: A New Way to Look at Polarimetric
    Radar Data. JTECH vol. 33 pp 551-562

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
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
        angle = dscfg.get('angle', 0)
        ang_tol = dscfg.get('ang_tol', 1.)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 50.)
        avg_type = dscfg.get('avg_type', 'mean')
        nvalid_min = dscfg.get('nvalid_min', 30)
        interp_kind = dscfg.get('interp_kind', 'none')

        if avg_type != 'mean' and avg_type != 'median':
            warn('Unsuported statistics '+avg_type)
            return None, None

        radar_aux = deepcopy(radar)
        # transform radar into ppi over the required elevation
        if radar_aux.scan_type == 'rhi':
            radar_aux = pyart.util.cross_section_rhi(
                radar_aux, [angle], el_tol=ang_tol)
        elif radar_aux.scan_type == 'ppi':
            radar_aux = radar_aux.extract_sweeps([int(angle)])
        else:
            warn('Error: unsupported scan type.')
            return None, None

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
            global_dict.update({'radar_out': qvp_aux})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1

        # modify metadata
        qvp = dscfg['global_data']['radar_out']

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
        values, _ = compute_directional_stats(
            radar_aux.fields[field_name]['data'], avg_type=avg_type,
            nvalid_min=nvalid_min, axis=0)

        # Project to vertical grid:
        qvp_data = project_to_vertical(
            values, radar_aux.gate_altitude['data'][0, :], qvp.range['data'],
            interp_kind=interp_kind)

        # Put data in radar object
        if np.size(qvp.fields[field_name]['data']) == 0:
            qvp.fields[field_name]['data'] = qvp_data.reshape(1, qvp.ngates)
        else:
            qvp.fields[field_name]['data'] = np.ma.concatenate(
                (qvp.fields[field_name]['data'],
                 qvp_data.reshape(1, qvp.ngates)))

        dscfg['global_data']['radar_out'] = qvp

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        qvp = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad


def process_rqvp(procstatus, dscfg, radar_list=None):
    """
    Computes range defined quasi vertical profiles, by averaging over height
    levels PPI data.

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
            The height resolution [m]. Default 2.
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to accept average. Default 30.
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'nearest'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation
        rmax : float
            ground range up to which the data is intended for use [m].
            Default 50000.
        weight_power : float
            Power p of the weighting function 1/abs(grng-(rmax-1))**p given to
            the data outside the desired range. -1 will set the weight to 0.
            Default 2.

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the QVP and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    Reference
    ---------
    Ryzhkov A., Zhang P., Reeves H., Kumjian M., Tschallener T., Trömel S.,
    Simmer C. 2016: Quasi-Vertical Profiles: A New Way to Look at Polarimetric
    Radar Data. JTECH vol. 33 pp 551-562

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
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
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 2.)
        avg_type = dscfg.get('avg_type', 'mean')
        nvalid_min = dscfg.get('nvalid_min', 30)
        interp_kind = dscfg.get('interp_kind', 'nearest')
        rmax = dscfg.get('rmax', 50000.)/1000.  # [Km]
        weight_power = dscfg.get('weight_power', 2.)

        if avg_type != 'mean' and avg_type != 'median':
            warn('Unsuported statistics '+avg_type)
            return None, None

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

        radar_aux = radar_ppi.extract_sweeps([0])

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
            qvp_aux.scan_type = 'rqvp'
            qvp_aux.sweep_mode['data'] = np.array(['qvp'])
            qvp_aux.sweep_start_ray_index['data'] = np.array(
                [0], dtype='int32')
            qvp_aux.fixed_angle['data'] = np.array([90.], dtype='float64')

            # ray dependent radar objects parameters
            qvp_aux.sweep_end_ray_index['data'] = np.array([-1], dtype='int32')
            qvp_aux.rays_per_sweep = np.array([0], dtype='int32')
            qvp_aux.azimuth['data'] = np.array([], dtype='float64')
            qvp_aux.elevation['data'] = np.array([], dtype='float64')
            qvp_aux.nrays = 0

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_out': qvp_aux})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1

        # modify metadata
        qvp = dscfg['global_data']['radar_out']

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

        # compute rQVP data
        val_interp = np.ma.masked_all((radar_ppi.nsweeps, qvp.ngates))
        grng_interp = np.ma.masked_all((radar_ppi.nsweeps, qvp.ngates))
        for sweep in range(radar_ppi.nsweeps):
            radar_aux = deepcopy(radar_ppi)
            radar_aux = radar_aux.extract_sweeps([sweep])

            # Compute QVP for this sweep
            values, _ = compute_directional_stats(
                radar_aux.fields[field_name]['data'], avg_type=avg_type,
                nvalid_min=nvalid_min, axis=0)

            height = radar_aux.gate_altitude['data'][0, :]

            # Project to grid
            val_interp[sweep, :] = project_to_vertical(
                values, height, qvp.range['data'], interp_kind=interp_kind)

            # compute ground range [Km]
            grng = np.sqrt(
                np.power(radar_aux.gate_x['data'][0, :], 2.) +
                np.power(radar_aux.gate_y['data'][0, :], 2.))/1000.

            # Project ground range to grid
            f = interp1d(
                height, grng, kind=interp_kind, bounds_error=False,
                fill_value='extrapolate')
            grng_interp[sweep, :] = f(qvp.range['data'])

        # Compute weight
        weight = np.ma.abs(grng_interp-(rmax-1.))
        weight[grng_interp <= rmax-1.] = 1./np.power(
            weight[grng_interp <= rmax-1.], 0.)

        if weight_power == -1:
            weight[grng_interp > rmax-1.] = 0.
        else:
            weight[grng_interp > rmax-1.] = 1./np.power(
                weight[grng_interp > rmax-1.], weight_power)

        # mask weights where there is no data
        mask = np.ma.getmaskarray(val_interp)
        weight = np.ma.masked_where(mask, weight)

        # Weighted average
        qvp_data = (
            np.ma.sum(val_interp*weight, axis=0)/np.ma.sum(weight, axis=0))

        if np.size(qvp.fields[field_name]['data']) == 0:
            qvp.fields[field_name]['data'] = qvp_data.reshape(1, qvp.ngates)
        else:
            qvp.fields[field_name]['data'] = np.ma.concatenate(
                (qvp.fields[field_name]['data'],
                 qvp_data.reshape(1, qvp.ngates)))

        dscfg['global_data']['radar_out'] = qvp

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        qvp = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': qvp})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad


def process_evp(procstatus, dscfg, radar_list=None):
    """
    Computes enhanced vertical profiles, by averaging over height levels
    PPI data.

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
        delta_rng, delta_azi : float
            maximum range distance [m] and azimuth distance [degree] from the
            central point of the evp containing data to average. Default 5000.
            and 10.
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 250.
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to consider the data valid when
            performing the averaging. Default 1
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'.
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the EVP and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
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
        lon = dscfg['lon']
        lat = dscfg['lat']
        latlon_tol = dscfg.get('latlon_tol', 0.0005)
        delta_rng = dscfg.get('delta_rng', 15000.)
        delta_azi = dscfg.get('delta_azi', 10.)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 250.)
        avg_type = dscfg.get('avg_type', 'mean')
        nvalid_min = dscfg.get('nvalid_min', 1)
        interp_kind = dscfg.get('interp_kind', 'none')
        if avg_type != 'mean' and avg_type != 'median':
            warn('Unsuported statistics '+avg_type)
            return None, None

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

        radar_aux = radar_ppi.extract_sweeps([0])

        # initialize dataset
        if dscfg['initialized'] == 0:
            evp_aux = deepcopy(radar_aux)
            # prepare space for field
            evp_aux.fields = dict()
            evp_aux.add_field(
                field_name, deepcopy(radar_aux.fields[field_name]))
            evp_aux.fields[field_name]['data'] = np.array([], dtype='float64')

            # fixed radar objects parameters
            evp_aux.range['data'] = np.arange(hmax/hres)*hres+hres/2.
            evp_aux.ngates = len(evp_aux.range['data'])

            evp_aux.time['units'] = pyart.io.make_time_unit_str(
                dscfg['timeinfo'])
            evp_aux.time['data'] = np.array([], dtype='float64')
            evp_aux.scan_type = 'evp'
            evp_aux.sweep_mode['data'] = np.array(['evp'])
            evp_aux.sweep_start_ray_index['data'] = np.array(
                [0], dtype='int32')
            evp_aux.fixed_angle['data'] = np.array([90.], dtype='float64')

            # ray dependent radar objects parameters
            evp_aux.sweep_end_ray_index['data'] = np.array([-1], dtype='int32')
            evp_aux.rays_per_sweep = np.array([0], dtype='int32')
            evp_aux.azimuth['data'] = np.array([], dtype='float64')
            evp_aux.elevation['data'] = np.array([], dtype='float64')
            evp_aux.nrays = 0

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_out': evp_aux})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1

        # modify metadata
        evp = dscfg['global_data']['radar_out']

        start_time = num2date(0, evp.time['units'], evp.time['calendar'])
        evp.time['data'] = np.append(
            evp.time['data'], (dscfg['timeinfo'] - start_time).total_seconds())
        evp.sweep_end_ray_index['data'][0] += 1
        evp.rays_per_sweep[0] += 1
        evp.nrays += 1

        evp.azimuth['data'] = np.ones((evp.nrays, ), dtype='float64')*0.
        evp.elevation['data'] = np.ones((evp.nrays, ), dtype='float64')*90.

        evp.gate_longitude['data'] = (
            np.ones((evp.nrays, evp.ngates), dtype='float64') * lon)
        evp.gate_latitude['data'] = (
            np.ones((evp.nrays, evp.ngates), dtype='float64') * lat)
        evp.gate_altitude['data'] = np.broadcast_to(
            evp.range['data'], (evp.nrays, evp.ngates))

        values = np.ma.array([], dtype=float)
        height = np.array([], dtype=float)
        for sweep in range(radar_ppi.nsweeps):
            radar_aux = deepcopy(radar_ppi)
            radar_aux = radar_aux.extract_sweeps([sweep])

            # find nearest gate to lat lon point
            ind_ray, _, azi, rng = find_nearest_gate(
                radar_aux, lat, lon, latlon_tol=latlon_tol)

            if ind_ray is None:
                continue

            # find neighbouring gates to be selected
            inds_ray, inds_rng = find_neighbour_gates(
                radar_aux, azi, rng, delta_azi=delta_azi, delta_rng=delta_rng)

            # keep only data we are interested in
            field = radar_aux.fields[field_name]['data'][:, inds_rng]
            field = field[inds_ray, :]

            vals, _ = compute_directional_stats(
                field, avg_type=avg_type, nvalid_min=nvalid_min, axis=0)
            values = np.ma.append(values, vals)

            height = np.append(
                height,
                radar_aux.gate_altitude['data'][ind_ray, inds_rng])

        # Project to vertical grid:
        evp_data = project_to_vertical(
            values, height, evp.range['data'], interp_kind=interp_kind)

        # Put data in radar object
        if np.size(evp.fields[field_name]['data']) == 0:
            evp.fields[field_name]['data'] = evp_data.reshape(1, evp.ngates)
        else:
            evp.fields[field_name]['data'] = np.ma.concatenate(
                (evp.fields[field_name]['data'],
                 evp_data.reshape(1, evp.ngates)))

        dscfg['global_data']['radar_out'] = evp

        new_dataset = dict()
        new_dataset.update({'radar_out': evp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        evp = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': evp})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad


def process_svp(procstatus, dscfg, radar_list=None):
    """
    Computes slanted vertical profiles, by averaging over height levels
    PPI data.

    Parameters
    ----------
    procstatus : int
        Processing status: 0 initializing, 1 processing volume,
        2 post-processing
    dscfg : dictionary of dictionaries
        data set configuration. Accepted Configuration Keywords::

        datatype : string. Dataset keyword
            The data type where we want to extract the point measurement
        angle : int or float
            If the radar object contains a PPI volume, the sweep number to
            use, if it contains an RHI volume the elevation angle.
            Default 0.
        ang_tol : float
            If the radar object contains an RHI volume, the tolerance in the
            elevation angle for the conversion into PPI
        lat, lon : float
            latitude and longitude of the point of interest [deg]
        latlon_tol : float
            tolerance in latitude and longitude in deg. Default 0.0005
        delta_rng, delta_azi : float
            maximum range distance [m] and azimuth distance [degree] from the
            central point of the svp containing data to average. Default 5000.
            and 10.
        hmax : float
            The maximum height to plot [m]. Default 10000.
        hres : float
            The height resolution [m]. Default 250.
        avg_type : str
            The type of averaging to perform. Can be either "mean" or "median"
            Default "mean"
        nvalid_min : int
            Minimum number of valid points to consider the data valid when
            performing the averaging. Default 1
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

    radar_list : list of Radar objects
        Optional. list of radar objects

    Returns
    -------
    new_dataset : dict
        dictionary containing the svp and a keyboard stating whether the
        processing has finished or not.
    ind_rad : int
        radar index

    """
    if procstatus == 0:
        return None, None

    if procstatus == 1:
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
        angle = dscfg.get('angle', 0)
        ang_tol = dscfg.get('ang_tol', 1.)
        lon = dscfg['lon']
        lat = dscfg['lat']
        latlon_tol = dscfg.get('latlon_tol', 0.0005)
        delta_rng = dscfg.get('delta_rng', 30000.)
        delta_azi = dscfg.get('delta_azi', 10.)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 250.)
        avg_type = dscfg.get('avg_type', 'mean')
        nvalid_min = dscfg.get('nvalid_min', 1)
        interp_kind = dscfg.get('interp_kind', 'none')
        if avg_type != 'mean' and avg_type != 'median':
            warn('Unsuported statistics '+avg_type)
            return None, None

        radar_aux = deepcopy(radar)
        # transform radar into ppi over the required elevation
        if radar_aux.scan_type == 'rhi':
            radar_aux = pyart.util.cross_section_rhi(
                radar_aux, [angle], el_tol=ang_tol)
        elif radar_aux.scan_type == 'ppi':
            radar_aux = radar_aux.extract_sweeps([int(angle)])
        else:
            warn('Error: unsupported scan type.')
            return None, None

        # initialize dataset
        if dscfg['initialized'] == 0:
            svp_aux = deepcopy(radar_aux)
            # prepare space for field
            svp_aux.fields = dict()
            svp_aux.add_field(
                field_name, deepcopy(radar_aux.fields[field_name]))
            svp_aux.fields[field_name]['data'] = np.array([], dtype='float64')

            # fixed radar objects parameters
            svp_aux.range['data'] = np.arange(hmax/hres)*hres+hres/2.
            svp_aux.ngates = len(svp_aux.range['data'])

            svp_aux.time['units'] = pyart.io.make_time_unit_str(
                dscfg['timeinfo'])
            svp_aux.time['data'] = np.array([], dtype='float64')
            svp_aux.scan_type = 'svp'
            svp_aux.sweep_mode['data'] = np.array(['svp'])
            svp_aux.sweep_start_ray_index['data'] = np.array(
                [0], dtype='int32')

            # ray dependent radar objects parameters
            svp_aux.sweep_end_ray_index['data'] = np.array([-1], dtype='int32')
            svp_aux.rays_per_sweep = np.array([0], dtype='int32')
            svp_aux.azimuth['data'] = np.array([], dtype='float64')
            svp_aux.elevation['data'] = np.array([], dtype='float64')
            svp_aux.nrays = 0

            global_dict = dict()
            global_dict.update({'start_time': dscfg['timeinfo']})
            global_dict.update({'radar_out': svp_aux})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1

        # modify metadata
        svp = dscfg['global_data']['radar_out']

        start_time = num2date(0, svp.time['units'], svp.time['calendar'])
        svp.time['data'] = np.append(
            svp.time['data'], (dscfg['timeinfo'] - start_time).total_seconds())
        svp.sweep_end_ray_index['data'][0] += 1
        svp.rays_per_sweep[0] += 1
        svp.nrays += 1

        svp.azimuth['data'] = np.ones((svp.nrays, ), dtype='float64')*0.
        svp.elevation['data'] = np.ones((svp.nrays, ), dtype='float64')*90.

        svp.gate_longitude['data'] = (
            np.ones((svp.nrays, svp.ngates), dtype='float64') * lon)
        svp.gate_latitude['data'] = (
            np.ones((svp.nrays, svp.ngates), dtype='float64') * lat)
        svp.gate_altitude['data'] = np.broadcast_to(
            svp.range['data'], (svp.nrays, svp.ngates))

        # find nearest gate to lat lon point
        ind_ray, _, azi, rng = find_nearest_gate(
            radar_aux, lat, lon, latlon_tol=latlon_tol)

        if ind_ray is None:
            values = np.ma.array([], dtype=float)
            height = np.array([], dtype=float)
        else:
            # find neighbouring gates to be selected
            inds_ray, inds_rng = find_neighbour_gates(
                radar_aux, azi, rng, delta_azi=delta_azi, delta_rng=delta_rng)

            # keep only data we are interested in
            field = radar_aux.fields[field_name]['data'][:, inds_rng]
            field = field[inds_ray, :]

            # compute values
            values, _ = compute_directional_stats(
                field, avg_type=avg_type, nvalid_min=nvalid_min, axis=0)

            height = radar_aux.gate_altitude['data'][ind_ray, inds_rng]

        # Project to vertical grid:
        svp_data = project_to_vertical(
            values, height, svp.range['data'], interp_kind=interp_kind)

        # Put data in radar object
        if np.size(svp.fields[field_name]['data']) == 0:
            svp.fields[field_name]['data'] = svp_data.reshape(1, svp.ngates)
        else:
            svp.fields[field_name]['data'] = np.ma.concatenate(
                (svp.fields[field_name]['data'],
                 svp_data.reshape(1, svp.ngates)))

        dscfg['global_data']['radar_out'] = svp

        new_dataset = dict()
        new_dataset.update({'radar_out': svp})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        svp = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': svp})
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
        interp_kind : str
            type of interpolation when projecting to vertical grid: 'none',
            or 'nearest', etc. Default 'none'
            'none' will select from all data points within the regular grid
            height bin the closest to the center of the bin.
            'nearest' will select the closest data point to the center of the
            height bin regardless if it is within the height bin or not.
            Data points can be masked values
            If another type of interpolation is selected masked values will be
            eliminated from the data points before the interpolation

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
        lon = dscfg['lon']
        lat = dscfg['lat']
        latlon_tol = dscfg.get('latlon_tol', 0.0005)
        hmax = dscfg.get('hmax', 10000.)
        hres = dscfg.get('hres', 50.)
        interp_kind = dscfg.get('interp_kind', 'none')

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
            global_dict.update({'radar_out': th_aux})
            dscfg['global_data'] = global_dict
            dscfg['initialized'] = 1

        # modify metadata
        th = dscfg['global_data']['radar_out']

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

        values = np.ma.array([], dtype=float)
        height = np.array([], dtype=float)
        for sweep in range(radar_ppi.nsweeps):
            radar_aux = deepcopy(radar_ppi.extract_sweeps([sweep]))

            # find nearest gate to lat lon point
            ind_ray, ind_rng, _, _ = find_nearest_gate(
                radar_aux, lat, lon, latlon_tol=latlon_tol)

            if ind_ray is None:
                continue

            values = np.ma.append(
                values,
                radar_aux.fields[field_name]['data'][ind_ray, ind_rng])
            height = np.append(
                height,
                radar_aux.gate_altitude['data'][ind_ray, ind_rng])

        # Project to vertical grid:
        th_data = project_to_vertical(
            values, height, th.range['data'], interp_kind=interp_kind)

        # Put data in radar object
        if np.size(th.fields[field_name]['data']) == 0:
            th.fields[field_name]['data'] = th_data.reshape(1, th.ngates)
        else:
            th.fields[field_name]['data'] = np.ma.concatenate(
                (th.fields[field_name]['data'],
                 th_data.reshape(1, th.ngates)))

        dscfg['global_data']['radar_out'] = th

        new_dataset = dict()
        new_dataset.update({'radar_out': th})
        new_dataset.update({'radar_type': 'temporal'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad

    if procstatus == 2:
        for datatypedescr in dscfg['datatype']:
            radarnr, _, datatype, _, _ = get_datatype_fields(datatypedescr)
            break

        ind_rad = int(radarnr[5:8])-1

        th = dscfg['global_data']['radar_out']

        new_dataset = dict()
        new_dataset.update({'radar_out': th})
        new_dataset.update({'radar_type': 'final'})
        new_dataset.update({'start_time': dscfg['global_data']['start_time']})

        return new_dataset, ind_rad
