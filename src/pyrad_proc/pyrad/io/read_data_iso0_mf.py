"""
pyrad.io.read_data_iso0_mf
==========================

Functions for reading HZT data

.. autosummary::
    :toctree: generated/

    iso2radar_data
    get_iso0_ref
    read_iso0_mf_data

"""
from warnings import warn
import datetime
import numpy as np

from pyart.core import geographic_to_cartesian_aeqd
from pyart.config import get_metadata


def iso2radar_data(radar, iso0_data, time_info, iso0_statistic='avg_by_dist',
                   field_name='height_over_iso0'):
    """
    get the iso0 value corresponding to each radar gate

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    iso0_data : dict
        dictionary containing the iso0 data and metadata from the model
    time_info : datetime object
        reference time of current radar volume
    iso0_statistic : str
        The statistic used to weight the iso0 points. Can be avg_by_dist,
        avg, min, max
    field_name : str
        name of HZT fields to convert (default height_over_iso0)

    Returns
    -------
    field_dict : dict
        dictionary containing the iso0 data in radar coordinates and the
        metadata

    """
    # get the relevant time indices for interpolation in time
    time_index = np.argmin(abs(iso0_data['fcst_time']-time_info))
    if time_info > iso0_data['fcst_time'][time_index]:
        time_index_future = time_index+1
        time_index_past = time_index
    else:
        time_index_future = time_index
        time_index_past = time_index-1

    # interpolate the iso0 ref in time
    if time_index_past == -1:
        # no interpolation: use data from time_index_future
        iso0_ref = get_iso0_ref(
            radar, iso0_data, time_index_future, statistic=iso0_statistic)
    elif time_index_future > iso0_data['fcst_time'].size:
        # no interpolation: use data from time_index_past
        iso0_ref = get_iso0_ref(
            radar, iso0_data, time_index_past, statistic=iso0_statistic)
    else:
        # interpolate between two time steps
        iso0_ref_past = get_iso0_ref(
            radar, iso0_data, time_index_past, statistic=iso0_statistic)
        iso0_ref_future = get_iso0_ref(
            radar, iso0_data, time_index_future, statistic=iso0_statistic)

        # put time in seconds from past forecast
        time_info_s = (
            time_info-iso0_data['fcst_time'][time_index_past]).total_seconds()
        fcst_time_s = (
            iso0_data['fcst_time'][time_index_future] -
            iso0_data['fcst_time'][time_index_past]).total_seconds()
        iso0_ref = np.interp(
            time_info_s, [0, fcst_time_s], [iso0_ref_past, iso0_ref_future])

    print('iso0_ref:', iso0_ref)

    # put field
    field_dict = get_metadata(field_name)
    field_dict['data'] = radar.gate_altitude['data']-iso0_ref

    return field_dict


def get_iso0_ref(radar, iso0_data, time_index, statistic='avg_by_dist'):
    """
    gets the is0 that has to be taken as reference. The iso0 is an average of
    the iso0 points available weighted by the distance to the radar. If no
    points are available is assumed that the iso0 is at ground level and the
    altitude of the radar is taken as reference.

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    iso0_data : dict
        dictionary containing the is0 data and metadata
    time_index : int
        the time index
    statistic : str
        The statistic used to weight the iso0 points. Can be avg_by_dist,
        avg, min, max

    Returns
    -------
    iso0_ref : float
        reference iso0 [m]

    """
    npoints_iso0 = iso0_data['npoints_iso0'][time_index]
    if npoints_iso0 == 0:
        # iso 0 below ground level: we take the radar altitude as reference
        warn('iso-0 below ground level. '
             'The radar altitude will be used as reference')
        return radar.altitude['data'][0]

    if statistic == 'avg':
        return np.mean(iso0_data['iso0_points'][time_index])
    if statistic == 'min':
        return np.min(iso0_data['iso0_points'][time_index])
    if statistic == 'max':
        return np.max(iso0_data['iso0_points'][time_index])

    if statistic != 'avg_by_dist':
        warn('unknown statistic '+statistic +
             '. Default avg_by_dist will be applied')
    # mean weighted by the distance from the iso0 point to the radar
    x_iso, y_iso = geographic_to_cartesian_aeqd(
        iso0_data['lon_points'][time_index],
        iso0_data['lat_points'][time_index], radar.longitude['data'][0],
        radar.latitude['data'][0])
    dist_iso = np.sqrt(x_iso*x_iso+y_iso*y_iso)
    return (
        np.sum(iso0_data['iso0_points'][time_index]*dist_iso) /
        np.sum(dist_iso))


def read_iso0_mf_data(fname):
    """
    Reads iso-0 degree data from a MF NWP in text format

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    iso0_data : dictionary
        dictionary with the data and metadata

    """
    try:
        with open(fname, 'r') as iso0file:

            run_time = []
            fcst_time = []
            npoints_iso0 = []
            lon_points = []
            lat_points = []
            iso0_points = []

            # skip the header
            next(iso0file)

            # read file contents
            fileend = 0
            while fileend == 0:
                # remove leading and trailing whitespace
                line = iso0file.readline()
                if line:
                    line = line.strip()

                    # ignore white lines
                    if not line:
                        continue

                    # ignore comments
                    if line.startswith('#'):
                        continue

                    line = line.partition('#')[0]  # Remove comments
                    line = line.strip()

                    vals = line.split()

                    if vals[0] == 'NIVEAU':
                        run_time_aux = datetime.datetime.strptime(
                            vals[6], '%Y%m%d%H%M%S')
                        npoinst_iso0_aux = int(vals[8])
                        run_time.append(run_time_aux)
                        fcst_time.append(
                            run_time_aux +
                            datetime.timedelta(hours=float(vals[4])))
                        npoints_iso0.append(npoinst_iso0_aux)

                        lon_points_aux = []
                        lat_points_aux = []
                        iso0_points_aux = []
                        for _ in range(npoinst_iso0_aux):
                            line = iso0file.readline()
                            line = line.strip()
                            vals2 = line.split()

                            if vals2[0] == 'LONGITUDE':
                                lon_points_aux.append(float(vals2[1])/1000.)
                                lat_points_aux.append(float(vals2[3])/1000.)
                                iso0_points_aux.append(float(vals2[5]))
                        lon_points.append(lon_points_aux)
                        lat_points.append(lat_points_aux)
                        iso0_points.append(iso0_points_aux)

                else:
                    fileend = 1

            iso0file.close()

    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None

    iso0_data = {
        'run_time': np.array(run_time),
        'fcst_time': np.array(fcst_time),
        'npoints_iso0': np.array(npoints_iso0),
        'lon_points': lon_points,
        'lat_points': lat_points,
        'iso0_points': iso0_points
    }

    return iso0_data
