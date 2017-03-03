"""
pyrad.io.read_data_cosmo
========================

Functions for reading COSMO data

.. autosummary::
    :toctree: generated/

    cosmo2radar_data
    cosmo2radar_coord
    get_cosmo_field
    read_cosmo_temp
    read_cosmo_coord
    _ncvar_to_dict
    _prepare_for_interpolation
    _put_radar_in_swiss_coord


"""

import netCDF4
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial import cKDTree

from pyart.core import wgs84_to_swissCH1903
from pyart.config import get_metadata, get_field_name

import time


def cosmo2radar_data(radar, cosmo_coord, cosmo_data, time_index=0,
                     slice_xy=True, slice_z=False, field_name=None):
    """
    get the COSMO value corresponding to each radar gate using nearest
    neighbour interpolation

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    cosmo_coord : dict
        dictionary containing the COSMO coordinates
    cosmo_data : dict
        dictionary containing the COSMO data
    time_index : int
        index of the forecasted data
    slice_xy : boolean
        if true the horizontal plane of the COSMO field is cut to the
        dimensions of the radar field
    slice_z : boolean
        if true the vertical plane of the COSMO field is cut to the dimensions
        of the radar field
    field_name : str
        name of COSMO parameter (default temperature)

    Returns
    -------
    cosmo_field : dict
        dictionary with the COSMO field and metadata

    """
    # debugging
    # start_time = time.time()

    # parse the field parameters
    if field_name is None:
        field_name = get_field_name('temperature')

    x_radar, y_radar, z_radar = _put_radar_in_swiss_coord(radar)

    (x_cosmo, y_cosmo, z_cosmo, nx, ny, nz, ind_xmin, ind_ymin, ind_zmin,
     values) = _prepare_for_interpolation(
        x_radar, y_radar, z_radar, cosmo_coord, slice_xy=slice_xy,
        slice_z=slice_z,
        cosmo_data=cosmo_data['variable']['data'][time_index, :, :, :])

    # find interpolation function
    interp_func = NearestNDInterpolator((z_cosmo, y_cosmo, x_cosmo), values)

    # debugging
    # print(" generating interpolation function takes %s seconds " %
    #      (time.time() - start_time))

    # interpolate
    # debugging
    # start_time = time.time()
    data_interp = interp_func((z_radar, y_radar, x_radar))
    # print(" interpolating takes %s seconds " % (time.time() - start_time))

    # put field
    cosmo_field = get_metadata(field_name)
    cosmo_field['data'] = data_interp

    return cosmo_field


def cosmo2radar_coord(radar, cosmo_coord, slice_xy=True, slice_z=False,
                      field_name=None):
    """
    Given the radar coordinates find the nearest COSMO model pixel

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    cosmo_coord : dict
        dictionary containing the COSMO coordinates
    slice_xy : boolean
        if true the horizontal plane of the COSMO field is cut to the
        dimensions of the radar field
    slice_z : boolean
        if true the vertical plane of the COSMO field is cut to the dimensions
        of the radar field

    Returns
    -------
    cosmo_ind_field : dict
        dictionary containing a field of COSMO indices and metadata

    """
    # debugging
    # start_time = time.time()

    # parse the field parameters
    if field_name is None:
        field_name = get_field_name('cosmo_index')

    x_radar, y_radar, z_radar = _put_radar_in_swiss_coord(radar)

    (x_cosmo, y_cosmo, z_cosmo, nx, ny, nz, ind_xmin, ind_ymin, ind_zmin,
     values) = _prepare_for_interpolation(
        x_radar, y_radar, z_radar, cosmo_coord, slice_xy=slice_xy,
        slice_z=slice_z, cosmo_data=None)

    tree = cKDTree(np.transpose((z_cosmo, y_cosmo, x_cosmo)))
    dd_vec, ind_vec = tree.query(np.transpose(
        (z_radar.flatten(), y_radar.flatten(), x_radar.flatten())), k=1)

    # put the index in the original cosmo coordinates
    nx_cosmo = len(cosmo_coord['x']['data'])
    ny_cosmo = len(cosmo_coord['y']['data'])
    nz_cosmo = len(cosmo_coord['z']['data'])

    ind_z = (ind_vec/(nx*ny)).astype(int)+ind_zmin
    ind_y = ((ind_vec-nx*ny*ind_z)/nx).astype(int)+ind_ymin
    ind_x = ((ind_vec-nx*ny*ind_z) % nx).astype(int)+ind_xmin
    ind_cosmo = (ind_x+nx_cosmo*ind_y+nx_cosmo*ny_cosmo*ind_z).astype(int)

    cosmo_ind_field = get_metadata(field_name)
    cosmo_ind_field['data'] = ind_cosmo.reshape(radar.nrays, radar.ngates)

    # debugging
    # print(" generating COSMO indices takes %s seconds " %
    #      (time.time() - start_time))

    return cosmo_ind_field


def get_cosmo_field(cosmo_data, cosmo_ind, time_index=0, field_name=None):
    """
    Get the COSMO data corresponding to each radar gate
    using a precomputed look up table of the nearest neighbour

    Parameters
    ----------
    cosmo_data : dict
        dictionary containing the COSMO data and metadata
    cosmo_limits : dict
        dictionary containing the x, y, z indices limiting the COSMO field
    cosmo_ind : dict
        dictionary containing a field of COSMO indices and metadata
    time_index : int
        index of the forecasted data
    field_name : str
        name of COSMO parameter (default temperature)

    Returns
    -------
    cosmo_field : dict
        dictionary with the COSMO field and metadata

    """
    # parse the field parameters
    if field_name is None:
        field_name = get_field_name('temperature')

    nrays, ngates = np.shape(cosmo_ind['data'])
    values = cosmo_data['variable']['data'][time_index, :, :, :].flatten()

    cosmo_field = get_metadata(field_name)
    cosmo_field['data'] = (
        values[cosmo_ind['data'].flatten()].reshape(nrays, ngates))

    return cosmo_field


def read_cosmo_temp(fname, celsius=False):
    """
    Reads COSMO temperature from a netcdf file

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    cosmo_temp : dictionary
        dictionary with the data and metadata

    """
    # read the data
    ncobj = netCDF4.Dataset(fname)
    ncvars = ncobj.variables

    # 4.1 Global attribute -> move to metadata dictionary
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])

    # 4.2 put variables in dictionary
    x_1 = _ncvar_to_dict(ncvars['x_1'])
    y_1 = _ncvar_to_dict(ncvars['y_1'])
    lon_1 = _ncvar_to_dict(ncvars['lon_1'])
    lat_1 = _ncvar_to_dict(ncvars['lat_1'])
    z_1 = _ncvar_to_dict(ncvars['z_1'])
    z_bnds_1 = _ncvar_to_dict(ncvars['z_bnds_1'])
    temp = _ncvar_to_dict(ncvars['T'])

    if celsius:
        temp['data'] -= 273.15
        temp['units'] = 'C'

    time = _ncvar_to_dict(ncvars['time'])

    # close object
    ncobj.close()

    cosmo_temp = {
        'metadata': metadata,
        'time': time,
        'x': x_1,
        'y': y_1,
        'z': z_1,
        'z_bnds': z_bnds_1,
        'lon': lon_1,
        'lat': lat_1,
        'variable': temp
    }

    return cosmo_temp


def read_cosmo_coord(fname):
    """
    Reads COSMO coordinates from a netcdf file

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    cosmo_coord : dictionary
        dictionary with the data and metadata

    """
    # read the data
    ncobj = netCDF4.Dataset(fname)
    ncvars = ncobj.variables

    # 4.1 Global attribute -> move to metadata dictionary
    metadata = dict([(k, getattr(ncobj, k)) for k in ncobj.ncattrs()])

    # 4.2 put variables in dictionary
    x_1 = _ncvar_to_dict(ncvars['x_1'])
    y_1 = _ncvar_to_dict(ncvars['y_1'])
    lon_1 = _ncvar_to_dict(ncvars['lon_1'])
    lat_1 = _ncvar_to_dict(ncvars['lat_1'])
    z_1 = _ncvar_to_dict(ncvars['z_1'])
    z_bnds_1 = _ncvar_to_dict(ncvars['z_bnds_1'])
    hfl = _ncvar_to_dict(ncvars['HFL'])
    hsurf = _ncvar_to_dict(ncvars['HSURF'])
    fr_land = _ncvar_to_dict(ncvars['FR_LAND'])

    # close object
    ncobj.close()

    cosmo_coord = {
        'metadata': metadata,
        'x': x_1,
        'y': y_1,
        'z': z_1,
        'z_bnds': z_bnds_1,
        'lon': lon_1,
        'lat': lat_1,
        'hfl': hfl,
        'hsurf': hsurf,
        'fr_land': fr_land,
    }

    return cosmo_coord


def _ncvar_to_dict(ncvar):
    """ Convert a NetCDF Dataset variable to a dictionary. """
    # copy all attributes
    d = dict((k, getattr(ncvar, k)) for k in ncvar.ncattrs())
    d.update({'data': ncvar[:]})
    if '_FillValue' in d:
        d['data'] = np.ma.asarray(d['data'])
        d['data'] = np.ma.masked_values(d['data'], float(d['_FillValue']))
    else:
        d['data'] = np.asarray(d['data'])

    return d


def _prepare_for_interpolation(x_radar, y_radar, z_radar, cosmo_coord,
                               slice_xy=True, slice_z=False, cosmo_data=None):
    """
    prepares the COSMO 3D volume for interpolation:
        1. if set slices the cosmo data to the area (or volume)
    covered by the radar
        2. creates the x, y, z grid for the interpolation

    Parameters
    ----------
    x_radar, y_radar, z_radar : arrays
        The Swiss coordinates of the radar
    cosmo_coord : dict
        dictionary containing the COSMO coordinates
    slice_xy : boolean
        if true the horizontal plane of the COSMO field is cut to the
        dimensions of the radar field
    slice_z : boolean
        if true the vertical plane of the COSMO field is cut to the dimensions
        of the radar field
    cosmo_data : array
        array containing the COSMO data for a particular time step

    Returns
    -------
    x_cosmo, y_cosmo, z_cosmo : 1D arrays
        arrays containing the flatten swiss coordinates of the COSMO data in
        the area of interest
    nx, ny, nz : ints
        the dimensions of the cut COSMO 3D volume
    values : 1D array
        array containing the flatten COSMO data in the area of interest

    """
    nx_cosmo = len(cosmo_coord['x']['data'])
    ny_cosmo = len(cosmo_coord['y']['data'])
    nz_cosmo = len(cosmo_coord['z']['data'])

    if slice_xy:
        # get the COSMO data within the radar range
        xmin = np.min(x_radar)
        xmax = np.max(x_radar)
        ymin = np.min(y_radar)
        ymax = np.max(y_radar)

        ind_xmin = np.where(cosmo_coord['x']['data'] < xmin)[0]
        if len(ind_xmin) == 0:
            ind_xmin = 0
        else:
            ind_xmin = ind_xmin[-1]

        ind_xmax = np.where(cosmo_coord['x']['data'] > xmax)[0]
        if len(ind_xmax) == 0:
            ind_xmax = nx_cosmo-1
        else:
            ind_xmax = ind_xmax[0]

        ind_ymin = np.where(cosmo_coord['y']['data'] < ymin)[0]
        if len(ind_ymin) == 0:
            ind_ymin = 0
        else:
            ind_ymin = ind_ymin[-1]

        ind_ymax = np.where(cosmo_coord['y']['data'] > ymax)[0]
        if len(ind_ymax) == 0:
            ind_ymax = ny_cosmo-1
        else:
            ind_ymax = ind_ymax[0]
    else:
        ind_xmin = 0
        ind_xmax = nx_cosmo-1
        ind_ymin = 0
        ind_ymax = ny_cosmo-1

    if slice_z:
        zmin = np.min(z_radar)
        zmax = np.max(z_radar)

        ind_z, ind_y, ind_x = np.where(cosmo_coord['hfl']['data'] < zmin)
        if len(ind_z) == 0:
            ind_zmin = 0
        else:
            ind_zmin = np.min(ind_z)
        ind_z, ind_y, ind_x = np.where(cosmo_coord['hfl']['data'] > zmax)
        if len(ind_z) == 0:
            ind_zmax = nz_cosmo-1
        else:
            ind_zmax = np.max(ind_z)
    else:
        ind_zmin = 0
        ind_zmax = nz_cosmo-1

    nx = ind_xmax-ind_xmin+1
    ny = ind_ymax-ind_ymin+1
    nz = ind_zmax-ind_zmin+1

    x_cosmo = cosmo_coord['x']['data'][ind_xmin:ind_xmax+1]
    y_cosmo = cosmo_coord['y']['data'][ind_ymin:ind_ymax+1]
    z_cosmo = cosmo_coord['hfl']['data'][
        ind_zmin:ind_zmax+1, ind_ymin:ind_ymax+1, ind_xmin:ind_xmax+1]

    x_cosmo = (
        np.broadcast_to(x_cosmo.reshape(1, 1, nx), (nz, ny, nx))).flatten()
    y_cosmo = (
        np.broadcast_to(y_cosmo.reshape(1, ny, 1), (nz, ny, nx))).flatten()
    z_cosmo = z_cosmo.flatten()

    values = None
    if cosmo_data is not None:
        values = cosmo_data[
            ind_zmin:ind_zmax+1, ind_ymin:ind_ymax+1, ind_xmin:ind_xmax+1]
        values = values.flatten()

    return (x_cosmo, y_cosmo, z_cosmo, nx, ny, nz, ind_xmin, ind_ymin,
            ind_zmin, values)


def _put_radar_in_swiss_coord(radar):
    """
    puts the Cartesian grid of the radar coordinates in Swiss coordinates

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates

    Returns
    -------
    x_radar, y_radar, z_radar : 2D arrays
        arrays containing swiss coordinates of the radar

    """
    x0, y0, z0 = wgs84_to_swissCH1903(
        radar.longitude['data'][0], radar.latitude['data'][0],
        radar.altitude['data'][0], no_altitude_transform=True)

    x_radar = radar.gate_x['data']+x0
    y_radar = radar.gate_y['data']+y0
    z_radar = radar.gate_altitude['data']

    return x_radar, y_radar, z_radar
