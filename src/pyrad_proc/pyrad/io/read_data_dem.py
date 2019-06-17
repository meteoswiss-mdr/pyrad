"""
pyrad.io.read_data_dem
========================

Functions for reading data derived from Digital Elevation Models (DEM)

.. autosummary::
    :toctree: generated/

    dem2radar_data
    read_idrisi_data
    read_idrisi_metadata
    _prepare_for_interpolation


"""

from warnings import warn
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from osgeo import gdal

from pyart.config import get_metadata
from pyart.aux_io import convert_data
from ..io.read_data_cosmo import _put_radar_in_swiss_coord

# from memory_profiler import profile

# import time


def dem2radar_data(radar, dem_data, slice_xy=True, field_name='visibility'):
    """
    get the DEM value corresponding to each radar gate using nearest
    neighbour interpolation

    Parameters
    ----------
    radar : Radar
        the radar object containing the information on the position of the
        radar gates
    dem_data : dict
        dictionary containing the DEM data
    slice_xy : boolean
        if true the horizontal plane of the DEM field is cut to the
        dimensions of the radar field
    field_names : str
        names of DEM fields to convert

    Returns
    -------
    dem_field : dict
        Dictionary with the DEM fields and metadata

    """
    # debugging
    # start_time = time.time()

    x_radar, y_radar, _ = _put_radar_in_swiss_coord(radar)

    (x_dem, y_dem, ind_xmin, ind_ymin, ind_xmax, ind_ymax) = (
        _prepare_for_interpolation(
            x_radar, y_radar, dem_data, slice_xy=slice_xy))

    if field_name not in dem_data:
        warn('DEM field '+field_name+' data not available')
        return None

    values = dem_data[field_name]['data'][
        ind_xmin:ind_xmax+1, ind_ymin:ind_ymax+1].flatten()
    # find interpolation function
    tree_options = {
        'compact_nodes': False,
        'balanced_tree': False
    }
    interp_func = NearestNDInterpolator(
        (x_dem, y_dem), values, tree_options=tree_options)

    del values

    # interpolate
    data_interp = interp_func((x_radar, y_radar))

    # put field
    field_dict = get_metadata(field_name)
    field_dict['data'] = data_interp.astype(float)

    del data_interp

    return field_dict


# @profile
def read_idrisi_data(fname, field_name, fill_value=-99.):
    """
    Reads DEM data from an IDRISI .rst file

    Parameters
    ----------
    fname : str
        name of the file to read
    field_name : str
        name of the readed variable
    fill_value : float
        The fill value

    Returns
    -------
    dem_data : dictionary
        dictionary with the data and metadata

    """
    # read the data
    try:
        raster = gdal.Open(fname)
        raster_array = raster.ReadAsArray()
        raster_array = np.ma.masked_equal(raster_array, fill_value)

        metadata = read_idrisi_metadata(fname)

        if metadata is None:
            return None

        field_dict = get_metadata(field_name)
        field_dict['data'] = np.transpose(raster_array)[:, ::-1]
        field_dict['units'] = metadata['value units']

        x = get_metadata('x')
        y = get_metadata('y')
        x['data'] = (
            np.arange(raster.RasterXSize)*metadata['resolution'] +
            metadata['resolution']/2.+metadata['min. X'])

        y['data'] = (
            np.arange(raster.RasterYSize)*metadata['resolution'] +
            metadata['resolution']/2.+metadata['min. Y'])

        dem_data = {
            'metadata': metadata,
            'x': x,
            'y': y,
            field_name: field_dict
        }

        return dem_data
    except EnvironmentError as ee:
        warn(str(ee))
        warn('Unable to read file '+fname)
        return None


def read_idrisi_metadata(fname):
    """
    Reads DEM metadata from a IDRISI .rdc file

    Parameters
    ----------
    fname : str
        name of the file to read

    Returns
    -------
    metadata : dictionary
        dictionary with the metadata

    """
    # read the data
    fname_rdc = fname.replace('.rst', '.rdc')

    try:
        metadata = dict()
        with open(fname_rdc, 'r', newline='') as txtfile:
            for line in txtfile:
                strs = line.split(':')
                metadata.update({
                    strs[0].strip(): convert_data(strs[1].strip())})

        return metadata
    except EnvironmentError:
        warn('Unable to read file '+fname_rdc)
        return None


def _prepare_for_interpolation(x_radar, y_radar, dem_coord, slice_xy=True):
    """
    prepares the DEM 2D volume for interpolation:
        1. if set slices the DEM data to the area
    covered by the radar
        2. creates the x, y grid for the interpolation

    Parameters
    ----------
    x_radar, y_radar : arrays
        The Swiss coordinates of the radar
    dem_coord : dict
        dictionary containing the DEM coordinates
    slice_xy : boolean
        if true the horizontal plane of the DEM field is cut to the
        dimensions of the radar field

    Returns
    -------
    x_dem, y_dem : 1D arrays
        arrays containing the flatten swiss coordinates of the DEM data in
        the area of interest
    ind_xmin, ind_ymin, ind_xmax, ind_ymax : ints
        the minimum and maximum indices of each dimension

    """
    nx_dem = len(dem_coord['x']['data'])
    ny_dem = len(dem_coord['y']['data'])

    if slice_xy:
        # get the D data within the radar range
        xmin = np.min(x_radar)
        xmax = np.max(x_radar)
        ymin = np.min(y_radar)
        ymax = np.max(y_radar)

        ind_xmin = np.where(dem_coord['x']['data'] < xmin)[0]
        if ind_xmin.size == 0:
            ind_xmin = 0
        else:
            ind_xmin = ind_xmin[-1]

        ind_xmax = np.where(dem_coord['x']['data'] > xmax)[0]
        if ind_xmax.size == 0:
            ind_xmax = nx_dem-1
        else:
            ind_xmax = ind_xmax[0]

        ind_ymin = np.where(dem_coord['y']['data'] < ymin)[0]
        if ind_ymin.size == 0:
            ind_ymin = 0
        else:
            ind_ymin = ind_ymin[-1]

        ind_ymax = np.where(dem_coord['y']['data'] > ymax)[0]
        if ind_ymax.size == 0:
            ind_ymax = ny_dem-1
        else:
            ind_ymax = ind_ymax[0]
    else:
        ind_xmin = 0
        ind_xmax = nx_dem-1
        ind_ymin = 0
        ind_ymax = ny_dem-1

    nx = ind_xmax-ind_xmin+1
    ny = ind_ymax-ind_ymin+1

    x_dem = dem_coord['x']['data'][ind_xmin:ind_xmax+1]
    y_dem = dem_coord['y']['data'][ind_ymin:ind_ymax+1]

    x_dem = (
        np.broadcast_to(x_dem.reshape(nx, 1), (nx, ny))).flatten()
    y_dem = (
        np.broadcast_to(y_dem.reshape(1, ny), (nx, ny))).flatten()

    return (x_dem, y_dem, ind_xmin, ind_ymin, ind_xmax, ind_ymax)
