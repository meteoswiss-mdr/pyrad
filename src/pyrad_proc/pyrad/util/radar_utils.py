"""
pyrad.util.radar_utils
======================

Miscellaneous functions dealing with radar data

.. autosummary::
    :toctree: generated/

    create_sun_hits_field
    create_sun_retrieval_field
    compute_quantiles_sweep
    compute_histogram_sweep

"""
from warnings import warn

import numpy as np

import pyart


def create_sun_hits_field(rad_el, rad_az, sun_el, sun_az, data, imgcfg):
    """
    creates a sun hits field from the position and power of the sun hits

    Parameters
    ----------
    rad_el, rad_az, sun_el, sun_az : ndarray 1D
        azimuth and elevation of the radar and the sun respectively in degree
    data : masked ndarray 1D
        the sun hit data
    imgcfg: dict
        a dictionary specifying the ranges and resolution of the field to
        create

    Returns
    -------
    field : masked ndarray 2D
        the sun hit field

    """
    if len(data.compressed()) == 0:
        warn('No valid sun hits to plot.')
        return None

    azmin = imgcfg['azmin']
    azmax = imgcfg['azmax']
    elmin = imgcfg['elmin']
    elmax = imgcfg['elmax']
    azres = imgcfg['azres']
    elres = imgcfg['elres']

    mask = np.ma.getmaskarray(data)
    rad_el = rad_el[~mask]
    rad_az = rad_az[~mask]
    sun_el = sun_el[~mask]
    sun_az = sun_az[~mask]
    data = data[~mask]

    d_el = rad_el-sun_el
    d_az = (rad_az-sun_az)*np.cos(sun_el*np.pi/180.)

    npix_az = int((azmax-azmin)/azres)
    npix_el = int((elmax-elmin)/elres)
    field = np.ma.zeros((npix_az, npix_el))
    field[:] = np.ma.masked

    ind_az = ((d_az+azmin)/azres).astype(int)
    ind_el = ((d_el+elmin)/elres).astype(int)

    field[ind_az, ind_el] = data

    return field


def create_sun_retrieval_field(par, imgcfg):
    """
    creates a sun retrieval field from the retrieval parameters

    Parameters
    ----------
    par : ndarray 1D
        the 5 retrieval parameters
    imgcfg: dict
        a dictionary specifying the ranges and resolution of the field to
        create

    Returns
    -------
    field : masked ndarray 2D
        the sun retrieval field

    """
    azmin = imgcfg['azmin']
    azmax = imgcfg['azmax']
    elmin = imgcfg['elmin']
    elmax = imgcfg['elmax']
    azres = imgcfg['azres']
    elres = imgcfg['elres']

    npix_az = int((azmax-azmin)/azres)
    npix_el = int((elmax-elmin)/elres)

    field = np.ma.zeros((npix_az, npix_el))
    field[:] = np.ma.masked

    d_az = np.array(np.array(range(npix_az))*azres+azmin)
    d_el = np.array(np.array(range(npix_el))*elres+elmin)

    d_az_mat = np.broadcast_to(d_az.reshape(npix_az, 1), (npix_az, npix_el))
    d_el_mat = np.broadcast_to(d_el.reshape(1, npix_el), (npix_az, npix_el))

    field = (par[0]+par[1]*d_az_mat+par[2]*d_el_mat+par[3]*d_az_mat*d_az_mat +
             par[4]*d_el_mat*d_el_mat)

    return field


def compute_quantiles_sweep(field, ray_start, ray_end, quantiles=None):
    """
    computes quantiles of a particular sweep

    Parameters
    ----------
    field : ndarray 2D
        the radar field
    ray_start, ray_end : int
        starting and ending ray indexes
    quantiles: float array
        list of quantiles to compute

    Returns
    -------
    quantiles : float array
        list of quantiles
    values : float array
        values at each quantile

    """
    if quantiles is None:
        quantiles = [10., 20., 30., 40., 50., 60., 70., 80., 90.]
        warn('No quantiles have been defined. Default ' + str(quantiles) +
             ' will be used')
    nquantiles = len(quantiles)
    values = np.ma.zeros(nquantiles)
    values[:] = np.ma.masked

    data_valid = field[ray_start:ray_end, :].compressed()
    for i in range(nquantiles):
        values[i] = np.percentile(data_valid, quantiles[i])

    return quantiles, values


def compute_histogram_sweep(field, ray_start, ray_end, field_name, step=None):
    """
    computes histogram of the data in a particular sweep

    Parameters
    ----------
    field : ndarray 2D
        the radar field
    ray_start, ray_end : int
        starting and ending ray indexes
    field_name: str
        name of the field
    step : float
        size of bin

    Returns
    -------
    bins : float array
        interval of each bin
    values : float array
        values at each bin

    """
    vmin, vmax = pyart.config.get_field_limits(field_name)
    if step is None:
        step = (vmax-vmin)/50.
        warn('No step has been defined. Default '+str(step)+' will be used')
    bins = np.linspace(vmin, vmax, num=int((vmax-vmin)/step))
    values = field[ray_start:ray_end, :].compressed()

    return bins, values
