"""
pyrad.graph.plots_grid
======================

Functions to plot data in a Cartesian grid format

.. autosummary::
    :toctree: generated/

    plot_surface
    plot_latitude_slice
    plot_longitude_slice
    plot_latlon_slice

"""

import numpy as np

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.pyplot as plt

import pyart

from .plots_aux import get_norm


def plot_surface(grid, field_name, level, prdcfg, fname_list):
    """
    plots a surface from gridded data

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    level : int
        level index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    History
    --------
    ????.??.?? created
    2017.08.?? -fvj- added option controlling dpi

    """
    dpi = 72
    if 'dpi' in prdcfg['ppiImageConfig']:
        dpi = prdcfg['ppiImageConfig']['dpi']

    norm, ticks, ticklabs = get_norm(field_name)

    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']
    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    fig.add_subplot(111, aspect='equal')
    lon_lines = np.arange(np.floor(prdcfg['ppiMapImageConfig']['lonmin']),
                          np.ceil(prdcfg['ppiMapImageConfig']['lonmax'])+1,
                          0.5)
    lat_lines = np.arange(np.floor(prdcfg['ppiMapImageConfig']['latmin']),
                          np.ceil(prdcfg['ppiMapImageConfig']['latmax'])+1,
                          0.5)
    display = pyart.graph.GridMapDisplay(grid)
    display.plot_basemap(lat_lines=lat_lines, lon_lines=lon_lines)
    display.plot_grid(field_name, level=level, norm=norm, ticks=ticks,
                      ticklabs=ticklabs)
    # display.plot_crosshairs(lon=lon, lat=lat)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()


def plot_latitude_slice(grid, field_name, lon, lat, prdcfg, fname_list):
    """
    plots a latitude slice from gridded data

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    lon, lat : float
        coordinates of the slice to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    History
    --------
    ????.??.?? created
    2017.08.?? -fvj- added option controlling dpi

    """
    dpi = 72
    if 'dpi' in prdcfg['rhiImageConfig']:
        dpi = prdcfg['rhiImageConfig']['dpi']

    norm, ticks, ticklabs = get_norm(field_name)

    xsize = prdcfg['rhiImageConfig']['xsize']
    ysize = prdcfg['rhiImageConfig']['ysize']
    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    display = pyart.graph.GridMapDisplay(grid)
    display.plot_latitude_slice(
        field_name, lon=lon, lat=lat, norm=norm, colorbar_orient='horizontal',
        ticks=ticks, ticklabs=ticklabs)
    ax.set_xlim(
        [prdcfg['rhiImageConfig']['xmin'], prdcfg['rhiImageConfig']['xmax']])
    ax.set_ylim(
        [prdcfg['rhiImageConfig']['ymin'], prdcfg['rhiImageConfig']['ymax']])

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()


def plot_longitude_slice(grid, field_name, lon, lat, prdcfg, fname_list):
    """
    plots a longitude slice from gridded data

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    lon, lat : float
        coordinates of the slice to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    History
    --------
    ????.??.?? created
    2017.08.?? -fvj- added option controlling dpi

    """
    dpi = 72
    if 'dpi' in prdcfg['rhiImageConfig']:
        dpi = prdcfg['rhiImageConfig']['dpi']

    norm, ticks, ticklabs = get_norm(field_name)

    xsize = prdcfg['rhiImageConfig']['xsize']
    ysize = prdcfg['rhiImageConfig']['ysize']
    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    display = pyart.graph.GridMapDisplay(grid)
    display.plot_longitude_slice(
        field_name, lon=lon, lat=lat, norm=norm, colorbar_orient='horizontal',
        ticks=ticks, ticklabs=ticklabs)
    ax.set_xlim(
        [prdcfg['rhiImageConfig']['xmin'], prdcfg['rhiImageConfig']['xmax']])
    ax.set_ylim(
        [prdcfg['rhiImageConfig']['ymin'], prdcfg['rhiImageConfig']['ymax']])

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()


def plot_latlon_slice(grid, field_name, coord1, coord2, prdcfg, fname_list):
    """
    plots a croos section crossing two points in the grid

    Parameters
    ----------
    grid : Grid object
        object containing the gridded data to plot
    field_name : str
        name of the radar field to plot
    coord1 : tupple of floats
        lat, lon of the first point
    coord2 : tupple of floats
        lat, lon of the second point
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    History
    --------
    ????.??.?? created
    2017.08.?? -fvj- added option controlling dpi

    """
    dpi = 72
    if 'dpi' in prdcfg['rhiImageConfig']:
        dpi = prdcfg['rhiImageConfig']['dpi']

    norm, ticks, ticklabs = get_norm(field_name)

    xsize = prdcfg['rhiImageConfig']['xsize']
    ysize = prdcfg['rhiImageConfig']['ysize']
    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    display = pyart.graph.GridMapDisplay(grid)
    display.plot_latlon_slice(
        field_name, coord1=coord1, coord2=coord2, norm=norm,
        colorbar_orient='vertical', ticks=ticks, ticklabs=ticklabs, fig=fig,
        ax=ax)
    # ax.set_ylim(
    #    [prdcfg['rhiImageConfig']['ymin'], prdcfg['rhiImageConfig']['ymax']])

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()
