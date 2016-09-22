"""
pyrad.graph.plots
=================

Functions to plot Pyrad datasets

.. autosummary::
    :toctree: generated/

    plot_ppi
    plot_rhi
    plot_cappi
    plot_bscope
    plot_timeseries
    plot_timeseries_comp

"""

from warnings import warn

import matplotlib.pyplot as plt
import numpy as np

import pyart


def plot_ppi(radar, field_name, ind_el, prdcfg, fname):
    """
    plots a PPI

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot

    field_name : str
        name of the radar field to plot

    ind_el : int
        sweep index to plot

    prdcfg : dict
        dictionary containing the product configuration

    fname : str
        name of the file where to store the plot

    Returns
    -------
    no return

    """
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=72)
    display.plot_ppi(field_name, sweep=ind_el)
    display.set_limits(
        ylim=[prdcfg['ppiImageConfig']['ymin'],
              prdcfg['ppiImageConfig']['ymax']],
        xlim=[prdcfg['ppiImageConfig']['xmin'],
              prdcfg['ppiImageConfig']['xmax']])
    display.plot_range_rings([10, 20, 30, 40])
    display.plot_cross_hair(5.)

    fig.savefig(fname)

    plt.close()


def plot_rhi(radar, field_name, ind_az, prdcfg, fname):
    """
    plots an RHI

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot

    field_name : str
        name of the radar field to plot

    ind_az : int
        sweep index to plot

    prdcfg : dict
        dictionary containing the product configuration

    fname : str
        name of the file where to store the plot

    Returns
    -------
    no return

    """
    display = pyart.graph.RadarDisplay(radar)
    fig = plt.figure(figsize=[prdcfg['rhiImageConfig']['xsize'],
                              prdcfg['rhiImageConfig']['ysize']],
                     dpi=72)
    display.plot_rhi(field_name, sweep=ind_az, reverse_xaxis=False)
    display.set_limits(
        ylim=[prdcfg['rhiImageConfig']['ymin'],
              prdcfg['rhiImageConfig']['ymax']],
        xlim=[prdcfg['rhiImageConfig']['xmin'],
              prdcfg['rhiImageConfig']['xmax']])
    display.plot_cross_hair(5.)

    fig.savefig(fname)
    plt.close()


def plot_bscope(radar, field_name, ind_sweep, prdcfg, fname):
    """
    plots a B-Scope (angle-range representation)

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot

    field_name : str
        name of the radar field to plot

    ind_sweep : int
        sweep index to plot

    prdcfg : dict
        dictionary containing the product configuration

    fname : str
        name of the file where to store the plot

    Returns
    -------
    no return

    """
    radar_aux = radar.extract_sweeps([ind_sweep])
    if radar_aux.scan_type == 'ppi':
        ang = np.sort(radar_aux.azimuth['data'])
        ind_ang = np.argsort(radar_aux.azimuth['data'])
        ang_min = np.min(radar_aux.azimuth['data'])
        ang_max = np.max(radar_aux.azimuth['data'])
        field = radar_aux.fields[field_name]['data'][ind_ang, :]
        labely = 'azimuth angle (degrees)'
    elif radar_aux.scan_type == 'rhi':
        ang = np.sort(radar_aux.elevation['data'])
        ind_ang = np.argsort(radar_aux.elevation['data'])
        ang_min = np.min(radar_aux.elevation['data'])
        ang_max = np.max(radar_aux.elevation['data'])
        field = radar_aux.fields[field_name]['data'][ind_ang, :]
        labely = 'elevation angle (degrees)'
    else:
        field = radar_aux.fields[field_name]['data']
        ang = np.array(range(radar_aux.nrays))
        ang_min = 0
        ang_max = radar_aux.nrays-1
        labely = 'ray number'

    # display data
    titl = pyart.graph.common.generate_title(radar_aux, field_name, ind_sweep)

    # set the colorbar/yaxis label.
    if 'standard_name' in radar_aux.fields[field_name]:
        standard_name = radar_aux.fields[field_name]['standard_name']
    elif 'long_name' in radar_aux.fields[field_name]:
        standard_name = radar_aux.fields[field_name]['long_name']
    else:
        standard_name = field_name

    if 'units' in radar_aux.fields[field_name]:
        units = radar_aux.fields[field_name]['units']
    else:
        units = '?'

    label = pyart.graph.common.generate_colorbar_label(standard_name, units)

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=72)
    ax = fig.add_subplot(111)
    if radar_aux.ngates == 1:
        plt.plot(ang, field, 'bx')
        plt.xlabel(labely)
        plt.ylabel(label)
        plt.title(titl)
    else:
        cmap = pyart.config.get_field_colormap(field_name)
        vmin, vmax = pyart.config.get_field_limits(field_name)

        rmin = radar_aux.range['data'][0]/1000.
        rmax = radar_aux.range['data'][-1]/1000.
        cax = ax.imshow(
            field, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax,
            extent=(rmin, rmax, ang_min, ang_max), aspect='auto')
        plt.xlabel('Range (km)')
        plt.ylabel(labely)
        plt.title(titl)

        cb = fig.colorbar(cax)
        cb.set_label(label)

    fig.savefig(fname)
    plt.close()


def plot_cappi(radar, field_name, altitude, prdcfg, fname):
    """
    plots a Constant Altitude Plan Position Indicator CAPPI

    Parameters
    ----------
    radar : Radar object
        object containing the radar data to plot

    field_name : str
        name of the radar field to plot

    altitude : float
        the altitude [m MSL] to be plotted

    prdcfg : dict
        dictionary containing the product configuration

    fname : str
        name of the file where to store the plot

    Returns
    -------
    no return

    """
    xmin = prdcfg['ppiImageConfig']['xmin']
    xmax = prdcfg['ppiImageConfig']['xmax']
    ymin = prdcfg['ppiImageConfig']['ymin']
    ymax = prdcfg['ppiImageConfig']['ymax']

    # cartesian mapping
    grid = pyart.map.grid_from_radars(
        (radar,), grid_shape=(1, 241, 241),
        grid_limits=((altitude, altitude), (ymin*1000., ymax*1000.),
                     (xmin*1000., xmax*1000.)),
        fields=[field_name])

    # display data
    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=72)
    ax = fig.add_subplot(111)
    cmap = pyart.config.get_field_colormap(field_name)
    vmin, vmax = pyart.config.get_field_limits(field_name)
    titl = pyart.graph.common.generate_grid_title(grid, field_name, 0)

    cax = ax.imshow(
        grid.fields[field_name]['data'][0], extent=(xmin, xmax, ymin, ymax),
        origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
    plt.xlabel('East West distance from radar(km)')
    plt.ylabel('North South distance from radar(km)')
    plt.title(titl)

    # plot the colorbar and set the label.
    if 'standard_name' in grid.fields[field_name]:
        standard_name = grid.fields[field_name]['standard_name']
    elif 'long_name' in grid.fields[field_name]:
        standard_name = grid.fields[field_name]['long_name']
    else:
        standard_name = field_name

    if 'units' in grid.fields[field_name]:
        units = grid.fields[field_name]['units']
    else:
        units = '?'

    label = pyart.graph.common.generate_colorbar_label(standard_name, units)

    cb = fig.colorbar(cax)
    cb.set_label(label)

    fig.savefig(fname)
    plt.close()


def plot_timeseries(date, value, fname, labelx='Time [UTC]', labely='Value',
                    label1='Sensor', titl='Time Series', period=0):
    """
    plots a time series

    Parameters
    ----------
    date : datetime object
        time of the time series
    value : float array
        values of the time series
    fname : str
        name of the file where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    label1 : str
        The label of the legend
    titl : str
        The figure title
    period : float
        measurement period in seconds used to compute accumulation. If 0 no
        accumulation is computed

    Returns
    -------
    no return

    """
    if period > 0:
        value *= (period/3600.)
        value = np.ma.cumsum(value)

    fig = plt.figure(figsize=[10, 6])
    plt.plot(date, value, label=label1)
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    fig.savefig(fname)
    plt.close()


def plot_timeseries_comp(date1, value1, date2, value2, fname,
                         labelx='Time [UTC]', labely='Value',
                         label1='Sensor 1', label2='Sensor 2',
                         titl='Time Series Comparison', period1=0, period2=0):
    """
    plots 2 time series in the same graph

    Parameters
    ----------
    date1 : datetime object
        time of the first time series

    value1 : float array
        values of the first time series

    date2 : datetime object
        time of the second time series

    value2 : float array
        values of the second time series

    fname : str
        name of the file where to store the plot

    labelx : str
        The label of the X axis

    labely : str
        The label of the Y axis

    label1, label2 : str
        legend label for each time series

    titl : str
        The figure title

     period1, period2 : float
        measurement period in seconds used to compute accumulation. If 0 no
        accumulation is computed

    Returns
    -------
    no return

    """
    if (period1 > 0) and (period2 > 0):
        value1 *= (period1/3600.)
        value1 = np.ma.cumsum(value1)

        value2 *= (period2/3600.)
        value2 = np.ma.cumsum(value2)

    fig = plt.figure(figsize=[10, 6])
    plt.plot(date1, value1, 'b', label=label1)
    plt.plot(date2, value2, 'r', label=label2)
    plt.legend(loc='best')
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    fig.savefig(fname)
    plt.close()
