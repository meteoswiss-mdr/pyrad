"""
pyrad.graph.plots
=================

Functions to plot Pyrad datasets

.. autosummary::
    :toctree: generated/

    plot_ppi
    plot_rhi
    plot_bscope
    plot_cappi
    plot_quantiles
    plot_histogram
    plot_timeseries
    plot_timeseries_comp
    plot_sun_hits
    plot_sun_retrieval_ts
    get_colobar_label
    get_field_name

"""

from warnings import warn

import matplotlib.pyplot as plt
import numpy as np

import pyart

from ..util.radar_utils import compute_quantiles_sweep
from ..util.radar_utils import compute_histogram_sweep


def plot_ppi(radar, field_name, ind_el, prdcfg, fname, plot_type='PPI',
             step=None, quantiles=None):
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
    plot_type : str
        type of plot (PPI, QUANTILES or HISTOGRAM)
    step : float
        step for histogram plotting
    quantiles : float array
        quantiles to plot

    Returns
    -------
    fname : str
        the name of the created plot file

    """
    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=72)

    if plot_type == 'PPI':
        display = pyart.graph.RadarDisplay(radar)
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

    elif plot_type == 'QUANTILES':
        quantiles, values = compute_quantiles_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_el],
            radar.sweep_end_ray_index['data'][ind_el], quantiles=quantiles)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_el)
        labely = get_colobar_label(radar.fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname, labelx='quantile',
                       labely=labely, titl=titl)

    elif plot_type == 'HISTOGRAM':
        bins, values = compute_histogram_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_el],
            radar.sweep_end_ray_index['data'][ind_el], field_name, step=step)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_el)
        labelx = get_colobar_label(radar.fields[field_name], field_name)

        plot_histogram(bins, values, fname, labelx=labelx,
                       labely='Number of Samples', titl=titl)

    return fname


def plot_rhi(radar, field_name, ind_az, prdcfg, fname, plot_type='PPI',
             step=None, quantiles=None):
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
    plot_type : str
        type of plot (PPI, QUANTILES or HISTOGRAM)
    step : float
        step for histogram plotting
    quantiles : float array
        quantiles to plot

    Returns
    -------
    fname : str
        the name of the created plot file

    """
    fig = plt.figure(figsize=[prdcfg['rhiImageConfig']['xsize'],
                              prdcfg['rhiImageConfig']['ysize']],
                     dpi=72)

    if plot_type == 'RHI':
        display = pyart.graph.RadarDisplay(radar)
        display.plot_rhi(field_name, sweep=ind_az, reverse_xaxis=False)
        display.set_limits(
            ylim=[prdcfg['rhiImageConfig']['ymin'],
                  prdcfg['rhiImageConfig']['ymax']],
            xlim=[prdcfg['rhiImageConfig']['xmin'],
                  prdcfg['rhiImageConfig']['xmax']])
        display.plot_cross_hair(5.)

        fig.savefig(fname)
        plt.close()
    elif plot_type == 'QUANTILES':
        quantiles, values = compute_quantiles_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_az],
            radar.sweep_end_ray_index['data'][ind_az], quantiles=quantiles)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_az)
        labely = get_colobar_label(radar.fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname, labelx='quantile',
                       labely=labely, titl=titl)

    elif plot_type == 'HISTOGRAM':
        bins, values = compute_histogram_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_az],
            radar.sweep_end_ray_index['data'][ind_az], field_name, step=step)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_az)
        labelx = get_colobar_label(radar.fields[field_name], field_name)

        plot_histogram(bins, values, fname, labelx=labelx,
                       labely='Number of Samples', titl=titl)

    return fname


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
    fname : str
        the name of the created plot file

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
    label = get_colobar_label(radar_aux.fields[field_name], field_name)

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
            extent=(rmin, rmax, ang_min, ang_max), aspect='auto',
            interpolation='none')
        plt.xlabel('Range (km)')
        plt.ylabel(labely)
        plt.title(titl)

        cb = fig.colorbar(cax)
        cb.set_label(label)

    fig.savefig(fname)
    plt.close()

    return fname


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
    fname : str
        the name of the created plot file

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
        origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='none')
    plt.xlabel('East West distance from radar(km)')
    plt.ylabel('North South distance from radar(km)')
    plt.title(titl)

    # plot the colorbar and set the label.
    label = get_colobar_label(grid.fields[field_name], field_name)
    cb = fig.colorbar(cax)
    cb.set_label(label)

    fig.savefig(fname)
    plt.close()

    return fname


def plot_quantiles(quant, value, fname, labelx='quantile', labely='value',
                   titl='quantile'):
    """
    plots quantiles

    Parameters
    ----------
    quant : array
        quantiles to be plotted
    value : array
        values of each quantie
    fname : str
        name of the file where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title

    Returns
    -------
    fname : str
        the name of the created plot file

    """
    fig = plt.figure(figsize=[10, 6])
    plt.plot(quant, value, 'bx-')
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    fig.savefig(fname)
    plt.close()

    return fname


def plot_histogram(bins, values, fname, labelx='bins',
                   labely='Number of Samples', titl='histogram'):
    """
    plots histogram

    Parameters
    ----------
    quant : array
        quantiles to be plotted
    value : array
        values of each quantie
    fname : str
        name of the file where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title

    Returns
    -------
    fname : str
        the name of the created plot file

    """
    fig = plt.figure(figsize=[10, 6])
    plt.hist(values, bins=bins)
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    fig.savefig(fname)
    plt.close()

    return fname


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
    fname : str
        the name of the created plot file

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

    return fname


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
    fname : str
        the name of the created plot file

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

    return fname


def plot_sun_hits(field, field_name, fname, prdcfg):
    """
    plots the sun hits

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
    fname : str
        the name of the created plot file

    """
    azmin = prdcfg['sunhitsImageConfig']['azmin']
    azmax = prdcfg['sunhitsImageConfig']['azmax']
    elmin = prdcfg['sunhitsImageConfig']['elmin']
    elmax = prdcfg['sunhitsImageConfig']['elmax']

    field_dict = pyart.config.get_metadata(field_name)

    # display data
    fig = plt.figure(figsize=[prdcfg['sunhitsImageConfig']['xsize'],
                              prdcfg['sunhitsImageConfig']['ysize']],
                     dpi=72)
    ax = fig.add_subplot(111)
    cmap = pyart.config.get_field_colormap(field_name)
    vmin, vmax = pyart.config.get_field_limits(field_name)
    titl = (prdcfg['timeinfo'].strftime('%Y-%m-%d') + '\n' +
            get_field_name(field_dict, field_name))

    cax = ax.imshow(
        field, extent=(azmin, azmax, elmin, elmax),
        origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='none')
    plt.xlabel('rad_az-sun_az (deg)')
    plt.ylabel('rad_el-sun_el (deg)')
    plt.title(titl)

    # plot the colorbar and set the label.
    label = get_colobar_label(field_dict, field_name)
    cb = fig.colorbar(cax)
    cb.set_label(label)

    fig.savefig(fname)
    plt.close()

    return fname


def plot_sun_retrieval_ts(sun_retrieval, data_type, fname):
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
    fname : str
        the name of the created plot file

    """
    labelx = 'Date'
    titl = 'Sun retrieval Time Series'

    value_std = None
    sun_ref = None
    date = sun_retrieval[0]
    if data_type == 'nhits_h':
        value = sun_retrieval[1]
        labely = 'Number of sun hits H channel'
        vmin = 0
        vmax = 30
    elif data_type == 'el_width_h':
        value = sun_retrieval[2]
        labely = 'Elevation beamwidth H channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'az_width_h':
        value = sun_retrieval[3]
        labely = 'Azimuth beamwidth H channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'el_bias_h':
        value = sun_retrieval[4]
        labely = 'Elevation bias H channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'az_bias_h':
        value = sun_retrieval[5]
        labely = 'Azimuth bias H channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'dBm_sun_est':
        value = sun_retrieval[6]
        value_std = sun_retrieval[7]
        labely = 'Sun Power H channel (dBm)'
        vmin = -110.
        vmax = -90.
    elif data_type == 'nhits_v':
        value = sun_retrieval[8]
        labely = 'Number of sun hits V channel'
        vmin = 0
        vmax = 30
    elif data_type == 'el_width_v':
        value = sun_retrieval[9]
        labely = 'Elevation beamwidth V channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'az_width_v':
        value = sun_retrieval[10]
        labely = 'Azimuth beamwidth V channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'el_bias_v':
        value = sun_retrieval[11]
        labely = 'Elevation bias V channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'az_bias_v':
        value = sun_retrieval[12]
        labely = 'Azimuth bias V channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'dBmv_sun_est':
        value = sun_retrieval[13]
        value_std = sun_retrieval[14]
        labely = 'Sun Power V channel (dBm)'
        vmin = -110.
        vmax = -90.
    elif data_type == 'nhits_zdr':
        value = sun_retrieval[15]
        labely = 'Number of sun hits ZDR'
        vmin = 0
        vmax = 30
    elif data_type == 'ZDR_sun_est':
        value = sun_retrieval[16]
        value_std = sun_retrieval[17]
        labely = 'Sun ZDR (dB)'
        vmin = -2.
        vmax = 2.

    fig = plt.figure(figsize=[10, 6])
    plt.plot(date, value)
    if value_std is not None:
        plt.plot(date, value+value_std, 'r')
        plt.plot(date, value-value_std, 'r')
    if sun_ref is not None:
        plt.plot(date, sun_ref, 'r')
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    axes = plt.gca()
    axes.set_ylim([vmin, vmax])

    fig.savefig(fname)
    plt.close()

    return fname


def get_colobar_label(field_dict, field_name):
    """
    creates the colorbar label using field metadata

    Parameters
    ----------
    field_dict : dict
        dictionary containing field metadata
    field_name : str
        name of the field

    Returns
    -------
    label : str
        colorbar label

    """
    if 'standard_name' in field_dict:
        standard_name = field_dict['standard_name']
    elif 'long_name' in field_dict:
        standard_name = field_dict['long_name']
    else:
        standard_name = field_name

    if 'units' in field_dict:
        units = field_dict['units']
    else:
        units = '?'

    return pyart.graph.common.generate_colorbar_label(standard_name, units)


def get_field_name(field_dict, field):
    """
    Return a nice field name for a particular field

    Parameters
    ----------
    field_dict : dict
        dictionary containing field metadata
    field : str
        name of the field

    Returns
    -------
    field_name : str
        the field name

    """
    if 'standard_name' in field_dict:
        field_name = field_dict['standard_name']
    elif 'long_name' in field_dict:
        field_name = field_dict['long_name']
    else:
        field_name = str(field)
    field_name = field_name.replace('_', ' ')
    field_name = field_name[0].upper() + field_name[1:]

    return field_name
