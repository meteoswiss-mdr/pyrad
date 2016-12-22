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
    plot_density
    plot_scatter
    plot_quantiles
    plot_histogram
    plot_histogram2
    plot_timeseries
    plot_timeseries_comp
    plot_monitoring_ts
    plot_sun_hits
    plot_sun_retrieval_ts
    get_colobar_label
    get_field_name

"""

from copy import deepcopy
from warnings import warn

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

import pyart

from ..util.radar_utils import compute_quantiles_sweep
from ..util.radar_utils import compute_quantiles_from_hist
from ..util.radar_utils import compute_histogram_sweep


def plot_ppi(radar, field_name, ind_el, prdcfg, fname_list, plot_type='PPI',
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
    fname_list : list of str
        list of names of the files where to store the plot
    plot_type : str
        type of plot (PPI, QUANTILES or HISTOGRAM)
    step : float
        step for histogram plotting
    quantiles : float array
        quantiles to plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

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

        for i in range(len(fname_list)):
            fig.savefig(fname_list[i])
        plt.close()

    elif plot_type == 'QUANTILES':
        quantiles, values = compute_quantiles_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_el],
            radar.sweep_end_ray_index['data'][ind_el], quantiles=quantiles)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_el)
        labely = get_colobar_label(radar.fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname_list, labelx='quantile',
                       labely=labely, titl=titl)

    elif plot_type == 'HISTOGRAM':
        bins, values = compute_histogram_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_el],
            radar.sweep_end_ray_index['data'][ind_el], field_name, step=step)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_el)
        labelx = get_colobar_label(radar.fields[field_name], field_name)

        plot_histogram(bins, values, fname_list, labelx=labelx,
                       labely='Number of Samples', titl=titl)

    return fname_list


def plot_rhi(radar, field_name, ind_az, prdcfg, fname_list, plot_type='PPI',
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
    fname_list : list of str
        list of names of the files where to store the plot
    plot_type : str
        type of plot (PPI, QUANTILES or HISTOGRAM)
    step : float
        step for histogram plotting
    quantiles : float array
        quantiles to plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

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

        for i in range(len(fname_list)):
            fig.savefig(fname_list[i])
        plt.close()
    elif plot_type == 'QUANTILES':
        quantiles, values = compute_quantiles_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_az],
            radar.sweep_end_ray_index['data'][ind_az], quantiles=quantiles)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_az)
        labely = get_colobar_label(radar.fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname_list, labelx='quantile',
                       labely=labely, titl=titl)

    elif plot_type == 'HISTOGRAM':
        bins, values = compute_histogram_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_az],
            radar.sweep_end_ray_index['data'][ind_az], field_name, step=step)

        titl = pyart.graph.common.generate_title(radar, field_name, ind_az)
        labelx = get_colobar_label(radar.fields[field_name], field_name)

        plot_histogram(bins, values, fname_list, labelx=labelx,
                       labely='Number of Samples', titl=titl)

    return fname_list


def plot_bscope(radar, field_name, ind_sweep, prdcfg, fname_list):
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
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

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

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_cappi(radar, field_name, altitude, prdcfg, fname_list):
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
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

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

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_density(hist_obj, hist_type, field_name, ind_sweep, prdcfg,
                 fname_list, quantiles=[25., 50., 75.], ref_value=0.):
    """
    density plot (angle-values representation)

    Parameters
    ----------
    hist_obj : histogram object
        object containing the histogram data to plot
    hist_type : str
        type of histogram (instantaneous data or cumulative)
    field_name : str
        name of the radar field to plot
    ind_sweep : int
        sweep index to plot
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    quantiles : array
        the quantile lines to plot
    ref_value : float
        the reference value

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    hist_obj_aux = hist_obj.extract_sweeps([ind_sweep])
    if hist_obj_aux.scan_type == 'ppi':
        ang = np.sort(hist_obj_aux.azimuth['data'])
        ind_ang = np.argsort(hist_obj_aux.azimuth['data'])
        ang_min = np.min(hist_obj_aux.azimuth['data'])
        ang_max = np.max(hist_obj_aux.azimuth['data'])
        field = hist_obj_aux.fields[field_name]['data'][ind_ang, :]
        labelx = 'azimuth angle (degrees)'
    elif hist_obj_aux.scan_type == 'rhi':
        ang = np.sort(hist_obj_aux.elevation['data'])
        ind_ang = np.argsort(hist_obj_aux.elevation['data'])
        ang_min = np.min(hist_obj_aux.elevation['data'])
        ang_max = np.max(hist_obj_aux.elevation['data'])
        field = hist_obj_aux.fields[field_name]['data'][ind_ang, :]
        labelx = 'elevation angle (degrees)'
    else:
        field = hist_obj_aux.fields[field_name]['data']
        ang = np.array(range(hist_obj_aux.nrays))
        ang_min = 0
        ang_max = hist_obj_aux.nrays-1
        labelx = 'ray number'

    # compute percentiles of the histogram
    az_percentile_ref = np.ma.empty(len(ang))
    az_percentile_ref[:] = np.ma.masked
    az_percentile_low = deepcopy(az_percentile_ref)
    az_percentile_high = deepcopy(az_percentile_ref)
    for ray in range(len(ang)):
        quantiles, values_ray = compute_quantiles_from_hist(
            hist_obj.range['data'], field[ray, :], quantiles=quantiles)

        az_percentile_low[ray] = values_ray[0]
        az_percentile_ref[ray] = values_ray[1]
        az_percentile_high[ray] = values_ray[2]

    quantiles, values_sweep = compute_quantiles_from_hist(
        hist_obj.range['data'], np.ma.sum(field, axis=0),
        quantiles=quantiles)

    # mask 0 data
    field = np.ma.masked_where(field == 0, field)

    # display data
    if hist_type == 'instant':
        titl = pyart.graph.common.generate_title(
            hist_obj_aux, field_name, ind_sweep)
    else:
        titl = (
                '{:.1f}'.format(hist_obj_aux.fixed_angle['data'][0])+' Deg. ' +
                pyart.graph.common.generate_radar_time_begin(
                    hist_obj).strftime('%Y-%m-%d') + '\n' +
                get_field_name(hist_obj.fields[field_name], field_name))
    label = 'Number of Points'
    labely = get_colobar_label(hist_obj_aux.fields[field_name], field_name)

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=72)
    ax = fig.add_subplot(111)

    cmap = pyart.config.get_field_colormap(field_name)
    vmin, vmax = pyart.config.get_field_limits(field_name)

    rmin = hist_obj_aux.range['data'][0]/1000.
    rmax = hist_obj_aux.range['data'][-1]/1000.
    cax = ax.imshow(
        np.ma.transpose(field), origin='lower', cmap=cmap, vmin=0.,
        vmax=np.max(field), extent=(ang_min, ang_max, vmin, vmax),
        aspect='auto', interpolation='none')

    # plot reference
    plt.plot(ang, np.zeros(len(ang))+ref_value, 'k--')

    # plot quantiles
    plt.plot(ang, np.zeros(len(ang))+values_sweep[1], 'r')
    plt.plot(ang, np.zeros(len(ang))+values_sweep[0], 'r--')
    plt.plot(ang, np.zeros(len(ang))+values_sweep[2], 'r--')

    plt.plot(ang, az_percentile_ref, 'k')
    plt.plot(ang, az_percentile_low, 'k--')
    plt.plot(ang, az_percentile_high, 'k--')

    plt.autoscale(enable=True, axis='both', tight=True)

    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    cb = fig.colorbar(cax)
    cb.set_label(label)

    metadata = (
            'npoints: '+str(np.ma.sum(field))+'\n' +
            str(quantiles[1]) + ' quant: '+str(values_sweep[1])+'\n' +
            str(quantiles[0]) + ' quant: '+str(values_sweep[0])+'\n' +
            str(quantiles[2]) + ' quant: '+str(values_sweep[2])+'\n')

    plt.text(0.05, 0.95, metadata, horizontalalignment='left',
             verticalalignment='top', transform=ax.transAxes)

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_scatter(bins1, bins2, hist_2d, field_name1, field_name2, fname_list,
                 prdcfg, metadata=None):
    """
    2D histogram

    Parameters
    ----------
    bins1, bins2 : float array2
        the bins of each field
    hist_2d : ndarray 2D
        the 2D histogram
    field_name1, field_name2 : str
        the names of each field
    fname_list : list of str
        list of names of the files where to store the plot
    prdcfg : dict
        product configuration dictionary
    metadata : str
        a string with metadata to write in the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    # mask 0 data
    hist_2d = np.ma.masked_where(hist_2d == 0, hist_2d)

    # display data
    titl = 'colocated radar gates'
    label = 'Number of Points'
    labelx = field_name1
    labely = field_name2

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=72)
    ax = fig.add_subplot(111)

    cmap = pyart.config.get_field_colormap(field_name1)

    cax = ax.imshow(
        hist_2d, origin='lower', cmap=cmap, vmin=0., vmax=np.max(hist_2d),
        extent=(bins1[0], bins1[-1], bins2[0], bins2[-1]),
        aspect='auto', interpolation='none')

    # plot reference
    plt.plot(bins1, bins2, 'k--')
    plt.autoscale(enable=True, axis='both', tight=True)

    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    cb = fig.colorbar(cax)
    cb.set_label(label)

    if metadata is not None:
        plt.text(0.05, 0.95, metadata, horizontalalignment='left',
                 verticalalignment='top', transform=ax.transAxes)

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_quantiles(quant, value, fname_list, labelx='quantile', labely='value',
                   titl='quantile'):
    """
    plots quantiles

    Parameters
    ----------
    quant : array
        quantiles to be plotted
    value : array
        values of each quantile
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig = plt.figure(figsize=[10, 6])
    plt.plot(quant, value, 'bx-')
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_histogram(bins, values, fname_list, labelx='bins',
                   labely='Number of Samples', titl='histogram'):
    """
    computes and plots histogram

    Parameters
    ----------
    bins : array
        histogram bins
    values : array
        data values
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig = plt.figure(figsize=[10, 6])
    plt.hist(values, bins=bins)
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_histogram2(bins, hist, fname_list, labelx='bins',
                    labely='Number of Samples', titl='histogram'):
    """
    plots histogram

    Parameters
    ----------
    quant : array
        histogram bins
    hist : array
        values for each bin
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig = plt.figure(figsize=[10, 6])
    plt.bar(bins, hist, width=bins[1]-bins[0])
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_timeseries(date, value, fname_list, labelx='Time [UTC]',
                    labely='Value', label1='Sensor', titl='Time Series',
                    period=0, timeformat=None):
    """
    plots a time series

    Parameters
    ----------
    date : datetime object
        time of the time series
    value : float array
        values of the time series
    fname_list : list of str
        list of names of the files where to store the plot
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
    timeformat : str
        Specifies the date and time format on the x axis

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    if period > 0:
        value *= (period/3600.)
        value = np.ma.cumsum(value)

    fig, ax = plt.subplots(figsize=[10, 6])
    ax.plot(date, value, label=label1)
    ax.set_title(titl)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)

    if (timeformat is not None):
        ax.xaxis.set_major_formatter(mdates.DateFormatter(timeformat))

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_timeseries_comp(date1, value1, date2, value2, fname_list,
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
    fname_list : list of str
        list of names of the files where to store the plot
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
    fname_list : list of str
        list of names of the created plots

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

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_monitoring_ts(date, np_t, cquant, lquant, hquant, field_name,
                       fname_list, ref_value=None, labelx='Time [UTC]',
                       labely='Value', titl='Time Series'):
    """
    plots a time series of monitoring data

    Parameters
    ----------
    date : datetime object
        time of the time series
    cquant, lquant, hquant : float array
        values of the central, low and high quantiles
    field_name : str
        name of the field
    fname_list : list of str
        list of names of the files where to store the plot
    ref_value : float
        the reference value
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    vmin, vmax = pyart.config.get_field_limits(field_name)

    fig = plt.figure(figsize=[10, 6])

    ax = fig.add_subplot(2, 1, 1)
    plt.plot(date, cquant)
    plt.plot(date, lquant, 'r')
    plt.plot(date, hquant, 'r')
    if ref_value is not None:
        plt.plot(date, np.zeros(len(date))+ref_value, 'k--')
    plt.ylabel(labely)
    plt.title(titl)

    axes = plt.gca()
    axes.set_ylim([vmin, vmax])

    ax = fig.add_subplot(2, 1, 2)
    plt.plot(date, np_t)

    plt.ylabel('Number of Samples')
    plt.xlabel(labelx)

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_sun_hits(field, field_name, fname_list, prdcfg):
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
    fname_list : list of str
        list of names of the files where to store the plot

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

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

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


def plot_sun_retrieval_ts(sun_retrieval, data_type, fname_list):
    """
    plots a time series

    Parameters
    ----------
    date : datetime object
        time of the time series
    value : float array
        values of the time series
    fname_list : list of str
        list of names of the files where to store the plot
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
    fname_list : list of str
        list of names of the created plots

    """
    labelx = 'Date'
    titl = 'Sun retrieval Time Series'

    value_std = None
    ref = None
    date = sun_retrieval[1]
    if data_type == 'nhits_h':
        value = sun_retrieval[2]
        labely = 'Number of sun hits H channel'
        vmin = 0
        vmax = np.max(sun_retrieval[2])+1
    elif data_type == 'el_width_h':
        value = sun_retrieval[3]
        labely = 'Elevation beamwidth H channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'az_width_h':
        value = sun_retrieval[4]
        labely = 'Azimuth beamwidth H channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'el_bias_h':
        value = sun_retrieval[5]
        ref = np.zeros(len(value))
        labely = 'Elevation pointing bias H channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'az_bias_h':
        value = sun_retrieval[6]
        ref = np.zeros(len(value))
        labely = 'Azimuth pointing bias H channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'dBm_sun_est':
        value = sun_retrieval[7]
        value_std = sun_retrieval[8]
        ref = sun_retrieval[19]
        labely = 'Sun Power H channel (dBm)'
        vmin = -110.
        vmax = -90.
    elif data_type == 'rx_bias_h':
        value = sun_retrieval[7]-sun_retrieval[19]
        value_std = sun_retrieval[8]
        ref = np.zeros(len(value))
        labely = 'Receiver bias H channel (dB)'
        vmin = -5.
        vmax = 5.
    elif data_type == 'nhits_v':
        value = sun_retrieval[9]
        labely = 'Number of sun hits V channel'
        vmin = 0
        vmax = np.max(sun_retrieval[9])+1
    elif data_type == 'el_width_v':
        value = sun_retrieval[10]
        labely = 'Elevation beamwidth V channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'az_width_v':
        value = sun_retrieval[11]
        labely = 'Azimuth beamwidth V channel (Deg)'
        vmin = 0.
        vmax = 4.
    elif data_type == 'el_bias_v':
        value = sun_retrieval[12]
        ref = np.zeros(len(value))
        labely = 'Elevation pointing bias V channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'az_bias_v':
        value = sun_retrieval[13]
        ref = np.zeros(len(value))
        labely = 'Azimuth pointing bias V channel (Deg)'
        vmin = -2.
        vmax = 2.
    elif data_type == 'dBmv_sun_est':
        value = sun_retrieval[14]
        value_std = sun_retrieval[15]
        ref = sun_retrieval[19]
        labely = 'Sun Power V channel (dBm)'
        vmin = -110.
        vmax = -90.
    elif data_type == 'rx_bias_v':
        value = sun_retrieval[14]-sun_retrieval[19]
        value_std = sun_retrieval[15]
        ref = np.zeros(len(value))
        labely = 'Receiver bias V channel (dB)'
        vmin = -5.
        vmax = 5.
    elif data_type == 'nhits_zdr':
        value = sun_retrieval[16]
        labely = 'Number of sun hits ZDR'
        vmin = 0
        vmax = np.max(sun_retrieval[16])+1
    elif data_type == 'ZDR_sun_est':
        value = sun_retrieval[17]
        value_std = sun_retrieval[18]
        ref = np.zeros(len(value))
        labely = 'Sun ZDR (dB)'
        vmin = -2.
        vmax = 2.

    mask = np.ma.getmaskarray(value)
    if mask.all():
        warn('Unable to create figure '+fname_list+'. No valid data')
        return None

    fig = plt.figure(figsize=[10, 6])
    plt.plot(date, value)
    if value_std is not None:
        plt.plot(date, value+value_std, 'r')
        plt.plot(date, value-value_std, 'r')
    if ref is not None:
        plt.plot(date, ref, 'k--')
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    axes = plt.gca()
    axes.set_ylim([vmin, vmax])

    for i in range(len(fname_list)):
        fig.savefig(fname_list[i])
    plt.close()

    return fname_list


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
