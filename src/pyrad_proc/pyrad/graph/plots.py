"""
pyrad.graph.plots
=================

Functions to plot Pyrad datasets

.. autosummary::
    :toctree: generated/

    plot_density
    plot_scatter
    plot_quantiles
    plot_histogram
    plot_histogram2
    plot_antenna_pattern
    plot_scatter_comp
    plot_sun_hits

"""

from copy import deepcopy

import numpy as np

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.pyplot as plt

import pyart

from .plots_aux import get_colobar_label, get_field_name

from ..util.radar_utils import compute_quantiles_from_hist


def plot_density(hist_obj, hist_type, field_name, ind_sweep, prdcfg,
                 fname_list, quantiles=[25., 50., 75.], ref_value=0.,
                 vmin=None, vmax=None):
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
    vmin, vmax : float
        Minim and maximum extend of the vertical axis

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
        ang_step = ang[1]-ang[0]
        field = hist_obj_aux.fields[field_name]['data'][ind_ang, :]
        labelx = 'azimuth angle (degrees)'
    elif hist_obj_aux.scan_type == 'rhi':
        ang = np.sort(hist_obj_aux.elevation['data'])
        ind_ang = np.argsort(hist_obj_aux.elevation['data'])
        ang_min = np.min(hist_obj_aux.elevation['data'])
        ang_max = np.max(hist_obj_aux.elevation['data'])
        ang_step = ang[1]-ang[0]
        field = hist_obj_aux.fields[field_name]['data'][ind_ang, :]
        labelx = 'elevation angle (degrees)'
    else:
        field = hist_obj_aux.fields[field_name]['data']
        ang = np.array(range(hist_obj_aux.nrays))
        ang_min = 0
        ang_max = hist_obj_aux.nrays-1
        ang_step = ang[1]-ang[0]
        labelx = 'ray number'

    # compute percentiles of the histogram
    az_percentile_ref = np.ma.masked_all(len(ang))
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

    dpi = 72
    if 'dpi' in prdcfg['ppiImageConfig']:
        dpi = prdcfg['ppiImageConfig']['dpi']

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=dpi)
    ax = fig.add_subplot(111)

    # define limits of field and color map
    cmap = pyart.config.get_field_colormap(field_name)
    step = hist_obj.range['data'][1]-hist_obj.range['data'][0]
    xmin = ang_min-ang_step/2.
    xmax = ang_max+ang_step/2.
    ymin = hist_obj.range['data'][0]-step/2.
    ymax = hist_obj.range['data'][-1]+step/2.

    cax = ax.imshow(
        np.ma.transpose(field), origin='lower', cmap=cmap, vmin=0.,
        vmax=np.max(field), extent=(xmin, xmax, ymin, ymax),
        aspect='auto', interpolation='none')
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.autoscale(False)

    # plot reference
    ax.plot(ang, np.zeros(len(ang))+ref_value, 'k--')

    # plot quantiles
    ax.plot(ang, np.zeros(len(ang))+values_sweep[1], 'r')
    ax.plot(ang, np.zeros(len(ang))+values_sweep[0], 'r--')
    ax.plot(ang, np.zeros(len(ang))+values_sweep[2], 'r--')

    ax.plot(ang, az_percentile_ref, 'k')
    ax.plot(ang, az_percentile_low, 'k--')
    ax.plot(ang, az_percentile_high, 'k--')

    # ax.autoscale(enable=True, axis='both', tight=True)

    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    cb = fig.colorbar(cax)
    cb.set_label(label)

    val_quant0_str = '--'
    if values_sweep[0] is not np.ma.masked:
        val_quant0_str = '{:.3f}'.format(values_sweep[0])
    val_quant1_str = '--'
    if values_sweep[1] is not np.ma.masked:
        val_quant1_str = '{:.3f}'.format(values_sweep[1])
    val_quant2_str = '--'
    if values_sweep[2] is not np.ma.masked:
        val_quant2_str = '{:.3f}'.format(values_sweep[2])

    metadata = (
        'npoints: '+str(np.ma.sum(field))+'\n' +
        str(quantiles[1])+' quant: '+val_quant1_str+'\n' +
        str(quantiles[0])+' quant: '+val_quant0_str+'\n' +
        str(quantiles[2])+' quant: '+val_quant2_str+'\n')

    ax.text(0.05, 0.05, metadata, horizontalalignment='left',
            verticalalignment='bottom', transform=ax.transAxes)

    # Turn on the grid
    ax.grid()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list


def plot_scatter(bin_edges1, bin_edges2, hist_2d, field_name1, field_name2, fname_list,
                 prdcfg, metadata=None, lin_regr=None, lin_regr_slope1=None,
                 rad1_name='RADAR001', rad2_name='RADAR002'):
    """
    2D histogram

    Parameters
    ----------
    bin_edges1, bin_edges2 : float array2
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
    lin_regr : tupple with 2 values
        the coefficients for a linear regression
    lin_regr_slope1 : float
        the intercep point of a linear regression of slope 1
    rad1_name, rad2_name : str
        name of the radars which data is used

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
    labelx = rad1_name+' '+field_name1
    labely = rad2_name+' '+field_name2

    dpi = 72
    if 'dpi' in prdcfg['ppiImageConfig']:
        dpi = prdcfg['ppiImageConfig']['dpi']

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=dpi)
    ax = fig.add_subplot(111)

    cmap = pyart.config.get_field_colormap(field_name1)

    cax = ax.imshow(
        np.ma.transpose(hist_2d), origin='lower', cmap=cmap, vmin=0.,
        vmax=np.max(hist_2d),
        extent=(bin_edges1[0], bin_edges1[-1], bin_edges2[0], bin_edges2[-1]),
        aspect='auto', interpolation='none')
    ax.autoscale(False)

    # plot reference
    step1 = bin_edges1[1]-bin_edges1[0]
    bin_centers1 = bin_edges1[0:-1]+step1/2.

    step2 = bin_edges2[1]-bin_edges2[0]
    bin_centers2 = bin_edges2[0:-1]+step2/2.

    ax.plot(bin_centers1, bin_centers2, 'k--')

    # plot linear regression
    if lin_regr is not None:
        ax.plot(bin_centers1, lin_regr[0]*bin_centers1+lin_regr[1], 'r')
    if lin_regr_slope1 is not None:
        ax.plot(bin_centers1, bin_centers1+lin_regr_slope1, 'g')

    # ax.autoscale(enable=True, axis='both', tight=True)

    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    cb = fig.colorbar(cax)
    cb.set_label(label)

    if metadata is not None:
        plt.text(0.05, 0.95, metadata, horizontalalignment='left',
                 verticalalignment='top', transform=ax.transAxes)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list


def plot_quantiles(quant, value, fname_list, labelx='quantile', labely='value',
                   titl='quantile', vmin=None, vmax=None, dpi=72):
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
    vmin, vmax: float
        Lower/Upper limit of data values
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)
    ax.plot(quant, value, 'bx-')
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_title(titl)

    # Turn on the grid
    ax.grid()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list


def plot_histogram(bin_edges, values, fname_list, labelx='bins',
                   labely='Number of Samples', titl='histogram', dpi=72):
    """
    computes and plots histogram

    Parameters
    ----------
    bin_edges : array
        histogram bin edges
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
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig = plt.figure(figsize=[10, 6], dpi=dpi)
    plt.hist(values, bins=bin_edges)
    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list


def plot_histogram2(bin_centers, hist, fname_list, labelx='bins',
                    labely='Number of Samples', titl='histogram', dpi=72,
                    ax=None, fig=None, save_fig=True, color=None, alpha=None,
                    invert_xaxis=False):
    """
    plots histogram

    Parameters
    ----------
    bin_centers : array
        histogram bin centers
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
    dpi : int
        dots per inch
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure
    color : str
        color of the bars
    alpha : float
        parameter controling the transparency
    invert_xaxis : bool
        If true inverts the x axis

    Returns
    -------
    fname_list or fig, ax: list of str
        list of names of the created plots

    """
    if fig is None:
        fig = plt.figure(figsize=[10, 6], dpi=dpi)
        ax = fig.add_subplot(111)
    else:
        ax.autoscale(False)

    ax.bar(
        bin_centers, hist, width=bin_centers[1]-bin_centers[0], color=color,
        alpha=alpha)
    if invert_xaxis:
        ax.invert_xaxis()

    plt.xlabel(labelx)
    plt.ylabel(labely)
    plt.title(titl)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close()

        return fname_list

    return (fig, ax)


def plot_antenna_pattern(antpattern, fname_list, labelx='Angle [Deg]',
                         linear=False, twoway=False, title='Antenna Pattern',
                         ymin=None, ymax=None, dpi=72):
    """
    plots an antenna pattern

    Parameters
    ----------
    antpattern : dict
        dictionary with the angle and the attenuation
    value : float array
        values of the time series
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    linear : boolean
        if true data is in linear units
    linear : boolean
        if true data represents the two way attenuation
    titl : str
        The figure title
    ymin, ymax: float
        Lower/Upper limit of y axis
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    waystr = 'One-way '
    if twoway:
        waystr = 'Two-way '

    linstr = 'Att [dB]'
    if linear:
        linstr = 'Att [-]'
        mainbeam = (
            antpattern['angle'][antpattern['attenuation'] >= 10.**(-0.3)])
    else:
        mainbeam = antpattern['angle'][antpattern['attenuation'] >= -3.]
    beamwidth = np.abs(np.max(mainbeam)-np.min(mainbeam))

    labely = waystr+linstr

    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    ax.plot(antpattern['angle'], antpattern['attenuation'])
    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=ymin, top=ymax)

    metadata = '3-dB beamwidth: '+'{:.2f}'.format(float(beamwidth))
    ax.text(0.05, 0.95, metadata, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list


def plot_scatter_comp(value1, value2, fname_list, labelx='Sensor 1',
                      labely='Sensor 2', titl='Scatter', axis=None,
                      metadata=None, dpi=72):
    """
    plots the scatter between two time series

    Parameters
    ----------
    value1 : float array
        values of the first time series
    value2 : float array
        values of the second time series
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    titl : str
        The figure title
    axis : str
        type of axis
    metadata : string
        a string containing metadata
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    max_value = np.ma.max([np.max(value1), np.max(value2)])

    fig, ax = plt.subplots(figsize=[7, 7], dpi=dpi)

    ax.plot(value1, value2, 'bx')
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_title(titl)

    if axis == 'equal':
        ax.axis([0, max_value, 0, max_value])
        ax.plot([0, max_value], [0, max_value], 'k--')
        ax.set(adjustable='box-forced', aspect='equal')

    if metadata is not None:
        plt.text(0.05, 0.95, metadata, horizontalalignment='left',
                 verticalalignment='top', transform=ax.transAxes)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
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
    dpi = 72
    if 'dpi' in prdcfg['sunhitsImageConfig']:
        dpi = prdcfg['sunhitsImageConfig']['dpi']

    fig = plt.figure(figsize=[prdcfg['sunhitsImageConfig']['xsize'],
                              prdcfg['sunhitsImageConfig']['ysize']],
                     dpi=dpi)
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

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list
