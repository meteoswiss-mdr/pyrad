"""
pyrad.graph.plots_spectra
=========================

Functions to plot spectral data

.. autosummary::
    :toctree: generated/

    plot_range_Doppler
    plot_Doppler
    plot_complex_range_Doppler
    plot_amp_phase_range_Doppler
    plot_complex_Doppler
    plot_amp_phase_Doppler

"""

import numpy as np

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.pyplot as plt

import pyart

from .plots_aux import get_colobar_label, get_norm
from .plots_aux import generate_complex_range_Doppler_title
from .plots_aux import generate_complex_Doppler_title


def plot_range_Doppler(spectra, field_name, ray, prdcfg, fname_list,
                       xaxis_info='Doppler_velocity', titl=None,
                       clabel=None, vmin=None, vmax=None):
    """
    Makes a range-Doppler plot

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra to plot
    field_name : str
        name of the field to plot
    ray : int
        ray index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity' or 'Doppler_frequency'
    titl : str or None
        The plot title
    clabel : str or None
        The color bar label
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'
    xres = np.abs(xaxis[1]-xaxis[0])

    yaxis = spectra.range['data']/1000.
    yres = np.abs(yaxis[1]-yaxis[0])
    ylabel = 'range (km)'

    xaxis_lim = np.append(xaxis-xres/2, xaxis[-1]+xres/2)
    yaxis_lim = np.append(yaxis-yres/2, yaxis[-1]+yres/2)

    field_2D = spectra.fields[field_name]['data'][ray, :, 0:xaxis.size]

    X, Y = np.meshgrid(xaxis_lim, yaxis_lim)

    # display data
    norm = None
    cmap = None
    ticks = None
    ticklabs = None
    if field_name is not None:
        field_dict = pyart.config.get_metadata(field_name)
        if clabel is None:
            clabel = get_colobar_label(field_dict, field_name)

        cmap = pyart.config.get_field_colormap(field_name)

        norm, ticks, ticklabs = get_norm(field_name)
        if vmin is None or vmax is None:
            vmin = vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                vmin, vmax = pyart.config.get_field_limits(field_name)
        else:
            norm = None
    else:
        if clabel is None:
            clabel = 'value'
        if vmin is None:
            vmin = np.ma.min(field_2D)
        if vmax is None:
            vmax = np.ma.max(field_2D)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_range_Doppler_title(
            spectra, field_name, ray, datetime_format=None)

    ax = fig.add_subplot(111)
    cax = ax.pcolormesh(
        X, Y, field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titl)

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_Doppler(spectra, field_name, ray, rng, prdcfg, fname_list,
                 xaxis_info='Doppler_velocity', ylabel=None,
                 titl=None, vmin=None, vmax=None):
    """
    Makes a Doppler plot

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra to plot
    field_name : str
        name of the field to plot
    ray, rng : int
        ray and rng index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity' or 'Doppler_frequency'
    ylabel : str or None
        The label of the y-axis
    titl : str or None
        The plot title
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'

    field = spectra.fields[field_name]['data'][ray, rng, 0:xaxis.size]

    # display data
    if field_name is not None:
        if ylabel is None:
            field_dict = pyart.config.get_metadata(field_name)
            ylabel = get_colobar_label(field_dict, field_name)
        if vmin is None or vmax is None:
            vmin, vmax = pyart.config.get_field_limits(field_name)
    else:
        if vmin is None:
            vmin = np.ma.min(field)
        if vmax is None:
            vmax = np.ma.max(field)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, ray, rng, datetime_format=None)

    ax = fig.add_subplot(111)
    ax.plot(xaxis, field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title(titl)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_complex_range_Doppler(spectra, field_name, ray, prdcfg, fname_list,
                               xaxis_info='Doppler_velocity', titl=None,
                               clabel=None, vmin=None, vmax=None):
    """
    Makes a complex range-Doppler plot. Plotting separately the real and the
    imaginary part

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra to plot
    field_name : str
        name of the field to plot
    ray : int
        ray index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity' or 'Doppler_frequency'
    titl : str or None
        The plot title
    clabel : str or None
        The label of color bar
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'
    xres = np.abs(xaxis[1]-xaxis[0])

    yaxis = spectra.range['data']/1000.
    yres = np.abs(yaxis[1]-yaxis[0])
    ylabel = 'range (km)'

    xaxis_lim = np.append(xaxis-xres/2, xaxis[-1]+xres/2)
    yaxis_lim = np.append(yaxis-yres/2, yaxis[-1]+yres/2)

    field_2D = spectra.fields[field_name]['data'][ray, :, 0:xaxis.size]

    re_field_2D = np.real(field_2D)
    im_field_2D = np.imag(field_2D)

    X, Y = np.meshgrid(xaxis_lim, yaxis_lim)

    # display data
    norm = None
    cmap = None
    ticks = None
    ticklabs = None
    if field_name is not None:
        field_dict = pyart.config.get_metadata(field_name)
        if clabel is None:
            clabel = get_colobar_label(field_dict, field_name)

        cmap = pyart.config.get_field_colormap(field_name)

        norm, ticks, ticklabs = get_norm(field_name)
        if vmin is None or vmax is None:
            vmin = vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                vmin, vmax = pyart.config.get_field_limits(field_name)
        else:
            norm = None
    else:
        if clabel is None:
            clabel = 'value'
        if vmin is None:
            vmin = np.ma.min([np.ma.min(re_field_2D), np.ma.min(im_field_2D)])
        if vmax is None:
            vmax = np.ma.max([np.ma.max(re_field_2D), np.ma.max(im_field_2D)])

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_range_Doppler_title(
            spectra, field_name, ray, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolormesh(
        X, Y, re_field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Real part')

    ax = fig.add_subplot(122)
    cax = ax.pcolormesh(
        X, Y, im_field_2D, cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Imaginary part')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_amp_phase_range_Doppler(spectra, field_name, ray, prdcfg, fname_list,
                                 xaxis_info='Doppler_velocity', titl=None,
                                 ampli_vmin=None, ampli_vmax=None,
                                 phase_vmin=None, phase_vmax=None,
                                 clabel=None):
    """
    Makes a complex range-Doppler plot plotting separately the module and the
    phase of the signal

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra to plot
    field_name : str
        name of the field to plot
    ray : int
        ray index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity' or 'Doppler_frequency'
    titl : str or None
        The plot title
    ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
        The value limits
    clabel : str or None
        The label of color bar

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'
    xres = np.abs(xaxis[1]-xaxis[0])

    yaxis = spectra.range['data']/1000.
    yres = np.abs(yaxis[1]-yaxis[0])
    ylabel = 'range (km)'

    xaxis_lim = np.append(xaxis-xres/2, xaxis[-1]+xres/2)
    yaxis_lim = np.append(yaxis-yres/2, yaxis[-1]+yres/2)

    field_2D = spectra.fields[field_name]['data'][ray, :, 0:xaxis.size]

    ampli_field_2D = np.ma.abs(field_2D)
    phase_field_2D = np.ma.angle(field_2D, deg=True)

    X, Y = np.meshgrid(xaxis_lim, yaxis_lim)

    # display data
    norm = None
    cmap = None
    ticks = None
    ticklabs = None
    if field_name is not None:
        field_dict = pyart.config.get_metadata(field_name)
        if clabel is None:
            clabel = get_colobar_label(field_dict, field_name)

        cmap = pyart.config.get_field_colormap(field_name)

        norm, ticks, ticklabs = get_norm(field_name)
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin = ampli_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                ampli_vmin, ampli_vmax = pyart.config.get_field_limits(
                    field_name)
        else:
            norm = None

        if phase_vmin is None or phase_vmax is None:
            phase_vmin = phase_vmax = None
            if norm is None:  # if norm is set do not override with vmin/vmax
                phase_vmin, phase_vmax = (-180., 180.)
        else:
            norm = None
    else:
        if clabel is None:
            clabel = 'value'
        if ampli_vmin is None:
            ampli_vmin = np.ma.min(ampli_field_2D)
        if ampli_vmax is None:
            ampli_vmax = np.ma.max(ampli_field_2D)
        if phase_vmin is None:
            phase_vmin = np.ma.min(phase_field_2D)
        if phase_vmax is None:
            phase_vmax = np.ma.max(phase_field_2D)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_range_Doppler_title(
            spectra, field_name, ray, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    cax = ax.pcolormesh(
        X, Y, ampli_field_2D, cmap=cmap, vmin=ampli_vmin, vmax=ampli_vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Amplitude')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    ax = fig.add_subplot(122)
    cax = ax.pcolormesh(
        X, Y, phase_field_2D, cmap=cmap, vmin=phase_vmin, vmax=phase_vmax,
        norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Phase')

    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    cb.set_label(clabel)

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_complex_Doppler(spectra, field_name, ray, rng, prdcfg, fname_list,
                         xaxis_info='Doppler_velocity', ylabel=None,
                         titl=None, vmin=None, vmax=None):
    """
    Makes a complex Doppler plot plotting separately the real and the
    imaginary parts

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra to plot
    field_name : str
        name of the field to plot
    ray, rng : int
        ray and range index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity' or 'Doppler_frequency'
    ylabel : str or None
        The label of the y-axis
    titl : str or None
        The plot title
    vmin, vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'

    field = spectra.fields[field_name]['data'][ray, rng, 0:xaxis.size]

    re_field = np.real(field)
    im_field = np.imag(field)

    # display data
    if field_name is not None:
        if ylabel is None:
            field_dict = pyart.config.get_metadata(field_name)
            ylabel = get_colobar_label(field_dict, field_name)
        if vmin is None or vmax is None:
            vmin, vmax = pyart.config.get_field_limits(field_name)
    else:
        if vmin is None:
            vmin = np.ma.min([np.ma.min(re_field), np.ma.min(im_field)])
        if vmax is None:
            vmax = np.ma.max([np.ma.max(re_field), np.ma.max(im_field)])

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, ray, rng, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    ax.plot(xaxis, re_field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title('Real part')

    ax = fig.add_subplot(122)
    ax.plot(xaxis, im_field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=vmin, top=vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title('Imaginary part')

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list


def plot_amp_phase_Doppler(spectra, field_name, ray, rng, prdcfg, fname_list,
                           xaxis_info='Doppler_velocity', ylabel=None,
                           titl=None, ampli_vmin=None, ampli_vmax=None,
                           phase_vmin=None, phase_vmax=None):
    """
    Makes a complex Doppler plot plotting separately the module and the phase
    of the signal

    Parameters
    ----------
    spectra : radar spectra object
        object containing the spectra to plot
    field_name : str
        name of the field to plot
    ray, rng : int
        ray and range index
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    xaxis_info : str
        Type of x-axis. Can be 'Doppler_velocity' or 'Doppler_frequency'
    ylabel : str or None
        The label of the y-axis
    titl : str or None
        The plot title
    ampli_vmin, ampli_vmax, phase_vmin, phase_vmax : float or None
        The value limits

    Returns
    -------
    fname_list : list of str
        list of names of the saved plots

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']

    if xaxis_info == 'Doppler_velocity':
        xaxis = spectra.Doppler_velocity['data'][ray, :].compressed()
        xlabel = 'Doppler velocity (m/s)'
    elif xaxis_info == 'Doppler_frequency':
        xaxis = spectra.Doppler_frequency['data'][ray, :].compressed()
        xlabel = 'Doppler frequency (Hz)'

    field = spectra.fields[field_name]['data'][ray, rng, 0:xaxis.size]

    ampli_field = np.ma.abs(field)
    phase_field = np.ma.angle(field, deg=True)

    # display data
    if field_name is not None:
        if ylabel is None:
            field_dict = pyart.config.get_metadata(field_name)
            ylabel = get_colobar_label(field_dict, field_name)
        if ampli_vmin is None or ampli_vmax is None:
            ampli_vmin, ampli_vmax = pyart.config.get_field_limits(
                field_name)
        if phase_vmin is None or phase_vmax is None:
            phase_vmin, phase_vmax = (-180., 180.)
    else:
        if ampli_vmin is None:
            ampli_vmin = np.ma.min(ampli_field)
        if ampli_vmax is None:
            ampli_vmax = np.ma.max(ampli_field)
        if phase_vmin is None:
            phase_vmin = np.ma.min(phase_field)
        if phase_vmax is None:
            phase_vmax = np.ma.max(phase_field)

    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    if titl is None:
        titl = generate_complex_Doppler_title(
            spectra, field_name, ray, rng, datetime_format=None)

    plt.suptitle(titl)

    ax = fig.add_subplot(121)
    ax.plot(xaxis, ampli_field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=ampli_vmin, top=ampli_vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title('Amplitude')

    ax = fig.add_subplot(122)
    ax.plot(xaxis, phase_field, marker='x')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=phase_vmin, top=phase_vmax)
    ax.set_xlim([xaxis[0], xaxis[-1]])
    ax.set_title('Phase')

    # Make a tight layout
    fig.tight_layout()

    plt.subplots_adjust(top=0.8)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close(fig)

    return fname_list
