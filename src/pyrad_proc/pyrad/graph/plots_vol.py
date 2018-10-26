"""
pyrad.graph.plots
=================

Functions to plot radar volume data

.. autosummary::
    :toctree: generated/

    plot_ppi
    plot_ppi_map
    plot_rhi
    plot_bscope
    plot_time_range
    plot_cappi
    plot_traj
    plot_rhi_contour
    plot_ppi_contour
    plot_pos
    plot_rhi_profile
    plot_along_coord
    plot_field_coverage
    _plot_time_range

"""
from warnings import warn

import numpy as np

import matplotlib as mpl
mpl.use('Agg')

# Increase a bit font size
mpl.rcParams.update({'font.size': 16})
mpl.rcParams.update({'font.family':  "sans-serif"})

import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable

import pyart

from .plots_aux import get_colobar_label, get_norm
from .plots import plot_quantiles, plot_histogram

from ..util.radar_utils import compute_quantiles_sweep
from ..util.radar_utils import compute_histogram_sweep


def plot_ppi(radar, field_name, ind_el, prdcfg, fname_list, plot_type='PPI',
             titl=None, step=None, quantiles=None, save_fig=True):
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
    titl : str
        Plot title
    step : float
        step for histogram plotting
    quantiles : float array
        quantiles to plot
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure
        

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    History
    --------
    ????.??.?? created
    2017.08.?? -fvj- added option controlling dpi
    2017-08.23 -jgr- minor graphical changes

    """
    if plot_type == 'PPI':
        dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

        norm, ticks, ticklabs = get_norm(field_name)

        xsize = prdcfg['ppiImageConfig']['xsize']
        ysize = prdcfg['ppiImageConfig']['ysize']
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)

        ax = fig.add_subplot(111, aspect='equal')

        display = pyart.graph.RadarDisplay(radar)
        display.plot_ppi(
            field_name, title=titl, sweep=ind_el, norm=norm, ticks=ticks,
            ticklabs=ticklabs, fig=fig)

        display.set_limits(
            ylim=[prdcfg['ppiImageConfig']['ymin'],
                  prdcfg['ppiImageConfig']['ymax']],
            xlim=[prdcfg['ppiImageConfig']['xmin'],
                  prdcfg['ppiImageConfig']['xmax']])
        if 'rngRing' in prdcfg['ppiImageConfig']:
            if prdcfg['ppiImageConfig']['rngRing'] > 0:
                display.plot_range_rings(np.arange(
                    0., radar.range['data'][-1]/1000.,
                    prdcfg['ppiImageConfig']['rngRing']))
        display.plot_cross_hair(5.)

        # Turn on the grid
        ax.grid()

        # Make a tight layout
        fig.tight_layout()

        if save_fig:
            for fname in fname_list:
                fig.savefig(fname, dpi=dpi)
            plt.close()

            return fname_list

        return (fig, ax)

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


def plot_ppi_map(radar, field_name, ind_el, prdcfg, fname_list):
    """
    plots a PPI on a geographic map

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

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    History
    --------
    ????.??.?? created
    2017.08.?? -fvj- added option controlling dpi
    2017-08.22 -jgr- changed colortable behavior: created now here
                     instead than on pyart

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    norm, ticks, ticklabs = get_norm(field_name)

    xsize = prdcfg['ppiImageConfig']['xsize']
    ysize = prdcfg['ppiImageConfig']['ysize']
    fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    lon_lines = np.arange(np.floor(prdcfg['ppiMapImageConfig']['lonmin']),
                          np.ceil(prdcfg['ppiMapImageConfig']['lonmax'])+1,
                          0.5)
    lat_lines = np.arange(np.floor(prdcfg['ppiMapImageConfig']['latmin']),
                          np.ceil(prdcfg['ppiMapImageConfig']['latmax'])+1,
                          0.5)

    display_map = pyart.graph.RadarMapDisplayCartopy(radar)
    display_map.plot_ppi_map(
        field_name, sweep=ind_el, norm=norm, ticks=ticks,
        ticklabs=ticklabs, min_lon=prdcfg['ppiMapImageConfig']['lonmin'],
        max_lon=prdcfg['ppiMapImageConfig']['lonmax'],
        min_lat=prdcfg['ppiMapImageConfig']['latmin'],
        max_lat=prdcfg['ppiMapImageConfig']['latmax'],
        resolution=prdcfg['ppiMapImageConfig']['mapres'],
        lat_lines=lat_lines, lon_lines=lon_lines,
        maps_list=prdcfg['ppiMapImageConfig']['maps'],
        colorbar_flag=False)

    if 'rngRing' in prdcfg['ppiMapImageConfig']:
        if prdcfg['ppiMapImageConfig']['rngRing'] > 0:
            rng_rings = np.arange(
                0., radar.range['data'][-1]/1000.,
                prdcfg['ppiMapImageConfig']['rngRing'])
            for rng_ring in rng_rings:
                display_map.plot_range_ring(rng_ring)

    # Adapt the axes of the colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)

    display_map.plot_colorbar(mappable=display_map.plots[0],
                              cax=cax, field=field_name, ax=ax)

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi, bbox_inches='tight')

    plt.close()

    return fname_list


def plot_rhi(radar, field_name, ind_az, prdcfg, fname_list, plot_type='RHI',
             titl=None, step=None, quantiles=None, save_fig=True):
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
    titl : str
        Plot title
    step : float
        step for histogram plotting
    quantiles : float array
        quantiles to plot
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str
        list of names of the created plots
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

     History
    --------
    ????.??.?? created
    2017.08.?? -fvj- added option controlling dpi
    2017-08.23 -jgr- minor graphical changes

    """
    if plot_type == 'RHI':
        dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

        norm, ticks, ticklabs = get_norm(field_name)

        xsize = prdcfg['rhiImageConfig']['xsize']
        ysize = prdcfg['rhiImageConfig']['ysize']
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')
        display = pyart.graph.RadarDisplay(radar)
        display.plot_rhi(
            field_name, title=titl, sweep=ind_az, norm=norm, ticks=ticks,
            ticklabs=ticklabs, colorbar_orient='horizontal',
            reverse_xaxis=False)
        display.set_limits(
            ylim=[prdcfg['rhiImageConfig']['ymin'],
                  prdcfg['rhiImageConfig']['ymax']],
            xlim=[prdcfg['rhiImageConfig']['xmin'],
                  prdcfg['rhiImageConfig']['xmax']])
        display.plot_cross_hair(5.)

        # Turn on the grid
        ax.grid()

        # Make a tight layout
        fig.tight_layout()

        if save_fig:
            for fname in fname_list:
                fig.savefig(fname, dpi=dpi)
            plt.close()

            return fname_list

        return (fig, ax)

    elif plot_type == 'QUANTILES':
        quantiles, values = compute_quantiles_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_az],
            radar.sweep_end_ray_index['data'][ind_az], quantiles=quantiles)

        if titl is None:
            titl = pyart.graph.common.generate_title(
                radar, field_name, ind_az)
        labely = get_colobar_label(radar.fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname_list, labelx='quantile',
                       labely=labely, titl=titl)

    elif plot_type == 'HISTOGRAM':
        bins, values = compute_histogram_sweep(
            radar.fields[field_name]['data'],
            radar.sweep_start_ray_index['data'][ind_az],
            radar.sweep_end_ray_index['data'][ind_az], field_name, step=step)

        if titl is None:
            titl = pyart.graph.common.generate_title(
                radar, field_name, ind_az)
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
    norm, ticks, ticklabs = get_norm(field_name)

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

    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=dpi)
    ax = fig.add_subplot(111)
    if radar_aux.ngates == 1:
        plt.plot(ang, field, 'bx')
        plt.xlabel(labely)
        plt.ylabel(label)
        plt.title(titl)
    else:
        cmap = pyart.config.get_field_colormap(field_name)

        vmin = vmax = None
        if norm is None:  # if norm is set do not override with vmin/vmax
            vmin, vmax = pyart.config.get_field_limits(field_name)

        rmin = radar_aux.range['data'][0]/1000.
        rmax = radar_aux.range['data'][-1]/1000.
        cax = ax.imshow(
            field, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, norm=norm,
            extent=(rmin, rmax, ang_min, ang_max), aspect='auto',
            interpolation='none')
        plt.xlabel('Range (km)')
        plt.ylabel(labely)
        plt.title(titl)

        cb = fig.colorbar(cax)
        if ticks is not None:
            cb.set_ticks(ticks)
        if ticklabs:
            cb.set_ticklabels(ticklabs)
        cb.set_label(label)

    # Make a tight layout
    fig.tight_layout()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list


def plot_time_range(radar, field_name, ind_sweep, prdcfg, fname_list):
    """
    plots a time-range plot

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
    field = radar_aux.fields[field_name]['data']

    # display data
    titl = pyart.graph.common.generate_title(radar_aux, field_name, ind_sweep)
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    xsize = prdcfg['ppiImageConfig'].get('xsize', 10)
    ysize = prdcfg['ppiImageConfig'].get('ysize', 8)

    rng_aux = radar_aux.range['data']/1000.
    rng_res = rng_aux[1]-rng_aux[0]
    rng_aux = np.append(rng_aux-rng_res/2., rng_aux[-1]+rng_res/2.)

    time_res = np.mean(radar_aux.time['data'][1:]-radar_aux.time['data'][0:-1])
    time_aux = np.append(radar_aux.time['data'], radar_aux.time['data'][-1]+time_res)
    return _plot_time_range(
        time_aux, rng_aux, field, field_name, fname_list, titl=titl,
        figsize=[xsize, ysize], dpi=dpi)


def plot_cappi(radar, field_name, altitude, prdcfg, fname_list,
               save_fig=True):
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
    save_fig : bool
        if true save the figure. If false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    norm, ticks, ticklabs = get_norm(field_name)

    xmin = prdcfg['ppiImageConfig']['xmin']
    xmax = prdcfg['ppiImageConfig']['xmax']
    ymin = prdcfg['ppiImageConfig']['ymin']
    ymax = prdcfg['ppiImageConfig']['ymax']

    wfunc = prdcfg.get('wfunc', 'NEAREST_NEIGHBOUR')
    cappi_res = prdcfg.get('res', 500.)

    # number of grid points in cappi
    ny = int((ymax-ymin)*1000./cappi_res)+1
    nx = int((xmax-xmin)*1000./cappi_res)+1

    # parameters to determine the gates to use for each grid point
    beamwidth = 1.
    beam_spacing = 1.
    if 'radar_beam_width_h' in radar.instrument_parameters:
        beamwidth = radar.instrument_parameters[
            'radar_beam_width_h']['data'][0]

    if radar.ray_angle_res is not None:
        beam_spacing = radar.ray_angle_res['data'][0]

    lat = float(radar.latitude['data'])
    lon = float(radar.longitude['data'])
    alt = 0.

    # cartesian mapping
    grid = pyart.map.grid_from_radars(
        (radar,), gridding_algo='map_to_grid', weighting_function=wfunc,
        roi_func='dist_beam', h_factor=1.0, nb=beamwidth, bsp=beam_spacing,
        min_radius=cappi_res/2.,
        grid_shape=(1, ny, nx),
        grid_limits=((altitude, altitude), (ymin*1000., ymax*1000.),
                     (xmin*1000., xmax*1000.)),
        grid_origin=(lat, lon), grid_origin_alt=alt,
        fields=[field_name])

    # display data
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                              prdcfg['ppiImageConfig']['ysize']],
                     dpi=dpi)
    ax = fig.add_subplot(111, aspect='equal')
    cmap = pyart.config.get_field_colormap(field_name)

    vmin = vmax = None
    if norm is None:  # if norm is set do not override with vmin/vmax
        vmin, vmax = pyart.config.get_field_limits(field_name)

    titl = pyart.graph.common.generate_grid_title(grid, field_name, 0)

    cax = ax.imshow(
        grid.fields[field_name]['data'][0], extent=(xmin, xmax, ymin, ymax),
        origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, norm=norm,
        interpolation='none')
    plt.xlabel('East West distance from radar(km)')
    plt.ylabel('North South distance from radar(km)')
    plt.title(titl)

    # plot the colorbar and set the label.
    cb = fig.colorbar(cax)
    if ticks is not None:
        cb.set_ticks(ticks)
    if ticklabs:
        cb.set_ticklabels(ticklabs)
    label = get_colobar_label(grid.fields[field_name], field_name)
    cb.set_label(label)

    # Make a tight layout
    fig.tight_layout()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close()

        return fname_list

    return (fig, ax)


def plot_traj(rng_traj, azi_traj, ele_traj, time_traj, prdcfg, fname_list,
              rad_alt=None, rad_tstart=None, ax=None, fig=None,
              save_fig=True):
    """
    plots a trajectory on a Cartesian surface

    Parameters
    ----------
    rng_traj, azi_traj, ele_traj : float array
        antenna coordinates of the trajectory [m and deg]
    time_traj : datetime array
        trajectory time
    prdcfg : dict
        dictionary containing the product configuration
    fname_list : list of str
        list of names of the files where to store the plot
    rad_alt : float or None
        radar altitude [m MSL]
    rad_tstart : datetime object or None
        start time of the radar scan
    surface_alt : float
        surface altitude [m MSL]
    color_ref : str
        What the color code represents. Can be 'None', 'rel_altitude',
        'altitude' or 'time'
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    color_ref = prdcfg.get('color_ref', 'None')
    if 'altitude' not in prdcfg and color_ref == 'rel_altitude':
        warn('Unable to plot trajectory relative to surface altitude. ' +
             'Unknown surface altitude.')
        color_ref = 'None'
    if rad_tstart is None and color_ref == 'time':
        warn('Unable to plot trajectory relative to radar scan start time. ' +
             'Unknown radar scan start time.')
        color_ref = 'None'
    if (rad_alt is None and
            (color_ref == 'rel_altitude' or color_ref == 'altitude')):
        warn('Unable to plot trajectory altitude. ' +
             'Unknown radar altitude.')
        color_ref = 'None'

    x, y, z = pyart.core.antenna_to_cartesian(
        rng_traj/1000., azi_traj, ele_traj)

    if color_ref == 'rel_altitude':
        h = z+rad_alt
        h_rel = h-prdcfg['altitude']

        marker = 'x'
        col = h_rel
        cmap = 'coolwarm'
        norm = plt.Normalize(-2000., 2000.)
        cb_label = 'Altitude relative to CAPPI [m]'
        plot_cb = True
    elif color_ref == 'altitude':
        h = z+rad_alt

        marker = 'x'
        col = h
        cmap = 'Greys'
        norm = plt.Normalize(h.min(), h.max())
        cb_label = 'Altitude [m MSL]'
        plot_cb = True
    elif color_ref == 'time':
        td_vec = time_traj-rad_tstart
        tt_s = []
        for td in td_vec:
            tt_s.append(td.total_seconds())
        tt_s = np.asarray(tt_s)

        marker = 'x'
        col = tt_s
        cmap = 'Greys'
        norm = plt.Normalize(tt_s.min(), tt_s.max())
        cb_label = 'Time from start of radar scan [s]'
        plot_cb = True
    else:
        col = 'k'
        marker = 'x'
        cmap = None
        plot_cb = False

    # display data
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)
    if fig is None:
        fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                                  prdcfg['ppiImageConfig']['ysize']],
                         dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')
    else:
        ax.autoscale(False)

    cax = ax.scatter(
        x/1000., y/1000., c=col, marker=marker, alpha=0.5, cmap=cmap,
        norm=norm)

    # plot colorbar
    if plot_cb:
        cb = fig.colorbar(cax, orientation='horizontal')
        cb.set_label(cb_label)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close()

        return fname_list

    return (fig, ax)


def plot_rhi_contour(radar, field_name, ind_az, prdcfg, fname_list,
                     contour_values=None, linewidths=1.5, ax=None, fig=None,
                     save_fig=True):
    """
    plots contour data on an RHI

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
    contour_values : float array
        list of contours to plot
    linewidths : float
        width of the contour lines
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    # get contour intervals
    if contour_values is None:
        field_dict = pyart.config.get_metadata(field_name)
        if 'boundaries' in field_dict:
            vmin = field_dict['boundaries'][0]
            vmax = field_dict['boundaries'][-1]
            num = len(field_dict['boundaries'])
        else:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            num = 10

        contour_values = np.linspace(vmin, vmax, num=num)

    # get data and position
    display = pyart.graph.RadarDisplay(radar)
    data = display._get_data(field_name, ind_az, None, True, None)

    x_edges, y_edges, z_edges = display._get_x_y_z(ind_az, True, True)
    delta_x = x_edges[1:, 1:]-x_edges[:-1, :-1]
    delta_y = y_edges[1:, 1:]-y_edges[:-1, :-1]
    delta_z = z_edges[1:, 1:]-z_edges[:-1, :-1]

    x = x_edges[:-1, :-1]+delta_x/2.
    y = y_edges[:-1, :-1]+delta_y/2.
    z = z_edges[:-1, :-1]+delta_z/2.

    R = np.sqrt(x ** 2 + y ** 2) * np.sign(x)

    # display data
    if fig is None:
        xsize = prdcfg['rhiImageConfig']['xsize']
        ysize = prdcfg['rhiImageConfig']['ysize']
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')

        ax.contour(R, z, data, contour_values, colors='k',
                   linewidths=linewidths)

        display._set_title(field_name, ind_az, None, ax)
        display._label_axes_rhi((None, None), ax)

        display.set_limits(
            ylim=[prdcfg['rhiImageConfig']['ymin'],
                  prdcfg['rhiImageConfig']['ymax']],
            xlim=[prdcfg['rhiImageConfig']['xmin'],
                  prdcfg['rhiImageConfig']['xmax']])
        display.plot_cross_hair(5.)

        # Turn on the grid
        ax.grid()

        # Make a tight layout
        fig.tight_layout()
    else:
        ax.autoscale(False)
        ax.contour(R, z, data, contour_values, colors='k',
                   linewidths=linewidths)



    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close()

        return fname_list

    return (fig, ax)


def plot_ppi_contour(radar, field_name, ind_el, prdcfg, fname_list,
                     contour_values=None, linewidths=1.5, ax=None, fig=None,
                     save_fig=True):
    """
    plots contour data on a PPI

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
    contour_values : float array
        list of contours to plot
    linewidths : float
        width of the contour lines
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    dpi = prdcfg['ppiImageConfig'].get('dpi', 72)

    # get contour intervals
    if contour_values is None:
        field_dict = pyart.config.get_metadata(field_name)
        if 'boundaries' in field_dict:
            vmin = field_dict['boundaries'][0]
            vmax = field_dict['boundaries'][-1]
            num = len(field_dict['boundaries'])
        else:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            num = 10

        contour_values = np.linspace(vmin, vmax, num=num)

    # get data and position
    display = pyart.graph.RadarDisplay(radar)
    data = display._get_data(field_name, ind_el, None, True, None)

    x_edges, y_edges = display._get_x_y(ind_el, True, True)
    delta_x = x_edges[1:, 1:]-x_edges[:-1, :-1]
    delta_y = y_edges[1:, 1:]-y_edges[:-1, :-1]

    x = x_edges[:-1, :-1]+delta_x/2.
    y = y_edges[:-1, :-1]+delta_y/2.

    # display data
    if fig is None:
        xsize = prdcfg['ppiImageConfig']['xsize']
        ysize = prdcfg['ppiImageConfig']['ysize']
        fig = plt.figure(figsize=[xsize, ysize], dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')

        ax.contour(x, y, data, contour_values, colors='k',
                   linewidths=linewidths)

        display._set_title(field_name, ind_el, None, ax)
        display._label_axes_ppi((None, None), ax)

        display.set_limits(
            ylim=[prdcfg['ppiImageConfig']['ymin'],
                  prdcfg['ppiImageConfig']['ymax']],
            xlim=[prdcfg['ppiImageConfig']['xmin'],
                  prdcfg['ppiImageConfig']['xmax']])
        display.plot_cross_hair(5.)

        # Turn on the grid
        ax.grid()

        # Make a tight layout
        fig.tight_layout()
    else:
        ax.autoscale(False)
        ax.contour(x, y, data, contour_values, colors='k',
                   linewidths=linewidths)

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close()

        return fname_list

    return (fig, ax)


def plot_pos(lat, lon, alt, fname_list, ax=None, fig=None, save_fig=True,
             sort_altitude='No', dpi=72, alpha=1., cb_label='height [m MSL]',
             titl='Position', xlabel='Lon [Deg]', ylabel='Lat [Deg]',
             limits=None, vmin=None, vmax=None):
    """
    plots a trajectory on a Cartesian surface

    Parameters
    ----------
    lat, lon, alt : float array
        Points coordinates
    fname_list : list of str
        list of names of the files where to store the plot
    fig : Figure
        Figure to add the colorbar to. If none a new figure will be created
    ax : Axis
        Axis to plot on. if fig is None a new axis will be created
    save_fig : bool
        if true save the figure if false it does not close the plot and
        returns the handle to the figure
    sort_altitude : str
        String indicating whether to sort the altitude data. Can be 'No',
        'Lowest_on_top' or 'Highest_on_top'
    dpi : int
        Pixel density
    alpha : float
        Transparency
    cb_label : str
        Color bar label
    titl : str
        Plot title
    limits : tupple or None
        The limits of the field to plot

    Returns
    -------
    fname_list : list of str or
    fig, ax : tupple
        list of names of the saved plots or handle of the figure an axes

    """
    if sort_altitude == 'Lowest_on_top' or sort_altitude == 'Highest_on_top':
        ind = np.argsort(alt)
        if sort_altitude == 'Lowest_on_top':
            ind = ind[::-1]
        lat = lat[ind]
        lon = lon[ind]
        alt = alt[ind]

    if vmin is None:
        vmin = alt.min()
    if vmax is None:
        vmax = alt.max()
    marker = 'x'
    col = alt
    cmap = 'viridis'
    norm = plt.Normalize(vmin, vmax)

    if fig is None:
        fig = plt.figure(figsize=[10, 8], dpi=dpi)
        ax = fig.add_subplot(111, aspect='equal')
    else:
        ax.autoscale(False)

    cax = ax.scatter(
        lon, lat, c=col, marker=marker, alpha=alpha, cmap=cmap, norm=norm)

    if limits is not None:
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])

    # plot colorbar
    cb = fig.colorbar(cax, orientation='horizontal')
    cb.set_label(cb_label)

    plt.title(titl)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Turn on the grid
    ax.grid()

    if save_fig:
        for fname in fname_list:
            fig.savefig(fname, dpi=dpi)
        plt.close()

        return fname_list

    return (fig, ax)


def plot_rhi_profile(data_list, hvec, fname_list, labelx='Value',
                     labely='Height (m MSL)', labels=['Mean'],
                     title='RHI profile', colors=None, linestyles=None,
                     vmin=None, vmax=None, hmin=None, hmax=None, dpi=72):
    """
    plots an RHI profile

    Parameters
    ----------
    data_list : list of float array
        values of the profile
    hvec : float array
        height points of the profile
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    labels : array of str
        The label of the legend
    title : str
        The figure title
    colors : array of str
        Specifies the colors of each line
    linestyles : array of str
        Specifies the line style of each line
    vmin, vmax: float
        Lower/Upper limit of data values
    hmin, hmax: float
        Lower/Upper limit of altitude
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    lab = None
    col = None
    lstyle = None
    for i, data in enumerate(data_list):
        if labels is not None:
            lab = labels[i]
        if colors is not None:
            col = colors[i]
        if linestyles is not None:
            lstyle = linestyles[i]
        ax.plot(
            data, hvec, label=lab, color=col, linestyle=lstyle, marker='x')

    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_xlim(left=vmin, right=vmax)
    ax.set_ylim(bottom=hmin, top=hmax)
    ax.legend(loc='best')

    # Turn on the grid
    ax.grid()

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list


def plot_along_coord(xval_list, yval_list, fname_list, labelx='coord',
                     labely='Value', labels=None, title='Plot along coordinate',
                     colors=None, linestyles=None, ymin=None, ymax=None,
                     dpi=72):
    """
    plots a time series

    Parameters
    ----------
    xval_list : list of float arrays
        the x values, range, azimuth or elevation
    yval_list : list of float arrays
        the y values. Parameter to plot
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    labels : array of str
        The label of the legend
    title : str
        The figure title
    colors : array of str
        Specifies the colors of each line
    linestyles : array of str
        Specifies the line style of each line
    ymin, ymax: float
        Lower/Upper limit of y axis
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    lab = None
    col = None
    lstyle = None

    for i, xval in enumerate(xval_list):
        yval = yval_list[i]
        if labels is not None:
            lab = labels[i]
        if colors is not None:
            col = colors[i]
        if linestyles is not None:
            lstyle = linestyles[i]
        ax.plot(xval, yval, label=lab, color=col, linestyle=lstyle)

    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=ymin, top=ymax)
    ax.legend(loc='best')

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list


def plot_field_coverage(xval_list, yval_list, fname_list,
                        labelx='Azimuth (deg)', labely='Range extension [m]',
                        labels=None, title='Field coverage', ymin=None,
                        ymax=None, xmeanval=None, ymeanval=None,
                        labelmeanval=None, dpi=72):
    """
    plots a time series

    Parameters
    ----------
    xval_list : list of float arrays
        the x values, azimuth
    yval_list : list of float arrays
        the y values. Range extension
    fname_list : list of str
        list of names of the files where to store the plot
    labelx : str
        The label of the X axis
    labely : str
        The label of the Y axis
    labels : array of str
        The label of the legend
    title : str
        The figure title
    ymin, ymax : float
        Lower/Upper limit of y axis
    xmeanval, ymeanval : float array
        the x and y values of a mean along elevation
    labelmeanval : str
        the label of the mean
    dpi : int
        dots per inch

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
    fig, ax = plt.subplots(figsize=[10, 6], dpi=dpi)

    lab = None

    for i, xval in enumerate(xval_list):
        yval = yval_list[i]
        if labels is not None:
            lab = labels[i]
        ax.plot(xval, yval, label=lab, linestyle='None', marker='o',
                fillstyle='full')

    if xmeanval is not None and ymeanval is not None:
        ax.plot(xmeanval, ymeanval, label=labelmeanval, linestyle='-',
                color='r', marker='x')

    ax.set_title(title)
    ax.set_xlabel(labelx)
    ax.set_ylabel(labely)
    ax.set_ylim(bottom=ymin, top=ymax)
    if labels is not None:
        ax.legend(loc='best')

    for fname in fname_list:
        fig.savefig(fname, dpi=dpi)
    plt.close()

    return fname_list


def _plot_time_range(rad_time, rad_range, rad_data, field_name, fname_list,
                     titl='Time-Range plot',
                     xlabel='time (s from start time)', ylabel='range (Km)',
                     clabel=None,
                     vmin=None, vmax=None, figsize=[10, 8], dpi=72):
    """
    plots a time-range plot

    Parameters
    ----------
    rad_time : Radar object
        object containing the radar data to plot
    rad_range : str
        name of the radar field to plot
    rad_data : int
        sweep index to plot
    field_name : str or None
        field name. Used to define plot characteristics
    fname_list : list of str
        list of names of the files where to store the plot
    titl : str
        Plot title
    xlabel, ylabel : str
        x- and y-axis labels
    clabel : str or None
        colorbar label
    vmin, vmax : float
        min and max values of the color bar
    figsize : list
        figure size [xsize, ysize]
    dpi : int
        dpi

    Returns
    -------
    fname_list : list of str
        list of names of the created plots

    """
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
            vmin = np.ma.min(rad_data)
        if vmax is None:
            vmax = np.ma.max(rad_data)

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111)

    T, R = np.meshgrid(rad_time, rad_range)
    cax = ax.pcolormesh(
        T, R, np.ma.transpose(rad_data), cmap=cmap, vmin=vmin, vmax=vmax,
        norm=norm)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(titl)

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
    plt.close()

    return fname_list
