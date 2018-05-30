"""
pyrad.prod.process_vol_products
===============================

Functions for obtaining Pyrad products from a radar volume dataset

.. autosummary::
    :toctree: generated/

    generate_vol_products

"""

from copy import deepcopy
from warnings import warn

import numpy as np
from netCDF4 import num2date

import pyart

from ..io.io_aux import get_save_dir, make_filename, get_fieldname_pyart

from ..io.write_data import write_cdf, write_rhi_profile, write_field_coverage
from ..io.write_data import write_last_state, write_histogram, write_quantiles

from ..graph.plots_vol import plot_ppi, plot_ppi_map, plot_rhi, plot_cappi
from ..graph.plots_vol import plot_bscope, plot_rhi_profile, plot_along_coord
from ..graph.plots_vol import plot_field_coverage, plot_time_range
from ..graph.plots import plot_quantiles, plot_histogram
from ..graph.plots_aux import get_colobar_label, get_field_name

from ..util.radar_utils import get_ROI, compute_profile_stats
from ..util.radar_utils import compute_histogram, compute_quantiles
from ..util.stat_utils import quantiles_weighted


def generate_vol_products(dataset, prdcfg):
    """
    Generates radar volume products.

    Parameters
    ----------
    dataset : Radar
        radar object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'PPI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(dataset.fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(dataset.fixed_angle['data'] == el)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'ppi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        step = None
        quantiles = None
        plot_type = 'PPI'
        if 'plot_type' in prdcfg:
            plot_type = prdcfg['plot_type']
        if 'step' in prdcfg:
            step = prdcfg['step']
        if 'quantiles' in prdcfg:
            quantiles = prdcfg['quantiles']

        plot_ppi(dataset, field_name, ind_el, prdcfg, fname_list,
                 plot_type=plot_type, step=step, quantiles=quantiles)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'PPI_MAP':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(dataset.fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(dataset.fixed_angle['data'] == el)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'ppi_map', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_ppi_map(dataset, field_name, ind_el, prdcfg, fname_list)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'RHI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        az_vec = np.sort(dataset.fixed_angle['data'])
        az = az_vec[prdcfg['anglenr']]
        ind_az = np.where(dataset.fixed_angle['data'] == az)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'rhi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='az'+'{:.1f}'.format(az),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        step = None
        quantiles = None
        plot_type = 'RHI'
        if 'plot_type' in prdcfg:
            plot_type = prdcfg['plot_type']
        if 'step' in prdcfg:
            step = prdcfg['step']
        if 'quantiles' in prdcfg:
            quantiles = prdcfg['quantiles']

        plot_rhi(dataset, field_name, ind_az, prdcfg, fname_list,
                 plot_type=plot_type, step=step, quantiles=quantiles)

        print('----- save to '+' '.join(fname_list))

        return fname_list


    elif prdcfg['type'] == 'PSEUDOPPI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        try:
            xsect = pyart.util.cross_section_rhi(
                dataset, [prdcfg['angle']], el_tol=prdcfg['EleTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'ppi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='el'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            step = None
            quantiles = None
            plot_type = 'PPI'
            if 'plot_type' in prdcfg:
                plot_type = prdcfg['plot_type']
            if 'step' in prdcfg:
                step = prdcfg['step']
            if 'quantiles' in prdcfg:
                quantiles = prdcfg['quantiles']

            plot_ppi(xsect, field_name, 0, prdcfg, fname_list,
                     plot_type=plot_type, step=step, quantiles=quantiles)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                'No data found at elevation ' + str(prdcfg['angle']) +
                '. Skipping product ' + prdcfg['type'])

            return None

    elif prdcfg['type'] == 'RHI_PROFILE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined parameters
        rangeStart = prdcfg.get('rangeStart', 0.)
        rangeStop = prdcfg.get('rangeStop', 25000.)
        heightResolution = prdcfg.get('heightResolution', 500.)
        hmin_user = prdcfg.get('heightMin', None)
        hmax_user = prdcfg.get('heightMax', None)
        quantity = prdcfg.get('quantity', 'quantiles')
        quantiles = prdcfg.get('quantiles', np.array([25., 50., 75.]))
        nvalid_min = prdcfg.get('nvalid_min', 4)

        fixed_span = prdcfg.get('fixed_span', 1)
        vmin = None
        vmax = None
        if fixed_span:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            if 'vmin' in prdcfg:
                vmin = prdcfg['vmin']
            if 'vmax' in prdcfg:
                vmax = prdcfg['vmax']

        # create new radar object with only data for the given rhi and range
        az_vec = np.sort(dataset.fixed_angle['data'])
        az = az_vec[prdcfg['anglenr']]
        ind_az = np.where(dataset.fixed_angle['data'] == az)[0][0]

        new_dataset = dataset.extract_sweeps([ind_az])
        field = new_dataset.fields[field_name]
        rng_mask = np.logical_and(new_dataset.range['data'] >= rangeStart,
                                  new_dataset.range['data'] <= rangeStop)
        field['data'] = field['data'][:, rng_mask]
        new_dataset.range['data'] = new_dataset.range['data'][rng_mask]
        new_dataset.ngates = len(new_dataset.range['data'])
        new_dataset.init_gate_x_y_z()
        new_dataset.init_gate_longitude_latitude()
        new_dataset.init_gate_altitude()

        new_dataset.fields = dict()
        new_dataset.add_field(field_name, field)

        # compute quantities
        if hmin_user is None:
            minheight = (
                round(np.min(dataset.gate_altitude['data']) /
                      heightResolution)*heightResolution-heightResolution)
        else:
            minheight = hmin_user
        if hmax_user is None:
            maxheight = (
                round(np.max(dataset.gate_altitude['data']) /
                      heightResolution)*heightResolution+heightResolution)
        else:
            maxheight = hmax_user
        nlevels = int((maxheight-minheight)/heightResolution)

        h_vec = minheight+np.arange(nlevels)*heightResolution+heightResolution/2.
        vals, val_valid = compute_profile_stats(
            field['data'], new_dataset.gate_altitude['data'], h_vec,
            heightResolution, quantity=quantity, quantiles=quantiles/100.,
            nvalid_min=nvalid_min)

        # plot data
        if quantity == 'mean':
            data = [vals[:, 0], vals[:, 1], vals[:, 2]]
            labels = ['Mean', 'Min', 'Max']
            colors = ['b', 'k', 'k']
            linestyles = ['-', '--', '--']
        elif quantity == 'mode':
            data = [vals[:, 0], vals[:, 2], vals[:, 4]]
            labels = ['Mode', '2nd most common', '3rd most common']
            colors = ['b', 'k', 'r']
            linestyles = ['-', '--', '--']
        else:
            data = [vals[:, 1], vals[:, 0], vals[:, 2]]
            labels = [
                str(quantiles[1])+'-percentile',
                str(quantiles[0])+'-percentile',
                str(quantiles[2])+'-percentile']
            colors = ['b', 'k', 'k']
            linestyles = ['-', '--', '--']

        labelx = get_colobar_label(dataset.fields[field_name], field_name)
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'az'+'{:.1f}'.format(az)+'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'rhi_profile', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles)

        print('----- save to '+' '.join(fname_list))

        fname = make_filename(
            'rhi_profile', prdcfg['dstype'], prdcfg['voltype'],
            ['csv'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'])[0]

        fname = savedir+fname

        if quantity == 'mode':
            data.append(vals[:, 1])
            labels.append('% points mode')
            data.append(vals[:, 3])
            labels.append('% points 2nd most common')
            data.append(vals[:, 5])
            labels.append('% points 3rd most common')

        sector = {
            'rmin': rangeStart,
            'rmax': rangeStop,
            'az': az
        }
        write_rhi_profile(
            h_vec, data, val_valid, labels, fname, datatype=labelx,
            timeinfo=prdcfg['timeinfo'], sector=sector)

        print('----- save to '+fname)

        # TODO: add Cartesian interpolation option

        return fname

    elif prdcfg['type'] == 'PROFILE_STATS':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # mask unclassified data
        field = deepcopy(dataset.fields[field_name]['data'])
        if prdcfg['voltype'] == 'hydro':
            field = np.ma.masked_equal(field, 0)

        # user defined parameters
        heightResolution = prdcfg.get('heightResolution', 100.)
        hmin_user = prdcfg.get('heightMin', None)
        hmax_user = prdcfg.get('heightMax', None)
        quantity = prdcfg.get('quantity', 'quantiles')
        quantiles = prdcfg.get('quantiles', np.array([25., 50., 75.]))
        nvalid_min = prdcfg.get('nvalid_min', 4)

        fixed_span = prdcfg.get('fixed_span', 1)
        vmin = None
        vmax = None
        if fixed_span:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            if 'vmin' in prdcfg:
                vmin = prdcfg['vmin']
            if 'vmax' in prdcfg:
                vmax = prdcfg['vmax']

        # compute quantities
        if hmin_user is None:
            minheight = (round(
                np.min(dataset.gate_altitude['data']) /
                heightResolution)*heightResolution-heightResolution)
        else:
            minheight = hmin_user
        if hmax_user is None:
            maxheight = (round(
                np.max(dataset.gate_altitude['data']) /
                heightResolution)*heightResolution+heightResolution)
        else:
            maxheight = hmax_user
        nlevels = int((maxheight-minheight)/heightResolution)

        h_vec = minheight+np.arange(nlevels)*heightResolution+heightResolution/2.
        vals, val_valid = compute_profile_stats(
            field, dataset.gate_altitude['data'], h_vec, heightResolution,
            quantity=quantity, quantiles=quantiles/100.,
            nvalid_min=nvalid_min)

        # plot data
        if quantity == 'mean':
            data = [vals[:, 0], vals[:, 1], vals[:, 2]]
            labels = ['Mean', 'Min', 'Max']
            colors = ['b', 'k', 'k']
            linestyles = ['-', '--', '--']
        elif quantity == 'mode':
            data = [vals[:, 0], vals[:, 2], vals[:, 4]]
            labels = ['Mode', '2nd most common', '3rd most common']
            colors = ['b', 'k', 'r']
            linestyles = ['-', '--', '--']
        else:
            data = [vals[:, 1], vals[:, 0], vals[:, 2]]
            labels = [
                str(quantiles[1])+'-percentile',
                str(quantiles[0])+'-percentile',
                str(quantiles[2])+'-percentile']
            colors = ['b', 'k', 'k']
            linestyles = ['-', '--', '--']

        labelx = get_colobar_label(dataset.fields[field_name], field_name)
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'rhi_profile', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        fname = make_filename(
            'rhi_profile', prdcfg['dstype'], prdcfg['voltype'],
            ['csv'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        if quantity == 'mode':
            data.append(vals[:, 1])
            labels.append('% points mode')
            data.append(vals[:, 3])
            labels.append('% points 2nd most common')
            data.append(vals[:, 5])
            labels.append('% points 3rd most common')

        write_rhi_profile(
            h_vec, data, val_valid, labels, fname, datatype=labelx,
            timeinfo=prdcfg['timeinfo'])

        print('----- save to '+fname)

        # TODO: add Cartesian interpolation option

        return fname

    elif prdcfg['type'] == 'WIND_PROFILE':
        # user defined parameters
        heightResolution = prdcfg.get('heightResolution', 100.)
        hmin_user = prdcfg.get('heightMin', None)
        hmax_user = prdcfg.get('heightMax', None)
        min_ele = prdcfg.get('min_ele', 5.)
        max_ele = prdcfg.get('max_ele', 85.)

        fixed_span = prdcfg.get('fixed_span', 1)
        vmin = None
        vmax = None
        if fixed_span:
            vmin, vmax = pyart.config.get_field_limits('eastward_wind_component')
            if 'vmin' in prdcfg:
                vmin = prdcfg['vmin']
            if 'vmax' in prdcfg:
                vmax = prdcfg['vmax']

        u_vel = deepcopy(dataset.fields['eastward_wind_component']['data'])
        v_vel = deepcopy(dataset.fields['northward_wind_component']['data'])
        w_vel = deepcopy(dataset.fields['vertical_wind_component']['data'])
        std_vel = deepcopy(dataset.fields['retrieved_velocity_std']['data'])
        diff_vel = deepcopy(dataset.fields['velocity_difference']['data'])

        # remove azimuth information
        u_vel_aux = np.ma.empty((dataset.nsweeps, dataset.ngates), dtype=float)
        u_vel_aux[:] = np.ma.masked
        v_vel_aux = np.ma.empty((dataset.nsweeps, dataset.ngates), dtype=float)
        v_vel_aux[:] = np.ma.masked
        w_vel_aux = np.ma.empty((dataset.nsweeps, dataset.ngates), dtype=float)
        w_vel_aux[:] = np.ma.masked
        std_vel_aux = np.ma.empty((dataset.nsweeps, dataset.ngates), dtype=float)
        std_vel_aux[:] = np.ma.masked
        ngates_aux = np.zeros((dataset.nsweeps, dataset.ngates), dtype=int)
        gate_altitude_aux = np.empty((dataset.nsweeps, dataset.ngates), dtype=float)
        for ind_sweep in range(dataset.nsweeps):
            ind_start = dataset.sweep_start_ray_index['data'][ind_sweep]
            ind_end = dataset.sweep_end_ray_index['data'][ind_sweep]
            for ind_rng in range(dataset.ngates):
                u_vel_aux[ind_sweep, ind_rng] = u_vel[ind_start, ind_rng]
                v_vel_aux[ind_sweep, ind_rng] = v_vel[ind_start, ind_rng]
                w_vel_aux[ind_sweep, ind_rng] = w_vel[ind_start, ind_rng]
                std_vel_aux[ind_sweep, ind_rng] = std_vel[ind_start, ind_rng]
                gate_altitude_aux[ind_sweep, ind_rng] = (
                    dataset.gate_altitude['data'][ind_start, ind_rng])
                ngates_aux[ind_sweep, ind_rng] = (
                    diff_vel[ind_start:ind_end, ind_rng].compressed().size)

        # exclude low elevations in the computation of vertical velocities
        std_w_vel_aux = deepcopy(std_vel_aux)
        ngates_w_aux = deepcopy(ngates_aux)
        ind = np.where(dataset.fixed_angle['data'] < min_ele)[0]
        if ind.size > 0:
            w_vel_aux[ind, :] = np.ma.masked
            std_w_vel_aux[ind, :] = np.ma.masked
            ngates_w_aux[ind, :] = 0

        # exclude hig elevations in the computation of horizontal velocities
        ind = np.where(dataset.fixed_angle['data'] > max_ele)[0]
        if ind.size > 0:
            u_vel_aux[ind, :] = np.ma.masked
            v_vel_aux[ind, :] = np.ma.masked
            std_vel_aux[ind, :] = np.ma.masked
            ngates_aux[ind, :] = 0

        # compute quantities
        if hmin_user is None:
            minheight = (round(
                np.min(dataset.gate_altitude['data']) /
                heightResolution)*heightResolution-heightResolution)
        else:
            minheight = hmin_user
        if hmax_user is None:
            maxheight = (round(
                np.max(dataset.gate_altitude['data']) /
                heightResolution)*heightResolution+heightResolution)
        else:
            maxheight = hmax_user
        nlevels = int((maxheight-minheight)/heightResolution)

        h_vec = minheight+np.arange(nlevels)*heightResolution+heightResolution/2.

        u_vals, val_valid = compute_profile_stats(
            u_vel_aux, gate_altitude_aux, h_vec, heightResolution,
            quantity='regression_mean', std_field=std_vel_aux,
            np_field=ngates_aux)
        v_vals, val_valid = compute_profile_stats(
            v_vel_aux, gate_altitude_aux, h_vec, heightResolution,
            quantity='regression_mean', std_field=std_vel_aux,
            np_field=ngates_aux)
        w_vals, w_val_valid = compute_profile_stats(
            w_vel_aux, gate_altitude_aux, h_vec, heightResolution,
            quantity='regression_mean', std_field=std_w_vel_aux,
            np_field=ngates_w_aux)

        # plot u wind data
        u_data = [u_vals[:, 0], u_vals[:, 0]+u_vals[:, 1],
                  u_vals[:, 0]-u_vals[:, 1]]
        labels = ['Regression mean', '+std', '-std']
        colors = ['b', 'k', 'k']
        linestyles = ['-', '--', '--']

        labelx = get_colobar_label(
            dataset.fields['eastward_wind_component'],
            'eastward_wind_component')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset.fields['eastward_wind_component'],
                'eastward_wind_component'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'wind_vel_h_u',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            u_data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        # plot v wind data
        v_data = [v_vals[:, 0], v_vals[:, 0]+v_vals[:, 1],
                  v_vals[:, 0]-v_vals[:, 1]]
        labels = ['Regression mean', '+std', '-std']
        colors = ['b', 'k', 'k']
        linestyles = ['-', '--', '--']

        labelx = get_colobar_label(
            dataset.fields['northward_wind_component'],
            'northward_wind_component')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset.fields['northward_wind_component'],
                'northward_wind_component'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'wind_vel_h_v',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            v_data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        # plot vertical wind data
        w_data = [w_vals[:, 0], w_vals[:, 0]+w_vals[:, 1],
                  w_vals[:, 0]-w_vals[:, 1]]
        labels = ['Regression mean', '+std', '-std']
        colors = ['b', 'k', 'k']
        linestyles = ['-', '--', '--']

        labelx = get_colobar_label(
            dataset.fields['vertical_wind_component'],
            'vertical_wind_component')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(
                dataset.fields['vertical_wind_component'],
                'vertical_wind_component'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'wind_vel_v',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            w_data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        # plot horizontal wind magnitude
        mag = np.ma.sqrt(
            np.ma.power(u_vals[:, 0], 2.)+np.ma.power(v_vals[:, 0], 2.))
        mag_data = [mag]
        labels = ['Regression mean']
        colors = ['b']
        linestyles = ['-']

        field_dict = pyart.config.get_metadata('wind_speed')
        labelx = get_colobar_label(field_dict, 'wind_speed')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(field_dict, 'wind_speed'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'WIND_SPEED',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            mag_data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=vmin, vmax=vmax,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        # plot horizontal wind direction
        wind_dir = 90.-np.ma.arctan2(u_vals[:, 0]/mag, v_vals[:, 0]/mag)*180./np.pi+180.
        wind_dir[wind_dir >= 360.] = wind_dir[wind_dir >= 360.]-360.
        dir_data = [wind_dir]
        labels = ['Regression mean']
        colors = ['b']
        linestyles = ['-']

        field_dict = pyart.config.get_metadata('wind_direction')
        labelx = get_colobar_label(field_dict, 'wind_direction')
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(field_dict, 'wind_direction'))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        prdcfginfo = 'hres'+str(int(heightResolution))
        fname_list = make_filename(
            'wind_profile', prdcfg['dstype'], 'WIND_DIRECTION',
            prdcfg['imgformat'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_rhi_profile(
            dir_data, h_vec, fname_list, labelx=labelx, labely='Height (m MSL)',
            labels=labels, title=titl, colors=colors,
            linestyles=linestyles, vmin=0., vmax=360.,
            hmin=minheight, hmax=maxheight)

        print('----- save to '+' '.join(fname_list))

        fname = make_filename(
            'wind_profile', prdcfg['dstype'], 'WIND',
            ['csv'], prdcfginfo=prdcfginfo,
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        data = [
            u_vals[:, 0], u_vals[:, 1], np.ma.asarray(val_valid),
            v_vals[:, 0], v_vals[:, 1], np.ma.asarray(val_valid),
            w_vals[:, 0], w_vals[:, 1], np.ma.asarray(w_val_valid),
            mag, wind_dir]
        labels = [
            'u_wind', 'std_u_wind', 'np_u_wind',
            'v_wind', 'std_v_wind', 'np_v_wind',
            'w_wind', 'std_w_wind', 'np_w_wind',
            'mag_h_wind', 'dir_h_wind']

        write_rhi_profile(
            h_vec, data, val_valid, labels, fname, datatype=None,
            timeinfo=prdcfg['timeinfo'])

        print('----- save to '+fname)

        return fname

    elif prdcfg['type'] == 'PSEUDOPPI_MAP':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        try:
            xsect = pyart.util.cross_section_rhi(
                dataset, [prdcfg['angle']], el_tol=prdcfg['EleTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'ppi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='el'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            plot_ppi_map(xsect, field_name, 0, prdcfg, fname_list)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                'No data found at elevation ' + str(prdcfg['angle']) +
                '. Skipping product ' + prdcfg['type'])

            return None

    elif prdcfg['type'] == 'PSEUDORHI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        try:
            xsect = pyart.util.cross_section_ppi(
                dataset, [prdcfg['angle']], az_tol=prdcfg['AziTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname_list = make_filename(
                'rhi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='az'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'])

            for i, fname in enumerate(fname_list):
                fname_list[i] = savedir+fname

            step = None
            quantiles = None
            plot_type = 'RHI'
            if 'plot_type' in prdcfg:
                plot_type = prdcfg['plot_type']
            if 'step' in prdcfg:
                step = prdcfg['step']
            if 'quantiles' in prdcfg:
                quantiles = prdcfg['quantiles']

            plot_rhi(xsect, field_name, 0, prdcfg, fname_list,
                     plot_type=plot_type, step=step, quantiles=quantiles)

            print('----- save to '+' '.join(fname_list))

            return fname_list
        except EnvironmentError:
            warn(
                ' No data found at azimuth ' +
                str(prdcfg['angle'])+'. Skipping product ' +
                prdcfg['type'])
            return None

    elif prdcfg['type'] == 'CAPPI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'cappi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='alt'+'{:.1f}'.format(prdcfg['altitude']),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_cappi(
            dataset, field_name, prdcfg['altitude'], prdcfg, fname_list)
        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'PLOT_ALONG_COORD':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        if dataset.scan_type != 'ppi' and dataset.scan_type != 'rhi':
            warn('This product is only available for PPI or RHI volumes')
            return None

        colors = None
        if 'colors' in prdcfg:
            colors = prdcfg['colors']

        if prdcfg['mode'] == 'ALONG_RNG':
            value_start = 0.
            if 'value_start' in prdcfg:
                value_start = prdcfg['value_start']
            value_stop = np.max(dataset.range['data'])
            if 'value_stop' in prdcfg:
                value_stop = prdcfg['value_stop']

            rng_mask = np.logical_and(dataset.range['data'] >= value_start,
                                      dataset.range['data'] <= value_stop)

            x = dataset.range['data'][rng_mask]

            xvals = []
            yvals = []
            valid_azi = []
            valid_ele = []
            if dataset.scan_type == 'ppi':
                for i in range(len(prdcfg['fix_elevations'])):
                    d_el = np.abs(dataset.fixed_angle['data'] -
                                  prdcfg['fix_elevations'][i])
                    min_d_el = np.min(d_el)
                    if min_d_el > prdcfg['AngTol']:
                        warn('No elevation angle found for fix_elevation ' +
                             str(prdcfg['fix_elevations'][i]))
                        continue
                    ind_sweep = np.argmin(d_el)
                    new_dataset = dataset.extract_sweeps([ind_sweep])

                    try:
                        dataset_line = pyart.util.cross_section_ppi(
                            new_dataset, [prdcfg['fix_azimuths'][i]],
                            az_tol=prdcfg['AngTol'])
                    except EnvironmentError:
                        warn(' No data found at azimuth ' +
                             prdcfg['fix_azimuths'][i]+' and elevation ' +
                             prdcfg['fix_elevations'][i])
                        continue
                    yvals.append(
                        dataset_line.fields[field_name]['data'][0, rng_mask])
                    xvals.append(x)
                    valid_azi.append(dataset_line.azimuth['data'][0])
                    valid_ele.append(dataset_line.elevation['data'][0])
            else:
                for i in range(len(prdcfg['fix_azimuths'])):
                    d_az = np.abs(dataset.fixed_angle['data'] -
                                  prdcfg['fix_azimuths'][i])
                    min_d_az = np.min(d_az)
                    if min_d_az > prdcfg['AngTol']:
                        warn('No azimuth angle found for fix_azimuth ' +
                             str(prdcfg['fix_azimuths'][i]))
                        continue
                    ind_sweep = np.argmin(d_az)
                    new_dataset = dataset.extract_sweeps([ind_sweep])

                    try:
                        dataset_line = pyart.util.cross_section_rhi(
                            new_dataset, [prdcfg['fix_elevations'][i]],
                            el_tol=prdcfg['AngTol'])
                    except EnvironmentError:
                        warn(
                            ' No data found at azimuth ' +
                            prdcfg['fix_azimuths'][i]+' and elevation ' +
                            prdcfg['fix_elevations'][i])
                        continue
                    yvals.append(
                        dataset_line.fields[field_name]['data'][0, rng_mask])
                    xvals.append(x)
                    valid_azi.append(dataset_line.azimuth['data'][0])
                    valid_ele.append(dataset_line.elevation['data'][0])

            if not yvals:
                warn('No data found')
                return None

            labelx = 'Range (m)'

            labels = list()
            for i, azi in enumerate(valid_azi):
                labels.append(
                    'azi '+'{:.1f}'.format(azi) +
                    ' ele '+'{:.1f}'.format(valid_ele[i]))

        elif prdcfg['mode'] == 'ALONG_AZI':
            value_start = np.min(dataset.azimuth['data'])
            if 'value_start' in prdcfg:
                value_start = prdcfg['value_start']
            value_stop = np.max(dataset.azimuth['data'])
            if 'value_stop' in prdcfg:
                value_stop = prdcfg['value_stop']

            yvals = []
            xvals = []
            valid_rng = []
            valid_ele = []
            for i in range(len(prdcfg['fix_ranges'])):
                d_rng = np.abs(dataset.range['data'] -
                               prdcfg['fix_ranges'][i])
                min_d_rng = np.min(d_rng)
                if min_d_rng > prdcfg['RngTol']:
                    warn('No range gate found for fix_range ' +
                         str(prdcfg['fix_ranges'][i]))
                    continue
                ind_rng = np.argmin(d_rng)

                if dataset.scan_type == 'ppi':
                    d_el = np.abs(dataset.fixed_angle['data'] -
                                  prdcfg['fix_elevations'][i])
                    min_d_el = np.min(d_el)
                    if min_d_el > prdcfg['AngTol']:
                        warn('No elevation angle found for fix_elevation ' +
                             str(prdcfg['fix_elevations'][i]))
                        continue
                    ind_sweep = np.argmin(d_el)
                    new_dataset = dataset.extract_sweeps([ind_sweep])
                else:
                    try:
                        new_dataset = pyart.util.cross_section_rhi(
                            dataset, [prdcfg['fix_elevations'][i]],
                            el_tol=prdcfg['AngTol'])
                    except EnvironmentError:
                        warn(
                            ' No data found at range ' +
                            prdcfg['fix_ranges'][i]+' and elevation ' +
                            prdcfg['fix_elevations'][i])
                        continue
                if value_start < value_stop:
                    azi_mask = np.logical_and(
                        new_dataset.azimuth['data'] >= value_start,
                        new_dataset.azimuth['data'] <= value_stop)
                else:
                    azi_mask = np.logical_or(
                        new_dataset.azimuth['data'] >= value_start,
                        new_dataset.azimuth['data'] <= value_stop)
                yvals.append(
                    new_dataset.fields[field_name]['data'][azi_mask, ind_rng])
                xvals.append(new_dataset.azimuth['data'][azi_mask])
                valid_rng.append(new_dataset.range['data'][ind_rng])
                valid_ele.append(new_dataset.elevation['data'][0])

            if not yvals:
                warn('No data found')
                return None

            labelx = 'Azimuth Angle (deg)'

            labels = list()
            for i, rng in enumerate(valid_rng):
                labels.append(
                    'rng '+'{:.1f}'.format(rng) +
                    ' ele '+'{:.1f}'.format(valid_ele[i]))

        elif prdcfg['mode'] == 'ALONG_ELE':
            value_start = np.min(dataset.elevation['data'])
            if 'value_start' in prdcfg:
                value_start = prdcfg['value_start']
            value_stop = np.max(dataset.elevation['data'])
            if 'value_stop' in prdcfg:
                value_stop = prdcfg['value_stop']

            yvals = []
            xvals = []
            valid_rng = []
            valid_azi = []
            for i in range(len(prdcfg['fix_ranges'])):
                d_rng = np.abs(dataset.range['data'] -
                               prdcfg['fix_ranges'][i])
                min_d_rng = np.min(d_rng)
                if min_d_rng > prdcfg['RngTol']:
                    warn('No range gate found for fix_range ' +
                         str(prdcfg['fix_ranges'][i]))
                    continue
                ind_rng = np.argmin(d_rng)

                if dataset.scan_type == 'ppi':
                    try:
                        new_dataset = pyart.util.cross_section_ppi(
                            dataset, [prdcfg['fix_azimuths'][i]],
                            az_tol=prdcfg['AngTol'])
                    except EnvironmentError:
                        warn(
                            ' No data found at range ' +
                            prdcfg['fix_ranges'][i]+' and elevation ' +
                            prdcfg['fix_azimuths'][i])
                        continue
                else:
                    d_az = np.abs(dataset.fixed_angle['data'] -
                                  prdcfg['fix_azimuths'][i])
                    min_d_az = np.min(d_az)
                    if min_d_az > prdcfg['AngTol']:
                        warn('No azimuth angle found for fix_azimuth ' +
                             str(prdcfg['fix_azimuths'][i]))
                        continue
                    ind_sweep = np.argmin(d_az)
                    new_dataset = dataset.extract_sweeps([ind_sweep])
                ele_mask = np.logical_and(
                    new_dataset.elevation['data'] >= value_start,
                    new_dataset.elevation['data'] <= value_stop)
                yvals.append(
                    new_dataset.fields[field_name]['data'][ele_mask, ind_rng])
                xvals.append(new_dataset.elevation['data'][ele_mask])
                valid_rng.append(new_dataset.range['data'][ind_rng])
                valid_azi.append(new_dataset.elevation['data'][0])
            if not yvals:
                warn('No data found')
                return None
            labelx = 'Elevation Angle (deg)'

            labels = list()
            for i, rng in enumerate(valid_rng):
                labels.append(
                    'rng '+'{:.1f}'.format(rng) +
                    ' azi '+'{:.1f}'.format(valid_azi[i]))
        else:
            warn('Unknown plotting mode '+prdcfg['mode'])
            return None

        labely = get_colobar_label(dataset.fields[field_name], field_name)
        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            prdcfg['mode'], prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_along_coord(
            xvals, yvals, fname_list, labelx=labelx, labely=labely,
            labels=labels, title=titl, colors=colors)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'BSCOPE_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        ang_vec = np.sort(dataset.fixed_angle['data'])
        ang = ang_vec[prdcfg['anglenr']]
        ind_ang = np.where(dataset.fixed_angle['data'] == ang)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'b-scope', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='ang'+'{:.1f}'.format(ang),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_bscope(dataset, field_name, ind_ang, prdcfg, fname_list)
        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'TIME_RANGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        ang_vec = np.sort(dataset.fixed_angle['data'])
        ang = ang_vec[prdcfg['anglenr']]
        ind_ang = np.where(dataset.fixed_angle['data'] == ang)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'time-range', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='ang'+'{:.1f}'.format(ang),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_time_range(dataset, field_name, ind_ang, prdcfg, fname_list)
        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'HISTOGRAM':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        step = prdcfg.get('step', None)
        write_data = prdcfg.get('write_data', 0)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'histogram', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        bin_edges, values = compute_histogram(
            dataset.fields[field_name]['data'], field_name, step=step)

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        labelx = get_colobar_label(dataset.fields[field_name], field_name)

        plot_histogram(bin_edges, values, fname_list, labelx=labelx,
                       labely='Number of Samples', titl=titl)

        print('----- save to '+' '.join(fname_list))

        if write_data:
            fname = savedir+make_filename(
                'histogram', prdcfg['dstype'], prdcfg['voltype'],
                ['csv'], timeinfo=prdcfg['timeinfo'],
                runinfo=prdcfg['runinfo'])[0]

            hist, bin_edges_aux = np.histogram(values, bins=bin_edges)
            write_histogram(
                bin_edges, hist, fname, datatype=prdcfg['voltype'], step=step)
            print('----- save to '+fname)

            return fname

        return fname_list

    elif prdcfg['type'] == 'QUANTILES':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # mask unclassified data
        field = deepcopy(dataset.fields[field_name]['data'])
        if prdcfg['voltype'] == 'hydro':
            field = np.ma.masked_equal(field, 0)

        # user defined variables
        quantiles = prdcfg.get('quantiles', None)
        write_data = prdcfg.get('write_data', 0)

        fixed_span = prdcfg.get('fixed_span', 1)
        vmin = None
        vmax = None
        if fixed_span:
            vmin, vmax = pyart.config.get_field_limits(field_name)
            if 'vmin' in prdcfg:
                vmin = prdcfg['vmin']
            if 'vmax' in prdcfg:
                vmax = prdcfg['vmax']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'quantiles', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        quantiles, values = compute_quantiles(field, quantiles=quantiles)

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        labely = get_colobar_label(dataset.fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname_list, labelx='quantile',
                       labely=labely, titl=titl, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        if write_data:
            fname = savedir+make_filename(
                'quantiles', prdcfg['dstype'], prdcfg['voltype'],
                ['csv'], timeinfo=prdcfg['timeinfo'],
                runinfo=prdcfg['runinfo'])[0]

            write_quantiles(
                quantiles, values, fname, datatype=prdcfg['voltype'])
            print('----- save to '+fname)

            return fname

        return fname_list

    elif prdcfg['type'] == 'FIELD_COVERAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        threshold = None
        if 'threshold' in prdcfg:
            threshold = prdcfg['threshold']
        nvalid_min = 5.
        if 'nvalid_min' in prdcfg:
            nvalid_min = prdcfg['nvalid_min']

        ele_res = 1.
        if 'ele_res' in prdcfg:
            ele_res = prdcfg['ele_res']
        azi_res = 2.
        if 'azi_res' in prdcfg:
            azi_res = prdcfg['azi_res']

        ele_min = 0.
        if 'ele_min' in prdcfg:
            ele_min = prdcfg['ele_min']
        ele_max = 30.
        if 'ele_max' in prdcfg:
            ele_max = prdcfg['ele_max']
        ele_step = 5.
        if 'ele_step' in prdcfg:
            ele_step = prdcfg['ele_step']

        ele_sect_start = None
        if 'ele_sect_start' in prdcfg:
            ele_sect_start = prdcfg['ele_sect_start']
        ele_sect_stop = None
        if 'ele_sect_stop' in prdcfg:
            ele_sect_stop = prdcfg['ele_sect_stop']
        quantiles = np.array([10., 20., 30., 40., 50., 60., 70., 80., 90.])
        if 'quantiles' in prdcfg:
            quantiles = np.array(prdcfg['quantiles'])

        # get coverage per ray
        field_coverage = np.ma.empty(dataset.nrays)
        field_coverage[:] = np.ma.masked

        for i in range(dataset.nrays):
            mask = np.ma.getmaskarray(
                dataset.fields[field_name]['data'][i, :])
            if threshold is not None:
                ind = np.where(np.logical_and(
                    ~mask,
                    dataset.fields[field_name]['data'][i, :] >= threshold))[0]
            else:
                ind = np.where(~mask)[0]
            if len(ind) > nvalid_min:
                field_coverage[i] = (dataset.range['data'][ind[-1]] -
                                     dataset.range['data'][ind[0]])

        # group coverage per elevation sectors
        nsteps = int((ele_max-ele_min)/ele_step)  # number of steps
        nele = int(ele_step/ele_res)  # number of elev per step
        ele_steps_vec = np.arange(nsteps)*ele_step+ele_min

        yval = []
        xval = []
        labels = []
        for i in range(nsteps-1):
            yval_aux = np.ma.array([])
            xval_aux = np.array([])
            for j in range(nele):
                ele_target = ele_steps_vec[i]+j*ele_res
                d_ele = np.abs(dataset.elevation['data']-ele_target)
                ind_ele = np.where(d_ele < prdcfg['AngTol'])[0]
                if ind_ele.size == 0:
                    continue
                yval_aux = np.ma.concatenate(
                    [yval_aux, field_coverage[ind_ele]])
                xval_aux = np.concatenate(
                    [xval_aux, dataset.azimuth['data'][ind_ele]])
            yval.append(yval_aux)
            xval.append(xval_aux)
            labels.append('ele '+'{:.1f}'.format(ele_steps_vec[i])+'-' +
                          '{:.1f}'.format(ele_steps_vec[i+1])+' deg')

        # get mean value per azimuth for a specified elevation sector
        xmeanval = None
        ymeanval = None
        quantval = None
        labelmeanval = None
        if ele_sect_start is not None and ele_sect_stop is not None:
            ind_ele = np.where(np.logical_and(
                dataset.elevation['data'] >= ele_sect_start,
                dataset.elevation['data'] <= ele_sect_stop))
            field_coverage_sector = field_coverage[ind_ele]
            azi_sector = dataset.azimuth['data'][ind_ele]
            nazi = int((np.max(dataset.azimuth['data']) -
                        np.min(dataset.azimuth['data']))/azi_res+1)

            xmeanval = np.arange(nazi)*azi_res+np.min(dataset.azimuth['data'])
            ymeanval = np.ma.empty(nazi)
            ymeanval[:] = np.ma.masked
            for i in range(nazi):
                d_azi = np.abs(azi_sector-xmeanval[i])
                ind_azi = np.where(d_azi < prdcfg['AngTol'])[0]
                if ind_azi.size == 0:
                    continue
                ymeanval[i] = np.ma.mean(field_coverage_sector[ind_azi])
            labelmeanval = ('ele '+'{:.1f}'.format(ele_sect_start)+'-' +
                            '{:.1f}'.format(ele_sect_stop)+' deg mean val')

            meanval, quantval, nvalid = quantiles_weighted(
                field_coverage_sector, quantiles=quantiles/100.)

        # plot field coverage
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'coverage', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        plot_field_coverage(
            xval, yval, fname_list, labels=labels, title=titl, ymin=0.,
            ymax=np.max(dataset.range['data'])+60000., xmeanval=xmeanval,
            ymeanval=ymeanval, labelmeanval=labelmeanval)

        print('----- save to '+' '.join(fname_list))

        fname = make_filename(
            'coverage', prdcfg['dstype'], prdcfg['voltype'],
            ['csv'], timeinfo=prdcfg['timeinfo'])[0]

        fname = savedir+fname

        if quantval is not None:
            data_type = get_colobar_label(
                dataset.fields[field_name], field_name)
            write_field_coverage(
                quantiles, quantval, ele_sect_start, ele_sect_stop,
                np.min(xmeanval), np.max(xmeanval), threshold, nvalid_min,
                data_type, prdcfg['timeinfo'], fname)

            print('----- save to '+fname)

        return fname

    elif prdcfg['type'] == 'CDF':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        quantiles = None
        if 'quantiles' in prdcfg:
            quantiles = prdcfg['quantiles']

        sector = {
            'rmin': None,
            'rmax': None,
            'azmin': None,
            'azmax': None,
            'elmin': None,
            'elmax': None,
            'hmin': None,
            'hmax': None}

        if 'sector' in prdcfg:
            if 'rmin' in prdcfg['sector']:
                sector['rmin'] = prdcfg['sector']['rmin']
            if 'rmax' in prdcfg['sector']:
                sector['rmax'] = prdcfg['sector']['rmax']
            if 'azmin' in prdcfg['sector']:
                sector['azmin'] = prdcfg['sector']['azmin']
            if 'azmax' in prdcfg['sector']:
                sector['azmax'] = prdcfg['sector']['azmax']
            if 'elmin' in prdcfg['sector']:
                sector['elmin'] = prdcfg['sector']['elmin']
            if 'elmax' in prdcfg['sector']:
                sector['elmax'] = prdcfg['sector']['elmax']
            if 'hmin' in prdcfg['sector']:
                sector['hmin'] = prdcfg['sector']['hmin']
            if 'hmax' in prdcfg['sector']:
                sector['hmax'] = prdcfg['sector']['hmax']

        vismin = None
        if 'vismin' in prdcfg:
            vismin = prdcfg['vismin']

        absolute = False
        if 'absolute' in prdcfg:
            absolute = prdcfg['absolute']

        use_nans = False
        nan_value = 0.
        if 'use_nans' in prdcfg:
            use_nans = prdcfg['use_nans']
            if 'nan_value' in prdcfg:
                nan_value = prdcfg['nan_value']

        filterclt = False
        if 'filterclt' in prdcfg:
            filterclt = prdcfg['filterclt']

        filterprec = np.array([], dtype=int)
        if 'filterprec' in prdcfg:
            filterprec = prdcfg['filterprec']

        data = deepcopy(dataset.fields[field_name]['data'])

        # define region of interest
        roi_flag = get_ROI(dataset, field_name, sector)
        data = data[roi_flag == 1]

        ntot = np.size(roi_flag[roi_flag == 1])

        if ntot == 0:
            warn('No radar gates found in sector')
            return None

        # get number of gates with clutter and mask them
        nclut = -1
        if filterclt:
            echoID_field = get_fieldname_pyart('echoID')
            if echoID_field in dataset.fields:
                echoID_ROI = dataset.fields[echoID_field]['data'][
                    roi_flag == 1]
                nclut = len(echoID_ROI[echoID_ROI == 2])
                data[echoID_ROI == 2] = np.ma.masked

        # get number of blocked gates and filter according to visibility
        nblocked = -1
        if vismin is not None:
            vis_field = get_fieldname_pyart('VIS')
            if vis_field in dataset.fields:
                vis_ROI = dataset.fields[vis_field]['data'][roi_flag == 1]
                nblocked = len(vis_ROI[vis_ROI < vismin])
                data[vis_ROI < vismin] = np.ma.masked

        # filter according to precip type
        nprec_filter = -1
        if filterprec.size > 0:
            hydro_field = get_fieldname_pyart('hydro')
            if hydro_field in dataset.fields:
                hydro_ROI = dataset.fields[hydro_field]['data'][roi_flag == 1]
                nprec_filter = 0
                for ind_hydro in filterprec:
                    nprec_filter += len(hydro_ROI[hydro_ROI == ind_hydro])
                    data[hydro_ROI == ind_hydro] = np.ma.masked

        if absolute:
            data = np.ma.abs(data)

        mask = np.ma.getmaskarray(data)
        nnan = np.count_nonzero(mask)

        if nnan == ntot:
            warn('No valid radar gates found in sector')
            return None

        if use_nans:
            data[mask] = nan_value

        # count and filter outliers
        quantiles_lim, values_lim = compute_quantiles(
            data, quantiles=[0.2, 99.8])
        if values_lim.mask[0] or values_lim.mask[1]:
            warn('No valid radar gates found in sector')
            return None

        nsmall = np.count_nonzero(data.compressed() < values_lim[0])
        nlarge = np.count_nonzero(data.compressed() > values_lim[1])
        noutliers = nlarge+nsmall
        data = data[np.logical_and(
            data >= values_lim[0], data <= values_lim[1])]

        # number of values used for cdf computation
        ncdf = np.size(data.compressed())

        quantiles, values = compute_quantiles(data, quantiles=quantiles)

        # plot CDF
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'cdf', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        labelx = get_colobar_label(dataset.fields[field_name], field_name)

        plot_quantiles(values, quantiles/100., fname_list, labelx=labelx,
                       labely='Cumulative probability', titl=titl)

        print('----- save to '+' '.join(fname_list))

        # store cdf values
        fname = make_filename(
            'cdf', prdcfg['dstype'], prdcfg['voltype'],
            ['txt'], timeinfo=prdcfg['timeinfo'])[0]

        fname = savedir+fname

        write_cdf(
            quantiles, values, ntot, nnan, nclut, nblocked, nprec_filter,
            noutliers, ncdf, fname, use_nans=use_nans, nan_value=nan_value,
            filterprec=filterprec, vismin=vismin, sector=sector,
            datatype=labelx, timeinfo=prdcfg['timeinfo'])

        print('----- save to '+fname)

        return fname

    elif prdcfg['type'] == 'SAVEVOL':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        new_dataset = deepcopy(dataset)
        new_dataset.fields = dict()
        new_dataset.add_field(field_name, dataset.fields[field_name])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], prdcfg['voltype'], ['nc'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        pyart.io.cfradial.write_cfradial(fname, new_dataset)
        print('saved file: '+fname)

        return fname

    elif prdcfg['type'] == 'SAVEALL':
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], 'all_fields', ['nc'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        pyart.io.cfradial.write_cfradial(fname, dataset)
        print('saved file: '+fname)

        return fname

    elif prdcfg['type'] == 'SAVESTATE':
        if prdcfg['lastStateFile'] is None:
            warn('Unable to save last state file. File name not specified')
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        max_time = np.max(dataset.time['data'])
        units = dataset.time['units']
        calendar = dataset.time['calendar']
        last_date = num2date(max_time, units, calendar)

        write_last_state(last_date, prdcfg['lastStateFile'])
        print('saved file: '+prdcfg['lastStateFile'])

        return prdcfg['lastStateFile']

    else:
        warn(' Unsupported product type: ' + prdcfg['type'])
        return None
