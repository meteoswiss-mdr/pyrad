"""
pyrad.prod.process_grid_products
================================

Functions for obtaining Pyrad products from gridded datasets

.. autosummary::
    :toctree: generated/

    generate_grid_products

"""

from warnings import warn

import pyart

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename

from ..graph.plots_grid import plot_surface
from ..graph.plots_grid import plot_longitude_slice, plot_latitude_slice
from ..graph.plots_grid import plot_latlon_slice


def generate_grid_products(dataset, prdcfg):
    """
    generates grid products

    Parameters
    ----------
    dataset : grid
        grid object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'SURFACE_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        level = 0
        if 'level' in prdcfg:
            level = prdcfg['level']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'surface', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='l'+str(level),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_surface(dataset, field_name, level, prdcfg, fname_list)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'LATITUDE_SLICE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        lon = dataset.origin_longitude['data'][0]
        lat = dataset.origin_latitude['data'][0]
        if 'lon' in prdcfg:
            lon = prdcfg['lon']
        if 'lat' in prdcfg:
            lat = prdcfg['lat']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'lat_slice', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='lat'+'{:.2f}'.format(lat),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_latitude_slice(dataset, field_name, lon, lat, prdcfg, fname_list)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'LONGITUDE_SLICE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        lon = dataset.origin_longitude['data'][0]
        lat = dataset.origin_latitude['data'][0]
        if 'lon' in prdcfg:
            lon = prdcfg['lon']
        if 'lat' in prdcfg:
            lat = prdcfg['lat']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'lon_slice', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='lon'+'{:.2f}'.format(lon),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_longitude_slice(
            dataset, field_name, lon, lat, prdcfg, fname_list)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'CROSS_SECTION':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        lon1 = dataset.point_longitude['data'][0, 0, 0]
        lat1 = dataset.point_latitude['data'][0, 0, 0]

        lon2 = dataset.point_longitude['data'][0, -1, -1]
        lat2 = dataset.point_latitude['data'][0, -1, -1]
        if 'coord1' in prdcfg:
            if 'lon' in prdcfg['coord1']:
                lon1 = prdcfg['coord1']['lon']
            if 'lat' in prdcfg['coord1']:
                lat1 = prdcfg['coord1']['lat']
        if 'coord2' in prdcfg:
            if 'lon' in prdcfg['coord2']:
                lon2 = prdcfg['coord2']['lon']
            if 'lat' in prdcfg['coord2']:
                lat2 = prdcfg['coord2']['lat']

        coord1 = (lon1, lat1)
        coord2 = (lon2, lat2)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'lonlat', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='lon-lat1_'+'{:.2f}'.format(lon1)+'-' +
            '{:.2f}'.format(lat1)+'_lon-lat2_' +
            '{:.2f}'.format(lon2)+'-'+'{:.2f}'.format(lat2),
            timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_latlon_slice(
            dataset, field_name, coord1, coord2, prdcfg, fname_list)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'SAVEVOL':
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

        fname = make_filename(
            'savevol', prdcfg['dstype'], prdcfg['voltype'], ['nc'],
            timeinfo=prdcfg['timeinfo'])[0]

        fname = savedir+fname

        pyart.io.write_grid(fname, dataset, write_point_x_y_z=True,
                            write_point_lon_lat_alt=True)
        print('saved file: '+fname)

        return fname

    else:
        warn(' Unsupported product type: ' + prdcfg['type'])
        return None
