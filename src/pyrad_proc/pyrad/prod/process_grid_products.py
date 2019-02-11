"""
pyrad.prod.process_grid_products
================================

Functions for obtaining Pyrad products from gridded datasets

.. autosummary::
    :toctree: generated/

    generate_sparse_grid_products
    generate_grid_products

"""

from warnings import warn

import pyart

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename

from ..graph.plots_grid import plot_surface
from ..graph.plots_grid import plot_longitude_slice, plot_latitude_slice
from ..graph.plots_grid import plot_latlon_slice
from ..graph.plots_vol import plot_pos
from ..graph.plots_aux import get_colobar_label, get_field_name


def generate_sparse_grid_products(dataset, prdcfg):
    """
    generates products defined by sparse points. Accepted product types:
        'SURFACE_IMAGE': Generates a surface image
            User defined parameters:
                'field_limits': list of floats
                    The limits of the surface to plot [deg]
                    lon0, lon1, lat0, lat1

    Parameters
    ----------
    dataset : dictionary containing the points and their values

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
        if field_name not in dataset['fields']:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'surface', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=prdcfg['timeinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        cb_label = get_colobar_label(
            dataset['fields'][field_name], field_name)

        titl = (prdcfg['timeinfo'].strftime('%Y-%m-%dT%H:%M%SZ')+'\n' +
                get_field_name(dataset['fields'][field_name], field_name))

        if 'field_limits' in prdcfg:
            field_limits = prdcfg['field_limits']
        else:
            field_limits = dataset['field_limits']

        # get colobar limits
        vmin, vmax = pyart.config.get_field_limits(field_name)
        plot_pos(
            dataset['lat'], dataset['lon'],
            dataset['fields'][field_name]['data'],
            fname_list, cb_label=cb_label,
            titl=titl, limits=field_limits, vmin=vmin,
            vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    warn(' Unsupported product type: ' + prdcfg['type'])
    return None

def generate_grid_products(dataset, prdcfg):
    """
    generates grid products. Accepted product types:
        'CROSS_SECTION': Plots a cross-section of gridded data
            User defined parameters:
                coord1, coord2: dict
                    The two lat-lon coordinates marking the limits. They have
                    the keywords 'lat' and 'lon' [degree]. The altitude limits
                    are defined by the parameters in 'rhiImageConfig' in the
                    'loc' configuration file
        'LATITUDE_SLICE': Plots a cross-section of gridded data over a
            constant latitude.
            User defined parameters:
                lon, lat: floats
                    The starting point of the cross-section. The ending point
                    is defined by the parameters in 'rhiImageConfig' in the
                    'loc' configuration file
        'LONGITUDE_SLICE': Plots a cross-ection of gridded data over a
            constant longitude.
            User defined parameters:
                lon, lat: floats
                    The starting point of the cross-section. The ending point
                    is defined by the parameters in 'rhiImageConfig' in the
                    'loc' configuration file
        'SAVEVOL': save the gridded data in a C/F radial file
        'SURFACE_IMAGE': Plots a surface image of gridded data.
            User defined parameters:
                level: int
                    The altitude level to plot. The rest of the parameters are
                    defined by the parameters in 'ppiImageConfig' and
                    'ppiMapImageConfig' in the 'loc' configuration file

    Parameters
    ----------
    dataset : grid
        grid object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    None or name of generated files

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

    if prdcfg['type'] == 'LATITUDE_SLICE':
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

    if prdcfg['type'] == 'LONGITUDE_SLICE':
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

    if prdcfg['type'] == 'CROSS_SECTION':
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

    if prdcfg['type'] == 'SAVEVOL':
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

    warn(' Unsupported product type: ' + prdcfg['type'])
    return None
