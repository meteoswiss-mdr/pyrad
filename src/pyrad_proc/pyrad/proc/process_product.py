"""
pyrad.proc.process_product
==========================

Functions for obtaining Pyrad products from the datasets

.. autosummary::
    :toctree: generated/

    get_product_type
    generate_vol_products
    generate_sun_hits_products
    generate_timeseries_products
    get_save_dir
    make_filename

"""

import os
from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.read_data import get_fieldname_rainbow, read_timeseries
from ..io.read_data import get_sensor_data

from ..io.write_data import write_timeseries, generate_field_name_str
from ..io.write_data import write_sun_hits, write_sun_retrieval

from ..graph.plots import plot_ppi, plot_rhi, plot_cappi, plot_bscope
from ..graph.plots import plot_timeseries, plot_timeseries_comp
from ..graph.plots import plot_quantiles, get_colobar_label, plot_sun_hits

from ..util.radar_utils import create_sun_hits_field
from ..util.radar_utils import create_sun_retrieval_field


def get_product_type(product_type):
    """
    maps the product type into its processing function

    Parameters
    ----------
    product_type : str
        product type, i.e. 'VOL', etc.

    Returns
    -------
    func_name : str
        pyrad function used to generate the product

    """
    if product_type == 'VOL':
        func_name = 'generate_vol_products'
    elif product_type == 'TIMESERIES':
        func_name = 'generate_timeseries_products'
    elif product_type == 'SUN_HITS':
        func_name = 'generate_sun_hits_products'
    else:
        print('ERROR: Unknown dataset type')

    return func_name


def generate_sun_hits_products(dataset, prdcfg):
    """
    generates sun hits products

    Parameters
    ----------
    dataset : tuple
        radar object and sun hits dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """
    if prdcfg['type'] == 'WRITE_SUN_HITS':
        if 'sun_hits' in dataset:
            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname = make_filename(
                'info', prdcfg['dstype'], 'detected', 'csv',
                timeinfo=prdcfg['timeinfo'], timeformat='%Y%m%d')

            write_sun_hits(dataset['sun_hits'], savedir+fname)

            print('saved sun hits file: '+savedir+fname)

    elif prdcfg['type'] == 'PLOT_SUN_HITS':
        if 'sun_hits_final' in dataset:
            field_name = get_fieldname_rainbow(prdcfg['voltype'])
            if prdcfg['voltype'] in dataset['sun_hits_final']:
                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                    prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

                fname = make_filename(
                    'detected', prdcfg['dstype'], prdcfg['voltype'],
                    prdcfg['convertformat'], timeinfo=prdcfg['timeinfo'],
                    timeformat='%Y%m%d')

                field = create_sun_hits_field(
                    dataset['sun_hits_final']['rad_el'],
                    dataset['sun_hits_final']['rad_az'],
                    dataset['sun_hits_final']['sun_el'],
                    dataset['sun_hits_final']['sun_az'],
                    dataset['sun_hits_final'][prdcfg['voltype']],
                    prdcfg['sunhitsImageConfig'])

                if field is not None:
                    plot_sun_hits(field, field_name, savedir+fname, prdcfg)

                    print('saved image: '+savedir+fname)
            else:
                warn(
                    ' Field type ' + prdcfg['voltype'] +
                    ' not available in data set. Skipping product ' +
                    prdcfg['type'])

    elif prdcfg['type'] == 'WRITE_SUN_RETRIEVAL':
        if 'sun_retrieval' in dataset:
            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                prdcfg['prdname'], timeinfo=None)

            fname = make_filename(
                'info', prdcfg['dstype'], 'retrieval', 'csv')

            write_sun_retrieval(dataset['sun_retrieval'], savedir+fname)

            print('saved sun retrieval file: '+savedir+fname)

    elif prdcfg['type'] == 'PLOT_SUN_RETRIEVAL':
        if 'sun_retrieval' in dataset:
            field_name = get_fieldname_rainbow(prdcfg['voltype'])
            par = None
            if field_name == 'sun_est_power_h':
                par = 'par_h'
            elif field_name == 'sun_est_power_v':
                par = 'par_v'
            elif field_name == 'sun_est_differential_reflectivity':
                par = 'par_zdr'

            if par in dataset['sun_retrieval'] and par is not None:
                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                    prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

                fname = make_filename(
                    'retrieval', prdcfg['dstype'], prdcfg['voltype'],
                    prdcfg['convertformat'], timeinfo=prdcfg['timeinfo'],
                    timeformat='%Y%m%d')

                field = create_sun_retrieval_field(
                    dataset['sun_retrieval'][par],
                    prdcfg['sunhitsImageConfig'])

                if field is not None:
                    plot_sun_hits(field, field_name, savedir+fname, prdcfg)

                    print('saved image: '+savedir+fname)
            else:
                warn(
                    ' Field type ' + prdcfg['voltype'] +
                    ' not available in data set. Skipping product ' +
                    prdcfg['type'])

    else:
        if 'radar' in dataset:
            generate_vol_products(dataset['radar'], prdcfg)


def generate_vol_products(dataset, prdcfg):
    """
    generates radar volume products

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
    if prdcfg['type'] == 'PPI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            el_vec = np.sort(dataset.fixed_angle['data'])
            el = el_vec[prdcfg['anglenr']]
            ind_el = np.where(dataset.fixed_angle['data'] == el)[0][0]

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname = make_filename(
                'ppi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['convertformat'], prdcfginfo='el'+'{:.1f}'.format(el),
                timeinfo=prdcfg['timeinfo'])

            step = None
            quantiles = None
            plot_type = 'PPI'
            if 'plot_type' in prdcfg:
                plot_type = prdcfg['plot_type']
            if 'step' in prdcfg:
                step = prdcfg['step']
            if 'quantiles' in prdcfg:
                quantiles = prdcfg['quantiles']

            plot_ppi(dataset, field_name, ind_el, prdcfg, savedir+fname,
                     plot_type=plot_type, step=step, quantiles=quantiles)

            print('saved figure: '+savedir+fname)
        else:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'RHI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            az_vec = np.sort(dataset.fixed_angle['data'])
            az = az_vec[prdcfg['anglenr']]
            ind_az = np.where(dataset.fixed_angle['data'] == az)[0][0]

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname = make_filename(
                'rhi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['convertformat'], prdcfginfo='az'+'{:.1f}'.format(az),
                timeinfo=prdcfg['timeinfo'])

            step = None
            quantiles = None
            plot_type = 'RHI'
            if 'plot_type' in prdcfg:
                plot_type = prdcfg['plot_type']
            if 'step' in prdcfg:
                step = prdcfg['step']
            if 'quantiles' in prdcfg:
                quantiles = prdcfg['quantiles']

            plot_rhi(dataset, field_name, ind_az, prdcfg, savedir+fname,
                     plot_type=plot_type, step=step, quantiles=quantiles)

            print('saved figure: '+savedir+fname)
        else:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'PSEUDOPPI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            try:
                xsect = pyart.util.cross_section_rhi(
                    dataset, [prdcfg['angle']], el_tol=prdcfg['EleTol'])

                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                    prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

                fname = make_filename(
                    'ppi', prdcfg['dstype'], prdcfg['voltype'],
                    prdcfg['convertformat'],
                    prdcfginfo='el'+'{:.1f}'.format(prdcfg['angle']),
                    timeinfo=prdcfg['timeinfo'])

                step = None
                quantiles = None
                plot_type = 'PPI'
                if 'plot_type' in prdcfg:
                    plot_type = prdcfg['plot_type']
                if 'step' in prdcfg:
                    step = prdcfg['step']
                if 'quantiles' in prdcfg:
                    quantiles = prdcfg['quantiles']

                plot_ppi(xsect, field_name, 0, prdcfg, savedir+fname,
                         plot_type=plot_type, step=step, quantiles=quantiles)

                print('saved figure: '+savedir+fname)
            except EnvironmentError:
                warn(
                    'No data found at elevation ' + str(prdcfg['angle']) +
                    '. Skipping product ' + prdcfg['type'])
        else:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'PSEUDORHI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            try:
                xsect = pyart.util.cross_section_ppi(
                    dataset, [prdcfg['angle']], az_tol=prdcfg['AziTol'])

                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                    prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

                fname = make_filename(
                    'rhi', prdcfg['dstype'], prdcfg['voltype'],
                    prdcfg['convertformat'],
                    prdcfginfo='az'+'{:.1f}'.format(prdcfg['angle']),
                    timeinfo=prdcfg['timeinfo'])

                step = None
                quantiles = None
                plot_type = 'RHI'
                if 'plot_type' in prdcfg:
                    plot_type = prdcfg['plot_type']
                if 'step' in prdcfg:
                    step = prdcfg['step']
                if 'quantiles' in prdcfg:
                    quantiles = prdcfg['quantiles']

                plot_rhi(xsect, field_name, 0, prdcfg, savedir+fname,
                         plot_type=plot_type, step=step, quantiles=quantiles)

                print('saved figure: '+savedir+fname)
            except EnvironmentError:
                warn(
                    ' No data found at azimuth ' +
                    str(prdcfg['angle'])+'. Skipping product ' +
                    prdcfg['type'])
        else:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'CAPPI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname = make_filename(
                'cappi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['convertformat'],
                prdcfginfo='alt'+'{:.1f}'.format(prdcfg['altitude']),
                timeinfo=prdcfg['timeinfo'])

            plot_cappi(dataset, field_name, prdcfg['altitude'], prdcfg,
                       savedir+fname)
            print('saved figure: '+savedir+fname)
        else:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    if prdcfg['type'] == 'BSCOPE_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            ang_vec = np.sort(dataset.fixed_angle['data'])
            ang = ang_vec[prdcfg['anglenr']]
            ind_ang = np.where(dataset.fixed_angle['data'] == ang)[0][0]

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname = make_filename(
                'b-scope', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['convertformat'],
                prdcfginfo='ang'+'{:.1f}'.format(ang),
                timeinfo=prdcfg['timeinfo'])

            plot_bscope(dataset, field_name, ind_ang, prdcfg, savedir+fname)
            print('saved figure: '+savedir+fname)
        else:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    if prdcfg['type'] == 'QUANTILES':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            ang_vec = np.sort(dataset.fixed_angle['data'])
            ang = ang_vec[prdcfg['anglenr']]
            ind_ang = np.where(dataset.fixed_angle['data'] == ang)[0][0]

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname = make_filename(
                'quantiles', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['convertformat'],
                prdcfginfo='ang'+'{:.1f}'.format(ang),
                timeinfo=prdcfg['timeinfo'])

            quant = dataset.range['data']
            value = dataset.fields[field_name]['data'][ind_ang, :]
            titl = pyart.graph.common.generate_title(
                dataset, field_name, ind_ang)
            labely = get_colobar_label(dataset.fields[field_name], field_name)

            plot_quantiles(quant, value, savedir+fname, labelx='quantile',
                           labely=labely, titl=titl)

            print('saved figure: '+savedir+fname)
        else:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'SAVEVOL':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            new_dataset = deepcopy(dataset)
            new_dataset.fields = dict()
            new_dataset.add_field(field_name, dataset.fields[field_name])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname = make_filename(
                'savevol', prdcfg['dstype'], prdcfg['voltype'], 'nc',
                timeinfo=prdcfg['timeinfo'])

            pyart.io.cfradial.write_cfradial(savedir+fname, new_dataset)
            print('saved file: '+savedir+fname)
        else:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'SAVEALL':
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], 'all_fields', 'nc',
            timeinfo=prdcfg['timeinfo'])

        pyart.io.cfradial.write_cfradial(savedir+fname, dataset)
        print('saved file: '+savedir+fname)

    else:
        warn(' Unsupported product type: ' + prdcfg['type'])


def generate_timeseries_products(dataset, prdcfg):
    """
    generates time series products

    Parameters
    ----------
    dataset : dictionary
        radar object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """
    if prdcfg['type'] == 'PLOT_AND_WRITE_POINT':
        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], 'csv',
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        write_timeseries(dataset, savedir+csvfname)
        print('saved CSV file: '+savedir+csvfname)

        date, value = read_timeseries(savedir+csvfname)

        if date is not None:
            figfname = make_filename(
                'ts', prdcfg['dstype'], dataset['datatype'],
                prdcfg['convertformat'], prdcfginfo=gateinfo,
                timeinfo=date[0], timeformat='%Y%m%d')

            label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
            titl = ('Time Series '+date[0].strftime('%Y-%m-%d'))

            labely = generate_field_name_str(dataset['datatype'])

            plot_timeseries(
                date, value, savedir+figfname, labelx='Time UTC',
                labely=labely, label1=label1, titl=titl)
            print('saved figure: '+savedir+figfname)

    elif prdcfg['type'] == 'PLOT_CUMULATIVE_POINT':
        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
            prdcfg['prdid'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], 'csv',
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        date, value = read_timeseries(savedir+csvfname)

        if date is not None:
            figfname = make_filename(
                'ts_cum', prdcfg['dstype'], dataset['datatype'],
                prdcfg['convertformat'], prdcfginfo=gateinfo,
                timeinfo=date[0], timeformat='%Y%m%d')

            label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
            titl = ('Time Series Acc. '+date[0].strftime('%Y-%m-%d'))

            labely = 'Radar estimated rainfall accumulation (mm)'

            plot_timeseries(
                date, value, savedir+figfname, labelx='Time UTC',
                labely=labely, label1=label1, titl=titl,
                period=prdcfg['ScanPeriod']*60.)
            print('saved figure: '+savedir+figfname)

    elif prdcfg['type'] == 'COMPARE_POINT':
        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
            prdcfg['prdid'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], 'csv',
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        radardate, radarvalue = read_timeseries(savedir_ts+csvfname)

        if radardate is not None:
            sensordate, sensorvalue, sensortype, period = get_sensor_data(
                radardate[0], dataset['datatype'], prdcfg)

            if sensordate is not None:
                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                    prdcfg['prdname'], timeinfo=radardate[0])

                figfname = make_filename(
                    'ts_comp', prdcfg['dstype'], dataset['datatype'],
                    prdcfg['convertformat'], prdcfginfo=gateinfo,
                    timeinfo=radardate[0], timeformat='%Y%m%d')

                label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
                label2 = sensortype+' '+prdcfg['sensorid']
                titl = ('Time Series Comp. ' +
                        radardate[0].strftime('%Y-%m-%d'))

                labely = generate_field_name_str(dataset['datatype'])

                plot_timeseries_comp(
                    radardate, radarvalue, sensordate, sensorvalue,
                    savedir+figfname, labelx='Time UTC', labely=labely,
                    label1=label1, label2=label2, titl=titl)
                print('saved figure: '+savedir+figfname)

    elif prdcfg['type'] == 'COMPARE_CUMULATIVE_POINT':
        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
            prdcfg['dsname'], prdcfg['prdid'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], 'csv',
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        radardate, radarvalue = read_timeseries(savedir_ts+csvfname)

        if radardate is not None:
            sensordate, sensorvalue, sensortype, period2 = get_sensor_data(
                radardate[0], dataset['datatype'], prdcfg)

            if sensordate is not None:
                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'], prdcfg['dsname'],
                    prdcfg['prdname'], timeinfo=radardate[0])

                figfname = make_filename(
                    'ts_cumcomp', prdcfg['dstype'], dataset['datatype'],
                    prdcfg['convertformat'], prdcfginfo=gateinfo,
                    timeinfo=radardate[0], timeformat='%Y%m%d')

                label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
                label2 = sensortype+' '+prdcfg['sensorid']
                titl = ('Time Series Acc. Comp. ' +
                        radardate[0].strftime('%Y-%m-%d'))
                labely = 'Rainfall accumulation (mm)'

                plot_timeseries_comp(
                    radardate, radarvalue, sensordate, sensorvalue,
                    savedir+figfname, labelx='Time UTC', labely=labely,
                    label1=label1, label2=label2, titl=titl,
                    period1=prdcfg['ScanPeriod']*60., period2=period2)
                print('saved figure: '+savedir+figfname)

    else:
        warn(' Unsupported product type: ' + prdcfg['type'])


def get_save_dir(basepath, procname, dsname, prdname, timeinfo=None,
                 timeformat='%Y-%m-%d'):
    """
    obtains the path to a product directory and eventually creates it

    Parameters
    ----------
    basepath : str
        product base path
    procname : str
        name of processing space
    dsname : str
        data set name
    prdname : str
        product name
    timeinfo : datetime
        time info to generate the date directory. If None there is no time
        format in the path
    timeformat : str
        Optional. The time format.

    Returns
    -------
    savedir : str
        path to product

    """
    if timeinfo is None:
        savedir = basepath+procname+'/'+dsname+'/'+prdname+'/'
    else:
        daydir = timeinfo.strftime('%Y-%m-%d')
        savedir = basepath+procname+'/'+daydir+'/'+dsname+'/'+prdname+'/'

    if os.path.isdir(savedir):
        pass
    else:
        os.makedirs(savedir)

    return savedir


def make_filename(prdtype, dstype, dsname, ext, prdcfginfo=None,
                  timeinfo=None, timeformat='%Y%m%d%H%M%S'):
    """
    creates a product file name

    Parameters
    ----------
    timeinfo : datetime
        time info to generate the date directory

    prdtype : str
        product type, i.e. 'ppi', etc.

    dstype : str
        data set type, i.e. 'raw', etc.

    dsname : str
        data set name

    ext : str
        file name extension, i.e. 'png'

    prdcfginfo : str
        Optional. string to add product configuration information, i.e. 'el0.4'

    timeformat : str
        Optional. The time format

    Returns
    -------
    fname : str
        file name

    """
    if timeinfo is None:
        timeinfostr = ''
    else:
        timeinfostr = timeinfo.strftime(timeformat)+'_'

    if prdcfginfo is None:
        cfgstr = ''
    else:
        cfgstr = '_'+prdcfginfo

    fname = timeinfostr+prdtype+'_'+dstype+'_'+dsname+cfgstr+'.'+ext

    return fname
