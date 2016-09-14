"""
pyrad.proc.process_product
==========================

Functions for obtaining Pyrad products from the datasets

.. autosummary::
    :toctree: generated/

    get_product_type
    generate_vol_products
    generate_timeseries_products
    get_save_dir
    make_filename

"""

import os
from copy import deepcopy

import numpy as np

import pyart

from ..io.read_data import get_fieldname_rainbow, read_timeseries
from ..io.read_data import get_sensor_data
from ..io.write_data import write_timeseries, generate_field_name_str
from ..graph.plots import plot_ppi, plot_rhi, plot_cappi 
from ..graph.plots import plot_timeseries, plot_timeseries_comp


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
    else:
        print('ERROR: Unknown dataset type')

    return func_name


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
                prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
                prdcfg['dsname'], prdcfg['prdname'])

            fname = make_filename(
                prdcfg['timeinfo'], 'ppi', prdcfg['dstype'],
                prdcfg['voltype'], prdcfg['convertformat'],
                prdcfginfo='el'+'{:.1f}'.format(el))

            plot_ppi(dataset, field_name, ind_el, prdcfg, savedir+fname)
            print('saved figure: '+savedir+fname)
        else:
            print(
                'WARNING: Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'RHI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            az_vec = np.sort(dataset.fixed_angle['data'])
            az = az_vec[prdcfg['anglenr']]
            ind_az = np.where(dataset.fixed_angle['data'] == az)[0][0]

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
                prdcfg['dsname'], prdcfg['prdname'])

            fname = make_filename(
                prdcfg['timeinfo'], 'rhi', prdcfg['dstype'],
                prdcfg['voltype'], prdcfg['convertformat'],
                prdcfginfo='az'+'{:.1f}'.format(az))

            plot_rhi(dataset, field_name, ind_az, prdcfg, savedir+fname)
            print('saved figure: '+savedir+fname)
        else:
            print(
                'WARNING: Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'PSEUDOPPI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            xsect = pyart.util.cross_section_rhi(dataset, [prdcfg['angle']],
                                                 el_tol=prdcfg['EleTol'])
            if xsect is not None:
                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'],
                    prdcfg['timeinfo'], prdcfg['dsname'], prdcfg['prdname'])
    
                fname = make_filename(
                    prdcfg['timeinfo'], 'ppi', prdcfg['dstype'],
                    prdcfg['voltype'], prdcfg['convertformat'],
                    prdcfginfo='el'+'{:.1f}'.format(prdcfg['angle']))
    
                plot_ppi(xsect, field_name, 0, prdcfg, savedir+fname)
                print('saved figure: '+savedir+fname)
            else:
                print(
                    'WARNING: No data found at elevation ' +
                    str(prdcfg['angle']) +'. Skipping product ' +
                    prdcfg['type'])
        else:
            print(
                'WARNING: Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'PSEUDORHI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            xsect = pyart.util.cross_section_ppi(dataset, [prdcfg['angle']],
                                                 az_tol=prdcfg['AziTol'])
            if xsect is not None:
                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'],
                    prdcfg['timeinfo'], prdcfg['dsname'], prdcfg['prdname'])

                fname = make_filename(
                    prdcfg['timeinfo'], 'rhi', prdcfg['dstype'],
                    prdcfg['voltype'], prdcfg['convertformat'],
                    prdcfginfo='az'+'{:.1f}'.format(prdcfg['angle']))

                plot_rhi(xsect, field_name, 0, prdcfg, savedir+fname)
                print('saved figure: '+savedir+fname)
            else:
                print(
                    'WARNING: No data found at azimuth ' +
                    str(prdcfg['angle'])+'. Skipping product ' +
                    prdcfg['type'])
        else:
            print(
                'WARNING: Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'CAPPI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
                prdcfg['dsname'], prdcfg['prdname'])

            fname = make_filename(
                prdcfg['timeinfo'], 'cappi', prdcfg['dstype'],
                prdcfg['voltype'], prdcfg['convertformat'],
                prdcfginfo='alt'+'{:.1f}'.format(prdcfg['altitude']))

            plot_cappi(dataset, field_name, prdcfg['altitude'], prdcfg,
                       savedir+fname)
            print('saved figure: '+savedir+fname)
        else:
            print(
                'WARNING: Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
                
    elif prdcfg['type'] == 'SAVEVOL':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            new_dataset = deepcopy(dataset)
            new_dataset.fields = dict()
            new_dataset.add_field(field_name, dataset.fields[field_name])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
                prdcfg['dsname'], prdcfg['prdname'])

            fname = make_filename(
                prdcfg['timeinfo'], 'savevol', prdcfg['dstype'],
                prdcfg['voltype'], 'nc')

            pyart.io.cfradial.write_cfradial(savedir+fname, new_dataset)
            print('saved file: '+savedir+fname)
        else:
            print(
                'WARNING: Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

    elif prdcfg['type'] == 'SAVEALL':
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
            prdcfg['dsname'], prdcfg['prdname'])

        fname = make_filename(
            prdcfg['timeinfo'], 'savevol', prdcfg['dstype'],
            'all_fields', 'nc')

        pyart.io.cfradial.write_cfradial(savedir+fname, dataset)
        print('saved file: '+savedir+fname)

    else:
        print('WARNING: Unsupported product type: ' + prdcfg['type'])


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
            prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
            prdcfg['dsname'], prdcfg['prdname'])

        csvfname = make_filename(
            prdcfg['timeinfo'], 'ts', prdcfg['dstype'], dataset['datatype'],
            'csv', prdcfginfo=gateinfo, timeformat='%Y%m%d')
        write_timeseries(dataset, savedir+csvfname)
        print('saved CSV file: '+savedir+csvfname)

        date, value = read_timeseries(savedir+csvfname)

        if date is not None:
            figfname = make_filename(
                date[0], 'ts', prdcfg['dstype'], dataset['datatype'],
                prdcfg['convertformat'], prdcfginfo=gateinfo,
                timeformat='%Y%m%d')

            titl = ('Time Series '+date[0].strftime('%Y-%m-%d') +
                    ' (az, el, r): ('+az+', '+el+', '+r+')')

            labely = generate_field_name_str(dataset['datatype'])

            plot_timeseries(
                date, value, savedir+figfname, labelx='Time UTC',
                labely=labely, titl=titl)
            print('saved figure: '+savedir+figfname)

    elif prdcfg['type'] == 'PLOT_CUMULATIVE_POINT':
        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
            prdcfg['dsname'], prdcfg['prdid'])

        csvfname = make_filename(
            prdcfg['timeinfo'], 'ts', prdcfg['dstype'], dataset['datatype'],
            'csv', prdcfginfo=gateinfo, timeformat='%Y%m%d')

        date, value = read_timeseries(savedir+csvfname)

        if date is not None:
            figfname = make_filename(
                date[0], 'ts_cum', prdcfg['dstype'], dataset['datatype'],
                prdcfg['convertformat'], prdcfginfo=gateinfo,
                timeformat='%Y%m%d')

            titl = ('Time Series Acc. '+date[0].strftime('%Y-%m-%d') +
                    ' (az, el, r): ('+az+', '+el+', '+r+')')

            labely = 'Radar estimated rainfall accumulation (mm)'

            plot_timeseries(
                date, value, savedir+figfname, labelx='Time UTC',
                labely=labely, titl=titl, period=prdcfg['ScanPeriod']*60.)
            print('saved figure: '+savedir+figfname)

    elif prdcfg['type'] == 'COMPARE_POINT':
        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
            prdcfg['dsname'], prdcfg['prdid'])

        csvfname = make_filename(
            prdcfg['timeinfo'], 'ts', prdcfg['dstype'], dataset['datatype'],
            'csv', prdcfginfo=gateinfo, timeformat='%Y%m%d')

        radardate, radarvalue = read_timeseries(savedir_ts+csvfname)

        if radardate is not None:
            sensordate, sensorvalue, label2, period = get_sensor_data(
                radardate[0], dataset['datatype'], prdcfg)

            if sensordate is not None:
                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'], radardate[0],
                    prdcfg['dsname'], prdcfg['prdname'])

                figfname = make_filename(
                    radardate[0], 'ts_comp', prdcfg['dstype'],
                    dataset['datatype'], prdcfg['convertformat'],
                    prdcfginfo=gateinfo, timeformat='%Y%m%d')

                titl = ('Time Series Comp. ' +
                        radardate[0].strftime('%Y-%m-%d') +
                        ' (az, el, r): ('+az+', '+el+', '+r+') - '+label2 +
                        ' '+prdcfg['sensorid'])

                labely = generate_field_name_str(dataset['datatype'])

                plot_timeseries_comp(
                    radardate, radarvalue, sensordate, sensorvalue,
                    savedir+figfname, labelx='Time UTC', labely=labely,
                    label1='Radar', label2=label2, titl=titl)
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
            prdcfg['timeinfo'], 'ts', prdcfg['dstype'], dataset['datatype'],
            'csv', prdcfginfo=gateinfo, timeformat='%Y%m%d')

        radardate, radarvalue = read_timeseries(savedir_ts+csvfname)

        if radardate is not None:
            sensordate, sensorvalue, label2, period2 = get_sensor_data(
                radardate[0], dataset['datatype'], prdcfg)

            if sensordate is not None:
                savedir = get_save_dir(
                    prdcfg['basepath'], prdcfg['procname'], radardate[0],
                    prdcfg['dsname'], prdcfg['prdname'])

                figfname = make_filename(
                    radardate[0], 'ts_cumcomp', prdcfg['dstype'],
                    dataset['datatype'], prdcfg['convertformat'],
                    prdcfginfo=gateinfo, timeformat='%Y%m%d')

                titl = ('Time Series Acc. Comp. ' +
                        radardate[0].strftime('%Y-%m-%d') +
                        ' (az, el, r): ('+az+', '+el+', '+r+') - '+label2 +
                        ' ' + prdcfg['sensorid'])

                labely = 'Rainfall accumulation (mm)'

                plot_timeseries_comp(
                    radardate, radarvalue, sensordate, sensorvalue,
                    savedir+figfname, labelx='Time UTC', labely=labely,
                    label1='Radar', label2=label2, titl=titl,
                    period1=prdcfg['ScanPeriod']*60., period2=period2)
                print('saved figure: '+savedir+figfname)

    else:
        print('WARNING: Unsupported product type: ' + prdcfg['type'])


def get_save_dir(basepath, procname, timeinfo, dsname, prdname):
    """
    obtains the path to a product directory and eventually creates it

    Parameters
    ----------
    basepath : str
        product base path

    procname : str
        name of processing space

    timeinfo : datetime
        time info to generate the date directory

    dsname : str
        data set name

    prdname : str
        product name

    Returns
    -------
    savedir : str
        path to product

    """
    daydir = timeinfo.strftime('%Y-%m-%d')
    savedir = basepath+procname+'/'+daydir+'/'+dsname+'/'+prdname+'/'
    if os.path.isdir(savedir):
        pass
    else:
        os.makedirs(savedir)

    return savedir


def make_filename(timeinfo, prdtype, dstype, dsname, ext, prdcfginfo=None,
                  timeformat='%Y%m%d%H%M%S'):
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
    timeinfostr = timeinfo.strftime(timeformat)

    if prdcfginfo is None:
        cfgstr = ''
    else:
        cfgstr = '_'+prdcfginfo

    fname = timeinfostr+'_'+prdtype+'_'+dstype+'_'+dsname+cfgstr+'.'+ext

    return fname
