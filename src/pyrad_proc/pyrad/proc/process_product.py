"""
pyrad.proc.process_product
==========================

Functions for obtaining Pyrad products from the datasets

    get_product_type
    generate_vol_products
    get_save_dir
    make_filename

"""

import os

import matplotlib.pyplot as plt
import numpy as np

import pyart

from ..io.read_data import get_fieldname_rainbow


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
    else:
        print('ERROR: Unknown dataset type')

    return func_name


def generate_vol_products(dataset, prdcfg):
    """
    maps the product type into its processing function

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

            display = pyart.graph.RadarDisplay(dataset)
            fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                                      prdcfg['ppiImageConfig']['ysize']],
                             dpi=72)
            ax = fig.add_subplot(111)
            display.plot_ppi(field_name, sweep=ind_el)
            display.set_limits(
                ylim=[prdcfg['ppiImageConfig']['ymin'],
                      prdcfg['ppiImageConfig']['ymax']],
                xlim=[prdcfg['ppiImageConfig']['xmin'],
                      prdcfg['ppiImageConfig']['xmax']])
            display.plot_range_rings([10, 20, 30, 40])
            display.plot_cross_hair(5.)

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
                prdcfg['dsname'], prdcfg['prdname'])

            fname = make_filename(
                prdcfg['timeinfo'], 'ppi', prdcfg['dstype'],
                prdcfg['voltype'], prdcfg['convertformat'],
                prdcfginfo='el'+str(el))
            fig.savefig(savedir+fname)
            plt.close()
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

            display = pyart.graph.RadarDisplay(dataset)
            fig = plt.figure(figsize=[prdcfg['rhiImageConfig']['xsize'],
                                      prdcfg['rhiImageConfig']['ysize']],
                             dpi=72)
            ax = fig.add_subplot(111)
            display.plot_rhi(field_name, sweep=ind_az, reverse_xaxis=False)
            display.set_limits(
                ylim=[prdcfg['rhiImageConfig']['ymin'],
                      prdcfg['rhiImageConfig']['ymax']],
                xlim=[prdcfg['rhiImageConfig']['xmin'],
                      prdcfg['rhiImageConfig']['xmax']])
            display.plot_cross_hair(5.)

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
                prdcfg['dsname'], prdcfg['prdname'])

            fname = make_filename(
                prdcfg['timeinfo'], 'rhi', prdcfg['dstype'],
                prdcfg['voltype'], prdcfg['convertformat'],
                prdcfginfo='az'+str(az))
            fig.savefig(savedir+fname)
            plt.close()
            print('saved figure: '+savedir+fname)
        else:
            print(
                'WARNING: Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
    elif prdcfg['type'] == 'PSEUDOPPI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            xsect = pyart.util.cross_section_rhi(dataset, [prdcfg['angle']])

            display = pyart.graph.RadarDisplay(xsect)
            fig = plt.figure(figsize=[prdcfg['ppiImageConfig']['xsize'],
                                      prdcfg['ppiImageConfig']['ysize']],
                             dpi=72)
            ax = fig.add_subplot(111)
            display.plot_ppi(field_name, sweep=0)
            display.set_limits(
                ylim=[prdcfg['ppiImageConfig']['ymin'],
                      prdcfg['ppiImageConfig']['ymax']],
                xlim=[prdcfg['ppiImageConfig']['xmin'],
                      prdcfg['ppiImageConfig']['xmax']])
            display.plot_range_rings([10, 20, 30, 40])
            display.plot_cross_hair(5.)

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
                prdcfg['dsname'], prdcfg['prdname'])

            fname = make_filename(
                prdcfg['timeinfo'], 'ppi', prdcfg['dstype'],
                prdcfg['voltype'], prdcfg['convertformat'],
                prdcfginfo='el'+str(prdcfg['angle']))
            fig.savefig(savedir+fname)
            plt.close()
            print('saved figure: '+savedir+fname)
        else:
            print(
                'WARNING: Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
    elif prdcfg['type'] == 'PSEUDORHI_IMAGE':
        field_name = get_fieldname_rainbow(prdcfg['voltype'])
        if field_name in dataset.fields:
            xsect = pyart.util.cross_section_ppi(dataset, [prdcfg['angle']])

            display = pyart.graph.RadarDisplay(xsect)
            fig = plt.figure(figsize=[prdcfg['rhiImageConfig']['xsize'],
                                      prdcfg['rhiImageConfig']['ysize']],
                             dpi=72)
            ax = fig.add_subplot(111)
            display.plot_rhi(field_name, sweep=0, reverse_xaxis=False)
            display.set_limits(
                ylim=[prdcfg['rhiImageConfig']['ymin'],
                      prdcfg['rhiImageConfig']['ymax']],
                xlim=[prdcfg['rhiImageConfig']['xmin'],
                      prdcfg['rhiImageConfig']['xmax']])
            display.plot_cross_hair(5.)

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], prdcfg['timeinfo'],
                prdcfg['dsname'], prdcfg['prdname'])

            fname = make_filename(
                prdcfg['timeinfo'], 'rhi', prdcfg['dstype'],
                prdcfg['voltype'], prdcfg['convertformat'],
                prdcfginfo='az'+str(prdcfg['angle']))
            fig.savefig(savedir+fname)
            plt.close()
            print('saved figure: '+savedir+fname)
        else:
            print(
                'WARNING: Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
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


def make_filename(timeinfo, prdtype, dstype, dsname, ext, prdcfginfo=None):
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

    Returns
    -------
    fname : str
        file name

    """
    timeinfostr = timeinfo.strftime('%Y%m%d%H%M%S')

    if prdcfginfo is None:
        cfgstr = ''
    else:
        cfgstr = '_'+prdcfginfo

    fname = timeinfostr+'_'+prdtype+'_'+dstype+'_'+dsname+cfgstr+'.'+ext

    return fname
