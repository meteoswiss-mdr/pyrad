"""
pyrad.prod.product_aux
==========================

Auxiliary functions to generate products

.. autosummary::
    :toctree: generated/

    get_product_type
    get_save_dir
    make_filename

"""

import os
from copy import deepcopy
from warnings import warn

import numpy as np

import pyart


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
