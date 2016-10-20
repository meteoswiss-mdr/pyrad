"""
pyrad.prod.product_aux
==========================

Auxiliary functions to generate products

.. autosummary::
    :toctree: generated/

    get_product_type

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
