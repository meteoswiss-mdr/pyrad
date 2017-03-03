"""
pyrad.prod.product_aux
==========================

Auxiliary functions to generate products

.. autosummary::
    :toctree: generated/

    get_prodgen_func

"""

import pyart

from .process_product import generate_vol_products
from .process_product import generate_timeseries_products
from .process_product import generate_sun_hits_products
from .process_product import generate_monitoring_products
from .process_product import generate_colocated_gates_products
from .process_product import generate_intercomp_products
from .process_product import generate_time_avg_products
from .process_product import generate_cosmo_coord_products

from .process_traj_products import generate_traj_product


def get_prodgen_func(dsformat, dsname, dstype):
    """
    maps the dataset format into its processing function

    Parameters
    ----------
    dsformat : str
        dataset group, i.e. 'VOL', etc.

    Returns
    -------
    func : function
        pyrad function used to generate the products

    """

    if dsformat == 'VOL':
        func = generate_vol_products
    elif dsformat == 'TIMESERIES':
        func = generate_timeseries_products
    elif dsformat == 'TIMEAVG':
        func = generate_time_avg_products
    elif dsformat == 'SUN_HITS':
        func = generate_sun_hits_products
    elif dsformat == 'MONITORING':
        func = generate_monitoring_products
    elif dsformat == 'INTERCOMP':
        func = generate_intercomp_products
    elif dsformat == 'TRAJ_ONLY':
        func = generate_traj_product
    elif dsformat == 'COLOCATED_GATES':
        func = generate_colocated_gates_products
    elif dsformat == 'COSMO_COORD':
        func = generate_cosmo_coord_products
    else:
        raise ValueError("ERROR: Unknown dataset format '%s' of dataset '%s'"
                         "(dataset type '%s')" % (dsformat, dsname, dstype))

    return func
