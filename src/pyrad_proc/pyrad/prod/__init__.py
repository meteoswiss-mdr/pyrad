"""
======================================================
Products generation (:mod:`pyrad.prod`)
======================================================

.. currentmodule:: pyrad.prod

Initiate the products generation.

Auxiliary functions
===================

.. autosummary::
    :toctree: generated/

    get_dsformat_func

Product generation
==================

.. autosummary::
    :toctree: generated/

    generate_cosmo_coord_products
    generate_sun_hits_products
    generate_intercomp_products
    generate_colocated_gates_products
    generate_time_avg_products
    generate_qvp_products
    generate_vol_products
    generate_timeseries_products
    generate_monitoring_products
    generate_grid_products

"""

from .product_aux import get_prodgen_func

from .process_product import generate_vol_products
from .process_product import generate_timeseries_products
from .process_product import generate_sun_hits_products
from .process_product import generate_monitoring_products
from .process_product import generate_cosmo_coord_products
from .process_product import generate_grid_products
from .process_product import generate_intercomp_products
from .process_product import generate_colocated_gates_products
from .process_product import generate_time_avg_products
from .process_product import generate_qvp_products


from .process_traj_products import generate_traj_product

__all__ = [s for s in dir() if not s.startswith('_')]
