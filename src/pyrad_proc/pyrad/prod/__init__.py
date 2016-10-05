"""
======================================================
Products generation (:mod:`pyrad.prod`)
======================================================

.. currentmodule:: pyrad.prod

Initiate the products generation.

Auxiliary functions
==================

.. autosummary::
    :toctree: generated/

    get_product_type
    get_save_dir
    make_filename

Product generation
==================

.. autosummary::
    :toctree: generated/

    generate_vol_products
    generate_timeseries_products
    generate_sun_hits_products

"""

from .product_aux import get_product_type, get_save_dir, make_filename

from .process_product import generate_vol_products
from .process_product import generate_timeseries_products
from .process_product import generate_sun_hits_products

__all__ = [s for s in dir() if not s.startswith('_')]
