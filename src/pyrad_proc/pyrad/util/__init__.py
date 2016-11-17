"""
==================================
Utilities (:mod:`pyrad.util`)
==================================

.. currentmodule:: pyrad.util

Functions to read and write data and configuration files.

Radar Utilities
===============

.. autosummary::
    :toctree: generated/

    time_avg_range
    get_closest_solar_flux
    create_sun_hits_field
    create_sun_retrieval_field
    compute_quantiles
    compute_quantiles_from_hist
    compute_quantiles_sweep
    compute_histogram
    compute_histogram_sweep

"""

from .radar_utils import time_avg_range, get_closest_solar_flux
from .radar_utils import create_sun_hits_field, create_sun_retrieval_field
from .radar_utils import compute_histogram, compute_histogram_sweep
from .radar_utils import compute_quantiles, compute_quantiles_sweep
from .radar_utils import compute_quantiles_from_hist

__all__ = [s for s in dir() if not s.startswith('_')]
