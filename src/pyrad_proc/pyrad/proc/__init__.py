"""
======================================================
Dataset and products processing (:mod:`pyrad.proc`)
======================================================

.. currentmodule:: pyrad.proc

Initiate the dataset and products processing.

Dataset processing
==================

.. autosummary::
    :toctree: generated/

    process_raw

Product processing
=================

.. autosummary::
    :toctree: generated/

    generate_vol_products

"""

from .process_dataset import get_process_type, process_raw, process_echo_id
from .process_dataset import process_snr, process_l, process_cdr
from .process_dataset import process_correct_noise_rhohv
from .process_dataset import process_echo_filter, process_attenuation
from .process_dataset import process_rainrate, process_hydroclass

from .process_product import get_product_type, generate_vol_products

__all__ = [s for s in dir() if not s.startswith('_')]
