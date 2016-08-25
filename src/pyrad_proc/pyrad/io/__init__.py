"""
==================================
Input and output (:mod:`pyrad.io`)
==================================

.. currentmodule:: pyrad.io

Functions to read and write data and configuration files.

Reading radar data
==================

.. autosummary::
    :toctree: generated/

    get_data

Reading configuration files
===========================

.. autosummary::
    :toctree: generated/

    read_config

"""

from .config import read_config
from .read_data import get_data, get_datatypefields, get_file_list
from .read_data import get_datetime, get_datasetfields

__all__ = [s for s in dir() if not s.startswith('_')]
