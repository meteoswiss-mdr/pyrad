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
    read_status

Reading configuration files
===========================

.. autosummary::
    :toctree: generated/

    read_config

Reading other data
==================

.. autosummary::
    :toctree: generated/

    read_timeseries

"""

from .config import read_config

from .read_data import get_data, read_status, read_timeseries
from .read_data import get_datetime, get_datasetfields, get_file_list
from .read_data import get_datatypefields

from .write_data import write_timeseries, generate_field_name_str

__all__ = [s for s in dir() if not s.startswith('_')]
