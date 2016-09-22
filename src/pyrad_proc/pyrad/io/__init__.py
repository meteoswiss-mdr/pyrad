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
    read_rad4alp_cosmo

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
    get_sensor_data
    read_smn
    read_disdro_scattering

"""

from .config import read_config

from .read_data import get_data, read_status, read_timeseries, read_smn
from .read_data import get_datetime, get_datasetfields, get_file_list
from .read_data import get_datatypefields, get_sensor_data, read_rad4alp_cosmo
from .read_data import read_selfconsistency

from .write_data import write_timeseries, generate_field_name_str

__all__ = [s for s in dir() if not s.startswith('_')]
