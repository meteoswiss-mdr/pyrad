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

Reading other data
==================

.. autosummary::
    :toctree: generated/

    read_status
    read_rad4alp_cosmo
    read_timeseries
    get_sensor_data
    read_smn
    read_disdro_scattering

Writing data
==================

.. autosummary::
    :toctree: generated/

    write_timeseries
    write_sun_hits
    generate_field_name_str

"""

from .config import read_config

from .read_data_radar import get_data, get_datetime, get_datasetfields
from .read_data_radar import get_file_list, get_datatypefields

from .read_data_aux import read_status, read_timeseries, read_smn
from .read_data_aux import get_sensor_data, read_rad4alp_cosmo, read_sun_hits
from .read_data_aux import read_selfconsistency

from .write_data import write_timeseries, generate_field_name_str
from .write_data import write_sun_hits, write_sun_retrieval

__all__ = [s for s in dir() if not s.startswith('_')]
