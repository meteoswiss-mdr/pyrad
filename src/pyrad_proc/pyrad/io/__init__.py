"""
==================================
Input and output (:mod:`pyrad.io`)
==================================

.. currentmodule:: pyrad.io

Functions to read and write data and configuration files.

Reading configuration files
===========================

.. autosummary::
    :toctree: generated/

    read_config

Reading radar data
==================

.. autosummary::
    :toctree: generated/

    get_data
    
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
    read_sun_hits
    read_sun_hits_multiple_days
    read_sun_retrieval
    read_selfconsistency

Writing data
==================

.. autosummary::
    :toctree: generated/

    write_timeseries
    write_sun_hits
    write_sun_retrieval


Auxiliary functions
===================

.. autosummary::
    :toctree: generated/

    get_save_dir
    make_filename
    get_datetime
    get_datasetfields
    get_file_list
    get_datatypefields    
    get_fieldname_rainbow    
    generate_field_name_str

"""

from .config import read_config

from .read_data_radar import get_data

from .read_data_other import read_status, read_rad4alp_cosmo, read_timeseries
from .read_data_other import get_sensor_data, read_smn, read_disdro_scattering
from .read_data_other import read_sun_hits, read_sun_hits_multiple_days
from .read_data_other import read_sun_retrieval
from .read_data_other import read_selfconsistency

from .write_data import write_timeseries
from .write_data import write_sun_hits, write_sun_retrieval

from .io_aux import get_save_dir, make_filename
from .io_aux import get_datetime, get_dataset_fields
from .io_aux import get_file_list, get_datatype_fields
from .io_aux import get_fieldname_pyart
from .io_aux import generate_field_name_str

__all__ = [s for s in dir() if not s.startswith('_')]
