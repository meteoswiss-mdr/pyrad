"""
===========================================
processing flow control (:mod:`pyrad.flow`)
===========================================

.. currentmodule:: pyrad.flow

Functions to control the Pyrad data processing flow

.. autosummary::
    :toctree: generated/

    main

"""

from .flow_control import main

__all__ = [s for s in dir() if not s.startswith('_')]
