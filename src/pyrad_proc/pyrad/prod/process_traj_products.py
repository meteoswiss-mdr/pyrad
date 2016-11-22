"""
pyrad.prod.process_product
==========================

Functions for obtaining Pyrad products from the datasets

.. autosummary::
    :toctree: generated/

    generate_traj_products

"""

from warnings import warn

import pyart


def generate_traj_products(dataset, prdcfg):
    """
    Generates trajectory products

    Parameters
    ----------
    dataset : Radar
        radar object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """

    if prdcfg['type'] == 'TRAJ_PLOT':

        print("%s: to be done!!!" % prdcfg['type'])  # XXX
        return None
    elif prdcfg['type'] == 'TRAJ_TEXT':

        print("%s: to be done!!!" % prdcfg['type'])  # XXX
        return None

    else:
        raise Exception("ERROR: Unsupported product type: '%s' of dataset '%s'"
                        % (prdcfg['type'], prdcfg['dsname']))
