"""
pyrad.prod.process_product
==========================

Functions for obtaining Pyrad products from the datasets

.. autosummary::
    :toctree: generated/

    generate_traj_product

"""

import pyart

from ..io.io_aux import get_save_dir, make_filename


def generate_traj_product(traj, prdcfg):
    """
    Generates trajectory products

    Parameters
    ----------
    traj : Trajectory object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    None

    """

    if prdcfg['type'] == 'TRAJ_PLOT':

        print("%s: to be done!!!" % prdcfg['type'])  # XXX
        return None
    elif prdcfg['type'] == 'TRAJ_TEXT':

        timeinfo = traj.timevector[0]

        savedir = get_save_dir(prdcfg['basepath'], prdcfg['procname'],
                               prdcfg['dsname'], prdcfg['prdname'],
                               timeinfo=timeinfo)

        fname = make_filename('ts', prdcfg['dstype'], 'TRAJ', ['csv'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        return None

    else:
        raise Exception("ERROR: Unsupported product type: '%s' of dataset '%s'"
                        % (prdcfg['type'], prdcfg['dsname']))
