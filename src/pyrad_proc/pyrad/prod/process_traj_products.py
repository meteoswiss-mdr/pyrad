"""
pyrad.prod.process_product
==========================

Functions for obtaining Pyrad products from the datasets

.. autosummary::
    :toctree: generated/

    generate_traj_product

"""

from ..io.io_aux import get_save_dir, make_filename
from ..io.write_data import write_timeseries
from ..io.timeseries import TimeSeries


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

        description = ["Description:",
                       "Time series of a plane trajectory in radar coordinates."]

        ts = TimeSeries(description, traj.timevector, timeformat="%Y-%m-%d %H:%M:%S.%f")
        ts.add_dataseries("Elevation", "deg", traj.radar_list[0].elevation_vec)
        ts.add_dataseries("Azimuth", "deg", traj.radar_list[0].azimuth_vec)
        ts.add_dataseries("Range", "m", traj.radar_list[0].range_vec)

        write_timeseries(ts, savedir + fname[0])

        return None

    else:
        raise Exception("ERROR: Unsupported product type: '%s' of dataset '%s'"
                        % (prdcfg['type'], prdcfg['dsname']))
