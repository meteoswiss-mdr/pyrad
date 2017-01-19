"""
pyrad.prod.process_product
==========================

Functions for obtaining Pyrad products from the datasets

.. autosummary::
    :toctree: generated/

    generate_traj_product

"""

from ..io.io_aux import get_save_dir, make_filename
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

    dssavedir = prdcfg['dsname']
    if ('dssavename' in prdcfg):
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'TRAJ_PLOT':

        timeinfo = traj.time_vector[0]

        savedir = get_save_dir(prdcfg['basepath'], prdcfg['procname'],
                               dssavedir, prdcfg['prdname'],
                               timeinfo=timeinfo)

        ts = TimeSeries("", traj.time_vector,
                        timeformat="%Y-%m-%d %H:%M:%S.%f")

        if (prdcfg['datatype'] == 'EL'):
            fname = make_filename('ts', prdcfg['dstype'], 'TRAJ',
                                  prdcfg['imgformat'],
                                  prdcfginfo="EL", timeinfo=timeinfo,
                                  timeformat='%Y%m%d%H%M%S',
                                  runinfo=prdcfg['runinfo'])

            ts.add_dataseries("Elevation", "Elevation", "deg",
                              traj.radar_list[0].elevation_vec)
            ts.plot(savedir + fname[0])

        elif (prdcfg['datatype'] == 'AZ'):
            fname = make_filename('ts', prdcfg['dstype'], 'TRAJ',
                                  prdcfg['imgformat'],
                                  prdcfginfo="AZ", timeinfo=timeinfo,
                                  timeformat='%Y%m%d%H%M%S',
                                  runinfo=prdcfg['runinfo'])

            ts.add_dataseries("Azimuth", "Azimuth", "deg", traj.radar_list[0].azimuth_vec)
            ts.plot(savedir + fname[0])

        elif (prdcfg['datatype'] == 'RANGE'):
            fname = make_filename('ts', prdcfg['dstype'], 'TRAJ',
                                  prdcfg['imgformat'],
                                  prdcfginfo="RANGE", timeinfo=timeinfo,
                                  timeformat='%Y%m%d%H%M%S',
                                  runinfo=prdcfg['runinfo'])

            ts.add_dataseries("Range", "Range", "m", traj.radar_list[0].range_vec)
            ts.plot(savedir + fname[0])

        else:
            raise Exception("ERROR: Unknown datatype '%s' (dataset: '%s')" %
                            (prdcfg['datatype'], prdcfg['dsname']))

        return None

    elif prdcfg['type'] == 'TRAJ_TEXT':

        timeinfo = traj.time_vector[0]

        savedir = get_save_dir(prdcfg['basepath'], prdcfg['procname'],
                               dssavedir, prdcfg['prdname'],
                               timeinfo=timeinfo)

        fname = make_filename('ts', prdcfg['dstype'], 'TRAJ', ['csv'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        description = ["Description:",
                       "Time series of a plane trajectory in radar "
                       "coordinates."]

        ts = TimeSeries(description, traj.time_vector,
                        timeformat="%Y-%m-%d %H:%M:%S.%f")
        ts.add_dataseries("Elevation", "deg", traj.radar_list[0].elevation_vec)
        ts.add_dataseries("Azimuth", "deg", traj.radar_list[0].azimuth_vec)
        ts.add_dataseries("Range", "m", traj.radar_list[0].range_vec)

        ts.add_dataseries("Absolute Speed", "m/s", traj.radar_list[0].v_abs)
        ts.add_dataseries("Radial Speed", "m/s", traj.radar_list[0].v_r)
        ts.add_dataseries("Elevation Speed", "deg/s", traj.radar_list[0].v_el)
        ts.add_dataseries("Azimuth Speed", "deg/s", traj.radar_list[0].v_az)

        ts.write(savedir + fname[0])

        return None

    else:
        raise Exception("ERROR: Unsupported product type: '%s' of dataset '%s'"
                        % (prdcfg['type'], prdcfg['dsname']))
