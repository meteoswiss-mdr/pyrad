"""
pyrad.prod.process_product
==========================

Functions for obtaining Pyrad products from the datasets

.. autosummary::
    :toctree: generated/

    generate_occurrence_products
    generate_cosmo_coord_products
    generate_sun_hits_products
    generate_qvp_products

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from .process_vol_products import generate_vol_products

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename

from ..io.read_data_other import read_sun_retrieval

from ..io.write_data import write_sun_hits, write_sun_retrieval
from ..io.write_data import write_excess_gates

from ..graph.plots import plot_sun_hits, plot_sun_retrieval_ts

from ..util.radar_utils import create_sun_hits_field
from ..util.radar_utils import create_sun_retrieval_field


def generate_occurrence_products(dataset, prdcfg):
    """
    generates occurrence products

    Parameters
    ----------
    dataset : tuple
        radar object and metadata dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    instant = False
    if 'instant' in prdcfg:
        instant = prdcfg['instant']

    if not instant and not dataset['occu_final']:
        return None

    if prdcfg['type'] == 'WRITE_EXCESS_GATES':
        if not dataset['occu_final']:
            return None

        radar = dataset['radar_obj']
        if (('frequency_of_occurrence' not in radar.fields) or
                ('occurrence' not in radar.fields) or
                ('number_of_samples' not in radar.fields)):
            warn('Unable to create quantile excess gates file. '
                 'Missing data')
            return None

        dssavedir = prdcfg['dsname']
        if 'dssavename' in prdcfg:
            dssavedir = prdcfg['dssavename']

        quant_min = 95.
        if 'quant_min' in prdcfg:
            quant_min = prdcfg['quant_min']

        # get index of gates exceeding quantile
        freq_occu = radar.fields['frequency_of_occurrence'][
            'data']
        ind_ray, ind_rng = np.where(freq_occu > quant_min)
        if not ind_ray:
            warn('No data exceeds the frequency of occurrence ' +
                 str(quant_min)+' %')
            return None

        excess_dict = {
            'starttime': dataset['starttime'],
            'endtime': dataset['endtime'],
            'quant_min': quant_min,
            'ray_ind': ind_ray,
            'rng_ind': ind_rng,
            'ele': radar.elevation['data'][ind_ray],
            'azi': radar.azimuth['data'][ind_ray],
            'rng': radar.range['data'][ind_rng],
            'nsamples': (
                radar.fields['number_of_samples']['data'][ind_ray, ind_rng]),
            'occurrence': (
                radar.fields['occurrence']['data'][ind_ray, ind_rng]),
            'freq_occu': freq_occu[ind_ray, ind_rng]
        }
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['endtime'])

        fname = make_filename(
            'excess_gates', prdcfg['dstype'], prdcfg['prdname'], ['csv'],
            prdcfginfo='quant'+'{:.1f}'.format(quant_min),
            timeinfo=dataset['endtime'])

        fname = savedir+fname[0]

        fname = write_excess_gates(excess_dict, fname)

        if fname is not None:
            print('saved excess gates file: '+fname)

        return fname
    else:
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if ((field_name == 'frequency_of_occurrence') and
                (not dataset['occu_final'])):
            return None
        if dataset['occu_final']:
            prdcfg['timeinfo'] = dataset['endtime']

        return generate_vol_products(dataset['radar_obj'], prdcfg)


def generate_cosmo_coord_products(dataset, prdcfg):
    """
    generates COSMO coordinates products

    Parameters
    ----------
    dataset : tuple
        radar object containing the COSMO coordinates

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    if prdcfg['type'] == 'SAVEVOL':
        radar_obj = dataset['radar_obj']
        ind_rad = dataset['ind_rad']

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in radar_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        new_dataset = deepcopy(radar_obj)
        new_dataset.fields = dict()
        new_dataset.add_field(field_name, radar_obj.fields[field_name])

        savedir = prdcfg['cosmopath'][ind_rad]+'rad2cosmo1/'
        fname = 'rad2cosmo_'+prdcfg['voltype']+'_'+prdcfg['procname']+'.nc'

        pyart.io.cfradial.write_cfradial(savedir+fname, new_dataset)
        print('saved file: '+savedir+fname)

        return fname

    else:
        warn(' Unsupported product type: ' + prdcfg['type'])
        return None


def generate_sun_hits_products(dataset, prdcfg):
    """
    generates sun hits products

    Parameters
    ----------
    dataset : tuple
        radar object and sun hits dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    prdcfg['timeinfo'] = dataset['timeinfo']

    if prdcfg['type'] == 'WRITE_SUN_HITS':
        if 'sun_hits' not in dataset:
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname = make_filename(
            'info', prdcfg['dstype'], 'detected', ['csv'],
            timeinfo=dataset['timeinfo'], timeformat='%Y%m%d')[0]

        fname = savedir+fname

        write_sun_hits(dataset['sun_hits'], fname)

        print('saved sun hits file: '+fname)

        return fname[0]

    elif prdcfg['type'] == 'PLOT_SUN_HITS':
        if 'sun_hits_final' not in dataset:
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])

        if prdcfg['voltype'] not in dataset['sun_hits_final']:
            warn(
                ' Field type ' + prdcfg['voltype'] +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname_list = make_filename(
            'detected', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=dataset['timeinfo'],
            timeformat='%Y%m%d')

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        field = create_sun_hits_field(
            dataset['sun_hits_final']['rad_el'],
            dataset['sun_hits_final']['rad_az'],
            dataset['sun_hits_final']['sun_el'],
            dataset['sun_hits_final']['sun_az'],
            dataset['sun_hits_final'][prdcfg['voltype']],
            prdcfg['sunhitsImageConfig'])

        if field is None:
            warn(
                'Unable to create field '+prdcfg['voltype'] +
                ' Skipping product ' + prdcfg['type'])
            return None

        plot_sun_hits(field, field_name, fname_list, prdcfg)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'WRITE_SUN_RETRIEVAL':
        if 'sun_retrieval' not in dataset:
            return None

        timeinfo = None
        timeformat = None
        if 'add_date_in_fname' in prdcfg:
            if prdcfg['add_date_in_fname']:
                timeinfo = dataset['timeinfo']
                timeformat = '%Y'

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=None)

        fname = make_filename(
            'info', prdcfg['dstype'], 'retrieval', ['csv'], timeinfo=timeinfo,
            timeformat=timeformat, runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname

        write_sun_retrieval(dataset['sun_retrieval'], fname)

        print('saved sun retrieval file: '+fname)

        return fname

    elif prdcfg['type'] == 'PLOT_SUN_RETRIEVAL':
        if 'sun_retrieval' not in dataset:
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        par = None
        if field_name == 'sun_est_power_h':
            par = 'par_h'
        elif field_name == 'sun_est_power_v':
            par = 'par_v'
        elif field_name == 'sun_est_differential_reflectivity':
            par = 'par_zdr'

        if par not in dataset['sun_retrieval']:
            warn(
                ' Field type ' + prdcfg['voltype'] +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname_list = make_filename(
            'retrieval', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=dataset['timeinfo'],
            timeformat='%Y%m%d')

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        if dataset['sun_retrieval'][par] is None:
            warn(
                ' Invalid retrieval parameters. Skipping product ' +
                prdcfg['type'])
            return None

        field = create_sun_retrieval_field(
            dataset['sun_retrieval'][par], field_name,
            prdcfg['sunhitsImageConfig'],
            lant=dataset['sun_retrieval']['lant'])

        if field is not None:
            plot_sun_hits(field, field_name, fname_list, prdcfg)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    elif prdcfg['type'] == 'PLOT_SUN_RETRIEVAL_TS':
        if 'sun_retrieval' not in dataset:
            return None

        dpi = 72
        if 'dpi' in prdcfg:
            dpi = prdcfg['dpi']

        timeinfo = None
        timeformat = None
        if 'add_date_in_fname' in prdcfg:
            if prdcfg['add_date_in_fname']:
                timeinfo = dataset['timeinfo']
                timeformat = '%Y'

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=None)

        fname = make_filename(
            'info', prdcfg['dstype'], 'retrieval', ['csv'], timeinfo=timeinfo,
            timeformat=timeformat, runinfo=prdcfg['runinfo'])

        fname = savedir + fname[0]

        sun_retrieval = read_sun_retrieval(fname)

        if sun_retrieval[0] is None:
            warn(
                'Unable to read sun retrieval file '+fname)
            return None

        if len(sun_retrieval[0]) < 2:
            warn(
                'Unable to plot sun retrieval time series. ' +
                'Not enough data points.')
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=None)

        fname_list = make_filename(
            'retrieval_ts', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=timeinfo,
            timeformat=timeformat, runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        titl = (prdcfg['runinfo']+' Sun Retrieval ' +
                sun_retrieval[1][0].strftime('%Y%m%d')+'-' +
                sun_retrieval[1][-1].strftime('%Y%m%d'))
        figfname = plot_sun_retrieval_ts(
            sun_retrieval, prdcfg['voltype'], fname_list, titl=titl, dpi=dpi)

        if figfname is None:
            return None

        print('----- save to '+' '.join(fname_list))
        return fname_list

    else:
        if 'radar' in dataset:
            generate_vol_products(dataset['radar'], prdcfg)


def generate_qvp_products(dataset, prdcfg):
    """
    Generates quasi vertical profile products. Quasi vertical profiles
    come from azimuthal averaging of polarimetric radar data.

    Parameters
    ----------
    dataset : dict
        dictionary containing the radar object and a keyword stating the
        status of the processing

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    qvp_type = 'final'
    if 'qvp_type' in prdcfg:
        qvp_type = prdcfg['qvp_type']

    if qvp_type == 'final' and dataset['radar_type'] != 'final':
        return None

    prdcfg['timeinfo'] = dataset['start_time']
    return generate_vol_products(dataset['radar_obj'], prdcfg)
