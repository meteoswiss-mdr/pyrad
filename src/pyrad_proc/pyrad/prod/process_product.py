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
    generate_ml_products

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from .process_vol_products import generate_vol_products

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename

from ..io.read_data_sun import read_sun_retrieval
from ..io.read_data_other import read_ml_ts

from ..io.write_data import write_sun_hits, write_sun_retrieval
from ..io.write_data import write_excess_gates, write_ts_ml

from ..graph.plots import plot_sun_hits
from ..graph.plots_timeseries import plot_sun_retrieval_ts, plot_ml_ts

from ..util.radar_utils import create_sun_hits_field
from ..util.radar_utils import create_sun_retrieval_field


def generate_occurrence_products(dataset, prdcfg):
    """
    generates occurrence products. Accepted product types:
        'WRITE_EXCESS_GATES': Write the data that identifies radar gates
            with clutter that has a frequency of occurrence above a certain
            threshold.
            User defined parameters:
                quant_min: float
                    Minimum frequency of occurrence in percentage to keep the
                    gate as valid. Default 95.
        All the products of the 'VOL' dataset group

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

        radar = dataset['radar_out']
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
        if ind_ray.size == 0:
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

    field_name = get_fieldname_pyart(prdcfg['voltype'])
    if ((field_name == 'frequency_of_occurrence') and
            (not dataset['occu_final'])):
        return None
    if dataset['occu_final']:
        prdcfg['timeinfo'] = dataset['endtime']

    return generate_vol_products(dataset, prdcfg)


def generate_cosmo_coord_products(dataset, prdcfg):
    """
    generates COSMO coordinates products. Accepted product types:
        'SAVEVOL': Save an object containing the index of the COSMO model grid
            that corresponds to each radar gate in a C/F radial file.

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
        radar_obj = dataset['radar_out']
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

    warn(' Unsupported product type: ' + prdcfg['type'])
    return None


def generate_sun_hits_products(dataset, prdcfg):
    """
    generates sun hits products. Accepted product types:
        'PLOT_SUN_HITS': Plots in a sun-radar azimuth difference-sun-radar
            elevation difference grid the values of all sun hits obtained
            during the processing period
        'PLOT_SUN_RETRIEVAL': Plots in a sun-radar azimuth difference-sun-
            radar elevation difference grid the retrieved sun pattern
        'PLOT_SUN_RETRIEVAL_TS': Plots time series of the retrieved sun
            pattern parameters
            User defined parameters:
                dpi: int
                    The pixel density of the plot. Default 72
                add_date_in_fname: Bool
                    If true the year is added in the plot file name
        'WRITE_SUN_HITS': Writes the information concerning possible sun hits
            in a csv file
        'WRITE_SUN_RETRIEVAL': Writes the retrieved sun pattern parameters in
            a csv file.
            User defined parameters:
                add_date_in_fname: Bool
                    If true the year is added in the csv file name
        All the products of the 'VOL' dataset group

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

    if prdcfg['type'] == 'PLOT_SUN_HITS':
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

    if prdcfg['type'] == 'WRITE_SUN_RETRIEVAL':
        if 'sun_retrieval' not in dataset:
            return None

        timeinfo = None
        timeformat = None
        if prdcfg.get('add_date_in_fname', False):
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

    if prdcfg['type'] == 'PLOT_SUN_RETRIEVAL':
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

    if prdcfg['type'] == 'PLOT_SUN_RETRIEVAL_TS':
        if 'sun_retrieval' not in dataset:
            return None

        dpi = prdcfg.get('dpi', 72)
        timeinfo = None
        timeformat = None
        if prdcfg.get('add_date_in_fname', False):
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

    if 'radar_out' in dataset:
        return generate_vol_products(dataset, prdcfg)

    return None


def generate_qvp_products(dataset, prdcfg):
    """
    Generates quasi vertical profile-like products. Quasi vertical profiles
    come from azimuthal averaging of polarimetric radar data. With the
    variable 'qvp_type' the user decides if the product has to be generated
    at the end of the processing period ('final') or instantaneously
    ('instant')
    Accepted product types:
        All the products of the 'VOL' dataset group

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
    return generate_vol_products(dataset, prdcfg)


def generate_ml_products(dataset, prdcfg):
    """
    Generates melting layer products. Accepted product types:
        'ML_TS': Plots and writes a time series of the melting layer, i.e.
            the evolution of the average and standard deviation of the melting
            layer top and thickness and the the number of rays used in the
            retrieval.
            User defined parameters:
                dpi: int
                    The pixel density of the plot. Default 72
        'SAVE_ML': Saves an object containing the melting layer retrieval
            information in a C/F radial file
        All the products of the 'VOL' dataset group

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

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'ML_TS':
        dpi = prdcfg.get('dpi', 72)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], 'ml', ['csv'],
            timeinfo=prdcfg['timeinfo'], timeformat='%Y%m%d')[0]

        csvfname = savedir+csvfname

        ml_bottom = dataset['ml_obj'].fields['melting_layer_height']['data'][:, 0]
        ml_top = dataset['ml_obj'].fields['melting_layer_height']['data'][:, 1]

        ml_top_avg = np.ma.asarray(np.ma.mean(ml_top))
        ml_top_std = np.ma.asarray(np.ma.std(ml_top))
        thick = ml_top-ml_bottom
        thick_avg = np.ma.asarray(np.ma.mean(thick))
        thick_std = np.ma.asarray(np.ma.std(thick))
        nrays_valid = thick.compressed().size
        nrays_total = thick.size

        write_ts_ml(
            prdcfg['timeinfo'], ml_top_avg, ml_top_std, thick_avg, thick_std,
            nrays_valid, nrays_total, csvfname)

        print('saved CSV file: '+csvfname)

        (dt_ml_arr, ml_top_avg_arr, ml_top_std_arr, thick_avg_arr,
         thick_std_arr, nrays_valid_arr, nrays_total_arr) = (
             read_ml_ts(csvfname))

        if dt_ml_arr is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        figfname_list = make_filename(
            'ts', prdcfg['dstype'], 'ml', prdcfg['imgformat'],
            timeinfo=dt_ml_arr[0], timeformat='%Y%m%d')

        for i, figfname in enumerate(figfname_list):
            figfname_list[i] = savedir+figfname

        titl = dt_ml_arr[0].strftime('%Y-%m-%d')+' melting layer time series'

        plot_ml_ts(
            dt_ml_arr, ml_top_avg_arr, ml_top_std_arr, thick_avg_arr,
            thick_std_arr, nrays_valid_arr, nrays_total_arr, figfname_list,
            labelx='Time UTC', titl=titl, dpi=dpi)
        print('----- save to '+' '.join(figfname_list))

        return figfname_list

    if prdcfg['type'] == 'SAVE_ML':
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'saveml', prdcfg['dstype'], 'ml_h', ['nc'],
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])[0]

        fname = savedir+fname
        pyart.io.cfradial.write_cfradial(fname, dataset['ml_obj'])
        print('saved file: '+fname)

        return fname

    return generate_vol_products(dataset, prdcfg)
