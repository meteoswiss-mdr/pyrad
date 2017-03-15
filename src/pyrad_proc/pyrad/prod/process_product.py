"""
pyrad.prod.process_product
==========================

Functions for obtaining Pyrad products from the datasets

.. autosummary::
    :toctree: generated/

    generate_cosmo_coord_products
    generate_sun_hits_products
    generate_intercomp_products
    generate_colocated_gates_products
    generate_time_avg_products
    generate_vol_products
    generate_timeseries_products
    generate_monitoring_products

"""

from copy import deepcopy
from warnings import warn

import numpy as np

import pyart

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename
from ..io.io_aux import generate_field_name_str

from ..io.read_data_other import get_sensor_data, read_timeseries
from ..io.read_data_other import read_sun_retrieval, read_monitoring_ts
from ..io.read_data_other import read_intercomp_scores_ts

from ..io.write_data import write_ts_polar_data, write_monitoring_ts
from ..io.write_data import write_sun_hits, write_sun_retrieval
from ..io.write_data import write_colocated_gates, write_colocated_data
from ..io.write_data import write_colocated_data_time_avg, write_cdf
from ..io.write_data import write_intercomp_scores_ts, write_ts_cum

from ..graph.plots import plot_ppi, plot_rhi, plot_cappi, plot_bscope
from ..graph.plots import plot_timeseries, plot_timeseries_comp
from ..graph.plots import plot_quantiles, get_colobar_label, plot_sun_hits
from ..graph.plots import plot_sun_retrieval_ts, plot_histogram
from ..graph.plots import plot_histogram2, plot_density, plot_monitoring_ts
from ..graph.plots import get_field_name, get_colobar_label, plot_scatter
from ..graph.plots import plot_intercomp_scores_ts, plot_scatter_comp

from ..util.radar_utils import create_sun_hits_field, rainfall_accumulation
from ..util.radar_utils import create_sun_retrieval_field, get_ROI
from ..util.radar_utils import compute_histogram, compute_quantiles
from ..util.radar_utils import compute_quantiles_from_hist, compute_2d_stats


def generate_cosmo_coord_products(dataset, prdcfg):
    """
    generates COSMO coordinates products

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
        print('saved file: '+fname)

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
    if ('dssavename' in prdcfg):
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
            timeinfo=dataset['timeinfo'], timeformat='%Y%m%d')

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        write_sun_hits(dataset['sun_hits'], fname[0])

        print('saved sun hits file: '+fname[0])

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

        fname = make_filename(
            'detected', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=dataset['timeinfo'],
            timeformat='%Y%m%d')

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

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

        plot_sun_hits(field, field_name, fname, prdcfg)

        print('----- save to '+' '.join(fname))

        return savedir+fname

    elif prdcfg['type'] == 'WRITE_SUN_RETRIEVAL':
        if 'sun_retrieval' not in dataset:
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=None)

        fname = make_filename(
            'info', prdcfg['dstype'], 'retrieval', ['csv'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        write_sun_retrieval(dataset['sun_retrieval'], fname[0])

        print('saved sun retrieval file: '+fname[0])

        return fname[0]

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

        fname = make_filename(
            'retrieval', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=dataset['timeinfo'],
            timeformat='%Y%m%d')

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        if dataset['sun_retrieval'][par] is None:
            warn(
                ' Invalid retrieval parameters. Skipping product ' +
                prdcfg['type'])
            return None

        field = create_sun_retrieval_field(
            dataset['sun_retrieval'][par], prdcfg['sunhitsImageConfig'])

        if field is not None:
            plot_sun_hits(field, field_name, fname, prdcfg)

        print('----- save to '+' '.join(fname))

        return fname

    elif prdcfg['type'] == 'PLOT_SUN_RETRIEVAL_TS':
        if 'sun_retrieval' not in dataset:
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=None)

        fname = make_filename(
            'info', prdcfg['dstype'], 'retrieval', ['csv'])

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

        fname = make_filename(
            'retrieval_ts', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        plot_sun_retrieval_ts(
            sun_retrieval, prdcfg['voltype'], fname)

        print('----- save to '+' '.join(fname))

        return fname

    else:
        if 'radar' in dataset:
            generate_vol_products(dataset['radar'], prdcfg)


def generate_intercomp_products(dataset, prdcfg):
    """
    generates radar intercomparison products

    Parameters
    ----------
    dataset : tuple
        values of colocated gates dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """

    dssavedir = prdcfg['dsname']
    if ('dssavename' in prdcfg):
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'WRITE_INTERCOMP':
        if dataset['final']:
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'colocated_data', prdcfg['dstype'], prdcfg['voltype'],
            ['csv'], timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        fname = savedir+fname[0]

        write_colocated_data(dataset['intercomp_dict'], fname)

        print('saved colocated data file: '+fname)

        return fname
    if prdcfg['type'] == 'WRITE_INTERCOMP_TIME_AVG':
        if dataset['final']:
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'colocated_data', prdcfg['dstype'], prdcfg['voltype'],
            ['csv'], timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        fname = savedir+fname[0]

        write_colocated_data_time_avg(dataset['intercomp_dict'], fname)

        print('saved colocated time averaged data file: '+fname)

        return fname
    elif prdcfg['type'] == 'PLOT_SCATTER_INTERCOMP':
        if not dataset['final']:
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname = make_filename(
            'scatter', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], timeinfo=dataset['timeinfo'],
            timeformat='%Y%m%d')

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        step = None
        if 'step' in prdcfg:
            step = prdcfg['step']

        hist_2d, bins1, bins2, stats = compute_2d_stats(
            np.ma.asarray(dataset['intercomp_dict']['rad1_val']),
            np.ma.asarray(dataset['intercomp_dict']['rad2_val']),
            field_name, field_name, step1=step, step2=step)
        if hist_2d is None:
            return None

        metadata = (
            'npoints: '+str(stats['npoints'])+'\n' +
            'mode bias: '+'{:.2f}'.format(float(stats['modebias']))+'\n' +
            'median bias: '+'{:.2f}'.format(float(stats['medianbias']))+'\n' +
            'mean bias: '+'{:.2f}'.format(float(stats['meanbias']))+'\n' +
            'intercep slope 1: '+'{:.2f}'.format(
                float(stats['intercep_slope_1']))+'\n' +
            'corr: '+'{:.2f}'.format(float(stats['corr']))+'\n' +
            'slope: '+'{:.2f}'.format(float(stats['slope']))+'\n' +
            'intercep: '+'{:.2f}'.format(float(stats['intercep']))+'\n')

        plot_scatter(bins1, bins2, np.ma.asarray(hist_2d), field_name,
                     field_name, fname, prdcfg, metadata=metadata,
                     lin_regr=[stats['slope'], stats['intercep']],
                     lin_regr_slope1=stats['intercep_slope_1'],
                     rad1_name=dataset['intercomp_dict']['rad1_name'],
                     rad2_name=dataset['intercomp_dict']['rad2_name'])

        print('----- save to '+' '.join(fname))

        return fname
    elif prdcfg['type'] == 'PLOT_AND_WRITE_INTERCOMP_TS':
        if not dataset['final']:
            return None

        field_name = get_fieldname_pyart(prdcfg['voltype'])
        step = None
        if 'step' in prdcfg:
            step = prdcfg['step']

        rad1_name = dataset['intercomp_dict']['rad1_name']
        rad2_name = dataset['intercomp_dict']['rad2_name']

        hist_2d, bins1, bins2, stats = compute_2d_stats(
            np.ma.asarray(dataset['intercomp_dict']['rad1_val']),
            np.ma.asarray(dataset['intercomp_dict']['rad2_val']),
            field_name, field_name, step1=step, step2=step)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=None)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'], ['csv'],
            timeinfo=None, timeformat='%Y%m%d')[0]

        csvfname = savedir+csvfname

        write_intercomp_scores_ts(
            dataset['timeinfo'], stats, field_name, csvfname,
            rad1_name=rad1_name, rad2_name=rad2_name)
        print('saved CSV file: '+csvfname)

        (date_vec, np_vec, meanbias_vec, medianbias_vec, modebias_vec,
         corr_vec, slope_vec, intercep_vec, intercep_slope1_vec) = (
            read_intercomp_scores_ts(csvfname))

        if date_vec is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        figfname = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=None, timeformat='%Y%m%d')

        for i in range(len(figfname)):
            figfname[i] = savedir+figfname[i]

        titl = rad1_name+'-'+rad2_name+' '+field_name+' intercomparison'
        plot_intercomp_scores_ts(
            date_vec, np_vec, meanbias_vec, medianbias_vec, modebias_vec,
            corr_vec, slope_vec, intercep_vec, intercep_slope1_vec, figfname,
            ref_value=0., labelx='Time UTC', titl=titl)
        print('----- save to '+' '.join(figfname))

        return figfname
    else:
        warn(' Unsupported product type: ' + prdcfg['type'])
        return None


def generate_colocated_gates_products(dataset, prdcfg):
    """
    generates colocated gates products

    Parameters
    ----------
    dataset : tuple
        radar objects and colocated gates dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """

    dssavedir = prdcfg['dsname']
    if ('dssavename' in prdcfg):
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'WRITE_COLOCATED_GATES':
        if prdcfg['radar'] not in dataset:
            return None
        if 'coloc_dict' not in dataset[prdcfg['radar']]:
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], 'colocated_gates',
            prdcfg['prdname'], timeinfo=None)

        fname = make_filename(
            'info', prdcfg['dstype'], prdcfg['prdname'], ['csv'],
            timeinfo=None)

        fname = savedir+fname[0]

        write_colocated_gates(
            dataset[prdcfg['radar']]['coloc_dict'], fname)

        print('saved colocated gates file: '+fname)

        return fname

    else:
        if prdcfg['radar'] not in dataset:
            return None
        if 'radar' not in dataset[prdcfg['radar']]:
            return None

        prdcfg['timeinfo'] = None
        generate_vol_products(dataset[prdcfg['radar']]['radar'], prdcfg)


def generate_time_avg_products(dataset, prdcfg):
    """
    generates time average products

    Parameters
    ----------
    dataset : tuple
        radar objects and colocated gates dictionary

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    filename : str
        the name of the file created. None otherwise

    """
    prdcfg['timeinfo'] = dataset['timeinfo']

    return generate_vol_products(dataset['radar_obj'], prdcfg)


def generate_vol_products(dataset, prdcfg):
    """
    generates radar volume products

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

    dssavedir = prdcfg['dsname']
    if ('dssavename' in prdcfg):
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'PPI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(dataset.fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(dataset.fixed_angle['data'] == el)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'ppi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=prdcfg['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        step = None
        quantiles = None
        plot_type = 'PPI'
        if 'plot_type' in prdcfg:
            plot_type = prdcfg['plot_type']
        if 'step' in prdcfg:
            step = prdcfg['step']
        if 'quantiles' in prdcfg:
            quantiles = prdcfg['quantiles']

        plot_ppi(dataset, field_name, ind_el, prdcfg, fname,
                 plot_type=plot_type, step=step, quantiles=quantiles)

        print('----- save to '+' '.join(fname))

        return fname

    elif prdcfg['type'] == 'RHI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        az_vec = np.sort(dataset.fixed_angle['data'])
        az = az_vec[prdcfg['anglenr']]
        ind_az = np.where(dataset.fixed_angle['data'] == az)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'rhi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='az'+'{:.1f}'.format(az),
            timeinfo=prdcfg['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        step = None
        quantiles = None
        plot_type = 'RHI'
        if 'plot_type' in prdcfg:
            plot_type = prdcfg['plot_type']
        if 'step' in prdcfg:
            step = prdcfg['step']
        if 'quantiles' in prdcfg:
            quantiles = prdcfg['quantiles']

        plot_rhi(dataset, field_name, ind_az, prdcfg, fname,
                 plot_type=plot_type, step=step, quantiles=quantiles)

        print('----- save to '+' '.join(fname))

        return fname

    elif prdcfg['type'] == 'PSEUDOPPI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        try:
            xsect = pyart.util.cross_section_rhi(
                dataset, [prdcfg['angle']], el_tol=prdcfg['EleTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname = make_filename(
                'ppi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='el'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'])

            for i in range(len(fname)):
                fname[i] = savedir+fname[i]

            step = None
            quantiles = None
            plot_type = 'PPI'
            if 'plot_type' in prdcfg:
                plot_type = prdcfg['plot_type']
            if 'step' in prdcfg:
                step = prdcfg['step']
            if 'quantiles' in prdcfg:
                quantiles = prdcfg['quantiles']

            plot_ppi(xsect, field_name, 0, prdcfg, fname,
                     plot_type=plot_type, step=step, quantiles=quantiles)

            print('----- save to '+' '.join(fname))

            return fname
        except EnvironmentError:
            warn(
                'No data found at elevation ' + str(prdcfg['angle']) +
                '. Skipping product ' + prdcfg['type'])

            return None

    elif prdcfg['type'] == 'PSEUDORHI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        try:
            xsect = pyart.util.cross_section_ppi(
                dataset, [prdcfg['angle']], az_tol=prdcfg['AziTol'])

            savedir = get_save_dir(
                prdcfg['basepath'], prdcfg['procname'], dssavedir,
                prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

            fname = make_filename(
                'rhi', prdcfg['dstype'], prdcfg['voltype'],
                prdcfg['imgformat'],
                prdcfginfo='az'+'{:.1f}'.format(prdcfg['angle']),
                timeinfo=prdcfg['timeinfo'])

            for i in range(len(fname)):
                fname[i] = savedir+fname[i]

            step = None
            quantiles = None
            plot_type = 'RHI'
            if 'plot_type' in prdcfg:
                plot_type = prdcfg['plot_type']
            if 'step' in prdcfg:
                step = prdcfg['step']
            if 'quantiles' in prdcfg:
                quantiles = prdcfg['quantiles']

            plot_rhi(xsect, field_name, 0, prdcfg, fname,
                     plot_type=plot_type, step=step, quantiles=quantiles)

            print('----- save to '+' '.join(fname))

            return fname
        except EnvironmentError:
            warn(
                ' No data found at azimuth ' +
                str(prdcfg['angle'])+'. Skipping product ' +
                prdcfg['type'])
            return None

    elif prdcfg['type'] == 'CAPPI_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'cappi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='alt'+'{:.1f}'.format(prdcfg['altitude']),
            timeinfo=prdcfg['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        plot_cappi(dataset, field_name, prdcfg['altitude'], prdcfg, fname)
        print('----- save to '+' '.join(fname))

        return fname

    elif prdcfg['type'] == 'BSCOPE_IMAGE':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        ang_vec = np.sort(dataset.fixed_angle['data'])
        ang = ang_vec[prdcfg['anglenr']]
        ind_ang = np.where(dataset.fixed_angle['data'] == ang)[0][0]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'b-scope', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            prdcfginfo='ang'+'{:.1f}'.format(ang),
            timeinfo=prdcfg['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        plot_bscope(dataset, field_name, ind_ang, prdcfg, fname)
        print('----- save to '+' '.join(fname))

        return fname

    elif prdcfg['type'] == 'HISTOGRAM':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        step = None
        if 'step' in prdcfg:
            step = prdcfg['step']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'histogram', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        bins, values = compute_histogram(
            dataset.fields[field_name]['data'], field_name, step=step)

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        labelx = get_colobar_label(dataset.fields[field_name], field_name)

        plot_histogram(bins, values, fname, labelx=labelx,
                       labely='Number of Samples', titl=titl)

        print('----- save to '+' '.join(fname))

        return fname

    elif prdcfg['type'] == 'QUANTILES':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        quantiles = None
        if 'quantiles' in prdcfg:
            quantiles = prdcfg['quantiles']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'quantiles', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        quantiles, values = compute_quantiles(
            dataset.fields[field_name]['data'], quantiles=quantiles)

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        labely = get_colobar_label(dataset.fields[field_name], field_name)

        plot_quantiles(quantiles, values, fname, labelx='quantile',
                       labely=labely, titl=titl)

        print('----- save to '+' '.join(fname))

        return fname

    elif prdcfg['type'] == 'CDF':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        quantiles = None
        if 'quantiles' in prdcfg:
            quantiles = prdcfg['quantiles']

        sector = {
            'rmin': None,
            'rmax': None,
            'azmin': None,
            'azmax': None,
            'elmin': None,
            'elmax': None,
            'hmin': None,
            'hmax': None}

        if 'sector' in prdcfg:
            if 'rmin' in prdcfg['sector']:
                sector['rmin'] = prdcfg['sector']['rmin']
            if 'rmax' in prdcfg['sector']:
                sector['rmax'] = prdcfg['sector']['rmax']
            if 'azmin' in prdcfg['sector']:
                sector['azmin'] = prdcfg['sector']['azmin']
            if 'azmax' in prdcfg['sector']:
                sector['azmax'] = prdcfg['sector']['azmax']
            if 'elmin' in prdcfg['sector']:
                sector['elmin'] = prdcfg['sector']['elmin']
            if 'elmax' in prdcfg['sector']:
                sector['elmax'] = prdcfg['sector']['elmax']
            if 'hmin' in prdcfg['sector']:
                sector['hmin'] = prdcfg['sector']['hmin']
            if 'hmax' in prdcfg['sector']:
                sector['hmax'] = prdcfg['sector']['hmax']

        vismin = None
        if 'vismin' in prdcfg:
            vismin = prdcfg['vismin']

        absolute = False
        if 'absolute' in prdcfg:
            absolute = prdcfg['absolute']

        use_nans = False
        nan_value = 0.
        if 'use_nans' in prdcfg:
            use_nans = prdcfg['use_nans']
            if 'nan_value' in prdcfg:
                nan_value = prdcfg['nan_value']

        filterclt = False
        if 'filterclt' in prdcfg:
            filterclt = prdcfg['filterclt']

        filterprec = []
        if 'filterprec' in prdcfg:
            filterprec = prdcfg['filterprec']

        data = deepcopy(dataset.fields[field_name]['data'])

        # define region of interest
        roi_flag = get_ROI(dataset, field_name, sector)
        data = data[roi_flag == 1]

        ntot = np.size(roi_flag[roi_flag == 1])

        if ntot == 0:
            warn('No radar gates found in sector')
            return None

        # get number of gates with clutter and mask them
        nclut = -1
        if filterclt:
            echoID_field = get_fieldname_pyart('echoID')
            if echoID_field in dataset.fields:
                echoID_ROI = dataset.fields[echoID_field]['data'][
                    roi_flag == 1]
                nclut = len(echoID_ROI[echoID_ROI == 2])
                data[echoID_ROI == 2] = np.ma.masked

        # get number of blocked gates and filter according to visibility
        nblocked = -1
        if vismin is not None:
            vis_field = get_fieldname_pyart('VIS')
            if vis_field in dataset.fields:
                vis_ROI = dataset.fields[vis_field]['data'][roi_flag == 1]
                nblocked = len(vis_ROI[vis_ROI < vismin])
                data[vis_ROI < vismin] = np.ma.masked

        # filter according to precip type
        nprec_filter = -1
        if len(filterprec) > 0:
            hydro_field = get_fieldname_pyart('hydro')
            if hydro_field in dataset.fields:
                hydro_ROI = dataset.fields[hydro_field]['data'][roi_flag == 1]
                nprec_filter = 0
                for i in range(len(filterprec)):
                    nprec_filter += len(hydro_ROI[hydro_ROI == filterprec[i]])
                    data[hydro_ROI == filterprec[i]] = np.ma.masked

        if absolute:
            data = np.ma.abs(data)

        mask = np.ma.getmaskarray(data)
        nnan = np.count_nonzero(mask)

        if nnan == ntot:
            warn('No valid radar gates found in sector')
            return None

        if use_nans:
            data[mask] = nan_value

        # count and filter outliers
        quantiles_lim, values_lim = compute_quantiles(
            data, quantiles=[0.2, 99.8])
        nsmall = np.count_nonzero(data.compressed() < values_lim[0])
        nlarge = np.count_nonzero(data.compressed() > values_lim[1])
        noutliers = nlarge+nsmall
        data = data[np.logical_and(
            data >= values_lim[0], data <= values_lim[1])]

        # number of values used for cdf computation
        ncdf = np.size(data.compressed())

        quantiles, values = compute_quantiles(data, quantiles=quantiles)

        # plot CDF
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'cdf', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=prdcfg['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        titl = (
            pyart.graph.common.generate_radar_time_begin(
                dataset).isoformat() + 'Z' + '\n' +
            get_field_name(dataset.fields[field_name], field_name))

        labelx = get_colobar_label(dataset.fields[field_name], field_name)

        plot_quantiles(values, quantiles/100., fname, labelx=labelx,
                       labely='Cumulative probability', titl=titl)

        print('----- save to '+' '.join(fname))

        # store cdf values
        fname = make_filename(
            'cdf', prdcfg['dstype'], prdcfg['voltype'],
            ['txt'], timeinfo=prdcfg['timeinfo'])[0]

        fname = savedir+fname

        write_cdf(
            quantiles, values, ntot, nnan, nclut, nblocked, nprec_filter,
            noutliers, ncdf, fname, use_nans=use_nans, nan_value=nan_value,
            filterprec=filterprec, vismin=vismin, sector=sector,
            datatype=labelx, timeinfo=prdcfg['timeinfo'])

        print('----- save to '+fname)

        return fname

    elif prdcfg['type'] == 'SAVEVOL':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        new_dataset = deepcopy(dataset)
        new_dataset.fields = dict()
        new_dataset.add_field(field_name, dataset.fields[field_name])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], prdcfg['voltype'], ['nc'],
            timeinfo=prdcfg['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        pyart.io.cfradial.write_cfradial(fname[0], new_dataset)
        print('saved file: '+fname[0])

        return fname[0]

    elif prdcfg['type'] == 'SAVEALL':
        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], 'all_fields', ['nc'],
            timeinfo=prdcfg['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        pyart.io.cfradial.write_cfradial(fname[0], dataset)
        print('saved file: '+fname[0])

        return fname[0]

    else:
        warn(' Unsupported product type: ' + prdcfg['type'])
        return None


def generate_timeseries_products(dataset, prdcfg):
    """
    generates time series products

    Parameters
    ----------
    dataset : dictionary
        radar object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    no return

    """

    dssavedir = prdcfg['dsname']
    if ('dssavename' in prdcfg):
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'PLOT_AND_WRITE_POINT':
        if dataset['final']:
            return None

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        for i in range(len(csvfname)):
            csvfname[i] = savedir+csvfname[i]

        write_ts_polar_data(dataset, csvfname[0])
        print('saved CSV file: '+csvfname[0])

        date, value = read_timeseries(csvfname[0])

        if date is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        figfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=date[0], timeformat='%Y%m%d')

        for i in range(len(figfname)):
            figfname[i] = savedir+figfname[i]

        label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        titl = ('Time Series '+date[0].strftime('%Y-%m-%d'))

        labely = generate_field_name_str(dataset['datatype'])

        plot_timeseries(
            date, [value], figfname, labelx='Time UTC',
            labely=labely, labels=[label1], title=titl)
        print('----- save to '+' '.join(figfname))

        return figfname

    elif prdcfg['type'] == 'PLOT_CUMULATIVE_POINT':
        if dataset['final']:
            return None

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        for i in range(len(csvfname)):
            csvfname[i] = savedir+csvfname[i]

        date, value = read_timeseries(csvfname[0])

        if date is None:
            warn(
                'Unable to plot accumulation time series. No valid data')
            return None

        figfname = make_filename(
            'ts_cum', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=date[0], timeformat='%Y%m%d')

        for i in range(len(figfname)):
            figfname[i] = savedir+figfname[i]

        label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        titl = ('Time Series Acc. '+date[0].strftime('%Y-%m-%d'))

        labely = 'Radar estimated rainfall accumulation (mm)'

        plot_timeseries(
            date, [value], figfname, labelx='Time UTC',
            labely=labely, labels=[label1], title=titl,
            period=prdcfg['ScanPeriod']*60.)
        print('----- save to '+' '.join(figfname))

        return figfname

    elif prdcfg['type'] == 'COMPARE_POINT':
        if dataset['final']:
            return None

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        for i in range(len(csvfname)):
            csvfname[i] = savedir_ts+csvfname[i]

        radardate, radarvalue = read_timeseries(csvfname[0])
        if radardate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid radar data')
            return None

        sensordate, sensorvalue, sensortype, period = get_sensor_data(
            radardate[0], dataset['datatype'], prdcfg)
        if sensordate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid sensor data')
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=radardate[0])

        figfname = make_filename(
            'ts_comp', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')

        for i in range(len(figfname)):
            figfname[i] = savedir+figfname[i]

        label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        label2 = sensortype+' '+prdcfg['sensorid']
        titl = 'Time Series Comp. '+radardate[0].strftime('%Y-%m-%d')
        labely = generate_field_name_str(dataset['datatype'])

        plot_timeseries_comp(
            radardate, radarvalue, sensordate, sensorvalue, figfname,
            labelx='Time UTC', labely=labely, label1=label1, label2=label2,
            titl=titl)
        print('----- save to '+' '.join(figfname))

        return figfname

    elif prdcfg['type'] == 'COMPARE_CUMULATIVE_POINT':
        if dataset['final']:
            return None

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=prdcfg['timeinfo'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=prdcfg['timeinfo'],
            timeformat='%Y%m%d')

        for i in range(len(csvfname)):
            csvfname[i] = savedir_ts+csvfname[i]

        radardate, radarvalue = read_timeseries(csvfname[0])
        if radardate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid radar data')
            return None

        sensordate, sensorvalue, sensortype, period2 = get_sensor_data(
            radardate[0], dataset['datatype'], prdcfg)
        if sensordate is None:
            warn(
                'Unable to plot sensor comparison at point of interest. ' +
                'No valid sensor data')
            return None

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=radardate[0])

        figfname = make_filename(
            'ts_cumcomp', prdcfg['dstype'], dataset['datatype'],
            prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')

        for i in range(len(figfname)):
            figfname[i] = savedir+figfname[i]

        label1 = 'Radar (az, el, r): ('+az+', '+el+', '+r+')'
        label2 = sensortype+' '+prdcfg['sensorid']
        titl = ('Time Series Acc. Comp. ' +
                radardate[0].strftime('%Y-%m-%d'))
        labely = 'Rainfall accumulation (mm)'

        plot_timeseries_comp(
            radardate, radarvalue, sensordate, sensorvalue,
            figfname, labelx='Time UTC', labely=labely,
            label1=label1, label2=label2, titl=titl,
            period1=prdcfg['ScanPeriod']*60., period2=period2)
        print('----- save to '+' '.join(figfname))

        return figfname

    elif prdcfg['type'] == 'COMPARE_TIME_AVG':
        if not dataset['final']:
            return None

        az = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][0])
        el = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][1])
        r = '{:.1f}'.format(dataset['antenna_coordinates_az_el_r'][2])
        gateinfo = ('az'+az+'r'+r+'el'+el)

        savedir_ts = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdid'], timeinfo=dataset['time'])

        csvfname = make_filename(
            'ts', prdcfg['dstype'], dataset['datatype'], ['csv'],
            prdcfginfo=gateinfo, timeinfo=dataset['time'],
            timeformat='%Y%m%d')[0]

        radardate, radarvalue = read_timeseries(savedir_ts+csvfname)
        if radardate is None:
            warn(
                'Unable to compared time averaged data at POI. ' +
                'No valid radar data')
            return None

        sensordate, sensorvalue, sensortype, period2 = get_sensor_data(
            radardate[0], dataset['datatype'], prdcfg)
        if sensordate is None:
            warn(
                'Unable to compared time averaged data at POI. ' +
                'No valid sensor data')
            return None

        cum_time = 3600
        if 'cum_time' in prdcfg:
            cum_time = prdcfg['cum_time']

        base_time = 0
        if 'base_time' in prdcfg:
            base_time = prdcfg['base_time']

        sensordate_cum, sensorvalue_cum, np_sensor_cum = rainfall_accumulation(
            sensordate, sensorvalue, cum_time=cum_time, base_time=base_time,
            dropnan=False)

        radardate_cum, radarvalue_cum, np_radar_cum = rainfall_accumulation(
            radardate, radarvalue, cum_time=cum_time, base_time=base_time,
            dropnan=False)

        # find common time stamps
        ind = np.where(np.in1d(radardate_cum, sensordate_cum))[0]
        if len(ind) == 0:
            warn('No sensor data for radar data time stamps')
        radardate_cum2 = radardate_cum[ind]
        radarvalue_cum2 = radarvalue_cum[ind]
        np_radar_cum2 = np_radar_cum[ind]

        ind = np.where(np.in1d(sensordate_cum, radardate_cum2))[0]
        sensorvalue_cum2 = sensorvalue_cum[ind]
        np_sensor_cum2 = np_sensor_cum[ind]

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=radardate[0])

        fname = make_filename(
            str(cum_time)+'s_acc_ts_comp', prdcfg['dstype'],
            dataset['datatype'], ['csv'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')[0]

        fname = savedir+fname

        new_dataset = deepcopy(dataset)
        new_dataset.update({
            'time': radardate_cum2,
            'sensor_value': sensorvalue_cum2,
            'np_sensor': np_sensor_cum2,
            'radar_value': radarvalue_cum2,
            'np_radar': np_radar_cum2,
            'cum_time': cum_time})
        new_dataset.update(prdcfg)

        write_ts_cum(new_dataset, fname)

        print('saved CSV file: '+fname)

        figfname = make_filename(
            str(cum_time)+'s_acc_ts_comp', prdcfg['dstype'],
            dataset['datatype'], prdcfg['imgformat'], prdcfginfo=gateinfo,
            timeinfo=radardate[0], timeformat='%Y%m%d')

        for i in range(len(figfname)):
            figfname[i] = savedir+figfname[i]

        labelx = sensortype+' '+prdcfg['sensorid']+' (mm)'
        labely = 'Radar (az, el, r): ('+az+', '+el+', '+r+') (mm)'
        titl = (str(cum_time)+'s Acc. Comp. ' +
                radardate_cum[0].strftime('%Y-%m-%d'))

        plot_scatter_comp(
            sensorvalue_cum2, radarvalue_cum2, figfname, labelx=labelx,
            labely=labely, titl=titl, axis='equal')

        print('----- save to '+' '.join(figfname))

        return figfname

    # ================================================================
    elif prdcfg['type'] == 'PLOT_AND_WRITE':

        timeinfo = dataset.time_vector[0]

        savedir = get_save_dir(prdcfg['basepath'], prdcfg['procname'],
                               dssavedir, prdcfg['prdname'],
                               timeinfo=timeinfo)

        dstype_str = prdcfg['dstype'].lower().replace('_', '')
        fname = make_filename('ts', dstype_str, dataset.datatype,
                              ['csv'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        dataset.write(savedir + fname[0])

        fname = make_filename('ts', dstype_str, dataset.datatype,
                              prdcfg['imgformat'],
                              prdcfginfo=None, timeinfo=timeinfo,
                              timeformat='%Y%m%d%H%M%S',
                              runinfo=prdcfg['runinfo'])

        ymin = None
        ymax = None
        if ('ymin' in prdcfg):
            ymin = prdcfg['ymin']
        if ('ymax' in prdcfg):
            ymax = prdcfg['ymax']

        dataset.plot(savedir + fname[0], ymin=ymin, ymax=ymax)

        return None

    # ================================================================
    else:
        raise Exception("ERROR: Unsupported product type: '%s' of dataset '%s'"
                        % (prdcfg['type'], prdcfg['dsname']))


def generate_monitoring_products(dataset, prdcfg):

    # check the type of dataset required
    hist_type = 'cumulative'
    if 'hist_type' in prdcfg:
        hist_type = prdcfg['hist_type']

    if dataset['hist_type'] != hist_type:
        return None

    hist_obj = dataset['hist_obj']

    dssavedir = prdcfg['dsname']
    if ('dssavename' in prdcfg):
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'VOL_HISTOGRAM':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        timeformat = '%Y%m%d'
        titl = (
                pyart.graph.common.generate_radar_time_begin(
                    hist_obj).strftime('%Y-%m-%d') + '\n' +
                get_field_name(hist_obj.fields[field_name], field_name))
        if hist_type == 'instant':
            timeformat = '%Y%m%d%H%M%S'
            titl = (
                pyart.graph.common.generate_radar_time_begin(
                    hist_obj).isoformat() + 'Z' + '\n' +
                get_field_name(hist_obj.fields[field_name], field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname = make_filename(
            'histogram', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=dataset['timeinfo'], timeformat=timeformat)

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        labelx = get_colobar_label(hist_obj.fields[field_name], field_name)

        plot_histogram2(
            hist_obj.range['data'],
            np.sum(hist_obj.fields[field_name]['data'], axis=0),
            fname, labelx=labelx, labely='Number of Samples',
            titl=titl)

        print('----- save to '+' '.join(fname))

        return fname

    if prdcfg['type'] == 'PPI_HISTOGRAM':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(hist_obj.fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(hist_obj.fixed_angle['data'] == el)[0][0]

        timeformat = '%Y%m%d'
        titl = (
                '{:.1f}'.format(el)+' Deg. ' +
                pyart.graph.common.generate_radar_time_begin(
                    hist_obj).strftime('%Y-%m-%d') + '\n' +
                get_field_name(hist_obj.fields[field_name], field_name))
        if hist_type == 'instant':
            timeformat = '%Y%m%d%H%M%S'
            titl = (
                '{:.1f}'.format(el)+' Deg. ' +
                pyart.graph.common.generate_radar_time_begin(
                    hist_obj).isoformat() + 'Z' + '\n' +
                get_field_name(hist_obj.fields[field_name], field_name))

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname = make_filename(
            'ppi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=dataset['timeinfo'], timeformat=timeformat)

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        labelx = get_colobar_label(hist_obj.fields[field_name], field_name)

        sweep_start = hist_obj.sweep_start_ray_index['data'][ind_el]
        sweep_end = hist_obj.sweep_end_ray_index['data'][ind_el]
        values = hist_obj.fields[field_name]['data'][sweep_start:sweep_end, :]
        plot_histogram2(
            hist_obj.range['data'], np.sum(values, axis=0),
            fname, labelx=labelx, labely='Number of Samples',
            titl=titl)

        print('----- save to '+' '.join(fname))

        return fname

    if prdcfg['type'] == 'ANGULAR_DENSITY':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        el_vec = np.sort(hist_obj.fixed_angle['data'])
        el = el_vec[prdcfg['anglenr']]
        ind_el = np.where(hist_obj.fixed_angle['data'] == el)[0][0]

        timeformat = '%Y%m%d'
        if hist_type == 'instant':
            timeformat = '%Y%m%d%H%M%S'

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname = make_filename(
            'ppi', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='el'+'{:.1f}'.format(el),
            timeinfo=dataset['timeinfo'], timeformat=timeformat)

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        quantiles = np.array([25., 50., 75.])
        ref_value = 0.
        if 'quantiles' in prdcfg:
            quantiles = prdcfg['quantiles']
        if 'ref_value' in prdcfg:
            ref_value = prdcfg['ref_value']

        plot_density(
            hist_obj, hist_type, field_name, ind_el, prdcfg, fname,
            quantiles=quantiles, ref_value=ref_value)

        print('----- save to '+' '.join(fname))

        return fname

    elif prdcfg['type'] == 'VOL_TS':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        csvtimeinfo = None
        if hist_type == 'instant':
            csvtimeinfo = prdcfg['timeinfo']

        quantiles = np.array([25., 50., 75.])
        ref_value = 0.
        if 'quantiles' in prdcfg:
            quantiles = prdcfg['quantiles']
        if 'ref_value' in prdcfg:
            ref_value = prdcfg['ref_value']

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=csvtimeinfo)

        csvfname = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'], ['csv'],
            timeinfo=csvtimeinfo, timeformat='%Y%m%d')[0]

        csvfname = savedir+csvfname

        quantiles, values = compute_quantiles_from_hist(
            hist_obj.range['data'],
            np.ma.sum(hist_obj.fields[field_name]['data'], axis=0),
            quantiles=quantiles)

        start_time = pyart.graph.common.generate_radar_time_begin(hist_obj)
        np_t = np.ma.sum(hist_obj.fields[field_name]['data'], dtype=int)
        if np.ma.getmaskarray(np_t):
            np_t = 0

        write_monitoring_ts(
            start_time, np_t, values, quantiles, prdcfg['voltype'],
            csvfname)
        print('saved CSV file: '+csvfname)

        date, np_t_vec, cquant_vec, lquant_vec, hquant_vec = (
            read_monitoring_ts(csvfname))

        if date is None:
            warn(
                'Unable to plot time series. No valid data')
            return None

        figtimeinfo = None
        titldate = ''
        if hist_type == 'instant':
            figtimeinfo = date[0]
            titldate = date[0].strftime('%Y-%m-%d')

        figfname = make_filename(
            'ts', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'],
            timeinfo=figtimeinfo, timeformat='%Y%m%d')

        for i in range(len(figfname)):
            figfname[i] = savedir+figfname[i]

        titl = ('Monitoring Time Series '+titldate)

        labely = generate_field_name_str(prdcfg['voltype'])

        plot_monitoring_ts(
            date, np_t_vec, cquant_vec, lquant_vec, hquant_vec, field_name,
            figfname, ref_value=ref_value, labelx='Time UTC',
            labely=labely, titl=titl)
        print('----- save to '+' '.join(figfname))

        return figfname

    elif prdcfg['type'] == 'SAVEVOL':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in hist_obj.fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        new_dataset = deepcopy(hist_obj)
        new_dataset.fields = dict()
        new_dataset.add_field(field_name, hist_obj.fields[field_name])

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=dataset['timeinfo'])

        fname = make_filename(
            'savevol', prdcfg['dstype'], prdcfg['voltype'], ['nc'],
            timeinfo=dataset['timeinfo'])

        for i in range(len(fname)):
            fname[i] = savedir+fname[i]

        pyart.io.cfradial.write_cfradial(fname[0], new_dataset)
        print('saved file: '+fname[0])

        return fname[0]

    else:
        warn(' Unsupported product type: ' + prdcfg['type'])
        return None
