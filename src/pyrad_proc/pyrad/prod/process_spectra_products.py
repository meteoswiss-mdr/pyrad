"""
pyrad.prod.process_grid_products
================================

Functions for obtaining Pyrad products from spectra datasets

.. autosummary::
    :toctree: generated/

    generate_spectra_products

"""

from warnings import warn

from ..io.io_aux import get_fieldname_pyart
from ..io.io_aux import get_save_dir, make_filename

from ..graph.plots_spectra import plot_range_Doppler, plot_Doppler
from ..graph.plots_spectra import plot_complex_range_Doppler
from ..graph.plots_spectra import plot_amp_phase_range_Doppler
from ..graph.plots_spectra import plot_complex_Doppler, plot_amp_phase_Doppler


def generate_spectra_products(dataset, prdcfg):
    """
    generates grid products. Accepted product types:
        'COMPLEX_RANGE_DOPPLER': Plots the complex spectra range-Doppler
            User defined parameters:
                coord1, coord2: dict
                    The two lat-lon coordinates marking the limits. They have
                    the keywords 'lat' and 'lon' [degree]. The altitude limits
                    are defined by the parameters in 'rhiImageConfig' in the
                    'loc' configuration file

    Parameters
    ----------
    dataset : spectra
        spectra object

    prdcfg : dictionary of dictionaries
        product configuration dictionary of dictionaries

    Returns
    -------
    None or name of generated files

    """

    dssavedir = prdcfg['dsname']
    if 'dssavename' in prdcfg:
        dssavedir = prdcfg['dssavename']

    if prdcfg['type'] == 'RANGE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        ray = prdcfg.get('ray', 0)
        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'surface', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='ray'+str(ray),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_range_Doppler(
            dataset['radar_out'], field_name, ray, prdcfg, fname_list,
            xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        ray = prdcfg.get('ray', 0)
        rng = prdcfg.get('rng', 0)
        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'surface', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='ray'+str(ray)+'rng'+str(rng),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_Doppler(
            dataset['radar_out'], field_name, ray, rng, prdcfg, fname_list,
            xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'COMPLEX_RANGE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        ray = prdcfg.get('ray', 0)
        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'surface', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='ray'+str(ray),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_complex_range_Doppler(
            dataset['radar_out'], field_name, ray, prdcfg, fname_list,
            xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'COMPLEX_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        ray = prdcfg.get('ray', 0)
        rng = prdcfg.get('rng', 0)
        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        vmin = prdcfg.get('vmin', None)
        vmax = prdcfg.get('vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'surface', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='ray'+str(ray)+'rng'+str(rng),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_complex_Doppler(
            dataset['radar_out'], field_name, ray, rng, prdcfg, fname_list,
            xaxis_info=xaxis_info, vmin=vmin, vmax=vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'AMPLITUDE_PHASE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        ray = prdcfg.get('ray', 0)
        rng = prdcfg.get('rng', 0)
        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        ampli_vmin = prdcfg.get('ampli_vmin', None)
        ampli_vmax = prdcfg.get('ampli_vmax', None)
        phase_vmin = prdcfg.get('phase_vmin', None)
        phase_vmax = prdcfg.get('phase_vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'surface', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='ray'+str(ray)+'rng'+str(rng),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_amp_phase_Doppler(
            dataset['radar_out'], field_name, ray, rng, prdcfg, fname_list,
            xaxis_info=xaxis_info, ampli_vmin=ampli_vmin,
            ampli_vmax=ampli_vmax, phase_vmin=phase_vmin,
            phase_vmax=phase_vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    if prdcfg['type'] == 'AMPLITUDE_PHASE_RANGE_DOPPLER':
        field_name = get_fieldname_pyart(prdcfg['voltype'])
        if field_name not in dataset['radar_out'].fields:
            warn(
                ' Field type ' + field_name +
                ' not available in data set. Skipping product ' +
                prdcfg['type'])
            return None

        # user defined values
        ray = prdcfg.get('ray', 0)
        xaxis_info = prdcfg.get('xaxis_info', 'Doppler_velocity')
        ampli_vmin = prdcfg.get('ampli_vmin', None)
        ampli_vmax = prdcfg.get('ampli_vmax', None)
        phase_vmin = prdcfg.get('phase_vmin', None)
        phase_vmax = prdcfg.get('phase_vmax', None)

        savedir = get_save_dir(
            prdcfg['basepath'], prdcfg['procname'], dssavedir,
            prdcfg['prdname'], timeinfo=prdcfg['timeinfo'])

        fname_list = make_filename(
            'surface', prdcfg['dstype'], prdcfg['voltype'],
            prdcfg['imgformat'], prdcfginfo='ray'+str(ray),
            timeinfo=prdcfg['timeinfo'], runinfo=prdcfg['runinfo'])

        for i, fname in enumerate(fname_list):
            fname_list[i] = savedir+fname

        plot_amp_phase_range_Doppler(
            dataset['radar_out'], field_name, ray, prdcfg, fname_list,
            xaxis_info=xaxis_info, ampli_vmin=ampli_vmin,
            ampli_vmax=ampli_vmax, phase_vmin=phase_vmin,
            phase_vmax=phase_vmax)

        print('----- save to '+' '.join(fname_list))

        return fname_list

    warn(' Unsupported product type: ' + prdcfg['type'])
    return None
