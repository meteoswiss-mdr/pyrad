#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
main_trt
================================================

This program processes TRT data

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import glob
import os
from shutil import copy

import numpy as np

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt

from pyrad.io import read_trt_scores, read_trt_cell_lightning
from pyrad.util import belongs_roi_indices
from pyrad.graph import plot_scatter_comp, plot_pos, plot_timeseries

print(__doc__)


def main():
    """
    """
    basepath = '/store/msrad/radar/trt/data_analysis_new/'

    roi = {
        'lon': [8.9000010, 9.2000000, 9.4999970, 9.4999970, 8.9000010],
        'lat': [47.0000030, 47.0000030, 47.0000030, 47.5999930, 47.5999930]
    }

    # plot_rpc_time_series(basepath, 'cell_rimed_particles_column.csv')
    # plot_rimed_column_EUCLID_pn(basepath, 'cell_rimed_particles_column.csv', 'Santis_cell_euclid_np_lightning.csv', roi)
    # plot_rimed_column_EUCLID(basepath, 'cell_rimed_particles_column.csv', 'Santis_cell_euclid_lightning.csv', roi)
    # plot_rimed_column_LMA(basepath, 'cell_rimed_particles_column.csv', 'cell_LMA_flashes.csv', roi)
    # plot_max_LMA_rank(basepath, 'cell_LMA_flashes.csv')
    # plot_time_diff_max_EUCLID_density_max_rank(basepath, 'Santis_cell_scores.csv')
    # plot_time_diff_max_LMA_max_rank(basepath, 'cell_LMA_flashes.csv', roi= roi)
    # plot_time_diff_max_EUCLID_strokes_max_rank(basepath, 'Santis_cell_euclid_lightning.csv')
    # read_general_scores(basepath, 'Santis_cell_scores.csv')
    # plot_rank_LMA_flashes(basepath, 'cell_LMA_flashes.csv', roi=roi)
    # plot_rank_LMA_sources(basepath, 'cell_LMA_sources.csv', roi=roi)
    plot_LMA_EUCLID(basepath, 'cell_LMA_flashes.csv', 'Santis_cell_euclid_lightning.csv', roi)


    return

    traj_ID_uniques = np.unique(traj_ID, return_index=False)

    max_rank_time = np.array([])
    max_nflashes_time = np.array([])
    max_flash_dens_time = np.array([])
    max_nflashes = np.array([])
    max_flash_dens = np.array([])

    max_nflashes_rank = np.array([])
    max_flash_dens_rank = np.array([])

    tdiff_dens = np.array([])
    tdiff_flash = np.array([])
    for traj_ID_unique in traj_ID_uniques:
        inds = np.where(traj_ID == traj_ID_unique)[0]

        flash_dens_aux = flash_dens_cell[inds]
        nflashes_aux = nflashes_cell[inds]
        time_aux = time_cell[inds]
        rank_aux = rank_cell[inds]

        if flash_dens_aux.max() > 0:
            continue

        max_rank = rank_aux.max()

        print('Cell '+str(traj_ID_unique)+' max rank: '+str(max_rank))



#        # keep only time steps within LMA domain
#        _, is_roi = belongs_roi_indices(lat_cell[inds], lon_cell[inds], roi)
#
#        if is_roi == 'None' or is_roi == 'Some':
#            continue
#
#        max_rank = rank_aux.max()
#        inds_max_rank = np.where(rank_aux == max_rank)[0]
#        if inds_max_rank.size > 1:
#            print(str(inds_max_rank.size)+' time steps with rank max in cell '+str(traj_ID_unique))
#
#        max_nflashes = nflashes_aux.max()
#        inds_max_nflashes = np.where(nflashes_aux == max_nflashes)[0]
#        if inds_max_nflashes.size > 1:
#            print(str(inds_max_nflashes.size)+' time steps with nflashes max in cell '+str(traj_ID_unique))
#
#        max_flash_dens = flash_dens_aux.max()
#        inds_max_flash_dens = np.where(flash_dens_aux == max_flash_dens)[0]
#        if inds_max_flash_dens.size > 1:
#            print(str(inds_max_flash_dens.size)+' time steps with flash density max in cell '+str(traj_ID_unique))
#
#        print('\n')
#
#
#
#        max_rank_time = np.append(max_rank_time, time_cell[inds][np.argmax(rank_cell[inds])])
#        max_nflashes_time = np.append(max_nflashes_time, time_cell[inds][np.argmax(nflashes_cell[inds])])
#        max_flash_dens_time = np.append(max_flash_dens_time, time_cell[inds][np.argmax(flash_dens_cell[inds])])
#        max_nflashes = np.append(max_nflashes, nflashes_cell[inds].max())
#        max_flash_dens = np.append(max_flash_dens, flash_dens_cell[inds].max())
#
#        max_nflashes_rank = np.append(max_nflashes_rank, rank_cell[inds][np.argmax(nflashes_cell[inds])])
#        max_flash_dens_rank = np.append(max_flash_dens_rank, rank_cell[inds][np.argmax(flash_dens_cell[inds])])
#
#        tdiff_dens = np.append(tdiff_dens, (time_cell[inds][np.argmax(rank_cell[inds])]- time_cell[inds][np.argmax(flash_dens_cell[inds])]).total_seconds())
#        tdiff_flash = np.append(tdiff_flash, (time_cell[inds][np.argmax(rank_cell[inds])]- time_cell[inds][np.argmax(nflashes_cell[inds])]).total_seconds())


    return


    print(np.size(tdiff_dens[tdiff_dens > 0.]))
    print(np.size(tdiff_dens[tdiff_dens == 0.]))
    print(np.size(tdiff_dens[tdiff_dens < 0.]))

    print(np.size(tdiff_flash[tdiff_flash > 0.]))
    print(np.size(tdiff_flash[tdiff_flash == 0.]))
    print(np.size(tdiff_flash[tdiff_flash < 0.]))



    return





#    fname = basepath+'cell_LMA_flashes.csv'
#    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
#     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(fname)
#
#    inds, is_roi = belongs_roi_indices(lat_cell, lon_cell, roi)
#
#    print(inds.size)
#
#    traj_ID = traj_ID[inds]
#    time_cell = time_cell[inds]
#    rank_cell = rank_cell[inds]
#    flash_dens_cell = flash_dens_cell[inds]
#    nflashes_cell = nflashes_cell[inds]
#
#    time_cell_unique = np.unique(time_cell, return_index=False)
#    traj_ID_unique = np.unique(traj_ID, return_index=False)
#    print(time_cell_unique.size)
#    print(traj_ID_unique.size)


    fname = basepath+'cell_euclid_lightning.csv'
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(fname)

    traj_ID_unique = np.unique(traj_ID, return_index=False)
    print('Number of TRT cells in analysis: '+str(traj_ID_unique.size))

    inds, is_roi = belongs_roi_indices(lat_cell, lon_cell, roi)

    print(inds.size)

    traj_ID = traj_ID[inds]
    time_cell = time_cell[inds]
    rank_cell = rank_cell[inds]
    flash_dens_cell = flash_dens_cell[inds]
    nflashes_cell = nflashes_cell[inds]

    traj_ID_unique = np.unique(traj_ID, return_index=False)
    print('Number of TRT cells in LMA domain: '+str(traj_ID_unique.size))

    time_cell_unique = np.unique(time_cell, return_index=False)
    traj_ID_unique = np.unique(traj_ID, return_index=False)
    print(time_cell_unique.size)
    print(traj_ID_unique.size)

#    plot_scatter_comp(
#        flash_dens_cell[inds], rank_cell[inds]/10., [basepath+'LMA_domain_LMA_flash_density'],
#        labelx='flash density [flashes/km2]', labely='rank',
#        titl='LMA flash density vs Rank', axis=None, metadata=None, dpi=72)
#
#    print("----- plot to '%s'" % basepath+'LMA_flash_density')

    fname = basepath+'cell_LMA_sources.csv'
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(fname)

    traj_ID_unique = np.unique(traj_ID, return_index=False)
    print('Number of TRT cells with LMA data: '+str(traj_ID_unique.size))

    inds, is_roi = belongs_roi_indices(lat_cell, lon_cell, roi)

    print(inds.size)

    traj_ID = traj_ID[inds]
    time_cell = time_cell[inds]
    rank_cell = rank_cell[inds]
    flash_dens_cell = flash_dens_cell[inds]
    nflashes_cell = nflashes_cell[inds]

    traj_ID_unique_LMA = np.unique(traj_ID, return_index=False)
    print('Number of TRT cells with LMA data in LMA domain: '+str(traj_ID_unique_LMA.size))

#
#    plot_scatter_comp(
#        flash_dens_cell[inds], rank_cell[inds]/10., [basepath+'LMA_domain_LMA_sources_density'],
#        labelx='Sources density [sources/km2]', labely='rank',
#        titl='LMA sources density vs Rank', axis=None, metadata=None, dpi=72)
#
#    print("----- plot to '%s'" % basepath+'LMA_sources_density')

    fname = basepath+'cell_rimmed_particles_column.csv'
    (traj_ID_rad, time_cell_rad, lon_cell_rad, lat_cell_rad, area_cell_rad, rank_cell_rad,
     rh_min, rh_max) = read_trt_cell_lightning(fname)

    traj_ID_unique = np.unique(traj_ID_rad, return_index=False)
    print('Number of TRT cells with radar data: '+str(traj_ID_unique.size))



    print(inds.size)

    traj_ID_rad = traj_ID_rad[inds]
    time_cell_rad = time_cell_rad[inds]
    rank_cell_rad = rank_cell_rad[inds]
    rh_min = rh_min[inds]
    rh_max = rh_max[inds]
    rh_thickness = rh_max-rh_min

    traj_ID_unique_rad = np.unique(traj_ID_rad, return_index=False)
    print('Number of TRT cells with radar data in LMA domain: '+str(traj_ID_unique_rad.size))

    print(traj_ID_unique_LMA[np.isin(traj_ID_unique_LMA, traj_ID_rad, invert=True)])

    flash_dens_common = np.asarray([])
    rh_thickness_common = np.asarray([])
    rh_min_common = np.asarray([])
    nflashes_cell_common = np.asarray([])
    traj_ID_common = np.asarray([])
    for i, time_cell_rad_el in enumerate(time_cell_rad):
        ind = np.ma.where(np.logical_and(
            time_cell == time_cell_rad_el,
            traj_ID == traj_ID_rad[i]))[0]
        if ind.size == 0:
            continue
        flash_dens_common = np.append(flash_dens_common, flash_dens_cell[ind])
        rh_thickness_common = np.append(rh_thickness_common, rh_thickness[i])
        rh_min_common = np.append(rh_min_common, rh_min[i])
        nflashes_cell_common = np.append(nflashes_cell_common, nflashes_cell[ind])
        traj_ID_common = np.append(traj_ID_common, traj_ID_rad[i])

    traj_ID_unique = np.unique(traj_ID_common, return_index=False)
    print('Number of TRT cells with radar and LMA data in LMA domain: '+str(traj_ID_unique.size))


#
#
#
#    plot_pos(
#        rh_min_common, rh_thickness_common, flash_dens_common, [basepath+'LMA_domain_Rimmed_column_height-thickness_vs_Euclid_flashes_density_no_zeros'], ax=None, fig=None, save_fig=True,
#        sort_altitude='No', dpi=72, alpha=1., cb_label='Euclid CG flash density [Flashes/Km2]',
#        titl='Rimmed column vs Euclid CG flash density',
#        xlabel='Rimmed column thicknes [m]', ylabel='Rimmed column height [m MSL]',
#        limits=None, vmin=None, vmax=None)
#
#    plot_pos(
#        rh_min_common, rh_thickness_common, nflashes_cell_common, [basepath+'LMA_domain_Rimmed_column_height-thickness_vs_Euclid_flashes_no_zeros'], ax=None, fig=None, save_fig=True,
#        sort_altitude='No', dpi=72, alpha=1., cb_label='Euclid CG flashes',
#        titl='Rimmed column vs Euclid CG flashes',
#        xlabel='Rimmed column thicknes [m]', ylabel='Rimmed column height [m MSL]',
#        limits=None, vmin=None, vmax=20)
#
#
#    print(rh_thickness_common.size)
#    plot_scatter_comp(
#        rh_thickness_common, flash_dens_common, [basepath+'LMA_domain_Rimmed_column_height_vs_Euclid_flashes_density'],
#        labelx='Rimmed column height [m]', labely='Euclid CG flash density [flashes/Km2]',
#        titl='Rimmed column height vs Euclid CG flash density', axis=None, metadata=None, dpi=72)
#
#    print("----- plot to '%s'" % basepath+'Rimmed_column_height')
#
#    plot_scatter_comp(
#        rh_thickness_common, nflashes_cell_common, [basepath+'LMA_domain_Rimmed_column_height_vs_Euclid_flashes'],
#        labelx='Rimmed column height [m]', labely='Euclid CG flashes',
#        titl='Rimmed column height vs Euclid CG flashes', axis=None, metadata=None, dpi=72)
#
#    print("----- plot to '%s'" % basepath+'Rimmed_column_height')
#
#
#    plot_scatter_comp(
#        rh_min_common, flash_dens_common, [basepath+'LMA_domain_Rimmed_column_min_height_vs_Euclid_flashes_density'],
#        labelx='Rimmed column min height [m MSL]', labely='Euclid CG flash density [flashes/Km2]',
#        titl='Rimmed column min height vs Euclid CG flash density', axis=None, metadata=None, dpi=72)
#
#    print("----- plot to '%s'" % basepath+'Rimmed_column_height')
#
#    plot_scatter_comp(
#        rh_min_common, nflashes_cell_common, [basepath+'LMA_domain_Rimmed_column_min_height_vs_Euclid_flashes'],
#        labelx='Rimmed column min height [m MSL]', labely='Euclid CG flashes',
#        titl='Rimmed column min height vs Euclid CG flashes', axis=None, metadata=None, dpi=72)
#
#    print("----- plot to '%s'" % basepath+'Rimmed_column_height')






#    plot_scatter_comp(
#        rh_max[inds]-rh_min[inds], rank_cell[inds]/10., [basepath+'LMA_domain_Rimmed_column_height'],
#        labelx='Rimmed column height [m]', labely='rank',
#        titl='Rimmed column height vs Rank', axis=None, metadata=None, dpi=72)
#
#    print("----- plot to '%s'" % basepath+'Rimmed_column_height')




#    (traj_ID, time_flash_density_max, flash_density_max_rank,
#     flash_density_max_nflashes, flash_density_max_area, flash_density_max,
#     time_rank_max, rank_max) = read_trt_scores(fname)
#
#    plot_scatter_comp(
#        flash_density_max, flash_density_max_rank/10., [basepath+'flash_density'],
#        labelx='flash density max [flashes/km2]', labely='rank',
#        titl='Max flash density vs Rank', axis=None, metadata=None, dpi=72)
#
#    plot_scatter_comp(
#        flash_density_max, rank_max/10., [basepath+'flash_density_rank'],
#        labelx='flash density max [flashes/km2]', labely='rank max',
#        titl='Max flash density vs Max rank', axis=None, metadata=None, dpi=72)
#
#    plot_scatter_comp(
#        flash_density_max, (rank_max-flash_density_max_rank)/10., [basepath+'flash_density_rank_diff'],
#        labelx='flash density max [flashes/km2]', labely='rank max-rank',
#        titl='Max flash density vs Max rank - rank at max flash time', axis=None, metadata=None, dpi=72)
#
#    tdiff = np.empty(rank_max.size)
#    for i, tmax_rank in enumerate(time_rank_max):
#        tdiff[i] = (tmax_rank-time_flash_density_max[i]).total_seconds()
#    plot_scatter_comp(
#        flash_density_max, tdiff, [basepath+'flash_density_rank_time_diff'],
#        labelx='flash density max [flashes/km2]', labely='Time rank max - Time rank [s]',
#        titl='Max flash density vs Max rank time - rank at max flash time', axis=None, metadata=None, dpi=72)
#
#    traj_ID_sorted = traj_ID[np.argsort(flash_density_max)]
#    print(traj_ID_sorted[::-1])


def plot_rank_EUCLID_CG(basepath, fname):
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(basepath+fname)

    print('Cell steps: ', rank_cell.size)

    # plot data
    fname = basepath+'rank-EUCLID_CG_number.png'
    plot_scatter_comp(
        rank_cell/10., nflashes_cell, [fname],
        labelx='rank', labely='N Strokes',
        titl='Cell rank vs EUCLID CG strokes', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'rank-EUCLID_CG_density.png'
    plot_scatter_comp(
        rank_cell/10., flash_dens_cell, [fname],
        labelx='rank', labely='Strokes density [strokes/km^2]',
        titl='Cell rank vs EUCLID CG strokes density', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def plot_rimed_column_LMA(basepath, fname_rimed, fname_LMA, roi):
    # read rimed particles file
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     rm_hmin, rm_hmax) = read_trt_cell_lightning(basepath+fname_rimed)

    print('Cell steps: ', rank_cell.size)

    # keep only time steps within LMA domain
    inds, is_roi = belongs_roi_indices(lat_cell, lon_cell, roi)

    traj_ID = traj_ID[inds]
    time_cell = time_cell[inds]
    rank_cell = rank_cell[inds]
    rm_hmin = rm_hmin[inds]
    rm_hmax = rm_hmax[inds]

    # put non-valid min and max to 0
    rm_hmin = rm_hmin.filled(0.)
    rm_hmax = rm_hmax.filled(0.)


#    # keep only valid time steps
#    valid = np.logical_not(np.logical_or(
#        np.ma.getmaskarray(rm_hmin), np.ma.getmaskarray(rm_hmax)))
#
#    traj_ID = traj_ID[valid]
#    time_cell = time_cell[valid]
#    rank_cell = rank_cell[valid]
#    rm_hmax = rm_hmax[valid]
#    rm_hmin = rm_hmin[valid]

    print('Cell steps within LMA domain: ', rank_cell.size)

    # read LMA particles column file
    (traj_ID_LMA, time_cell_LMA, lon_cell_LMA, lat_cell_LMA, area_cell_LMA, rank_cell_LMA,
     nflashes_cell_LMA, flash_dens_cell_LMA) = read_trt_cell_lightning(
        basepath+fname_LMA)

    print('Cell steps: ', rank_cell_LMA.size)

    # keep only time steps within LMA domain
    inds, is_roi = belongs_roi_indices(lat_cell_LMA, lon_cell_LMA, roi)

    traj_ID_LMA = traj_ID_LMA[inds]
    time_cell_LMA = time_cell_LMA[inds]
    rank_cell_LMA = rank_cell_LMA[inds]
    nflashes_cell_LMA = nflashes_cell_LMA[inds]
    flash_dens_cell_LMA = flash_dens_cell_LMA[inds]

    print('Cell steps within LMA domain: ', rank_cell_LMA.size)

    # Get data from common TRT cells and time steps
    rm_hmin_common = np.asarray([])
    rm_hmax_common = np.asarray([])

    flash_dens_LMA_common = np.ma.asarray([])
    nflashes_LMA_common = np.ma.asarray([])

    for i, time_cell_LMA_el in enumerate(time_cell_LMA):
        ind = np.ma.where(np.logical_and(
            time_cell == time_cell_LMA_el,
            traj_ID == traj_ID_LMA[i]))[0]
        if ind.size == 0:
            continue
        rm_hmin_common = np.append(rm_hmin_common, rm_hmin[ind])
        rm_hmax_common = np.append(rm_hmax_common, rm_hmax[ind])

        flash_dens_LMA_common = np.append(flash_dens_LMA_common, flash_dens_cell_LMA[i])
        nflashes_LMA_common = np.append(nflashes_LMA_common, nflashes_cell_LMA[i])


    print('Common cell steps: ', rm_hmin_common.size)

    # plot data
    fname = basepath+'LMA_domain_rpc_min_height_vs_LMA_flashes.png'
    plot_scatter_comp(
        rm_hmin_common, nflashes_LMA_common, [fname],
        labelx='RPC base height [m MSL]', labely='N flashes',
        titl='RPC base height vs LMA flashes', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'LMA_domain_rpc_min_height_vs_LMA_flashes_density.png'
    plot_scatter_comp(
        rm_hmin_common, flash_dens_LMA_common, [fname],
        labelx='RPC base height [m MSL]', labely='Flash density [flashes/km^2]',
        titl='RPC base height vs LMA flashes density', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    # plot data
    fname = basepath+'LMA_domain_rpc_thickness_vs_LMA_flashes.png'
    plot_scatter_comp(
        rm_hmax_common-rm_hmin_common, nflashes_LMA_common, [fname],
        labelx='RPC length [m]', labely='N flashes',
        titl='RPC length vs LMA flashes', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'LMA_domain_rpc_thickness_vs_LMA_flashes_density.png'
    plot_scatter_comp(
        rm_hmax_common-rm_hmin_common, flash_dens_LMA_common, [fname],
        labelx='RPC length [m]', labely='Flash density [flashes/km^2]',
        titl='RPC length vs LMA flashes density', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def plot_rpc_time_series(basepath, fname):
    # read rimed particles file
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     rm_hmin, rm_hmax) = read_trt_cell_lightning(basepath+fname)

    traj_ID_uniques = np.unique(traj_ID, return_index=False)

    for traj_ID_unique in traj_ID_uniques:
        inds = np.where(traj_ID == traj_ID_unique)[0]
        td = time_cell[inds]
        rpc_hmin = rm_hmin[inds]
        rpc_hmax = rm_hmax[inds]

        figfname = basepath+str(traj_ID_unique)+'_rpc_base_altitude.png'
        plot_timeseries(
            td, [rpc_hmin], [figfname], labelx='Time UTC',
            labely='Altitude [m MSL]', title=str(traj_ID_unique)+' RPC base altitude')

        figfname = basepath+str(traj_ID_unique)+'_rpc_height.png'
        plot_timeseries(
            td, [rpc_hmax-rpc_hmin], [figfname], labelx='Time UTC',
            labely='Height [m]', title=str(traj_ID_unique)+' RPC height')



def plot_rimed_column_EUCLID(basepath, fname_rimed, fname_EUCLID, roi):
    # read rimed particles file
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     rm_hmin, rm_hmax) = read_trt_cell_lightning(basepath+fname_rimed)

    print('Cell steps: ', rank_cell.size)

    # keep only time steps within EUCLID domain
    inds, is_roi = belongs_roi_indices(lat_cell, lon_cell, roi)

    traj_ID = traj_ID[inds]
    time_cell = time_cell[inds]
    rank_cell = rank_cell[inds]
    rm_hmin = rm_hmin[inds]
    rm_hmax = rm_hmax[inds]

    # put non-valid min and max to 0
    rm_hmin = rm_hmin.filled(0.)
    rm_hmax = rm_hmax.filled(0.)


#    # keep only valid time steps
#    valid = np.logical_not(np.logical_or(
#        np.ma.getmaskarray(rm_hmin), np.ma.getmaskarray(rm_hmax)))
#
#    traj_ID = traj_ID[valid]
#    time_cell = time_cell[valid]
#    rank_cell = rank_cell[valid]
#    rm_hmax = rm_hmax[valid]
#    rm_hmin = rm_hmin[valid]

    print('Cell steps within LMA domain: ', rank_cell.size)

    # read EUCLID particles column file
    (traj_ID_EUCLID, time_cell_EUCLID, lon_cell_EUCLID, lat_cell_EUCLID, area_cell_EUCLID, rank_cell_EUCLID,
     nflashes_cell_EUCLID, flash_dens_cell_EUCLID) = read_trt_cell_lightning(
        basepath+fname_EUCLID)

    print('Cell steps: ', rank_cell_EUCLID.size)

    # keep only time steps within EUCLID domain
    inds, is_roi = belongs_roi_indices(lat_cell_EUCLID, lon_cell_EUCLID, roi)

    traj_ID_EUCLID = traj_ID_EUCLID[inds]
    time_cell_EUCLID = time_cell_EUCLID[inds]
    rank_cell_EUCLID = rank_cell_EUCLID[inds]
    nflashes_cell_EUCLID = nflashes_cell_EUCLID[inds]
    flash_dens_cell_EUCLID = flash_dens_cell_EUCLID[inds]

    print('Cell steps within LMA domain: ', rank_cell_EUCLID.size)

    # Get data from common TRT cells and time steps
    rm_hmin_common = np.asarray([])
    rm_hmax_common = np.asarray([])

    flash_dens_EUCLID_common = np.ma.asarray([])
    nflashes_EUCLID_common = np.ma.asarray([])

    for i, time_cell_EUCLID_el in enumerate(time_cell_EUCLID):
        ind = np.ma.where(np.logical_and(
            time_cell == time_cell_EUCLID_el,
            traj_ID == traj_ID_EUCLID[i]))[0]
        if ind.size == 0:
            continue
        rm_hmin_common = np.append(rm_hmin_common, rm_hmin[ind])
        rm_hmax_common = np.append(rm_hmax_common, rm_hmax[ind])

        flash_dens_EUCLID_common = np.append(flash_dens_EUCLID_common, flash_dens_cell_EUCLID[i])
        nflashes_EUCLID_common = np.append(nflashes_EUCLID_common, nflashes_cell_EUCLID[i])


    print('Common cell steps: ', rm_hmin_common.size)

    # plot data
    fname = basepath+'LMA_domain_rpc_min_height_vs_EUCLID_flashes.png'
    plot_scatter_comp(
        rm_hmin_common, nflashes_EUCLID_common, [fname],
        labelx='RPC base height [m MSL]', labely='N strokes',
        titl='RPC base height vs EUCLID strokes', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'LMA_domain_rpc_min_height_vs_EUCLID_flashes_density.png'
    plot_scatter_comp(
        rm_hmin_common, flash_dens_EUCLID_common, [fname],
        labelx='RPC base height [m MSL]', labely='Stroke density [strokes/km^2]',
        titl='RPC base height vs EUCLID stroke density', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    # plot data
    fname = basepath+'LMA_domain_rpc_thickness_vs_EUCLID_flashes.png'
    plot_scatter_comp(
        rm_hmax_common-rm_hmin_common, nflashes_EUCLID_common, [fname],
        labelx='RPC length [m]', labely='N strokes',
        titl='RPC length vs EUCLID strokes', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'LMA_domain_rpc_thickness_vs_EUCLID_flashes_density.png'
    plot_scatter_comp(
        rm_hmax_common-rm_hmin_common, flash_dens_EUCLID_common, [fname],
        labelx='RPC length [m]', labely='Strokes density [strokes/km^2]',
        titl='RPC length vs EUCLID strokes density', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def plot_rimed_column_EUCLID_pn(basepath, fname_rimed, fname_EUCLID, roi):
    # read rimed particles file
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     rm_hmin, rm_hmax) = read_trt_cell_lightning(basepath+fname_rimed)

    print('Cell steps: ', rank_cell.size)

    # keep only time steps within EUCLID domain
    inds, is_roi = belongs_roi_indices(lat_cell, lon_cell, roi)

    traj_ID = traj_ID[inds]
    time_cell = time_cell[inds]
    rank_cell = rank_cell[inds]
    rm_hmin = rm_hmin[inds]
    rm_hmax = rm_hmax[inds]

    # put non-valid min and max to 0
    rm_hmin = rm_hmin.filled(0.)
    rm_hmax = rm_hmax.filled(0.)


#    # keep only valid time steps
#    valid = np.logical_not(np.logical_or(
#        np.ma.getmaskarray(rm_hmin), np.ma.getmaskarray(rm_hmax)))
#
#    traj_ID = traj_ID[valid]
#    time_cell = time_cell[valid]
#    rank_cell = rank_cell[valid]
#    rm_hmax = rm_hmax[valid]
#    rm_hmin = rm_hmin[valid]

    print('Cell steps within LMA domain: ', rank_cell.size)

    # read EUCLID particles column file
    (traj_ID_EUCLID, time_cell_EUCLID, lon_cell_EUCLID, lat_cell_EUCLID,
     area_cell_EUCLID, rank_cell_EUCLID, nflashes_p_cell_EUCLID,
     nflashes_n_cell_EUCLID) = read_trt_cell_lightning(basepath+fname_EUCLID)

    print('Cell steps: ', rank_cell_EUCLID.size)

    # keep only time steps within EUCLID domain
    inds, is_roi = belongs_roi_indices(lat_cell_EUCLID, lon_cell_EUCLID, roi)

    traj_ID_EUCLID = traj_ID_EUCLID[inds]
    time_cell_EUCLID = time_cell_EUCLID[inds]
    rank_cell_EUCLID = rank_cell_EUCLID[inds]
    nflashes_p_cell_EUCLID = nflashes_p_cell_EUCLID[inds]
    nflashes_n_cell_EUCLID = nflashes_n_cell_EUCLID[inds]
    area_cell_EUCLID = area_cell_EUCLID[inds]

    print('Cell steps within LMA domain: ', rank_cell_EUCLID.size)

    # Get data from common TRT cells and time steps
    rm_hmin_common = np.asarray([])
    rm_hmax_common = np.asarray([])

    nflashes_p_cell_EUCLID_common = np.ma.asarray([])
    nflashes_n_cell_EUCLID_common = np.ma.asarray([])
    area_cell_EUCLID_common = np.ma.asarray([])

    for i, time_cell_EUCLID_el in enumerate(time_cell_EUCLID):
        ind = np.ma.where(np.logical_and(
            time_cell == time_cell_EUCLID_el,
            traj_ID == traj_ID_EUCLID[i]))[0]
        if ind.size == 0:
            continue
        rm_hmin_common = np.append(rm_hmin_common, rm_hmin[ind])
        rm_hmax_common = np.append(rm_hmax_common, rm_hmax[ind])

        nflashes_p_cell_EUCLID_common = np.append(
            nflashes_p_cell_EUCLID_common, nflashes_p_cell_EUCLID[i])
        nflashes_n_cell_EUCLID_common = np.append(
            nflashes_n_cell_EUCLID_common, nflashes_n_cell_EUCLID[i])
        area_cell_EUCLID_common = np.append(
            area_cell_EUCLID_common, area_cell_EUCLID[i])


    print('Common cell steps: ', rm_hmin_common.size)

    # plot data
    fname = basepath+'LMA_domain_rpc_min_height_vs_EUCLID_flashes.png'

    fig, ax = plot_scatter_comp(
        rm_hmin_common, nflashes_p_cell_EUCLID_common+nflashes_n_cell_EUCLID_common, [fname],
        labelx='RPC base height [m MSL]', labely='N strokes',
        titl='RPC base height vs EUCLID strokes', axis=None, metadata=None,
        dpi=72, save_fig=False, point_format='go')

    fig, ax = plot_scatter_comp(
        rm_hmin_common, nflashes_p_cell_EUCLID_common, [fname],
        labelx='RPC base height [m MSL]', labely='N strokes',
        titl='RPC base height vs EUCLID strokes', axis=None, metadata=None,
        dpi=72, ax=ax, fig=fig, save_fig=False, point_format='bx')

    plot_scatter_comp(
        rm_hmin_common, nflashes_n_cell_EUCLID_common, [fname],
        labelx='RPC base height [m MSL]', labely='N strokes',
        titl='RPC base height vs EUCLID strokes', axis=None, metadata=None,
        dpi=72, ax=ax, fig=fig, save_fig=True, point_format='r+')

    print("----- plot to '%s'" % fname)

    fname = basepath+'LMA_domain_rpc_min_height_vs_EUCLID_flashes_density.png'
    fig, ax = plot_scatter_comp(
        rm_hmin_common, (nflashes_p_cell_EUCLID_common+nflashes_n_cell_EUCLID_common)/area_cell_EUCLID_common,
        [fname], labelx='RPC base height [m MSL]',
        labely='Stroke density [strokes/km^2]',
        titl='RPC base height vs EUCLID stroke density', axis=None,
        metadata=None, dpi=72, save_fig=False, point_format='go')

    fig, ax = plot_scatter_comp(
        rm_hmin_common, nflashes_p_cell_EUCLID_common/area_cell_EUCLID_common,
        [fname], labelx='RPC base height [m MSL]',
        labely='Stroke density [strokes/km^2]',
        titl='RPC base height vs EUCLID stroke density', axis=None,
        metadata=None, dpi=72, ax=ax, fig=fig, save_fig=False, point_format='bx')

    plot_scatter_comp(
        rm_hmin_common, nflashes_n_cell_EUCLID_common/area_cell_EUCLID_common,
        [fname], labelx='RPC base height [m MSL]',
        labely='Stroke density [strokes/km^2]',
        titl='RPC base height vs EUCLID stroke density', axis=None,
        metadata=None, dpi=72, ax=ax, fig=fig, save_fig=True, point_format='r+')

    print("----- plot to '%s'" % fname)

    # plot data
    fname = basepath+'LMA_domain_rpc_thickness_vs_EUCLID_flashes.png'
    fig, ax = plot_scatter_comp(
        rm_hmax_common-rm_hmin_common,
        nflashes_p_cell_EUCLID_common+nflashes_n_cell_EUCLID_common, [fname],
        labelx='RPC length [m]', labely='N strokes',
        titl='RPC length vs EUCLID strokes', axis=None, metadata=None, dpi=72,
        save_fig=False, point_format='go')

    fig, ax = plot_scatter_comp(
        rm_hmax_common-rm_hmin_common, nflashes_p_cell_EUCLID_common, [fname],
        labelx='RPC length [m]', labely='N strokes',
        titl='RPC length vs EUCLID strokes', axis=None, metadata=None, dpi=72,
        ax=ax, fig=fig, save_fig=False, point_format='bx')

    plot_scatter_comp(
        rm_hmax_common-rm_hmin_common, nflashes_n_cell_EUCLID_common, [fname],
        labelx='RPC length [m]', labely='N strokes',
        titl='RPC length vs EUCLID strokes', axis=None, metadata=None, dpi=72,
        ax=ax, fig=fig, save_fig=True, point_format='r+')

    print("----- plot to '%s'" % fname)

    fname = basepath+'LMA_domain_rpc_thickness_vs_EUCLID_flashes_density.png'
    fig, ax = plot_scatter_comp(
        rm_hmax_common-rm_hmin_common,
        (nflashes_p_cell_EUCLID_common+nflashes_n_cell_EUCLID_common)/area_cell_EUCLID_common, [fname],
        labelx='RPC length [m]', labely='Strokes density [strokes/km^2]',
        titl='RPC length vs EUCLID strokes density', axis=None, metadata=None,
        dpi=72, save_fig=False, point_format='go')

    fig, ax = plot_scatter_comp(
        rm_hmax_common-rm_hmin_common,
        nflashes_p_cell_EUCLID_common/area_cell_EUCLID_common, [fname],
        labelx='RPC length [m]', labely='Strokes density [strokes/km^2]',
        titl='RPC length vs EUCLID strokes density', axis=None, metadata=None,
        dpi=72, ax=ax, fig=fig, save_fig=False, point_format='bx')

    plot_scatter_comp(
        rm_hmax_common-rm_hmin_common,
        nflashes_n_cell_EUCLID_common/area_cell_EUCLID_common, [fname],
        labelx='RPC length [m]', labely='Strokes density [strokes/km^2]',
        titl='RPC length vs EUCLID strokes density', axis=None, metadata=None,
        dpi=72, ax=ax, fig=fig, save_fig=True, point_format='r+')

    print("----- plot to '%s'" % fname)


def plot_rank_LMA_flashes(basepath, fname, roi=None):
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(basepath+fname)

    print('Cell steps: ', rank_cell.size)

    # plot data
    fname = basepath+'rank-LMA_flashes_density'
    plot_scatter_comp(
        rank_cell/10., flash_dens_cell, [fname],
        labelx='rank', labely='flash density [flashes/km2]',
        titl='Cell rank vs LMA flash density', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'rank-LMA_flashes'
    plot_scatter_comp(
        rank_cell/10., nflashes_cell, [fname],
        labelx='rank', labely='Flashes',
        titl='Cell rank vs LMA flashes', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    if roi is None:
        return

    # keep only time steps within LMA domain
    inds, is_roi = belongs_roi_indices(lat_cell, lon_cell, roi)

    nflashes_cell = nflashes_cell[inds]
    flash_dens_cell = flash_dens_cell[inds]
    rank_cell = rank_cell[inds]

    print('Cell steps: ', rank_cell.size)

    # plot data
    fname = basepath+'LMA_domain_rank-LMA_flashes_density'
    plot_scatter_comp(
        rank_cell/10., flash_dens_cell, [fname],
        labelx='rank', labely='flash density [flashes/km2]',
        titl='Cell rank vs LMA flash density\nReduced LMA domain only', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'LMA_domain_rank-LMA_flashes'
    plot_scatter_comp(
        rank_cell/10., nflashes_cell, [fname],
        labelx='rank', labely='Flashes',
        titl='Cell rank vs LMA flashes\nReduced LMA domain only', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def plot_rank_LMA_sources(basepath, fname, roi=None):
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(basepath+fname)

    print('Cell steps: ', rank_cell.size)

    # plot data
    fname = basepath+'rank-LMA_sources_density'
    plot_scatter_comp(
        rank_cell/10., flash_dens_cell, [fname],
        labelx='rank', labely='sources density [sources/km2]',
        titl='Cell rank vs LMA sources density', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'rank-LMA_sources'
    plot_scatter_comp(
        rank_cell/10., nflashes_cell, [fname],
        labelx='rank', labely='Sources',
        titl='Cell rank vs LMA sources', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    if roi is None:
        return

    # keep only time steps within LMA domain
    inds, is_roi = belongs_roi_indices(lat_cell, lon_cell, roi)

    nflashes_cell = nflashes_cell[inds]
    flash_dens_cell = flash_dens_cell[inds]
    rank_cell = rank_cell[inds]

    print('Cell steps: ', rank_cell.size)

    # plot data
    fname = basepath+'LMA_domain_rank-LMA_sources_density'
    plot_scatter_comp(
        rank_cell/10., flash_dens_cell, [fname],
        labelx='rank', labely='sources density [sources/km2]',
        titl='Cell rank vs LMA sources density\nReduced LMA domain only', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'LMA_domain_rank-LMA_sources'
    plot_scatter_comp(
        rank_cell/10., nflashes_cell, [fname],
        labelx='rank', labely='sources',
        titl='Cell rank vs LMA sources\nReduced LMA domain only', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def read_general_scores(basepath, fname):
    (traj_ID, time_flash_density_max, flash_density_max_rank,
     flash_density_max_nflashes, flash_density_max_area,
     flash_density_max, time_rank_max, rank_max) = read_trt_scores(basepath+fname)

    print('Total number of cells: '+str(traj_ID.size))
    print('weak: '+str(rank_max[rank_max<12].size))
    print('developing: '+str(rank_max[np.logical_and(rank_max >=12, rank_max <15)].size))
    print('moderate: '+str(rank_max[np.logical_and(rank_max >=15, rank_max <25)].size))
    print('severe: '+str(rank_max[np.logical_and(rank_max >=25, rank_max <35)].size))
    print('very severe: '+str(rank_max[rank_max>=35].size))

    print('cells without EUCLID flashes '+str(rank_max[flash_density_max_nflashes == 0].size))
    print('Max rank without EUCLID flashes '+str(np.max(rank_max[flash_density_max_nflashes == 0])))

    # plot data
    fname = basepath+'max_Euclid_flashes_density-rank'
    plot_scatter_comp(
        flash_density_max, flash_density_max_rank/10., [fname],
        labelx='max strokes density [strokes/km^2]', labely='rank',
        titl='Rank when max Euclid CG strokes density', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'max_Euclid_flashes-rank'
    plot_scatter_comp(
        flash_density_max_nflashes, flash_density_max_rank/10., [fname],
        labelx='max N strokes', labely='rank',
        titl='Rank when max Euclid CG strokes', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def plot_time_diff_max_EUCLID_density_max_rank(basepath, fname):
    (traj_ID, time_flash_density_max, flash_density_max_rank,
     flash_density_max_nflashes, flash_density_max_area,
     flash_density_max, time_rank_max, rank_max) = read_trt_scores(basepath+fname)

    #filter to data with flashes
    time_rank_max = time_rank_max[flash_density_max > 0.]
    time_flash_density_max = time_flash_density_max[flash_density_max > 0.]
    flash_density_max = flash_density_max[flash_density_max > 0.]

    tdiff = np.empty(flash_density_max.size)
    for i, tmax_rank in enumerate(time_rank_max):
        tdiff[i] = (tmax_rank-time_flash_density_max[i]).total_seconds()

    print('Max rank after max stroke density', np.size(tdiff[tdiff > 0.]))
    print('Max rank together max stroke density', np.size(tdiff[tdiff == 0.]))
    print('Max rank before max stroke density', np.size(tdiff[tdiff < 0.]))

    fname = basepath+'max_EUCLID_flashes_density-time_diff_rank_max.png'
    plot_scatter_comp(
        flash_density_max, tdiff, [fname],
        labelx='max stroke density [strokes/km^2]', labely='Time diff [s]',
        titl='Max EUCLID CG stroke density vs\nTime max rank - Time max stroke density',
        axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def plot_time_diff_max_EUCLID_strokes_max_rank(basepath, fname):
    # read EUCLID flashes file
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(basepath+fname)

    traj_ID_uniques = np.unique(traj_ID, return_index=False)

    time_rank_max = np.array([])
    time_nflashes_max = np.array([])
    max_nflashes = np.array([])

    tdiff = np.array([])
    for traj_ID_unique in traj_ID_uniques:
        inds = np.where(traj_ID == traj_ID_unique)[0]

        time_aux = time_cell[inds]
        nflashes_aux = nflashes_cell[inds]
        rank_aux = rank_cell[inds]

        if nflashes_aux.max() == 0:
            continue

        time_rank_max_aux = time_aux[np.argmax(rank_aux)]
        time_nflashes_max_aux = time_aux[np.argmax(nflashes_aux)]

        time_rank_max = np.append(time_rank_max, time_rank_max_aux)
        time_nflashes_max = np.append(time_nflashes_max, time_nflashes_max_aux)
        max_nflashes = np.append(max_nflashes, nflashes_aux.max())
        tdiff = np.append(tdiff, (time_rank_max_aux-time_nflashes_max_aux).total_seconds())

    print('Max rank after max stroke number', np.size(tdiff[tdiff > 0.]))
    print('Max rank together max stroke number', np.size(tdiff[tdiff == 0.]))
    print('Max rank before max stroke number', np.size(tdiff[tdiff < 0.]))

    fname = basepath+'max_EUCLID_flashes-time_diff_rank_max.png'
    plot_scatter_comp(
        max_nflashes, tdiff, [fname],
        labelx='max number of strokes', labely='Time diff [s]',
        titl='Max EUCLID CG stroke number vs\nTime max rank - Time max stroke number',
        axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def plot_time_diff_max_LMA_max_rank(basepath, fname, roi=None):
    # read LMA flashes file
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(basepath+fname)

    traj_ID_uniques = np.unique(traj_ID, return_index=False)

    time_rank_max = np.array([])
    time_nflashes_max = np.array([])
    time_density_max = np.array([])

    max_nflashes = np.array([])
    max_density = np.array([])

    tdiff_flashes = np.array([])
    tdiff_density = np.array([])
    for traj_ID_unique in traj_ID_uniques:
        inds = np.where(traj_ID == traj_ID_unique)[0]

        time_aux = time_cell[inds]
        nflashes_aux = nflashes_cell[inds]
        density_aux = flash_dens_cell[inds]
        rank_aux = rank_cell[inds]

        if nflashes_aux.max() == 0:
            continue

        if roi is not None:
            lat_cell_aux = lat_cell[inds]
            lon_cell_aux = lon_cell[inds]
            _, is_roi = belongs_roi_indices(lat_cell_aux, lon_cell_aux, roi)

            if is_roi in ('None', 'Some'):
                continue

        time_rank_max_aux = time_aux[np.argmax(rank_aux)]
        time_nflashes_max_aux = time_aux[np.argmax(nflashes_aux)]
        time_density_max_aux = time_aux[np.argmax(density_aux)]

        time_rank_max = np.append(time_rank_max, time_rank_max_aux)
        time_nflashes_max = np.append(time_nflashes_max, time_nflashes_max_aux)
        time_density_max = np.append(time_density_max, time_density_max_aux)

        max_nflashes = np.append(max_nflashes, nflashes_aux.max())
        max_density = np.append(max_density, density_aux.max())

        tdiff_flashes = np.append(tdiff_flashes, (time_rank_max_aux-time_nflashes_max_aux).total_seconds())
        tdiff_density = np.append(tdiff_density, (time_rank_max_aux-time_density_max_aux).total_seconds())

    print('Max rank after max stroke number', np.size(tdiff_flashes[tdiff_flashes > 0.]))
    print('Max rank together max stroke number', np.size(tdiff_flashes[tdiff_flashes == 0.]))
    print('Max rank before max stroke number', np.size(tdiff_flashes[tdiff_flashes < 0.]))

    print('Max rank after max stroke density', np.size(tdiff_density[tdiff_density > 0.]))
    print('Max rank together max stroke density', np.size(tdiff_density[tdiff_density == 0.]))
    print('Max rank before max stroke density', np.size(tdiff_density[tdiff_density < 0.]))

    if roi is not None:
        fname = basepath+'LMA_domain_max_LMA_flashes-time_diff_rank_max.png'
    else:
        fname = basepath+'max_LMA_flashes-time_diff_rank_max.png'
    plot_scatter_comp(
        max_nflashes, tdiff_flashes, [fname],
        labelx='max number of flashes', labely='Time diff [s]',
        titl='Max LMA flashes vs\nTime max rank - Time max flashes',
        axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    if roi is not None:
        fname = basepath+'LMA_domain_max_LMA_flashes_density-time_diff_rank_max.png'
    else:
        fname = basepath+'max_LMA_flashes_density-time_diff_rank_max.png'
    plot_scatter_comp(
        max_density, tdiff_density, [fname],
        labelx='max flash density', labely='Time diff [s]',
        titl='Max LMA flashes vs\nTime max rank - Time max flash density',
        axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def plot_max_LMA_rank(basepath, fname):
    # read LMA flashes file
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(basepath+fname)

    traj_ID_uniques = np.unique(traj_ID, return_index=False)

    max_nflashes = np.array([])
    max_density = np.array([])
    max_nflashes_rank = np.array([])
    max_density_rank = np.array([])

    tdiff_flashes = np.array([])
    tdiff_density = np.array([])
    for traj_ID_unique in traj_ID_uniques:
        inds = np.where(traj_ID == traj_ID_unique)[0]

        nflashes_aux = nflashes_cell[inds]
        density_aux = flash_dens_cell[inds]
        rank_aux = rank_cell[inds]

        max_nflashes = np.append(max_nflashes, nflashes_aux.max())
        max_nflashes_rank = np.append(max_nflashes_rank, rank_aux[np.argmax(nflashes_aux)])

        max_density = np.append(max_density, density_aux.max())
        max_density_rank = np.append(max_density_rank, rank_aux[np.argmax(density_aux)])

    fname = basepath+'max_LMA_flashes-rank.png'
    plot_scatter_comp(
        max_nflashes, max_nflashes_rank/10., [fname],
        labelx='max number of flashes', labely='rank',
        titl='Rank when max LMA flashes reached',
        axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)

    fname = basepath+'max_LMA_flashes_density-rank.png'
    plot_scatter_comp(
        max_density, max_density_rank/10., [fname],
        labelx='max flash density', labely='rank',
        titl='Rank when max LMA flashes density reached',
        axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def plot_LMA_EUCLID(basepath, fname_LMA, fname_EUCLID, roi):
    # read EUCLID flashes file
    (traj_ID, time_cell, lon_cell, lat_cell, area_cell, rank_cell,
     nflashes_cell, flash_dens_cell) = read_trt_cell_lightning(basepath+fname_EUCLID)

    print('Cell steps: ', rank_cell.size)

    # keep only time steps within LMA domain
    inds, is_roi = belongs_roi_indices(lat_cell, lon_cell, roi)

    traj_ID = traj_ID[inds]
    time_cell = time_cell[inds]
    rank_cell = rank_cell[inds]
    nflashes_cell = nflashes_cell[inds]
    flash_dens_cell = flash_dens_cell[inds]

    print('Cell steps within LMA domain: ', rank_cell.size)


    # read LMA particles column file
    (traj_ID_LMA, time_cell_LMA, lon_cell_LMA, lat_cell_LMA, area_cell_LMA, rank_cell_LMA,
     nflashes_cell_LMA, flash_dens_cell_LMA) = read_trt_cell_lightning(basepath+fname_LMA)

    print('Cell steps: ', rank_cell_LMA.size)

    # keep only time steps within LMA domain
    inds, is_roi = belongs_roi_indices(lat_cell_LMA, lon_cell_LMA, roi)

    traj_ID_LMA = traj_ID_LMA[inds]
    time_cell_LMA = time_cell_LMA[inds]
    rank_cell_LMA = rank_cell_LMA[inds]
    nflashes_cell_LMA = nflashes_cell_LMA[inds]
    flash_dens_cell_LMA = flash_dens_cell_LMA[inds]

    print('Cell steps within LMA domain: ', rank_cell_LMA.size)

    # Get data from common TRT cells and time steps
    flash_dens_common = np.asarray([])
    nflashes_cell_common = np.asarray([])

    flash_dens_LMA_common = np.ma.asarray([])
    nflashes_LMA_common = np.ma.asarray([])

    for i, time_cell_LMA_el in enumerate(time_cell_LMA):
        ind = np.ma.where(np.logical_and(
            time_cell == time_cell_LMA_el,
            traj_ID == traj_ID_LMA[i]))[0]
        if ind.size == 0:
            continue
        flash_dens_common = np.append(flash_dens_common, flash_dens_cell[ind])
        nflashes_cell_common = np.append(nflashes_cell_common, nflashes_cell[ind])

        flash_dens_LMA_common = np.append(flash_dens_LMA_common, flash_dens_cell_LMA[i])
        nflashes_LMA_common = np.append(nflashes_LMA_common, nflashes_cell_LMA[i])


    print('Common cell steps: ', flash_dens_common.size)

    # plot data
    fname = basepath+'LMA_domain_LMA_flashes-Euclid_flashes'
    plot_scatter_comp(
        nflashes_LMA_common, nflashes_cell_common, [fname],
        labelx='LMA flashes', labely='Euclid CG strokes',
        titl='LMA flashes vs Euclid CG strokes\nReduced LMA domain only', axis=None, metadata=None, dpi=72)

    print("----- plot to '%s'" % fname)


def _print_end_msg(text):
    """
    prints end message

    Parameters
    ----------
    text : str
        the text to be printed

    Returns
    -------
    Nothing

    """
    print(text + datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
