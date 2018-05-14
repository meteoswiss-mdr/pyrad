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

from pyrad.io import read_trt_traj_data
from pyrad.util import belongs_roi_indices
from pyrad.graph import plot_timeseries

print(__doc__)


def main():
    """
    """
    database_path = '/data/TRT/'

    time_dir_list = [
        '2017-06-29', '2017-06-30', '2017-07-10', '2017-07-18']

    roi = {
        'lon': [8.9000010, 9.2000000, 9.4999970, 9.4999970, 8.9000010],
        'lat': [47.0000030, 47.0000030, 47.0000030, 47.5999930, 47.5999930]
    }

    print("====== TRT cell processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== comparison finished: ")

    for i, time_dir in enumerate(time_dir_list):
        data_input_path = database_path+time_dir+'/TRTC_cell/'
        data_output_base = database_path+time_dir+'/TRTC_cell_plots/'

        flist = glob.glob(data_input_path+'*.trt')

        for fname in flist:
            print('Reading TRT trajectory file '+fname)
            (traj_ID, yyyymmddHHMM, lon, lat, ell_L, ell_S, ell_or, area,
             vel_x, vel_y, det, RANKr, CG_n, CG_p, CG, CG_percent_p, ET45,
             ET45m, ET15, ET15m, VIL, maxH, maxHm, POH, RANK, Dvel_x,
             Dvel_y, cell_contour) = read_trt_traj_data(fname)

            inds, is_roi = belongs_roi_indices(lat, lon, roi)

            if is_roi == 'None':
                continue
            elif is_roi == 'Some' and len(lat[inds]) < 3:
                continue

            data_output_path = data_output_base+is_roi+'/'
            if not os.path.isdir(data_output_path):
                os.makedirs(data_output_path)

            # copy file
            copy(fname, data_output_path)

            figfname = data_output_path+str(traj_ID[0])+'_area.png'
            plot_timeseries(
                yyyymmddHHMM, [area], [figfname], labelx='Time UTC',
                labely='area [km2]', title=str(traj_ID[0])+' cell area')

            figfname = data_output_path+str(traj_ID[0])+'_vel.png'
            plot_timeseries(
                yyyymmddHHMM, [vel_x, vel_y], [figfname], labelx='Time UTC',
                labely='Velocity [km2]', labels=['x speed', 'y speed'],
                title=str(traj_ID[0])+' cell velocity')

            figfname = data_output_path+str(traj_ID[0])+'_det.png'
            plot_timeseries(
                yyyymmddHHMM, [det], [figfname], labelx='Time UTC',
                labely='Detection threshold [dBZ]',
                title=str(traj_ID[0])+' cell detection threshold')

            figfname = data_output_path+str(traj_ID[0])+'_rank.png'
            plot_timeseries(
                yyyymmddHHMM, [RANKr], [figfname], labelx='Time UTC',
                labely='Rank [-]', title=str(traj_ID[0])+' cell rank')

            figfname = data_output_path+str(traj_ID[0])+'_lightning.png'
            plot_timeseries(
                yyyymmddHHMM, [CG_n, CG_p, CG], [figfname], labelx='Time UTC',
                labely='N flash [-]', labels=['CG-', 'CG+', 'CG'],
                title=str(traj_ID[0])+' flashes in cell')

            figfname = data_output_path+str(traj_ID[0])+'_ET.png'
            plot_timeseries(
                yyyymmddHHMM, [ET45, ET45m, ET15, ET15m], [figfname],
                labelx='Time UTC', labely='Echo Top [km]',
                labels=['ET45', 'ET45m', 'ET15', 'ET15m'],
                title=str(traj_ID[0])+' Echo top')

            figfname = data_output_path+str(traj_ID[0])+'_VIL.png'
            plot_timeseries(
                yyyymmddHHMM, [VIL], [figfname], labelx='Time UTC',
                labely='VIL [Kg/m2]', labels=['VIL'],
                title=str(traj_ID[0])+' VIL')

            figfname = data_output_path+str(traj_ID[0])+'_maxH.png'
            plot_timeseries(
                yyyymmddHHMM, [maxH, maxHm], [figfname], labelx='Time UTC',
                labely='Max. Echo Height [Km]', labels=['maxH', 'maxHm'],
                title=str(traj_ID[0])+' Height of Max. Reflectivity')

            figfname = data_output_path+str(traj_ID[0])+'_POH.png'
            plot_timeseries(
                yyyymmddHHMM, [POH], [figfname], labelx='Time UTC',
                labely='POH [%]', labels=['POH'],
                title=str(traj_ID[0])+' Probability of Hail')

            # plot position
            figfname = data_output_path+str(traj_ID[0])+'_pos.png'

            td_vec = yyyymmddHHMM-yyyymmddHHMM[0]
            tt_s = []
            for td in td_vec:
                tt_s.append(td.total_seconds())
            tt_s = np.asarray(tt_s)

            marker = 'x'
            col = tt_s
            cmap = 'viridis'
            norm = plt.Normalize(tt_s.min(), tt_s.max())
            cb_label = (
                'Time since '+yyyymmddHHMM[0].strftime('%Y-%m-%d %H:%M') +
                ' [s]')

            fig = plt.figure(figsize=[10, 8], dpi=72)
            ax = fig.add_subplot(111, aspect='equal')

            cax = ax.scatter(
                lon, lat, c=col, marker=marker, alpha=0.5, cmap=cmap,
                norm=norm)
            cb = fig.colorbar(cax, orientation='horizontal')
            cb.set_label(cb_label)

            plt.title(str(traj_ID[0])+' Cell Position')
            plt.xlabel('Lon [Deg]')
            plt.ylabel('Lat [Deg]')

            # Turn on the grid
            ax.grid()

            fig.savefig(figfname, dpi=72)
            plt.close()


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
