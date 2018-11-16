#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
lightning_traj_postproc
================================================

This program post-processes lightning trajectories

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import glob
from warnings import warn
import os
from copy import deepcopy

import numpy as np

from pyrad.io import read_histogram
from pyrad.graph import plot_histogram2

print(__doc__)


def main():
    """
    """
    # basepath = '/data/pyrad_products/rad4alp_hydro_PHA/'
    basepath = '/store/msrad/radar/pyrad_products/rad4alp_hydro_PHA/data_analysis_min10sources/'

    print("====== Lightning post-processing started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== Lightning post-processing finished: ")

#    prefix = ['All', 'no_CG', 'CGt']
#    dir = ['all_data/', 'no_CG/', 'CGt/']

    prefix = ['CGt', 'CGn', 'CGp']
    dir = ['CGt/', 'CGn/', 'CGp/']

    sources = ['allsources', 'firstsource']
    titles = ['Value at VHF source location', 'Value at flash origin location']
    datatypes = [
        'dBZc', 'entropy', 'hydro', 'hydro_prop', 'KDPc', 'nhydro', 'RhoHVc',
        'TEMP', 'ZDRc']
    labels = [
        'horizontal reflectivity [dBZ]', 'entropy [-]',
        'radar echo classification [-]', 'proportion of hydrometeors [%]',
        'specific differential phase [deg/km]',
        'Number of hydrometeor in radar gate',
        'copolar correlation coefficient [-]', 'temperature [deg Celsius]',
        'differential reflectivity [dB]']

#    sources = ['', '_first_source']
#    titles = ['Value at VHF source location', 'Value at flash origin location']
#    datatypes = ['alt', 'dBm']
#    labels = [
#        'VHF source altitude [m MSL]', 'VHF source power [dBm]']

    for source, titl in zip(sources, titles):
        for datatype, labelx in zip(datatypes, labels):
            hist1, bin_edges1 = read_histogram(
                basepath+dir[0]+prefix[0]+'_'+source+'_ts_trajlightning_'+datatype+'.csv')
            hist2, bin_edges2 = read_histogram(
                basepath+dir[1]+prefix[1]+'_'+source+'_ts_trajlightning_'+datatype+'.csv')
            hist3, bin_edges3 = read_histogram(
                basepath+dir[2]+prefix[2]+'_'+source+'_ts_trajlightning_'+datatype+'.csv')

#            hist1, bin_edges1 = read_histogram(
#                basepath+dir[0]+prefix[0]+'_Santis_hist_'+datatype+source+'.csv')
#            hist2, bin_edges2 = read_histogram(
#                basepath+dir[1]+prefix[1]+'_Santis_hist_'+datatype+source+'.csv')
#            hist3, bin_edges3 = read_histogram(
#                basepath+dir[2]+prefix[2]+'_Santis_hist_'+datatype+source+'.csv')

            if (not np.array_equal(bin_edges1, bin_edges2) or
                    not np.array_equal(bin_edges1, bin_edges3)):
                warn('Bin edges should be identical to group histograms')
                continue

            if hist1 is None or hist2 is None or hist3 is None:
                warn('Dataset not available')
                continue

            invert_xaxis=False
            if datatype == 'TEMP':
                invert_xaxis=True

            bin_res = bin_edges1[1]-bin_edges1[0]
            bin_centers = bin_edges1[1:]-bin_res/2.
            fname = basepath+'group_Santis_CG_hist_'+source+'_'+datatype+'.png'

            fig, ax = plot_histogram2(
                bin_centers, hist1, [fname], labelx=labelx, titl=titl, alpha=0.25,
                save_fig=False, color='b', invert_xaxis=invert_xaxis)

            fig, ax = plot_histogram2(
                bin_centers, hist2, [fname], labelx=labelx, titl=titl,
                ax=ax, fig=fig, save_fig=False, color='g', alpha=0.25,
                invert_xaxis=invert_xaxis)

            fname_list = plot_histogram2(
                bin_centers, hist3, [fname], labelx=labelx, titl=titl,
                ax=ax, fig=fig, save_fig=True, color='r', alpha=0.25,
                invert_xaxis=invert_xaxis)

            # Total number of values
            n1 = np.ma.sum(hist1)
            n2 = np.ma.sum(hist2)
            n3 = np.ma.sum(hist3)

            print(n1)
            print(n2)
            print(n3)

            # Mode
            print(bin_centers[np.ma.argmax(hist1)])
            print(bin_centers[np.ma.argmax(hist2)])
            print(bin_centers[np.ma.argmax(hist3)])

            # Median
            freq1 = np.ma.cumsum(hist1)/n1
            freq2 = np.ma.cumsum(hist2)/n2
            freq3 = np.ma.cumsum(hist3)/n3

            ind1 = np.where(freq1 >= 0.5)[0][0]
            ind2 = np.where(freq2 >= 0.5)[0][0]
            ind3 = np.where(freq3 >= 0.5)[0][0]

            print(bin_centers[ind1])
            print(bin_centers[ind2])
            print(bin_centers[ind3])

            if datatype == 'hydro':
                print(hist1/n1*100.)
                print(hist2/n2*100.)
                print(hist3/n3*100.)


            print('plotted '+''.join(fname_list))



#    titles = ['Flash duration', 'Flash occurrence time', 'Flash area']
#    datatypes = ['duration', 'time', 'area']
#    labels = ['duration [ms]', 'time [h UTC]', 'area [km2]']
#
#    for datatype, titl, labelx in zip(datatypes, titles, labels):
#
#        hist1, bin_edges1 = read_histogram(
#            basepath+dir[0]+prefix[0]+'_Santis_hist_'+datatype+'.csv')
#        hist2, bin_edges2 = read_histogram(
#            basepath+dir[1]+prefix[1]+'_Santis_hist_'+datatype+'.csv')
#        hist3, bin_edges3 = read_histogram(
#            basepath+dir[2]+prefix[2]+'_Santis_hist_'+datatype+'.csv')
#
#        if (not np.array_equal(bin_edges1, bin_edges2) or
#                not np.array_equal(bin_edges1, bin_edges3)):
#            warn('Bin edges should be identical to group histograms')
#            continue
#
#        if hist1 is None or hist2 is None or hist3 is None:
#            warn('Dataset not available')
#            continue
#
#        invert_xaxis=False
#        if datatype == 'TEMP':
#            invert_xaxis=True
#
#        bin_res = bin_edges1[1]-bin_edges1[0]
#        bin_centers = bin_edges1[1:]-bin_res/2.
#        fname = basepath+'group_Santis_CG_hist_'+datatype+'.png'
#
#        fig, ax = plot_histogram2(
#            bin_centers, hist1, [fname], labelx=labelx, titl=titl, alpha=0.25,
#            save_fig=False, color='b', invert_xaxis=invert_xaxis)
#
#        fig, ax = plot_histogram2(
#            bin_centers, hist2, [fname], labelx=labelx, titl=titl,
#            ax=ax, fig=fig, save_fig=False, color='g', alpha=0.25,
#            invert_xaxis=invert_xaxis)
#
#        fname_list = plot_histogram2(
#            bin_centers, hist3, [fname], labelx=labelx, titl=titl,
#            ax=ax, fig=fig, save_fig=True, color='r', alpha=0.25,
#            invert_xaxis=invert_xaxis)
#
#        # Total number of values
#        n1 = np.ma.sum(hist1)
#        n2 = np.ma.sum(hist2)
#        n3 = np.ma.sum(hist3)
#
#        print(n1)
#        print(n2)
#        print(n3)
#
#        # Mode
#        print(bin_centers[np.ma.argmax(hist1)])
#        print(bin_centers[np.ma.argmax(hist2)])
#        print(bin_centers[np.ma.argmax(hist3)])
#
#        # Median
#        freq1 = np.ma.cumsum(hist1)/n1
#        freq2 = np.ma.cumsum(hist2)/n2
#        freq3 = np.ma.cumsum(hist3)/n3
#
#        ind1 = np.where(freq1 >= 0.5)[0][0]
#        ind2 = np.where(freq2 >= 0.5)[0][0]
#        ind3 = np.where(freq3 >= 0.5)[0][0]
#
#        print(bin_centers[ind1])
#        print(bin_centers[ind2])
#        print(bin_centers[ind3])
#
#        print('plotted '+''.join(fname_list))




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
