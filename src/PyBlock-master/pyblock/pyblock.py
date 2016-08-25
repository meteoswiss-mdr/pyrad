"""
Python Polarimetric Radar Beam Blockage Calculation (PyBlock)
PyBlock v1.3


Author
------
Timothy J. Lang
NASA MSFC
timothy.j.lang@nasa.gov
(256) 961-7861


Description
-----------
    Calculates beam blockage from polarimetric radar data using the specific
differential phase (KDP) and fully self-consistent (FSC) methods of
Lang et al. (2009). The core class is BeamBlockSingleVolume, which obtains
blockage-relevant data from individual radar volumes. This class is invoked
repeatedly by the BeamBlockMultiVolume class to compile data from a large
number of files. Helper classes like BlockStats and RawDataStorage enable
the saving, loading, and plotting of blockage information relevant to both
methods.
    If you have a set of radar files to analyze, the easiest way to start is to
import this module and start ingesting data via the BeamBlockMultiVolume class.
e.g., multiblock = pyblock.BeamBlockMultiVolume(list_of_files, **kwargs). Set
keywords to fit the characteristics of your dataset. See DEFAULT_KW for what
program expects and how you might change those parameters.


Last Updated
------------
v1.3 - 08/07/2015
v1.2 - 07/02/2015
v1.1 - 04/29/2015
v1.0 - 04/08/2015


References
----------
Giangrande, S. E., and A. V. Ryzhkov, 2005: Calibration of Dual-Polarization
    Radar in the Presence of Partial Beam Blockage. J. Atmos. Oceanic Technol.,
    22, 1156–1166. doi: http://dx.doi.org/10.1175/JTECH1766.1
Lang, T. J., S. W. Nesbitt, and L. D. Carey, 2009: On the correction of
    partial beam blockage in polarimetric radar data. J. Atmos. Oceanic
    Technol., 26, 943–957.


Dependencies
------------
numpy, pyart, csu_radartools, dualpol, warnings, os, __future__, matplotlib,
statsmodels, gzip, pickle, six


Change Log
----------
v1.3 Major Changes (08/07/2015):
1. Made code Python 3 compliant.

v1.2 Major Changes (07/02/2015):
1. Made all code pep8 compliant.

v1.1 Major Changes (04/29/2015):
1. Added PrintBlock and _PlotBlock objects to support text printout and simple
   plotting methods for KdpMethodAnalysis and SelfConistentAnalysis.
2. Added capability to manually edit corrections.
3. Fixed bug that was preventing an azimuthally varying
   BeamBlockSingleVolume.range_thresh attribute if the user specified that.

v1.0 Functionality (04/08/2015)
1. Will compute beam blockage for any arbitrary set of polarimetric radar
   volumes via the KDP method of Lang et al. (2009). These results can be
   saved to file, turned into statistics, and plotted.
2. Compares KDP method's Zh and Zdr data inside and outside of blocked regions,
   in order to help with the derivation of corrections that need to be applied.
3. Collates data relevant for calculating blockage via the FSC method. Fits
   multiple-linear regression to determine self-consistency relationship and
   does rainfall integration comparisons to estimate beam blockage.

"""
from __future__ import division
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from warnings import warn
import statsmodels.api as sm
import os
import gzip
import pickle
import pyart
import dualpol
from csu_radartools import csu_misc
import six

VERSION = '1.3'
DATA_DIR = os.sep.join([os.path.dirname(__file__), 'data'])+'/'
DEFAULT_SND = DATA_DIR + 'default_sounding.txt'
RANGE_MULT = 1000.0  # m per km
DEFAULT_RANGE = [20, 90]  # km
DEFAULT_RAIN = [1.5, 2]  # KDP in deg/km
DEFAULT_DRIZZ = [-0.1, 0.1]  # KDP in deg/km
BAD = -32768
DEFAULT_SDP = 12
DEFAULT_BINS = 1.0
DEFAULT_ZH_UNCERTAINTY = 1.0
DEFAULT_ZDR_UNCERTAINTY = 0.1
XTICKS = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360]
STATS_KEYS = ['N', 'low', 'median', 'high']
DEFAULT_DOMAIN = [-20, 60]
DEFAULT_ZMAX = 55
DEFAULT_DELZ = 1.0
DEFAULT_DELK = 0.1  # Resolution (dBZ) of I2 integrations (FSC)
DEFAULT_KMAX = 350  # Multiply by delk to get maximum block considered by FSC
MAX_CORRECTION = DEFAULT_DELK * DEFAULT_KMAX

#####################################

DEFAULT_KW = {'sweep': 0, 'dz': 'ZH', 'dr': 'DR', 'dp': 'DP', 'rh': 'RH',
              'kd': None, 'ld': None, 'sounding': DEFAULT_SND,
              'verbose': False, 'thresh_sdp': DEFAULT_SDP, 'fhc_T_factor': 1,
              'fhc_weights': dualpol.DEFAULT_WEIGHTS, 'fhc_name': 'FH',
              'band': 'S', 'fhc_method': 'hybrid', 'kdp_method': 'CSU',
              'bad': BAD, 'use_temp': True, 'rain_dp_thresh': 100,
              'drizzle_dz_thresh': 25, 'vr': 'VR', 'magnetron': False,
              'bin_width': DEFAULT_BINS, 'rain_kdp_thresh': DEFAULT_RAIN,
              'drizzle_kdp_thresh': DEFAULT_DRIZZ, 'rng_thresh': DEFAULT_RANGE,
              'dsd_flag': False, 'output': 100, 'rain_dz_thresh': 39,
              'liquid_ice_flag': False, 'precip_flag': False,
              'maxZ': DEFAULT_ZMAX, 'debug': False}

kwargs = np.copy(DEFAULT_KW)

"""
kwargs descriptions
-------------------
sweep = Sweep number to examine
dz = Name of reflectivity field in Py-ART radar object
dr = Name of differential reflectivity field in Py-ART radar object
dp = Name of differential phase field in Py-ART radar object
rh = Name of correlation coefficient field in Py-ART radar object
kd = Name of specific differential phase field in Py-ART radar obj. (if avail.)
ld = Name of linear depolarization ratio field in Py-ART radar obj. (if avail.)
vr = Name of Doppler velocity field in Py-ART radar object
verbose = Set to True for text notifications
thresh_sdp = Threshold for specific differential phase (can vary spatially)
fhc_T_factor = Extra weighting to be used for T in FHC calculations
fhc_weights = Weights used for each polarimetric & T field in FHC calculations
fhc_name = Name to give to newly created FHC field
band = Wavelength band of radar ('S' or 'C' supported)
fhc_method = Method to use in FHC calculations
kdp_method = Method to use in KDP calculations
bad = Bad data value
use_temp = Set to False to not use T in FHC
rain_dp_thresh = Differential phase threshold below which data will not be used
                 in KDP-method blockage calculation unless reflectivity exceeds
                 rain_dz_thresh
drizzle_dz_thresh = Reflectivity threshold above which data will not be used in
                    determination of ZDR blockage magnitude
magnetron = Set to True if transmitter was a magnetron and you can thus remove
            second-trip by filtering on Doppler velocity
bin_width = Width of each bin in azimuth degrees for determining blockage
rain_kdp_thresh = Two-element tuple denoting min/max KDP values to consider in
                  estimation of reflectivity blockage via KDP method
drizzle_kdp_thresh = Two-element tuple denoting min/max KDP values to consider
                     in estimation of ZDR blockage
rng_thresh = Two-element tuple or list of tuples (same size as number of bins
             in blockage calculation) indicating range (km) to consider for
             analysis of blockage
dsd_flag = Set to True to also retrieve DSD parameters via DualPol
output = Number of files to process at a time before outputting raw data for
         analysis
rain_dz_thresh = Reflectivity threshold below which data will not be used
                 in KDP-method blockage calculation unless differential phase
                 exceeds rain_dp_thresh
liquid_ice_flag = Set to True to also retrieve liquid/ice mass via DualPol
precip_flag = Set to True to also retrieve rainfall rate via DualPol
maxZ = Integer, highest Z value to consider in FSC meth. (to avoid ice contam.)
debug = Set to True to check the BeamBlockSingleVolume processing loop
"""

#####################################


class BeamBlockSingleVolume(object):

    """
    Core class that processes single volume of radar data and isolates all the
    individual gates that meet the user-specified criteria for rain (used for
    determining reflectivity blockage via KDP method) as well as for drizzle
    (used for determining differential reflectivity blockage for both the KDP
    and FSC methods). An additional mask is applied to obtain good rain/drizzle
    data for the FSC method.

    Many kwargs are passed to dualpol.DualPolRetrieval, which is used to
    calculate KDP (if necessary) and also do hydrometeor identification.

    Py-ART is used to ingest the individual radar files. The csu_radartools
    module is used to help filter bad data like insects.
    """

    def __init__(self, filename, **kwargs):
        """
        Must specify names of key polarimetric radar fields.
        KDP is optional - will calculate if needed.
        Will use default sounding provided with package if none provided.
        Expects UWYO format for soundings.
        """
        kwargs = dualpol.check_kwargs(kwargs, DEFAULT_KW)
        self.debug = kwargs['debug']
        self.verbose = kwargs['verbose']
        if kwargs['debug']:
            self.main_loop(filename, kwargs)
        else:
            try:
                self.main_loop(filename, kwargs)
            except:
                warn('Failure in reading or analyzing file, moving on ...')
                self.Fail = True
                print('BeamBlockSingleVolume Read/Processing Fail =',
                      self.Fail)

    def main_loop(self, filename, kwargs):
        """Broken off into separate method to simplify code debugging"""
        radar = pyart.io.read(filename)
        if kwargs['debug']:
            print('debug radar keys', radar.fields.keys())
        self.sweep = kwargs['sweep']
        self.maxZ = np.int32(kwargs['maxZ'])
        self.rain_total_pts = 0
        self.drizz_total_pts = 0
        self.radar = radar.extract_sweeps([self.sweep])
        if kwargs['fhc_name'] in self.radar.fields.keys():
            kwargs['fhc_flag'] = False
        self.retrieve = dualpol.DualPolRetrieval(self.radar, **kwargs)
        self.retrieve.name_vr = kwargs['vr']
        self.thresh_sdp = kwargs['thresh_sdp']
        self.get_bad_data_mask(magnetron=kwargs['magnetron'])
        self.get_2d_azimuth_and_range()
        self.get_bins(kwargs['bin_width'])
        self.get_range_mask(kwargs['rng_thresh'])
        # For KDP Method
        self.partition_rain_data(
            kwargs['rain_kdp_thresh'], kwargs['rain_dz_thresh'],
            kwargs['rain_dp_thresh'])
        self.group_rain_data()
        self.partition_drizzle_data(kwargs['drizzle_kdp_thresh'],
                                    kwargs['drizzle_dz_thresh'])
        self.group_drizzle_data()
        # For FSC Method
        self.partition_fsc_data()
        self.group_fsc_data()

    def get_2d_azimuth_and_range(self):
        """
        Two-dimensionalize 1-D azimuth and range to simplify masking.
        """
        az = self.radar.azimuth['data']
        rng = self.radar.range['data'] / RANGE_MULT
        self.range, self.azimuth = np.meshgrid(rng, az)

    def get_bins(self, bin_width):
        """
        Provided an azimuth bin width, develop the azimuth bin structure.
        Assumes azimuth stays between 0 and 360 deg.
        """
        self.bin_width = bin_width
        self.nbins = get_index(360.0, bin_width)
        self.azimuth_indices = get_index(self.azimuth, bin_width)
        cond = self.azimuth_indices == self.nbins
        self.azimuth_indices[cond] = 0
        if self.verbose:
            print('azimuth_indices =', self.azimuth_indices[0])

    def get_range_mask(self, rng_thresh):
        """
        For each azimuth bin, develop a mask that will allow us to filter out
        data from undesired ranges.
        """
        test = rng_thresh[0]
        if not hasattr(test, '__len__'):
            self.fix_rng_thresh(rng_thresh)
        else:
            self.range_thresh = rng_thresh
        cond = 0 * self.range
        for i in np.arange(np.shape(self.range)[0]):
            index_test = self.azimuth_indices[i][0]
            subrange = [self.range_thresh[index_test][0],
                        self.range_thresh[index_test][1]]
            if self.verbose:
                print('i, subrange =', i, subrange)
            try:
                cond[i] = np.logical_or(self.range[i] < subrange[0],
                                        self.range[i] > subrange[1])
            except ValueError:
                warn('Likely range_thresh is messed up')
                print('debug subrange, range', subrange, self.range)
        self.range_mask = cond.astype(bool)

    def fix_rng_thresh(self, rng_thresh):
        """
        If user did not supply range thresholds that varied by azimuth, then
        populate all azimuth bins with constant thresholds before developing
        the range mask.
        """
        dummy = []
        for i in np.arange(self.nbins):
            dummy.append(rng_thresh)
        self.range_thresh = dummy

    def partition_rain_data(self, rain_kdp_thresh, rain_dz_thresh,
                            rain_dp_thresh):
        """Produces mask for all useful rain data in volume"""
        self.fhc = self.retrieve.extract_unmasked_data(self.retrieve.name_fhc)
        self.dp = self.retrieve.extract_unmasked_data(self.retrieve.name_dp)
        self.kd = self.retrieve.extract_unmasked_data(self.retrieve.name_kd)
        self.rain_not_fhc = self.fhc != 2
        self.rain_not_kdp = np.logical_or(self.kd < rain_kdp_thresh[0],
                                          self.kd > rain_kdp_thresh[1])
        self.rain_not_dzdp = np.logical_and(self.dz < rain_dz_thresh,
                                            self.dp < rain_dp_thresh)
        self.rain_good_mask = consolidate_masks(
            [self.bad_data_mask, self.range_mask, self.rain_not_fhc,
             self.rain_not_kdp, self.rain_not_dzdp])

    def partition_drizzle_data(self, drizzle_kdp_thresh, drizzle_dz_thresh):
        """Produces mask for all useful drizzle data in volume"""
        self.drizzle_not_fhc = self.fhc != 1
        self.drizzle_not_kdp = np.logical_or(self.kd < drizzle_kdp_thresh[0],
                                             self.kd > drizzle_kdp_thresh[1])
        self.drizzle_not_dz = self.dz > drizzle_dz_thresh
        self.drizzle_good_mask = consolidate_masks(
            [self.bad_data_mask, self.range_mask, self.drizzle_not_fhc,
             self.drizzle_not_kdp, self.drizzle_not_dz])

    def partition_fsc_data(self):
        """Consolidate and invert the bad masks, then add more good masks."""
        self.fsc_not_fhc = np.logical_and(self.fhc != 1, self.fhc != 2)
        self.fsc_good_mask = consolidate_masks(
            [self.bad_data_mask, self.range_mask, self.fsc_not_fhc])
        self.kd_mask = self.kd > -100
        self.dr_mask = self.dr > -100
        self.dz_mask = np.logical_and(self.dz > -100, self.dz < self.maxZ)
        cond = np.logical_and(self.kd_mask, self.dr_mask)
        cond = np.logical_and(self.dz_mask, cond)
        self.fsc_good_mask = np.logical_and(self.fsc_good_mask, cond)

    def group_rain_data(self):
        """Applies rain mask to volume, and bins up by azimuth"""
        indices = self.azimuth_indices[self.rain_good_mask]
        dzvals = self.dz[self.rain_good_mask]
        self.rain_total_pts += len(dzvals)
        self.rain_binned_dz = []
        for i in np.arange(self.nbins):
            self.rain_binned_dz.append(dzvals[indices == i])

    def group_drizzle_data(self):
        """Applies drizzle mask to volume, and bins up by azimuth"""
        indices = self.azimuth_indices[self.drizzle_good_mask]
        drvals = self.dr[self.drizzle_good_mask]
        self.drizz_total_pts += len(drvals)
        self.drizzle_binned_dr = []
        for i in np.arange(self.nbins):
            self.drizzle_binned_dr.append(drvals[indices == i])

    def group_fsc_data(self):
        """Applies FSC mask to volume, and bins up by azimuth"""
        indices = self.azimuth_indices[self.fsc_good_mask]
        dzvals = self.dz[self.fsc_good_mask]
        drvals = self.dr[self.fsc_good_mask]
        kdvals = self.kd[self.fsc_good_mask]
        self.fsc_binned_data = {}
        self.fsc_binned_data['DZ'] = []
        self.fsc_binned_data['DR'] = []
        self.fsc_binned_data['KD'] = []
        for i in np.arange(self.nbins):
            self.fsc_binned_data['DZ'].append(dzvals[indices == i])
            self.fsc_binned_data['DR'].append(drvals[indices == i])
            self.fsc_binned_data['KD'].append(kdvals[indices == i])

    def get_bad_data_mask(self, magnetron=False):
        """Develops mask to identify insects, second trip, and noise"""
        self.dz = self.retrieve.extract_unmasked_data(self.retrieve.name_dz)
        self.dr = self.retrieve.extract_unmasked_data(self.retrieve.name_dr)
        self.insect_mask = csu_misc.insect_filter(self.dz, self.dr)
        if magnetron:
            vr_array = self.retrieve.extract_unmasked_data(
                self.retrieve.name_vr)
            self.trip2_mask = csu_misc.second_trip_filter_magnetron(vr_array)
            new_mask = np.logical_or(self.insect_mask, self.trip2_mask)
        else:
            new_mask = self.insect_mask
        if hasattr(self.retrieve, 'name_sdp'):
            if not hasattr(self.thresh_sdp, '__len__'):
                self.thresh_sdp = 0.0 * self.dz + self.thresh_sdp
            sdp_array = \
                self.retrieve.extract_unmasked_data(self.retrieve.name_sdp)
            self.sdp_mask = csu_misc.differential_phase_filter(
                sdp_array, self.thresh_sdp)
            new_mask = np.logical_or(new_mask, self.sdp_mask)
        self.bad_data_mask = new_mask

#####################################


class BlockStats(object):

    """
    Helper class that enables ingest and plotting of blockage statistics.
    """

    def __init__(self):
        """
        Purpose of this class is to provide useful common methods to other
        classes. Hence, there is no point to populate __init__().
        """
        pass

    def load_stats(self, filename):
        """
        Load KDP-method azimuth-based blockage stats from file.
        """
        dtype = {'names': ('Azimuth', 'N_r', 'low_r', 'median_r', 'high_r',
                           'N_d', 'low_d', 'median_d', 'high_d'),
                 'formats': ('float', 'float', 'float', 'float', 'float',
                             'float', 'float', 'float', 'float')}
        tmp_data = np.loadtxt(filename, dtype=dtype, skiprows=2, delimiter=',')
        self.Azimuth = tmp_data['Azimuth']
        self.rain_stats = {}
        self.rain_stats['N'] = tmp_data['N_r']
        self.rain_stats['low'] = tmp_data['low_r']
        self.rain_stats['median'] = tmp_data['median_r']
        self.rain_stats['high'] = tmp_data['high_r']
        self.drizzle_stats = {}
        self.drizzle_stats['N'] = tmp_data['N_d']
        self.drizzle_stats['low'] = tmp_data['low_d']
        self.drizzle_stats['median'] = tmp_data['median_d']
        self.drizzle_stats['high'] = tmp_data['high_d']

    def make_plots(self):
        """
        Makes a two-panel plot of reflectivity and differential reflectivity
        blockage as functions of azimuth.
        """
        # Make Rainfall Reflectivity Plots
        fig = plt.figure(figsize=(7, 9.5))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        fig, ax1 = self._plot_medians(
            fig, ax1, 'rain_stats', [0, 360], [0, 60],
            '(a) KDP Method - Rainfall', 'Azimuth (deg)', 'Reflectivity (dBZ)')
        fig, ax2 = self._plot_medians(
            fig, ax2, 'drizzle_stats', [0, 360], [-4, 4],
            '(b) KDP Method - Drizzle',
            'Azimuth (deg)', 'Differential Reflectivity (dB)')
        plt.tight_layout()
        plt.savefig(self.image_dir + 'block_kdp_method' + self.image_ext)
        plt.close()

    def _plot_medians(self, fig, ax, var, xlim, ylim, title, xlabel, ylabel):
        """
        Internal method that actually produces an individual panel in the
        two-panel plot.
        """
        var = getattr(self, var)
        for j in np.arange(len(self.Azimuth)):
            ax.plot([self.Azimuth[j], self.Azimuth[j]],
                    [var['low'][j], var['high'][j]], 'k-')
            ax.plot([self.Azimuth[j]], [var['median'][j]], 'bD', ms=3, mew=0)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xticks(XTICKS)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return fig, ax

#####################################


class BeamBlockMultiVolume(BlockStats):

    """
    This class facilitates the ingest of multiple radar volumes to compute beam
    blockage as functions of individual variables (Z, ZDR) and azimuth.

    Key Attributes
    --------------
    rain_refl - Reflectivity data meeting rain criteria implemented by
                BeamBlockSingleVolume, as a function of azimuth.
    drizz_zdr - Differential reflectivity data meeting drizzle criteria
                implemented by BeamBlockSingleVolume, as a function of azimuth.
    fsc_data - Dict containing reflectivity, differential reflectivity, and
               specific differential phase data for use with the FSC method.
               Keys are 'DZ', 'DR', and 'KD' and data are arranged by azimuth.
    Azimuth - Array of azimuth bins.
    rain_stats - Dict of rain statistics as functions of azimuth (KDP method).
    drizzle_stats - Same but for drizzle statistics (KDP/FSC method).
    """

    def __init__(self, files, image_dir='./', image_ext='.png', save=None,
                 stats_name='blockage_stats_kdp_method.txt', **kwargs):
        """
        Initializing method. This checks and populates keywords, & figures out
        whether a list of files to process was provided, or statistics file was
        provided. If the former, then the files will start being processed.
        If the latter, then the statistics will be ingested and plotted.
        """
        kwargs = dualpol.check_kwargs(kwargs, DEFAULT_KW)
        self.kwargs = kwargs
        self.check_file_list(files)
        self.image_dir = image_dir
        self.image_ext = image_ext
        self.save = save
        self.stats_name = stats_name
        # If just provided a stats file, ingest it and make a plot.
        if self.stats_flag:
            print('Read stats file, making a plot')
            self.make_plots()
        # If provided list of radar volume filenames, start processing!
        else:
            self.get_sounding_list()
            self.process_multiple_volumes()

    def process_multiple_volumes(self):
        """
        Main processing method. This loops through the list of radar volume
        files, iteratively calling the beam blockage correction routines
        supplied by the BeamBlockSingleVolume class, then compiling the data
        from each file into master arrays that are periodically saved to file,
        both as raw data and as statistics, as well as plotted up.
        """
        for i, filen in enumerate(self.file_list):
            print('i=', i, 'of', len(self.file_list)-1)
            print('radar file:', os.path.basename(filen), '& sounding:',
                  os.path.basename(self.sounding_list[i]))
            self.kwargs['sounding'] = self.sounding_list[i]
            bb = BeamBlockSingleVolume(filen, **self.kwargs)
            if not hasattr(bb, 'Fail'):  # This avoids crashing on a bum file
                # Initialize dicts
                if i == 0 or not hasattr(self, 'rain_total_pts'):
                    self.initialize_data_lists(bb)
                # Keep track of total points from each file
                self.rain_total_pts += bb.rain_total_pts
                self.drizz_total_pts += bb.drizz_total_pts
                print('Total KDP method points - rain, drizzle:',
                      self.rain_total_pts, self.drizz_total_pts)
                # Populate dicts
                for j in np.arange(len(self.Azimuth)):
                    if len(bb.rain_binned_dz[j]) > 0:
                        self.rain_refl[np.str(j)] = \
                            np.append(self.rain_refl[np.str(j)],
                                      bb.rain_binned_dz[j])
                    if len(bb.drizzle_binned_dr[j]) > 0:
                        self.drizz_zdr[np.str(j)] =\
                            np.append(self.drizz_zdr[np.str(j)],
                                      bb.drizzle_binned_dr[j])
                    if len(bb.fsc_binned_data['DZ'][j]) > 0:
                        self.fsc_data['DZ'][np.str(j)] = \
                            np.append(self.fsc_data['DZ'][np.str(j)],
                                      bb.fsc_binned_data['DZ'][j])
                        self.fsc_data['DR'][np.str(j)] = \
                            np.append(self.fsc_data['DR'][np.str(j)],
                                      bb.fsc_binned_data['DR'][j])
                        self.fsc_data['KD'][np.str(j)] = \
                            np.append(self.fsc_data['KD'][np.str(j)],
                                      bb.fsc_binned_data['KD'][j])
                # Periodic saving/plotting of data & statistics
                if i % 10 == 0 or i == len(self.file_list)-1:
                    print('Saving statistics to file and plotting an image')
                    self.get_statistics()
                    self.make_plots()
                    self.write_stats()
                if i % self.kwargs['output'] == 0 or \
                        i == len(self.file_list)-1:
                    if i != 0:  # Don't bother saving data if we just started
                        print('Writing raw data to file')
                        store = RawDataStorage(obj=self, filename=self.save)
            print('')

    def initialize_data_lists(self, bb):
        """
        Initialize data lists. Basic structure is a 1-D array of 1-D variable-
        length lists. The master array is the same size as the number azimuth
        analysis bins (e.g., BeamBlockSingleVolume.bin_width=1 means size 360).
        The secondary lists, which are unique to each az bin, can be of any
        length and that depends on the amount of data meeting the specified
        criteria in each bin. This will typically grow as the number of files
        processed increases.

        The save attribute basically flags whether the user will start from
        scratch with a list of radar volume files, or if initial arrays will be
        populated with blockage data from a previous run (e.g., on a different
        dataset from the same radar).
        """
        if self.save is not None:
            try:
                store = RawDataStorage(filename=self.save)
                self.Azimuth = store.Azimuth
                self.rain_refl = store.rain_refl
                self.drizz_zdr = store.drizz_zdr
                self.fsc_data = store.fsc_data
                self.rain_total_pts = 0
                self.drizz_total_pts = 0
                for key in self.rain_refl.keys():
                    self.rain_total_pts += len(self.rain_refl[key])
                    self.drizz_total_pts += len(self.drizz_zdr[key])
                if bb.nbins != len(self.Azimuth):
                    wstr = 'bb.nbins, self.Azimuth = ' + str(bb.nbins) + \
                           ', ' + str(self.Azimuth)
                    warn(
                         'Wrong number of azimuth bins, going to #FAIL - ' +
                         wstr)
            except:
                self._start_from_scratch(bb)
        else:
            self._start_from_scratch(bb)

    def get_statistics(self):
        """
        Method to compute blockage statistics as functions of azimuth. Focus
        is on median values along with 95% confidence intervals.
        """
        length = len(self.Azimuth)
        self.rain_stats = {}
        self.drizzle_stats = {}
        for key in STATS_KEYS:
            self.rain_stats[key] = np.zeros(length)
            self.drizzle_stats[key] = np.zeros(length)
        for j in np.arange(length):
            self.rain_stats['N'][j], self.rain_stats['low'][j], \
                self.rain_stats['median'][j], self.rain_stats['high'][j] = \
                calc_median_ci(self.rain_refl[np.str(j)])
            self.drizzle_stats['N'][j], self.drizzle_stats['low'][j], \
                self.drizzle_stats['median'][j],\
                self.drizzle_stats['high'][j] = \
                calc_median_ci(self.drizz_zdr[np.str(j)])

    def write_stats(self):
        """
        Method to write out a blockage statistics file.
        """
        fileobj = open(self.image_dir+self.stats_name, 'w')
        wstr = 'Azimuth,#Samples(Rain),CI_LO(Rain),Med.(Rain),CI_HI(Rain),' + \
               '#Samples(Drizz),CI_LO(Drizz),Median(Drizz),CI_HI(Drizz)'
        fileobj.write(wstr+'\n')
        fileobj.write('------------------------------------------------------')
        for j in np.arange(len(self.Azimuth)):
            wstr = '\n' + str(self.Azimuth[j]) + ',' + \
                str(self.rain_stats['N'][j]) + ',' + \
                str(self.rain_stats['low'][j]) + ',' + \
                str(self.rain_stats['median'][j]) + ',' + \
                str(self.rain_stats['high'][j]) + ',' + \
                str(self.drizzle_stats['N'][j]) + ',' + \
                str(self.drizzle_stats['low'][j]) + ',' + \
                str(self.drizzle_stats['median'][j]) + ',' + \
                str(self.drizzle_stats['high'][j])
            fileobj.write(wstr)
        fileobj.close()

    def check_file_list(self, files):
        """
        Checks first if stats file, then if not checks if argument is list or
        a single basestring (i.e., one file)
        """
        self.stats_flag = False
        if isinstance(files, six.string_types):
            try:
                self.load_stats(files)
                self.stats_flag = True
            except:
                self.file_list = [files]
        else:
            self.file_list = files

    def get_sounding_list(self):
        """
        User is responsible for inputing properly formatted string
        file names to program. Limited error checking is done.
        """
        sounding = self.kwargs['sounding']
        if isinstance(sounding, six.string_types):
            sndlist = []
            for i in np.arange(len(self.file_list)):
                sndlist.append(sounding)
        else:
            if len(sounding) != len(self.file_list):
                warn('sounding list not same length as file_list')
                sndlist = []
                for i in np.arange(len(self.file_list)):
                    sndlist.append(sounding[0])
            else:
                sndlist = sounding
        self.sounding_list = sndlist

    def _start_from_scratch(self, bb):
        self.rain_total_pts = 0
        self.drizz_total_pts = 0
        self.rain_refl = {}
        self.drizz_zdr = {}
        self.fsc_data = {}
        self.fsc_data['DZ'] = {}
        self.fsc_data['DR'] = {}
        self.fsc_data['KD'] = {}
        self.Azimuth = np.zeros(bb.nbins)
        for j in np.arange(len(self.Azimuth)):
            self.rain_refl[np.str(j)] = []
            self.drizz_zdr[np.str(j)] = []
            self.fsc_data['DZ'][np.str(j)] = []
            self.fsc_data['DR'][np.str(j)] = []
            self.fsc_data['KD'][np.str(j)] = []
            self.Azimuth[j] = bb.bin_width * j

#####################################


class RawDataStorage(object):

    """
    Class that facilitates saving and loading of azimuth-binned reflectivity
    and differential reflectivity data. Uses pickle module & stores the saved
    object as a binary file.
    """

    def __init__(self, obj=None, filename=None):
        """
        Determine what was just passed to class.
        """
        if obj is not None:
            if hasattr(obj, 'kwargs'):
                self.save_raw_data(obj, filename)
            else:
                warn('Not sure was passed proper BeamBlockMultiVolume object')
        else:
            self.load_raw_data(filename)

    def save_raw_data(self, obj, filename):
        """
        Saves data to pickled file.
        """
        self.populate_attributes(obj)
        if filename is None:
            filename = './temporary_block_data.dat'
        with gzip.open(filename+'.gz', 'wb') as f:
            pickle.dump(self, f)

    def load_raw_data(self, filename):
        """
        Loads data from pickled file.
        """
        if filename is None:
            filename = './temporary_block_data.dat'
        with gzip.open(filename+'.gz', 'rb') as f:
            loadobj = pickle.load(f)
        self.populate_attributes(loadobj)

    def populate_attributes(self, obj):
        """
        Only certain attributes are actually saved.
        """
        self.rain_refl = obj.rain_refl
        self.drizz_zdr = obj.drizz_zdr
        self.Azimuth = obj.Azimuth
        self.fsc_data = obj.fsc_data

#####################################


class MaskHelper(object):

    """Helper class to feed common methods for masking data to other classes"""

    def __init__(self):
        pass

    def get_azimuth_mask(self, azimuths=[(0, 360)]):
        """
        azimuths = list of tuples, each containing a span of unblocked azimuths
                   to consider in the determination of the median Zh to compare
                   against blocked azimuths.
        """
        for i, interval in enumerate(azimuths):
            new_cond = np.logical_and(self.Azimuth >= interval[0],
                                      self.Azimuth < interval[1])
            if i == 0:
                cond = new_cond
            else:
                cond = np.logical_or(cond, new_cond)
        self.azimuth_mask = cond

#####################################


class KdpMethodAnalysis(BlockStats, MaskHelper):

    def __init__(self, filename, azimuths=[(0, 360)]):
        """azimuths = list of 2-element tuples defining blocked azimuths"""
        self.load_stats(filename)
        self.get_azimuth_mask(azimuths)
        self.get_unblocked_medians()
        self.isolate_blocked_data()
        self.calc_blockage_correction()

    def get_unblocked_medians(self):
        """
        Find median values and their confidence intervals in unblocked regions.
        """
        self.unblocked_rain_stats = {}
        self.unblocked_drizz_stats = {}
        for key in STATS_KEYS:
            self.unblocked_rain_stats[key] = \
                self.rain_stats[key][self.azimuth_mask]
            self.unblocked_drizz_stats[key] = \
                self.drizzle_stats[key][self.azimuth_mask]
        self.rain_N, self.rain_low, self.rain_median, self.rain_high = \
            calc_median_ci(self.unblocked_rain_stats['median'])
        self.drizz_N, self.drizz_low, self.drizz_median, self.drizz_high = \
            calc_median_ci(self.unblocked_drizz_stats['median'])

    def isolate_blocked_data(self):
        """
        Isolate data within user-defined blocked regions.
        """
        self.blocked_rain_stats = {}
        self.blocked_drizz_stats = {}
        for key in STATS_KEYS:
            self.blocked_rain_stats[key] = \
                self.rain_stats[key][~self.azimuth_mask]
            self.blocked_drizz_stats[key] = \
                self.drizzle_stats[key][~self.azimuth_mask]
        self.blocked_azimuths = self.Azimuth[~self.azimuth_mask]

    def calc_blockage_correction(self):
        self.length = len(self.blocked_azimuths)
        self.zh_adjustments = Corrections(self, zh_flag=True)
        self.zdr_adjustments = Corrections(self, zh_flag=False)

    def print_out_data(self, filename='suggested_kdp_corrections.txt'):
        self.printout = PrintBlock(self, filename=filename, kdp_flag=True)

    def plot_corrections(self, save='kdp_corrections.png', return_flag=False,
                         method='standard', zh_range=[0, 35],
                         zdr_range=[-5, 5]):
        """
        method = Only relevant to KdpMethodAnalysis plots, can be either
                 'standard', 'strict', or 'loose'
        save = Full path and filename to save image file to
        zdr_range = Range of ZDR corrections displayed
        zh_range = Range of ZH corrections displayed
        return_flag = Set to True to return figure and axes objects
        """
        plot = _PlotBlock(self, save=save, method=method, zh_range=zh_range,
                          zdr_range=zdr_range)
        if return_flag:
            return plot.fig, plot.ax1, plot.ax2

    def edit_corrections(self, azimuths, offsets, var='DZ', method='standard'):
        """
        azimuths = list of azimuths to adjust
        offsets = list of offsets to change to, same size and order as azimuths
        var = change to 'DR', 'ZDR', or 'ZD' to adjust ZDR rather than ZH
        method = 'standard', 'strict', or 'loose'
        """
        if not hasattr(azimuths, '__len__'):
            azimuths = [azimuths]
        if not hasattr(offsets, '__len__'):
            offsets = [offsets]
        if len(azimuths) != len(offsets):
            warn('Arguments must be same size, fix and try again')
            return
        for i, azimuth in enumerate(azimuths):
            if var.upper() in ['DR', 'ZDR', 'ZD']:
                self.zdr_adjustments.suggested_corrections = \
                    self._correct(self.zdr_adjustments.suggested_corrections,
                                  method, offsets[i], azimuth)
            else:
                self.zh_adjustments.suggested_corrections = \
                    self._correct(self.zh_adjustments.suggested_corrections,
                                  method, offsets[i], azimuth)

    def _correct(self, adjust, method, offset, azimuth):
        diff = np.abs(azimuth - adjust['azimuth'])
        index = np.argmin(diff)
        try:
            adjust[method][index] = offset
        except KeyError:
            warn('Bad method used, not adjusting: '+method)
        return adjust

#####################################


class SelfConsistentAnalysis(MaskHelper):

    """
    Class to facilitate the diagnosis of partial beam blockage via the fully
    self-consistent (FSC) method as applied in Lang et al. (2009). The anchor
    reference for this technique is Giangrande and Ryzhkov (2005).

    Since this can be confusing, let's explain how we treat the data.

    *self.Azimuth, self.blocked_azimuths, and self.unblocked_azimuths are float
     arrays that give the actual azimuth associated with each bin.
    *self.azimuth_indices and its blocked and unblocked variants are integer
     arrays that provide the relevant bin number for each azimuth. This
     distinction with the above allows the user to set wider or narrower bins
     than 1 deg. Then string conversions of these indices are the dictionary
     keys for the zh_adjustments and other associated attributes.

    A simple code snippet to get the zh_adjustments in a more typical array
    format is as follows:

    zh_adjustments_array = []
    for index in self.blocked_azimuth_indices:
        key = str(index)
        zh_adjustments_array.append(self.zh_adjustments[key])

    Then you could, for example, easily plot zh_adjustments_array vs.
    self.blocked_azimuths.
    """

    def __init__(self, filename, azimuths=[(0, 360)]):
        """
        filename = Name of RawStorageData file
        azimuths = List of tuples, each containing a span of unblocked azimuths
                   to consider in the determination of the rainfall
                   self-consistency relationship between Z, Zdr, and Kdp.
        """
        data = RawDataStorage(filename=filename)
        self.fsc_data = data.fsc_data
        self.Azimuth = data.Azimuth
        self.get_azimuth_mask(azimuths)
        self.unblocked_azimuths = self.Azimuth[self.azimuth_mask]
        self.blocked_azimuths = self.Azimuth[~self.azimuth_mask]
        self.azimuth_indices = np.int32(np.arange(len(self.Azimuth)))
        self.unblocked_azimuth_indices = \
            self.azimuth_indices[self.azimuth_mask]
        self.blocked_azimuth_indices = self.azimuth_indices[~self.azimuth_mask]
        self.zh_range = np.int32(np.arange(DEFAULT_ZMAX))

    def regress_unblocked_medians(self, zmin=0):
        """
        Performs weighted multiple linear regression to obtain coefficients for
        self-consistency relationship:
        Zestimate = b * log10(Kdp) + c * Zdr + a (cf. Eq. 1, Lang et al., 2009)
        This is done using median values for unblocked data, with Zh being
        binned up and Kdp and Zdr medians within each bin being determined.

        zmin = Minimum Zh value to consider in regression. We control for
               Kdp <= 0 but just in case it is still noisy in low Zh, you
               can use zmin to filter that out.

        Key output
        ----------
        self.beta_hat = list of coefficients for above eq. (order = [b, c, a])
        self.a, self.b, self.c = more direct storage for these coefficients
        """
        self._populate_unblocked_data()
        self._get_medians()
        self._get_weights()
        mask = np.logical_and(self.unblocked_medians['KD'] > 0,
                              self.unblocked_medians['DR'] > -900)
        mask = np.logical_and(mask, self.unblocked_medians['DZ'] > zmin)
        self.logKdp = np.log10(self.unblocked_medians['KD'])
        self.beta_hat = \
            weighted_multiple_linear_regression(
                self.unblocked_medians['DZ'][mask],
                [self.logKdp[mask], self.unblocked_medians['DR'][mask]],
                self.weights[mask])
        self.define_coefficients(self.beta_hat[2], self.beta_hat[0],
                                 self.beta_hat[1])

    def define_coefficients(self, a, b, c):
        """
        Define a, b, and c in the self-consistency relationship:
        Z = a + b * log10(Kdp) + c * Zdr
        """
        self.a = a
        self.b = b
        self.c = c

    def plot_histogram_unblocked(self, save=None, domain=DEFAULT_DOMAIN,
                                 zmin=0):
        """
        Plots a simple 2-D histogram for visualizing the results & performance
        of the regress_unblocked_data() method.

        Keywords
        --------
        save = Name of image file to save figure to.
        domain = 2-element tuple to narrow down the display range of histogram.
        zmin = Minimum Zh to consider in regression (if redone).
        """
        if not hasattr(self, 'a'):
            self.regress_unblocked_medians(zmin=zmin)
        data = linear_calc(
            [self.b, self.c, self.a],
            [np.log10(self.unblocked_data['KD']), self.unblocked_data['DR']])
        lim = DEFAULT_DOMAIN
        fig = plt.figure(figsize=(6, 5))
        hist = plt.hist2d(data, self.unblocked_data['DZ'], bins=80,
                          cmap='Blues', range=[lim, lim], cmin=1)
        plt.plot(domain, domain, 'r--', lw=2)
        plt.xlim(domain)
        plt.ylim(domain)
        plt.xlabel('Predicted DZ (dBZ)')
        plt.ylabel('Actual DZ (dBZ)')
        plt.colorbar(label='Number of points')
        if save is not None:
            plt.savefig(save)

    def compute_I1_and_I2(self, kdp_analysis=None):
        """
        Compute the I1 and I2 integrals from Eq. 2-3 in Lang et al. (2009).
        Loops thru each of the blocked azimuth bins. First it will
        integrate the I1 side of the self-consistency relationship. Then it
        integrates the I2 side of the relationship and checks a variety of
        possible Zh offsets to see which one best matches I1. The resolution
        of the offset steps is delk, and kmax*delk is the maximum possible Zh
        offset examined. The zh_adjustments attribute then holds all of the Zh
        offsets to apply at each azimuth bin. If it is 0 then there is no
        offset to apply at that azimuth. If it is equal to delk*(kmax-1)
        then the solution likely did not converge and the block is too strong.
        Note that blocks ~30 dBZ or more are extremely difficult to correct,
        and even 20+ dBZ corrections are likely associated with greater error
        than weaker blocks.

        Assumes <ZDR> should match inherent offset in unblocked data when
        Z=0 dBZ. This is used to develop the zdr_offset attribute.

        Key Outputs
        -----------
        zh_adjustments = Add to Zh at each azimuth bin to correct for blockage.
        zdr_offsets = Difference between <ZDR(Z=0)> at each blocked azimuth and
                      the same parameter in unblocked data. Apply to the ZDR
                      data in blocked regions to force it to match unblocked
                      behavior. This is an alternative to correcting ZDR via
                      the KDP method.
        """

        self._intitialize_for_integrations()
        for i, index in enumerate(self.blocked_azimuth_indices):
            key = str(index)
            self._get_expected_vals(key)
            self.I1[key] = np.dot(self.expected_kdp[key],
                                  DEFAULT_DELZ*self.number_of_obs[key])
            self._match_I2(i, key, kdp_analysis=kdp_analysis)
        print('Zh adjustments determined & stored as zh_adjustments')
        print('ZDR offsets stored as zdr_offsets attribute')

    def print_out_data(self, filename='suggested_fsc_corrections.txt'):
        self.printout = PrintBlock(self, filename=filename, kdp_flag=False)

    def plot_corrections(self, save='fsc_corrections.png', return_flag=False,
                         zdr_range=[-5, 5], zh_range=[0, 35]):
        """
        save = Full path and filename to save image file to
        zdr_range = Range of ZDR corrections displayed
        zh_range = Range of ZH corrections displayed
        return_flag = Set to True to return figure and axes objects
        """
        plot = _PlotBlock(self, save=save, zdr_range=zdr_range,
                          zh_range=zh_range)
        if return_flag:
            return plot.fig, plot.ax1, plot.ax2

    def edit_corrections(self, azimuths, offsets, var='DZ'):
        """
        azimuths = list of azimuths to adjust
        offsets = list of offsets to change to, same size and order as azimuths
        var = change to 'DR', 'ZDR', or 'ZD' to adjust ZDR rather than ZH
        """
        if not hasattr(azimuths, '__len__'):
            azimuths = [azimuths]
        if not hasattr(offsets, '__len__'):
            offsets = [offsets]
        if len(azimuths) != len(offsets):
            warn('Arguments must be same size, fix and try again')
            return
        for i, azimuth in enumerate(azimuths):
            diff = np.abs(azimuth-self.blocked_azimuths)
            index = np.argmin(diff)
            key = str(self.blocked_azimuth_indices[index])
            if var.upper() == 'DR' or var.upper() == 'ZDR' or \
                    var.upper() == 'ZD':
                self.zdr_offsets[key] = offsets[i]
            else:
                self.zh_adjustments[key] = offsets[i]

    def _intitialize_for_integrations(self):
        if not hasattr(self, 'a'):
            self.regress_unblocked_medians()
        self.kmax = DEFAULT_KMAX
        self.delk = DEFAULT_DELK
        self.delz = DEFAULT_DELZ
        self.zh_adjustments = {}
        self.zdr_offsets = {}
        self.I1 = {}
        self.I2 = {}
        self.expected_zdr = {}
        self.expected_kdp = {}
        self.number_of_obs = {}

    def _populate_unblocked_data(self):
        """
        Populates the unblocked data dictionary by consolidating all data from
        unblocked azimuths. Performs a simple test to make sure there is no
        weirdness in terms of data structure.
        """
        self.unblocked_data = {}
        self.unblocked_data['DZ'] = []
        self.unblocked_data['DR'] = []
        self.unblocked_data['KD'] = []
        for index in self.unblocked_azimuth_indices:
            self.unblocked_data['DZ'] = np.append(
                self.unblocked_data['DZ'], self.fsc_data['DZ'][str(index)])
            self.unblocked_data['DR'] = np.append(
                self.unblocked_data['DR'], self.fsc_data['DR'][str(index)])
            self.unblocked_data['KD'] = np.append(
                self.unblocked_data['KD'], self.fsc_data['KD'][str(index)])
        self._check_lengths()

    def _check_lengths(self):
        """Simple check to make sure everything is the same size"""
        if (len(self.unblocked_data['KD']) != len(self.unblocked_data['DZ']) or
                len(self.unblocked_data['KD']) !=
                len(self.unblocked_data['DR']) or
                len(self.unblocked_data['DZ']) !=
                len(self.unblocked_data['DR'])):
            wstr = str(self.unblocked_data['DZ']) + ' ' + \
                str(self.unblocked_data['DR']) + ' ' + \
                str(self.unblocked_data['KD'])
            warn('DZ DR KD not equal length, going to #FAIL: ' + wstr)

    def _get_medians(self):
        """Determine median Zdr/Kdp as function of Z for unblocked regions"""
        self.unblocked_medians = {}
        self.unblocked_medians['DZ'] = self.zh_range + 0.5
        self.unblocked_medians['KD'] = np.zeros(len(self.zh_range)) - 999.0
        self.unblocked_medians['DR'] = np.zeros(len(self.zh_range)) - 999.0
        self.unblocked_medians['N'] = np.zeros(len(self.zh_range))
        self._loop_to_get_medians()

    def _loop_to_get_medians(self):
        indices = np.int32(np.floor(self.unblocked_data['DZ']))
        for index in self.zh_range:
            mask = indices == index
            self.unblocked_medians['N'][index] = \
                np.size(self.unblocked_data['DZ'][mask])
            if self.unblocked_medians['N'][index] > 1:
                self.unblocked_medians['KD'][index] = \
                    np.median(self.unblocked_data['KD'][mask])
                self.unblocked_medians['DR'][index] = \
                    np.median(self.unblocked_data['DR'][mask])

    def _get_weights(self):
        """Weight by number of Kdp/Zdr observations in each Z bin"""
        self.weights = self.unblocked_medians['N'] * \
            self.unblocked_medians['KD']
        self.weights = self.weights / self.weights.sum()

    def _match_I2(self, ii, key, kdp_analysis=None):
        """Loop thru all possible reflectivity offsets up to delk*kmax"""
        self._get_zdr_offset(ii, key, kdp_analysis)
        itest = np.zeros(self.kmax)
        for i in np.arange(len(self.zh_range)):
            for k in np.arange(self.kmax):
                dztest = k * self.delk + self.zh_range + 0.5
                expon = dztest[i] / self.b - self.a / self.b - self.c * \
                    (self.expected_zdr[key][i] + self.zdr_offsets[key]) / \
                    self.b
                itest[k] += self.delz * self.number_of_obs[key][i] * \
                    10.0**expon
        diff = np.abs(self.I1[key] - itest)
        kmin = np.argmin(diff)
        # print(key, self.I1[key], itest[kmin], diff[kmin], self.delk*kmin)
        self.I2[key] = itest[kmin]
        self.zh_adjustments[key] = self.delk * kmin

    def _get_zdr_offset(self, ii, key, kdp_analysis):
        # Make sure this doesn't blow up if Zdr data missing at low Z values
        mask = self.unblocked_medians['DR'] > -100
        unblocked_zdr = self.unblocked_medians['DR'][mask]
        # Get zdr_offset at this azimuth
        if kdp_analysis is None:
            # If missing Zdr data at low Z, this will just be 0
            # Leaving this for now, not sure if worth fixing ...
            self.zdr_offsets[key] = unblocked_zdr[0] - \
                self.expected_zdr[key][0]
        else:
            # Use kdp_analysis-determined Zdr offset instead
            azdiff = np.abs(
                self.blocked_azimuths[ii] -
                kdp_analysis.zdr_adjustments.suggested_corrections['azimuth'])
            azin = np.argmin(azdiff)
            testzdr = 1.0 * \
                kdp_analysis.zdr_adjustments.suggested_corrections[
                    'standard'][azin]
            if testzdr > -100:
                self.zdr_offsets[key] = testzdr
            else:
                self.zdr_offsets[key] = 0.0

    def _get_expected_vals(self, key):
        """Obtain expected values for Zdr/Kdp for each blocked azimuth bin"""
        tmp_dz = np.array(self.fsc_data['DZ'][key])
        tmp_kd = np.array(self.fsc_data['KD'][key])
        tmp_dr = np.array(self.fsc_data['DR'][key])
        self.expected_zdr[key] = np.zeros(DEFAULT_ZMAX)
        self.expected_kdp[key] = np.zeros(DEFAULT_ZMAX)
        self.number_of_obs[key] = np.zeros(DEFAULT_ZMAX)
        dz_int = np.int32(np.rint(tmp_dz))
        for index in self.zh_range:
            mask = dz_int == index
            self.number_of_obs[key][index] = len(tmp_dr[mask])
            if self.number_of_obs[key][index] > 1:
                self.expected_zdr[key][index] = np.median(tmp_dr[mask])
                self.expected_kdp[key][index] = np.median(tmp_kd[mask])

#####################################


class Corrections(object):

    """
    Calculate the suggested reflectivity/ZDR corrections. There are currently
    three approaches, which will tend to agree in well-behaved data (i.e.,
    confidence intervals are narrow). In data with wide confidence intervals,
    some corrections may not end up being suggested due to lack of certainty.
    standard: Difference between median reflectivity/ZDR in blocked azimuth
              & median of unblocked azimuths is > 1 dBZ or 0.1 dB (default) &
              the difference is greater than half the 95% confidence interval
              at that azimuth.
    loose: standard conditions apply plus the difference between the high
            value in the 95% interval at the blocked azimuth and the median in
            unblocked azimuths is still greater than half the unblocked
            confidence interval.
    strict: Difference between median reflectivity/ZDR in blocked azimuth
            and median of unblocked azimuths is > 1 dBZ or 0.1 dB (default) &
            the difference is greater than the the difference between the high
            value in the 95% interval at the blocked azimuth and the median
            in unblocked azimuths.
    """

    def __init__(self, kdp_analysis, zh_flag=True, verbose=False):
        """
        kdp_analysis = KdpMethodAnalysis object
        zh_flag = Set to false to look at Zdr instead of Zh
        verbose = Set to True for additional text output
        """
        self.blocked_azimuths = kdp_analysis.blocked_azimuths
        self.length = kdp_analysis.length
        self._check_var(kdp_analysis, zh_flag, verbose)
        self.suggested_corrections = {}
        self.suggested_corrections['azimuth'] = self.blocked_azimuths
        for key in ['standard', 'strict', 'loose']:
            self.suggested_corrections[key] = np.zeros(self.length)
        self.determine_differences()
        self.determine_corrections(verbose=verbose)

    def determine_differences(self):
        """
        Prep work for determining the corrections ...
        Figure out first the differences!
        """
        self.difference = self.blocked_stats['median'] - self.median
        self.conf_bl = self.blocked_stats['high'] - self.blocked_stats['low']
        self.diff_ci_hi = self.blocked_stats['high'] - self.median
        self.conf_unbl = np.max(self.unblocked_stats['high'] -
                                self.unblocked_stats['low'])

    def determine_corrections(self, verbose=False):
        """
        Performs standard, loose, and strict correction calculations
        """
        for i, az in enumerate(self.blocked_azimuths):
            if np.abs(self.difference[i]) >= 1.0 * self.uncertainty:
                # 'standard' correction
                if np.abs(self.difference[i]) > 0.5 * self.conf_bl[i]:
                    self.suggested_corrections['standard'][i] = \
                           -1.0 * self.difference[i]
                    # 'loose' correction - subset of 'standard'
                    if np.abs(self.diff_ci_hi[i]) > 0.5 * self.conf_unbl:
                        self.suggested_corrections['loose'][i] = \
                                -1.0 * self.difference[i]
                # 'strict' correction
                if np.abs(self.difference[i]) > np.abs(
                        self.diff_ci_hi[i]):
                    self.suggested_corrections[
                        'strict'][i] = -1.0 * self.difference[i]
            for key in self.suggested_corrections.keys():
                if key != 'azimuth':
                    mask = np.abs(
                        self.suggested_corrections[key]) > MAX_CORRECTION
                    self.suggested_corrections[key][mask] = BAD
            if verbose:
                print(az, ["%.2f" % np.round(
                    self.suggested_corrections[key][i],
                    decimals=2) for key in self.suggested_corrections.keys()
                    if not key == 'azimuth'])

    def _check_var(self, kdp_analysis, zh_flag, verbose):
        """Choose the right variables to consider depending on zh_flag"""
        if zh_flag:
            self.blocked_stats = kdp_analysis.blocked_rain_stats
            self.unblocked_stats = kdp_analysis.unblocked_rain_stats
            self.median = kdp_analysis.rain_median
            self.uncertainty = DEFAULT_ZH_UNCERTAINTY
        else:
            self.blocked_stats = kdp_analysis.blocked_drizz_stats
            self.unblocked_stats = kdp_analysis.unblocked_drizz_stats
            self.median = kdp_analysis.drizz_median
            self.uncertainty = DEFAULT_ZDR_UNCERTAINTY

#####################################


class PrintBlock(object):

    """
    Class to print out blockage correction data to file as well as screen.
    Powers the print_out_data() methods in KdpMethodAnalysis and
    SelfConsistentAnalysis classes.
    """

    def __init__(self, obj, kdp_flag=True,
                 filename='suggested_corrections.txt'):
        """
        obj = KdpMethodAnalysis or SelfConsistentAnalysis object
        kdp_flag = Switch between printing KdpMethodAnalysis or
                   SelfConsistentAnalysis corrections
        filename = Name of file to write to, including path
        """
        if kdp_flag:
            self.print_kdp_data(obj, filename)
        else:
            self.print_fsc_data(obj, filename)

    def print_kdp_data(self, obj, filename):
        """
        Prints out Zh & Zdr corrections for all three approaches:
        standard, strict, & loose (values in dBZ or dB)
        """
        shortZH = obj.zh_adjustments.suggested_corrections
        shortDR = obj.zdr_adjustments.suggested_corrections
        fileobj = open(filename, 'w')
        wstr = 'Azimuth, ZH(std), ZH(strict), ZH(loose), ' + \
            'Azimuth, DR(std), DR(strict), DR(loose)'
        fileobj.write(wstr+'\n')
        print(wstr)
        wstr = '--------------------------------------------------------'
        fileobj.write(wstr+'\n')
        print(wstr)
        for i, az in enumerate(shortZH['azimuth']):
            wstr = str(az) + ',' + str(shortZH['standard'][i]) + ',' + \
                   str(shortZH['strict'][i]) + ',' + str(shortZH['loose'][i]) \
                   + ',' + str(shortDR['azimuth'][i]) + ',' + \
                   str(shortDR['standard'][i]) + ',' + \
                   str(shortDR['strict'][i]) + ',' + str(shortDR['loose'][i])
            fileobj.write(wstr+'\n')
            print(wstr)
        fileobj.close()
        print('')
        print('Data also written to', filename)

    def print_fsc_data(self, obj, filename):
        """
        Prints out Zh & Zdr corrections for SelfConsistentAnalysis, along with
        I1 and matched I2 integrated values.
        """
        fileobj = open(filename, 'w')
        wstr = 'Azimuth, I1, I2, ZH(fsc), DR(fsc)'
        fileobj.write(wstr+'\n')
        print(wstr)
        wstr = '--------------------------------------------------------'
        fileobj.write(wstr+'\n')
        print(wstr)
        for i, index in enumerate(obj.blocked_azimuth_indices):
            key = str(index)
            wstr = str(obj.blocked_azimuths[i]) + ',' + \
                str(obj.I1[key]) + ',' + str(obj.I2[key]) + ',' + \
                str(obj.zh_adjustments[key]) + ',' + str(obj.zdr_offsets[key])
            fileobj.write(wstr+'\n')
            print(wstr)
        fileobj.close()
        print('')
        print('Data also written to', filename)

#####################################


class _PlotBlock(object):

    """
    Helper class to execute the plotting of suggested blockage corrections
    from the KdpMethodAnalysis and SelfConsistentAnalysis classes.
    """

    def __init__(self, obj, method='standard', save='corrections.png',
                 zh_range=[0, 35], zdr_range=[-5, 5]):
        """
        obj = KdpMethodAnalysis or SelfConsistentAnalysis object
        method = Only relevant to KdpMethodAnalysis plots, can be either
                 'standard', 'strict', or 'loose'
        save = Full path and filename to save image file to
        zdr_range = Range of ZDR corrections displayed
        zh_range = Range of ZH corrections displayed
        """
        self.object = obj
        self.zh_range = zh_range
        self.zdr_range = zdr_range
        if hasattr(obj, 'zdr_adjustments'):
            self._translate_kdp_data(method)
        elif hasattr(obj, 'I1'):
            self._translate_fsc_data()
        else:
            warn('Not a KdpMethodAnalysis or SelfConsistentAnalysis object')
            return
        self._plot_corrections(save)

    def _plot_corrections(self, save='corrections.png'):
        """Simple plotting of Zh and Zdr corrections via matplotlib"""
        self.fig = plt.figure(figsize=(13, 6))
        self.ax1 = self.fig.add_subplot(211)
        self.ax1.plot(self.azimuths_zh, self.zh, 'kD', ms=3, label=self.label)
        plt.xlim(0, 360)
        plt.ylim(self.zh_range)
        plt.xlabel('Azimuth (deg)')
        plt.ylabel('Suggested Zh Correction (dBZ)')
        plt.legend(numpoints=1, loc='upper left')
        plt.title('(a) Reflectivity Corrections, '+self.label)
        plt.xticks(30.0*np.arange(13))
        self.ax2 = self.fig.add_subplot(212)
        self.ax2.plot(self.azimuths_zdr, self.zdr, 'kD',
                      ms=3, label=self.label)
        self.ax2.plot([0, 360], [0, 0], 'k--')
        plt.xlim(0, 360)
        plt.ylim(self.zdr_range)
        plt.xlabel('Azimuth (deg)')
        plt.ylabel('Suggested Zdr Correction (dB)')
        plt.legend(numpoints=1, loc='upper left')
        plt.title('(b) Differential Reflectivity Corrections, '+self.label)
        plt.xticks(30.0*np.arange(13))
        plt.tight_layout()
        plt.savefig(save)

    def _translate_kdp_data(self, method):
        """Translates the KdpMethodAnalysis object to plottable arrays"""
        self.label = 'KDP Method'
        shortZH = self.object.zh_adjustments.suggested_corrections
        shortDR = self.object.zdr_adjustments.suggested_corrections
        self.azimuths_zh = shortZH['azimuth']
        self.zh = shortZH[method]
        self.azimuths_zdr = shortDR['azimuth']
        self.zdr = shortDR[method]

    def _translate_fsc_data(self):
        """Translates the SelfConsistentAnalysis object to plottable arrays"""
        self.label = 'FSC Method'
        self.azimuths_zh = self.object.blocked_azimuths
        self.azimuths_zdr = self.azimuths_zh
        self.zh = 0.0 * self.azimuths_zh
        self.zdr = 0.0 * self.azimuths_zdr
        for i, index in enumerate(self.object.blocked_azimuth_indices):
            key = str(index)
            self.zh[i] = self.object.zh_adjustments[key]
            self.zdr[i] = self.object.zdr_offsets[key]

#####################################


def consolidate_masks(bad_masks):
    """
    bad_masks = list of bad data masks to apply to array.
    Returns a mask that will provide all good data.
    """
    if not hasattr(bad_masks, '__len__'):
        return ~bad_masks
    else:
        for i, mask in enumerate(bad_masks):
            if i == 0:
                new_mask = mask
            else:
                new_mask = np.logical_or(mask, new_mask)
        return ~new_mask


def get_index(azimuth, bin_width):
    """
    Obtains az bin indices, provided a bin_width. Rounds to nearest integer.
    For example, for 1-deg width, the 30-deg azimuth bin spans 29.5-30.5 deg.
    """
    index = azimuth / bin_width
    return np.int32(np.rint(index))


def calc_median_ci(array):
    """Calculate 95% confidence interval for median"""
    Nint = len(array)
    if Nint > 1:
        N = np.float(Nint)
        median = np.median(array)
        m1 = np.floor(0.5*(N-1.0))
        new_array = sorted(array)
        interv = 1.96 * np.sqrt(N) / 2.0
        index = np.int32(np.rint(np.float(m1) - interv))
        if index < 0:
            index = 0
        lo = new_array[index]
        index = np.int32(np.rint(np.float(m1+1) + interv))
        if index >= Nint:
            index = Nint - 1
        hi = new_array[index]
        # sigma = np.std(array)
        # lo = median - 1.96 * 1.25 * sigma / np.sqrt(N)
        # hi = median + 1.96 * 1.25 * sigma / np.sqrt(N)
        return Nint, lo, median, hi
    else:
        return Nint, BAD, BAD, BAD


def multiple_linear_regression(dependent_var, independent_vars):
    """
    Simple unweighted least-squares multiple linear regression using
    numpy.linalg module.
    """
    y = dependent_var
    x = independent_vars
    X = np.column_stack(x + [[1] * len(x[0])])
    beta_hat = np.linalg.lstsq(X, y)[0]
    return beta_hat


def weighted_multiple_linear_regression(dependent_var, independent_vars,
                                        weights):
    """
    Simple weighted least-squares multiple linear regression using
    statsmodel module.
    """
    Y = dependent_var
    x = independent_vars
    X = np.column_stack(x+[[1]*len(x[0])])
    wls_model = sm.WLS(Y, X, weights=weights)
    results = wls_model.fit()
    return results.params


def linear_calc(coef, independent_vars):
    """
    Use np.dot() to quickly evaluate a linear expression to recover the
    dependent variable from any number of independent variables.
    """
    x = independent_vars
    X = np.column_stack(x+[[1]*len(x[0])])
    return np.dot(X, coef)

#####################################
