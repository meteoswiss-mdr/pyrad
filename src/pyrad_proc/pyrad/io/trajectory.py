"""
pyrad.io.trajectory
===================

Trajectory class implementation for reading trajectory file.
Converting to different coordinate systems.

.. autosummary::
    :toctree: generated/

    Trajectory

"""

import re
import datetime
import locale
import numpy as np
from warnings import warn

import pyart


class Trajectory(object):
    """
    A class for reading and handling trajectory data from a file.

    Attributes
    ----------
    filename : str
        Path and name of the trajectory definition file
    starttime : datetime
        Start time of trajectory processing.
    endtime : datetime
        End time of trajectory processing.
    timevector : Array of datetime objects
        Array containing the trajectory time samples
    wgs84_lat_deg : Array of floats
        WGS84 latitude samples in radian
    wgs84_lon_deg : Array of floats
        WGS84 longitude samples in radian
    wgs84_alt_m : Array of floats
        WGS84 altitude samples in m

    Methods:
    --------
    add_radar : Add a radar
    get_start_time : Return time of first trajectory sample
    get_end_time : Return time of last trajectory sample
    _read_traj : Read trajectory from file

    """

    def __init__(self, filename, starttime=None, endtime=None):
        """
        Initalize the object.

        Parameters
        ----------
        filename : str
            Filename containing the trajectory samples.
        starttime : datetime
            Start time of trajectory processing. If not given, use
            the time of the first trajectory sample.
        endtime : datetime
            End time of trajectory processing. If not given, use
            the time of the last trajectory sample.
        """

        self.filename = filename
        self.starttime = starttime
        self.endtime = endtime

        self.timevector = []
        self.wgs84_lat_deg = np.array([], dtype=float)
        self.wgs84_lon_deg = np.array([], dtype=float)
        self.wgs84_alt_m = np.array([], dtype=float)

        self.radar_list = []

        try:
            self._read_traj()
        except:
            raise

    def add_radar(self, radar):
        """
        Add the coordinates (WGS84 longitude, latitude and non WGS84 altitude)
        of a radar to the radar_list.

        Parameters
        ----------
        radar : pyart radar object
            containing the radar coordinates

        """

        # Check if radar location is already in the radar list
        for rad in self.radar_list:
            if (rad.location_is_equal(radar.latitude['data'][0],
                                      radar.longitude['data'][0],
                                      radar.altitude['data'][0])):
                warn("WARNING: Tried to add the same radar twice to the"
                     " radar list")
                return

        rad = _Radar_Trajectory(radar.latitude['data'][0],
                                radar.longitude['data'][0],
                                radar.altitude['data'][0],
                                len(self.timevector))
        self.radar_list.append(rad)

        proj = dict({'proj': 'pyart_aeqd',
                     'lon_0': radar.longitude['data'][0],
                     'lat_0': radar.latitude['data'][0]})

        xvec, yvec = pyart.core.geographic_to_cartesian(
            self.wgs84_lon_deg, self.wgs84_lat_deg, proj)

        rvec, azvec, elvec = pyart.core.cartesian_to_antenna(
            xvec, yvec, self.wgs84_alt_m - rad.altitude)

        return

    def get_start_time(self):
        """
        Get time of first trajectory sample.

        Return
        ------
        datetime object
        """

        return self.timevector[0]

    def get_end_time(self):
        """
        Get time of last trajectory sample.

        Return
        ------
        datetime object
        """

        return self.timevector[-1]

    def _read_traj(self):
        """
        Read trajectory from file

        File format
        -----------
        Comment symbol: '#'
        Columns:
        1. Day (UTC) format: DD-MMM-YYYY
        2. Seconds since midnight (UTC)
        3. WGS84 latitude [radians]
        4. WGS84 longitude [radians]
        5. WGS84 altitude [m]

        """

        # check if the file can be read
        try:
            tfile = open(self.filename, "r")
        except:
            raise Exception("ERROR: Could not find|open trajectory file '" +
                            self.filename+"'")

        repat = re.compile("(\d+\-[A-Za-z]+\-\d+)\s+([\d\.]+)\s+"
                           "([\-\d\.]+)\s+([\-\d\.]+)\s+([\-\d\.]+)")

        try:
            loc_set = False
            loc = locale.getlocale()  # get current locale
            if (loc[0] != 'en_US'):
                try:
                    locale.setlocale(locale.LC_ALL, ('en_US', 'UTF-8'))
                except Exception as ee:
                    raise Exception("ERROR: Cannot set local 'en_US': %s"
                                    % str(ee))
                loc_set = True

            recording_started = True
            if (self.starttime is not None):
                recording_started = False
            recording_check_stop = False
            if (self.endtime is not None):
                recording_check_stop = True

            for line in tfile:
                line = line.strip()

                if len(line) == 0:
                    continue

                # ignore comments
                if line.startswith('#'):
                    continue

                line = line.partition('#')[0]  # Remove comments
                line = line.strip()

                mm = repat.match(line)
                if not mm:
                    print("WARNING: Format error in trajectory file '%s'"
                          " on line '%s'" % (self.filename, line),
                          file=sys.stderr)
                    continue

                # Get time stamp
                try:
                    sday = datetime.datetime.strptime(mm.group(1), "%d-%b-%Y")
                except Exception as ee:
                    print(datetime.datetime.utcnow().strftime("%d-%b-%Y"))
                    raise Exception("ERROR: Format error in traj file '%s' "
                                    "on line '%s' (%s)"
                                    % (self.filename, line, str(ee)))

                sday += datetime.timedelta(seconds=float(mm.group(2)))

                if (not recording_started):
                    if (sday < self.starttime):
                        continue
                    else:
                        recording_started = True

                if (recording_check_stop):
                    if (sday > self.endtime):
                        break

                self.timevector.append(sday)

                self.wgs84_lat_deg = np.append(
                    self.wgs84_lat_deg, [float(mm.group(3)) * 180. / np.pi])
                self.wgs84_lon_deg = np.append(
                    self.wgs84_lon_deg, [float(mm.group(4)) * 180. / np.pi])
                self.wgs84_alt_m = np.append(
                    self.wgs84_alt_m, [float(mm.group(5))])
        except:
            raise
        finally:
            tfile.close()
            if loc_set:
                locale.setlocale(locale.LC_ALL, loc)  # restore saved locale


class _Radar_Trajectory:
    """
    A class for holding the trajectory data assigned to a radar.

    Attributes
    ----------
    latitude : float
       WGS84 latitude [deg]
    longitude : float
       WGS84 longitude [deg]
    altitude : float
       altitude [m] (non WGS84)
    elevation_vec : float list
       Elevation values of the trajectory samples
    azimuth_vec : float list
       Azimuth values of the trajectory samples
    range_vec : float list
       Range values of the trajectory samples

    Methods:
    --------
    location_is_equal
    """

    def __init__(self, lat, lon, alt, nsamps):
        """
        Initalize the object.

        Parameters
        ----------
        lat, lon , alt : radar location coordinates
        nsamps : number of samples
        """

        self.latitude = lat
        self.longitude = lon
        self.altitude = alt

        self.elevation_vec = []
        self.azimuth_vec = []
        self.range_vec = []

    def location_is_equal(self, lat, lon, alt):
        """
        Check if the given coordinates are the same.

        Parameters
        ----------
        lat, lon , alt : radar location coordinates

        Return
        ------
        True if the radar location is equal, False otherwise
        """

        lat_tol = 0.002  # [deg]
        lon_tol = 0.002  # [deg]
        alt_tol = 2.0    # [m]

        if ((np.abs(self.latitude-lat) > lat_tol) or
                (np.abs(self.longitude-lon) > lon_tol) or
                (np.abs(self.altitude-alt) > alt_tol)):
            return False
        else:
            return True

    def add_traj_sample(self, el, az, rr):
        """
        Append a new trajectory sample in radar coordinates to the
        trajectory vectors.

        Parameters
        ----------
        el, az, rr : elevation, azimuth and range
        """

        self.elevation_vec.append(el)
        self.azimuth_vec.append(az)
        self.range_vec.append(rr)
