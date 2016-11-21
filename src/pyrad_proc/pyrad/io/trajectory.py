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
    wgs84_lat_rad : Array of floats
        WGS84 latitude samples in radian
    wgs84_lon_rad : Array of floats
        WGS84 longitude samples in radian
    wgs84_alt_m : Array of floats
        WGS84 altitude samples in m

    Methods:
    --------
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
        self.wgs84_lat_rad = []
        self.wgs84_lon_rad = []
        self.wgs84_alt_m = []

        try:
            self._read_traj()
        except:
            raise

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
            raise Exception("ERROR: Could not find|open trajectory file '"+self.filename+"'")

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
                    print("WARNING: Format error in trajectory file '%s' on line '%s'"
                          % (self.filename, line))
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

                self.wgs84_lat_rad.append(float(mm.group(3)))
                self.wgs84_lon_rad.append(float(mm.group(4)))
                self.wgs84_alt_m.append(float(mm.group(5)))
        except:
            raise
        finally:
            tfile.close()
            if loc_set:
                locale.setlocale(locale.LC_ALL, loc)  # restore saved locale
