"""
pyrad.io.timeseries
===================

TimeSeries class implementation for holding timeseries data.

.. autosummary::
    :toctree: generated/

    TimeSeries

"""


from ..graph.plots import plot_timeseries


class TimeSeries(object):
    """
    Holding timeseries data and metadata.

    Attributes
    ----------
    description : array of str
        Description of the data of the time series.
    time_vector : array of datetime objects
    timeformat : how to print the time
    dataseries : List of _dataSeries object holding the
                 data

    Methods:
    --------
    add_dataseries : Add a data series to the object
    write : Write time series to a file
    plot : Plot a figure of a time series

    """

    def __init__(self, desc, timevec, timeformat=None):
        """
        Initalize the object.

        Parameters
        ----------
        desc : array of str
        tvec : array of datetime
        """

        self.description = desc
        self.time_vector = timevec
        self.timeformat = timeformat
        self.dataseries = []

    def add_dataseries(self, label, unit, dataseries):
        """
        Add a new data series to the timeseries object.
        The length of the data vector must be the same as the
        length of the time vector.
        """
        if (len(dataseries) != len(self.time_vector)):
            raise Exception("ERROR: Number of data series sample do "
                            "not correspond to time vector ('%s')" % label)

        ds = _DataSeries(label, unit, dataseries)
        self.dataseries.append(ds)

    def write(self, fname):
        """
        """

        print("----- write to '%s'" % fname)

        try:
            tsfile = open(fname, "w")
        except:
            raise Exception("ERROR: Could not create file '%s'" % fname)

        print("# Weather radar timeseries data file", file=tsfile)
        print("# Project: MALSplus", file=tsfile)
        print("# Start : %s UTC" %
              self.time_vector[0].strftime("%Y-%m-%d %H:%M:%S"),
              file=tsfile)
        print("# End   : %s UTC" %
              self.time_vector[-1].strftime("%Y-%m-%d %H:%M:%S"),
              file=tsfile)
        print("# Header lines with comments are preceded by '#'", file=tsfile)
        for line in self.description:
            print("# %s" % line, file=tsfile)
        print("#", file=tsfile)

        # Make raw header
        if (self.timeformat is None):
            print("# Date, UTC [seconds since midnight]", end="", file=tsfile)
        else:
            print("# Date [%s]" % self.timeformat, end="", file=tsfile)

        for ds in self.dataseries:
            print(", %s [%s]" % (ds.label, ds.unit), end="", file=tsfile)
        print("", file=tsfile)

        # Store the data
        nsample = len(self.time_vector)
        for kk in range(nsample):
            if (self.timeformat is None):
                dt = self.time_vector[kk]
                daystr = dt.strftime("%d-%b-%Y")
                secs = dt.hour*3600. + dt.minute*60. + dt.second + \
                    dt.microsecond/1000000.
                print("%s, %14.4f" % (daystr, secs), end="", file=tsfile)
            else:
                print(self.time_vector[kk].strftime(self.timeformat), end="",
                      file=tsfile)

            for ds in self.dataseries:
                print(", %14.4f" % (ds.data[kk]), end="", file=tsfile)

            print("", file=tsfile)

        tsfile.close()

    def plot(self, fname, label):
        """
        Make a figure of a time series
        """

        found = False
        for ds in self.dataseries:
            if (ds.label == label):
                found = True
                break

        if (not found):
            raise Exception("ERROR: Undefined time series '%s'" % label)

        print("----- plot to '%s'" % fname)

        title = "Trajectory Time Series %s" % \
            self.time_vector[0].strftime("%Y-%m-%d")
        labely = "%s [%s]" % (ds.label, ds.unit)
        plot_timeseries(self.time_vector, ds.data, [fname], titl=title,
                        labely=labely, timeformat="%H:%M")


class _DataSeries(object):
    """
    Hold a data vector and some meta information.
    """

    def __init__(self, label, unit, data):
        """
        Initalize the object.
        """

        self.label = label
        self.unit = unit
        self.data = data
