"""
pyrad.io.timeseries
===================

TimeSeries class implementation for holding timeseries data.

.. autosummary::
    :toctree: generated/

    TimeSeries

"""


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
