"""
pyrad.util.stat_utils
======================

Miscellaneous functions dealing with statistics

.. autosummary::
    :toctree: generated/

    quantiles_weighted

"""

import numpy as np


def quantiles_weighted(values, weight_vector=None, quantiles=np.array([0.5]),
                       fill_value=None, weight_threshold=None):
    """
    Given a set of values and weights, compute the weighted
    quantile(s).
    """

    if (weight_vector is not None):
        if (weight_vector.size != values.shape[0]):
            raise Exception("ERROR: Unexpected size of weight vector "
                            "(%d instead of %d)" % (weight_vector.size, values.shape[0]))
    else:
        weight_vector = np.ones(values.shape[0], dtype=float)

    if (len(values.shape) > 1):
        # repeat weight vec
        weight_vector = np.repeat(weight_vector, values.shape[1]) \
            .reshape(weight_vector.size, values.shape[1])

        values = values.reshape(-1)
        weight_vector = weight_vector.reshape(-1)

    # remove nans
    if (fill_value is not None):
        values_ma = np.ma.masked_values(values, fill_value)
        if (values_ma.mask.any()):
            values = values_ma[~values_ma.mask]
            weight_vector = weight_vector[~values_ma.mask]

    nvalid = weight_vector.size
    if (nvalid < 3):
        return (None, np.array([None] * quantiles.size), None)

    total_weight = np.sum(weight_vector)

    # Average
    avg = np.sum(np.multiply(values, weight_vector)) / total_weight

    sorter = np.argsort(values, axis=None)
    values = values[sorter]
    weight_vector = weight_vector[sorter]

    if (weight_threshold is not None):
        if (total_weight < weight_threshold):
            return (avg, np.array([None] * quantiles.size), nvalid)

    weighted_quantiles = np.cumsum(weight_vector) - 0.5 * weight_vector

    weighted_quantiles /= total_weight

    # As done by np.percentile():
    # weighted_quantiles -= weighted_quantiles[0]
    # weighted_quantiles /= weighted_quantiles[-1]

    # Note: Does not extrapolate
    quants = np.interp(quantiles, weighted_quantiles, values)

    return (avg, quants, nvalid)
