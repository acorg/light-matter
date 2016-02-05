from __future__ import division

from math import log, floor
from warnings import warn


def _pp_scaleLog(dist, base):
    """
    Use a log scale to adjust a distance.

    By logging distances (usually with a base like 1.1) we can map a range
    of distances onto a smaller range. So some different input distances
    will have the same scaled valued. E.g., with base 1.1, the distance
    range 0-250 is mapped onto the range 0-57.

    The mapping is not linear. Many larger distances will be mapped to the
    same scaled distance. E.g., with base 1.1, 200-207 all are mapped to
    55, 208-228 are all mapped to 56, and 229-250 are all mapped to 57. So
    at the upper end of the distance range 0-250 the 1.1 logarithm base
    allows for about 10% sloppiness in adjusted distances. E.g., distance
    230 is adjusted to the same value (57) as distance 245.  On the other
    hand, values 0-12 are all mapped to distinct values, which means
    distances are strictly interpreted among those small values.

    @param dist: An C{int} distance, usually a distance between two
        features in an AA sequence.
    @param base: The C{float} logarithmic base to use in scaling.
    @raise ValueError: If C{base} is zero.
    @return: An C{int} scaled distance.
    """
    if dist > 0:
        result = int(log(dist) / log(base)) if base != 1.0 else dist
        return min([result, dist])
    elif dist < 0:
        result = -1 * int(log(-dist) / log(base)) if base != 1.0 else dist
        return max([result, dist])
    else:
        return 0

try:
    from light._distance import lib
except ImportError as e:
    warn('Could not import cffi scale function (%s). Using pure Python '
         'version.' % e)
    # Use the pure Python version.
    scaleLog = _pp_scaleLog
else:
    # Use the C version (see ../src/distance.py).
    #
    # Note that this function will return 0.0 if base is 0.0. This is
    # undesirable but fixing it would require writing a Python function to
    # test for it happening and raising ValueError if so. That partly
    # defeats the point of using a C extension. Seeing as it's unlikely
    # that we'll ever use a distance base of zero (and not get tripped up
    # by the test for that in parameters.py), and that the results of doing
    # so would be so spectacularly weird (all distances would be scaled to
    # zero), I thought it better to just call the C function directly. See
    # also the tests in ../test/test_distance.py
    scaleLog = lib.scale


def scaleLinear(dist, divisor):
    """
    Linearly scale a distance.

    @param dist: An C{int} distance, usually a delta.
    @param divisor: The C{int} number by which dist is divided.
    @raise ZeroDivisionError: If C{divisor} is zero.
    @return: An C{int} scaled distance.
    """
    return int(floor(dist / divisor))
