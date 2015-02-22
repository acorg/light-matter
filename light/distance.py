from math import log


def scale(dist, base):
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
    @return: An C{int} scaled distance.
    """
    if dist > 0:
        return int(log(dist) / log(base)) if base != 1.0 else dist
    elif dist < 0:
        return -1 * int(log(-dist) / log(base)) if base != 1.0 else dist
    else:
        return 0
