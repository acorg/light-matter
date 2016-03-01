from math import sqrt
import numpy as np

from light.string import MultilineString


class Stats(list):
    """
    Compute basic statistics on a set of values.

    @param description: A C{str} description of the variable whose values are
        being analyzed. Used for printing a summary string.
    """

    # Calculating variance is actually non-trivial. See e.g.,
    # https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance for
    # some discussion of online algorithms. For more detail see e.g.,
    # "Algorithms for Computing the Sample Variance: Analysis and
    # Recommendations" by Chan, Golub and LeVeque in The American
    # Statistician, Vol. 37, No. 3, August 1983.
    #
    # This implementation uses a non-online / two-pass approach. That has
    # the (potentially major) downside of the requirement of storing all
    # values in memory. But it has the upside of allowing us to just use
    # numpy's algorithms for the sake of simplicity and of also making it
    # possible to compute medians.
    #
    # Note that we are not passing a value for 'ddof' to np.var, so we are
    # calculating the exact population variance (as opposed to an estimate
    # of the variance based on a sample, which would require passing
    # ddof=1). For details, see
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.var.html

    def __init__(self, description='Summary'):
        self._description = description

    def summary(self):
        """
        Compute summary statistics on the data values.

        @raise ValueError: If no values have been added.
        @return: A C{dict} with keys (see below) for the computed values.
        """
        if self:
            values = np.array(self)
            variance = np.var(values)
            return {
                'count': len(values),
                'max': np.max(values),
                'mean': np.mean(values),
                'median': np.median(values),
                'min': np.min(values),
                'sd': sqrt(variance),
                'sum': np.sum(values),
                'variance': variance,
            }
        else:
            raise ValueError('No values have been added to Stats instance.')

    def __str__(self):
        """
        Summarize the variable's statistics.

        @return: A C{str} summary of the variable's statistics.
        """
        summary = self.summary()
        result = MultilineString()

        result.append('%s:' % self._description)
        result.indent()

        result.extend([
            'Count: %d' % summary['count'],
            'Max: %s' % summary['max'],
            'Mean: %.4f' % summary['mean'],
            'Median: %.4f' % summary['median'],
            'Min: %s' % summary['min'],
            'SD: %.4f' % summary['sd'],
            'Variance: %.4f' % summary['variance'],
        ])

        return str(result)
