#!/usr/bin/env python

import unittest
import argparse

from light.performance.runner import PerformanceTestRunner


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Run light-matter performance tests')

    parser.add_argument(
        '--startDir', default='.',
        help='The directory in which to (recursively) search for tests.')

    parser.add_argument(
        '--verbosity', default=1, type=int,
        help='Verbosity level. Use 0 to suppress all test progress output.')

    parser.add_argument(
        '--description', default='<not given>',
        help='A description of the code being tested.')

    args = parser.parse_args()
    suite = unittest.defaultTestLoader.discover(args.startDir,
                                                pattern='perf_*.py')
    result = PerformanceTestRunner(verbosity=args.verbosity).run(suite)
    result.save(description=args.description)
