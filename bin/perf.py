#!/usr/bin/env python

import sys
import unittest
import argparse
from os.path import basename

from light.performance.runner import PerformanceTestRunner, filterTestSuite


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

    parser.add_argument(
        '--testIdPattern', default=None,
        help='A test id prefix. Tests whose ids do not contain this pattern '
        'will not be run. The pattern is case-sensitive.')

    args = parser.parse_args()
    suite = unittest.defaultTestLoader.discover(args.startDir,
                                                pattern='perf_*.py')
    if args.testIdPattern:
        suite = filterTestSuite(args.testIdPattern, suite)
        if suite.countTestCases() == 0:
            print >>sys.stderr, '%s: No test cases match %r.' % (
                basename(sys.argv[0]), args.testIdPattern)
            sys.exit(1)

    result = PerformanceTestRunner(verbosity=args.verbosity).run(suite)
    result.save(description=args.description)
