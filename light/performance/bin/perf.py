#!/usr/bin/env python

from __future__ import print_function

import sys
import unittest
import argparse
from os.path import basename, dirname, exists, isdir, join
from time import gmtime, strftime
from os import mkdir

from light.parameters import DatabaseParameters, FindParameters
import light.performance
from light.performance.runner import PerformanceTestRunner, filterTestSuite

# The overall JSON result will be written to a file with this name, in the
# directory given by --outputDir on the command line (or its default value,
# based on the current date/time, if not specified).
_RESULT_FILENAME = 'result.json'

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Run light-matter performance tests')

    parser.add_argument(
        '--interactive', action='store_true', default=False,
        help=('If True, performance tests will display interactive windows of '
              'results as they are computed. If False, just compute results '
              'and store them.'))

    parser.add_argument(
        '--startDir',
        default=join(dirname(light.performance.__file__), 'test'),
        help='The directory in which to (recursively) search for tests.')

    parser.add_argument(
        '--outputDir',
        default='/tmp/lm-performance-' + strftime('%Y%m%d-%H:%M:%S', gmtime()),
        help=('The directory where test output (if any) should be written. '
              'This directory will be created if it doesn\'t exist.'))

    parser.add_argument(
        '--hidePassFail', default=False, action='store_true',
        help=('If given, suppress all test progress output. The overall '
              'performance suite result is always written to standard output, '
              'this option just suppresses the intermediate output showing '
              'the pass/fail status of individual tests.'))

    parser.add_argument(
        '--description', default='<not given>',
        help='A description of the code being tested.')

    parser.add_argument(
        '--testIdPrefix', default=None,
        help=('A test id prefix. Tests whose ids do not contain this pattern '
              'will not be run. The pattern is case-sensitive.'))

    FindParameters.addArgsToParser(parser)
    DatabaseParameters.addArgsToParser(parser)

    args = parser.parse_args()

    # Add the find and database parameters to args for the convenience of
    # the tests that may want them. Sanity check that those attributes
    # don't already exist on args.
    assert not (hasattr(args, 'findParams') or hasattr(args, 'dbParams'))
    args.findParams = FindParameters.fromArgs(args)
    args.dbParams = DatabaseParameters.fromArgs(args)

    outputDir = args.outputDir
    if exists(outputDir):
        if not isdir(outputDir):
            print('The specified output directory %r already exists, but it '
                  'is a file!' % outputDir)
    else:
        mkdir(outputDir)

    # Make the command-line args available to the tests. There doesn't seem
    # to be another way to pass information like this when using a unittest
    # test suite.
    light.performance.testArgs = args

    print('Writing output to directory %r.' % outputDir, file=sys.stderr)

    suite = unittest.defaultTestLoader.discover(args.startDir,
                                                pattern='perf_*.py')
    if args.testIdPrefix:
        suite = filterTestSuite(args.testIdPrefix, suite)
        if suite.countTestCases() == 0:
            print('%s: No test cases match %r.' % (
                basename(sys.argv[0]), args.testIdPrefix), file=sys.stderr)
            sys.exit(1)

    result = PerformanceTestRunner(args).run(suite)

    with open(join(outputDir, _RESULT_FILENAME), 'w') as fp:
        result.save(args, fp=fp)
