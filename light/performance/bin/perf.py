#!/usr/bin/env python

from __future__ import print_function

import sys
import unittest
import argparse
from os.path import basename, dirname, join
from time import gmtime, strftime

from light.parameters import DatabaseParameters, FindParameters
import light.performance
from light.performance.parameters import PARAMETER_SETS
from light.performance.runner import (
    PerformanceTestRunner, filterTestSuite, suiteTestNames)
from light.performance.utils import makeDir

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
        '--listTests', action='store_true', default=False,
        help=('If True, the names of the performance tests that would be run '
              'will be printed to standard output. The tests will not be '
              'run.'))

    parser.add_argument(
        '--parameterSet', action='append', dest='parameterSets',
        choices=sorted(PARAMETER_SETS),
        help=('The name of a canned parameter set to use. This option may be '
              'repeated.'))

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

    # If no parameter sets were specified, use the ones on the command line.
    if not args.parameterSets:
        args.parameterSets = ['command-line']

    if 'command-line' in args.parameterSets:
        PARAMETER_SETS['command-line'] = {
            'dbParams': DatabaseParameters.fromArgs(args),
            'findParams': FindParameters.fromArgs(args),
        }

    outputDir = args.outputDir
    makeDir(outputDir)

    # Make the command-line args available to the tests. There doesn't seem
    # to be another way to pass information like this when using a unittest
    # test suite.
    light.performance.testArgs = args

    print('Writing output to directory %r.' % outputDir, file=sys.stderr)

    suite = unittest.defaultTestLoader.discover(args.startDir,
                                                pattern='perf_*.py')
    if args.testIdPrefix:
        suite = filterTestSuite(args.testIdPrefix, suite)

    if args.listTests:
        print('\n'.join(suiteTestNames(suite)))
        sys.exit(0)

    if suite.countTestCases() == 0:
        print('%s: No test cases match %r.' % (
            basename(sys.argv[0]), args.testIdPrefix), file=sys.stderr)
        sys.exit(1)

    result = PerformanceTestRunner(args).run(suite)

    with open(join(outputDir, _RESULT_FILENAME), 'w') as fp:
        result.save(args, fp=fp)
