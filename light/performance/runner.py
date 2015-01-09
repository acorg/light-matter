# The code in this file is adapted from the TextTestRunner and TextTestResult
# classes that are distributed with unittest. The original code can be found in
# unittest/runner.py wherever you have Python installed (to find that location
# do this
#
#   $ python -c 'from unittest import runner; print runner.__file__'
#
# There are no tests for this code. If we want to add them, which we should
# do at some point (especially if we start adding to it), we might want to
# move it under light.performance and put tests in the test directory etc.

import sys
from time import time, gmtime, strftime
from ujson import dump

from unittest import TestResult


class _WritelnDecorator(object):
    """
    Used to decorate file-like objects with a handy 'writeln' method.
    """
    def __init__(self, stream):
        self.stream = stream

    def __getattr__(self, attr):
        if attr in ('stream', '__getstate__'):
            raise AttributeError(attr)
        return getattr(self.stream, attr)

    def writeln(self, arg=None):
        if arg:
            self.write(arg)
        self.write('\n')  # text-mode streams translate to \r\n if needed


class PerformanceTestResult(TestResult):
    """
    A test result class that collects results and their timings and provides a
    method for saving these to a file as JSON.

    @param stream: A file-like object.
    @param verbosity: An C{int} used to control the amount of information
        printed. Unused.
    """

    def __init__(self, stream, verbosity):
        super(PerformanceTestResult, self).__init__(stream, False, verbosity)
        self.stream = stream
        self.showAll = verbosity > 0
        self.status = None
        self.results = {}

    def save(self, fp=sys.stdout, description=None):
        """
        Write test results to C{fp} in JSON format, if no errors or failures
        have occurred.

        @param fp: A file-like object.
        @param description: A C{str} description of the code whose performance
            is being tested.
        """
        if not self.errors and not self.failures:
            dump({
                'UTC': strftime('%F %T', gmtime(self.startTestRunTime)),
                'description': description,
                'results': self.results,
                'startTestRunTime': self.startTestRunTime,
                'elapsed': self.stopTestRunTime - self.startTestRunTime,
                'testCount': self.testsRun,
                }, fp)

    def stopTest(self, test):
        super(PerformanceTestResult, self).stopTest(test)
        try:
            # Do not record details for tests that ask to be ignored.
            if test.ignore:
                return
        except AttributeError:
            pass
        elapsed = time() - self.startTime
        if self.status:
            result = {
                'elapsed': elapsed,
                'status': self.status,
            }
            try:
                details = test.details
            except AttributeError:
                pass
            else:
                result['details'] = details
            self.results[test.id()] = result

    def startTest(self, test):
        super(PerformanceTestResult, self).startTest(test)
        self.startTime = time()
        self.currentTest = test
        if self.showAll:
            self.stream.write(str(test))
            self.stream.write(' ... ')
            self.stream.flush()

    def addSuccess(self, test):
        super(PerformanceTestResult, self).addSuccess(test)
        self.status = 'success'
        if self.showAll:
            self.stream.writeln('ok')

    def addError(self, test, err):
        super(PerformanceTestResult, self).addError(test, err)
        if self.showAll:
            self.stream.writeln('ERROR')

    def addFailure(self, test, err):
        super(PerformanceTestResult, self).addFailure(test, err)
        self.status = 'fail'
        if self.showAll:
            self.stream.writeln('FAIL')

    def addSkip(self, test, reason):
        super(PerformanceTestResult, self).addSkip(test, reason)
        self.status = 'skip'
        if self.showAll:
            self.stream.writeln('skipped {0!r}'.format(reason))

    def addExpectedFailure(self, test, err):
        super(PerformanceTestResult, self).addExpectedFailure(test, err)
        self.status = 'expected fail'
        if self.showAll:
            self.stream.writeln('expected failure')

    def addUnexpectedSuccess(self, test):
        super(PerformanceTestResult, self).addUnexpectedSuccess(test)
        self.status = 'unexpected success'
        if self.showAll:
            self.stream.writeln('unexpected success')

    def printErrors(self):
        if self.showAll:
            self.stream.writeln()
        self.printErrorList('ERROR', self.errors)
        self.printErrorList('FAIL', self.failures)

    def printErrorList(self, flavor, errors):
        separator = '=' * 70
        for test, err in errors:
            self.stream.writeln(separator)
            self.stream.writeln('%s: %s' % (flavor, str(test)))
            self.stream.writeln(separator)
            self.stream.writeln('%s' % err)

    def startTestRun(self):
        self.startTestRunTime = time()

    def stopTestRun(self):
        self.stopTestRunTime = time()


class PerformanceTestRunner(object):
    """
    A performance test runner.

    @param stream: A file-like object.
    @param verbosity: An C{int} used to control the amount of information
        printed by the result class.
    """

    def __init__(self, stream=sys.stderr, verbosity=1):
        self.stream = _WritelnDecorator(stream)
        self.verbosity = verbosity

    def run(self, test):
        """
        Run the given test case or test suite.

        @param test: A C{TestCase} or C{TestSuite} instance.
        @return: An instance of L{PerformanceTestResult}.
        """
        result = PerformanceTestResult(self.stream, self.verbosity)
        result.startTestRun()
        startTime = time()
        try:
            test(result)
        finally:
            result.stopTestRun()
        stopTime = time()
        timeTaken = stopTime - startTime
        result.printErrors()
        run = result.testsRun
        if self.verbosity:
            self.stream.writeln('Ran %d performance test%s in %.3fs' %
                                (run, run != 1 and 's' or '', timeTaken))
            self.stream.writeln()

        try:
            results = map(len, (result.expectedFailures,
                                result.unexpectedSuccesses,
                                result.skipped))
        except AttributeError:
            expectedFails = unexpectedSuccesses = skipped = 0
        else:
            expectedFails, unexpectedSuccesses, skipped = results

        infos = []
        if not result.wasSuccessful():
            self.stream.write('FAILED')
            failed, errored = map(len, (result.failures, result.errors))
            if failed:
                infos.append('failures=%d' % failed)
            if errored:
                infos.append('errors=%d' % errored)
        else:
            if self.verbosity:
                self.stream.write('OK')
        if skipped:
            infos.append('skipped=%d' % skipped)
        if expectedFails:
            infos.append('expected failures=%d' % expectedFails)
        if unexpectedSuccesses:
            infos.append('unexpected successes=%d' % unexpectedSuccesses)
        if self.verbosity:
            if infos:
                self.stream.writeln(' (%s)' % (', '.join(infos),))
            else:
                self.stream.writeln()
        return result
