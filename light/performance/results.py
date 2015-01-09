import bz2
from json import loads
from collections import defaultdict


class PerformanceResult(object):
    """
    A class to work with results from performance tests.

    @param resultFilenames: Either a single C{str} filename or a C{list} of
        C{str} file names containing our performance test output
        (possibly bzip2 compressed).
    """

    def __init__(self, resultFilenames):
        if type(resultFilenames) == str:
            resultFilenames = [resultFilenames]
        else:
            self.resultFilenames = resultFilenames

        self.result = defaultdict(list)

        for filename in resultFilenames:
            if filename.endswith('.bz2'):
                self._fp = bz2.BZ2File(filename)
            else:
                self._fp = open(filename)

            line = self._fp.readline()
            if not line:
                raise ValueError('Result JSON file %r was empty.' % filename)
            try:
                fileResult = loads(line)
            except ValueError:
                raise ValueError('Content of file %r could not be converted '
                                 'to JSON.' % filename)

            for test, res in fileResult['results'].items():
                self.result[test].append({fileResult['startTestRunTime']: res})

    def showAllTests(self):
        """
        Return a C{list} of names of all tests present in the resultFilenames.
        """
        tests = []
        for test in self.result:
            tests.append(test)
        return tests

    def returnResult(self, testName):
        """
        Return the results of a test in all files.

        @param testName: the C{str} name of the test of which results should be
            returned.
        """
        try:
            return self.result[testName]
        except ValueError:
            raise ValueError('No test with name: %s present in result.' % (
                             testName))
