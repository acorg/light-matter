import bz2
from json import loads


class PerformanceResults(object):
    """
    A class to work with results from performance tests.

    @param resultFilenames: Either a single C{str} filename or a C{list} of
        C{str} file names containing our performance test output
        (possibly bzip2 compressed).
    @raise ValueError: If any of the result files are empty or the content
        can't be converted to JSON.
    """
    def __init__(self, resultFilenames):
        if type(resultFilenames) == str:
            resultFilenames = [resultFilenames]
        else:
            self.resultFilenames = resultFilenames

        self.tests = []

        for filename in resultFilenames:
            if filename.endswith('.bz2'):
                fp = bz2.BZ2File(filename)
            else:
                fp = open(filename)

            data = fp.read()
            if not data:
                raise ValueError('Result JSON file %r was empty.' % filename)
            try:
                fileResult = loads(data)
            except ValueError as e:
                raise ValueError('Content of file %r could not be converted '
                                 'to JSON (%s).' % (filename, e))
            else:
                self.tests.append(fileResult)

    def testNames(self):
        """
        Get the names of the tests present in the result files.

        Return a C{set} of names of all tests present in any result file.
        """
        allTests = set()
        for test in self.tests:
            allTests.update(test['results'])
        return allTests

    def resultsForTest(self, testName):
        """
        Return the results of a test in all files.

        @param testName: the C{str} name of the test of which results should be
            returned.
        @return: A generator that yields either results for the named test for
            each results file, or C{None} when the test was not present in a
            result file.
        """
        return (test['results'].get(testName) for test in self.tests)
