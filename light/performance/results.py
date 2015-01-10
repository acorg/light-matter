import bz2
from json import loads


class PerformanceResult(object):
    """
    A class to work with results from performance tests.

    @param resultFilenames: Either a single C{str} filename or a C{list} of
        C{str} file names containing our performance test output
        (possibly bzip2 compressed).
    @raise ValueError: If any of the resultFiles are empty or the content can't
        be converted to JSON.
    """
    def __init__(self, resultFilenames):
        if type(resultFilenames) == str:
            resultFilenames = [resultFilenames]
        else:
            self.resultFilenames = resultFilenames

        self.result = []

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
            except ValueError:
                raise ValueError('Content of file %r could not be converted '
                                 'to JSON.' % filename)

            self.result.append(fileResult)

    def showAllTests(self):
        """
        Return a C{list} of names of all tests present in the resultFilenames.
        """
        allTests = set()
        for tests in self.result:
            allTests.update(tests['results'].iterkeys())
        return allTests

    def returnResult(self, testName):
        """
        Return the results of a test in all files.

        @param testName: the C{str} name of the test of which results should be
            returned.
        @raise KeyError: If testName is not present in self.result.
        """
        return (results['results'].get(testName, None) for
                results in self.result)
