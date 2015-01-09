import bz2
from json import loads
from collections import defaultdict


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

        self.result = defaultdict(list)

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

            for test, res in fileResult['results'].items():
                self.result[test].append({fileResult['startTestRunTime']: res})
            self.result['resultInfo'] = {fileResult['startTestRunTime']: {
                'UTC': fileResult['UTC'],
                'description': fileResult['description'],
                'elapsed': fileResult['elapsed'],
                'testCount': fileResult['testCount'],
            },
            }

    def showAllTests(self):
        """
        Return a C{list} of names of all tests present in the resultFilenames.
        """
        return [test for test in self.result.keys() if test != 'resultInfo']

    def returnResult(self, testName):
        """
        Return the results of a test in all files.

        @param testName: the C{str} name of the test of which results should be
            returned.
        @raise KeyError: If testName is not present in self.result.
        """
        if not self.result[testName]:
            raise KeyError('%s not present in self.result' % testName)
        else:
            return self.result[testName]
