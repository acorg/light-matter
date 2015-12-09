import six
from six.moves import builtins
from unittest import TestCase

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from json import dumps
import bz2

from light.performance.results import PerformanceResults
from .sample_data import RESULT1, RESULT2
from test.mocking import mockOpen


class BZ2(object):
    """
    A BZ2File mock.
    """
    def __init__(self, data):
        self._data = data
        self._index = 0

    def close(self):
        pass

    def read(self):
        self._index += 1
        return self._data[self._index - 1]

    def __iter__(self):
        index = self._index
        self._index = len(self._data)
        return iter(self._data[index:])


class TestPerformanceResults(TestCase):
    """
    Tests for the PerformanceResults class.
    """
    def testEmptyFileValueError(self):
        """
        If an empty resultFile is given, a C{ValueError} must be raised when
        trying to read it.
        """
        mockOpener = mockOpen()
        with patch.object(builtins, 'open', mockOpener):
            error = "^Result JSON file 'file\\.json' was empty\\.$"
            six.assertRaisesRegex(
                self, ValueError, error, PerformanceResults, ['file.json'])

    def testNonJSONInput(self):
        """
        When given a file whose contents are not JSON, attempting to
        read the light matter results from it must raise a C{ValueError}.
        """
        mockOpener = mockOpen(read_data='not JSON\n')
        with patch.object(builtins, 'open', mockOpener):
            if six.PY3:
                error = ("^Content of file 'file\.json' could not be "
                         "converted to JSON \(Expecting value: line 1 column "
                         "1 \(char 0\)\)\.$")
            else:
                error = ("^Content of file 'file\.json' could not be "
                         "converted to JSON \(No JSON object could be "
                         "decoded\)\.$")
            six.assertRaisesRegex(
                self, ValueError, error, PerformanceResults, ['file.json'])

    def testCorrectJSONDictOneFile(self):
        """
        If one resultFile is given, the correct JSON dictionary must be made.
        """
        test = 'performance.perf_database.TestDatabase.testCreation'
        mockOpener = mockOpen(read_data=dumps(RESULT1) + '\n')
        with patch.object(builtins, 'open', mockOpener):
            performance = PerformanceResults(['file.json'])
            self.assertEqual({'elapsed': 4.3869e-05,
                              'status': 'success',
                              }, performance.tests[0]['results'][test])

    def testCorrectJSONDictOneCompressedFile(self):
        """
        If one compressed resultFile is given, the correct JSON dictionary
        must be made.
        """
        result = BZ2([dumps(RESULT1) + '\n'])
        test = 'performance.perf_database.TestDatabase.testCreation'
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.return_value = result
            performance = PerformanceResults(['file.json.bz2'])
            self.assertEqual({'elapsed': 4.3869e-05,
                              'status': 'success',
                              }, performance.tests[0]['results'][test])

    def testCorrectJSONDictTwoCompressedFiles(self):
        """
        If two compressed result files are given, two results must be present
        on the PerformanceResults instance.
        """
        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename):
                if self.first:
                    self.first = False
                    return BZ2([dumps(RESULT1) + '\n'])
                else:
                    return BZ2([dumps(RESULT2) + '\n'])

        sideEffect = SideEffect()
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            performance = PerformanceResults(['f1.json.bz2', 'f2.json.bz2'])
            self.assertEqual(2, len(performance.tests))

    def testTestNames(self):
        """
        testNames must return all test names correctly.
        """
        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename):
                if self.first:
                    self.first = False
                    return BZ2([dumps(RESULT1) + '\n'])
                else:
                    return BZ2([dumps(RESULT2) + '\n'])

        sideEffect = SideEffect()

        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            performance = PerformanceResults(['f1.json.bz2', 'f2.json.bz2'])
            allTests = list(performance.testNames())
            self.assertEqual(sorted(
                [
                    'performance.perf_database.TestDatabase.testAddSubjects',
                    'performance.perf_database.TestDatabase.testChecksum',
                    'performance.perf_database.TestDatabase.testChecksumEmpty',
                    'performance.perf_database.TestDatabase.testCreation',
                    'performance.perf_database.TestDatabase.testSomething',
                    'performance.perf_findSelf.TestFindSelf.testFindIdentical',
                ]),
                sorted(allTests))

    def testResultsForTestNameNotPresent(self):
        """
        If resultsForTest is asked to return results from a non-existent test,
        C{None} must be returned for both result files.
        """
        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename):
                if self.first:
                    self.first = False
                    return BZ2([dumps(RESULT1) + '\n'])
                else:
                    return BZ2([dumps(RESULT2) + '\n'])

        sideEffect = SideEffect()

        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            performance = PerformanceResults(['f1.json.bz2', 'f2.json.bz2'])
            result = list(performance.resultsForTest('weirdName'))
            self.assertEqual([None, None], result)

    def testResultsForTestSomeTestsMatch(self):
        """
        resultsForTest must return all results correctly, if asked to return
        results from a test that is present in just the second result file.
        """
        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename):
                if self.first:
                    self.first = False
                    return BZ2([dumps(RESULT1) + '\n'])
                else:
                    return BZ2([dumps(RESULT2) + '\n'])

        sideEffect = SideEffect()

        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            performance = PerformanceResults(['f1.json.bz2', 'f2.json.bz2'])
            test = 'performance.perf_database.TestDatabase.testSomething'
            result = list(performance.resultsForTest(test))
            self.assertEqual([None,
                              {'status': 'success',
                               'details': 0.1157001019,
                               'elapsed': 0.1762839317}],
                             result)

    def testResultsForTestAllTestsMatch(self):
        """
        resultsForTest must return all results correctly, if asked to return
        results from a test that is present in all the results.
        """
        class SideEffect(object):
            def __init__(self):
                self.first = True

            def sideEffect(self, _ignoredFilename):
                if self.first:
                    self.first = False
                    return BZ2([dumps(RESULT1) + '\n'])
                else:
                    return BZ2([dumps(RESULT2) + '\n'])

        sideEffect = SideEffect()

        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            performance = PerformanceResults(['f1.json.bz2', 'f2.json.bz2'])
            test = 'performance.perf_database.TestDatabase.testCreation'
            result = list(performance.resultsForTest(test))
            self.assertEqual([{'status': 'success', 'elapsed': 4.3869e-05},
                              {'status': 'success', 'elapsed': 0.100043869}],
                             result)
