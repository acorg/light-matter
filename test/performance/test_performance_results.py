from unittest import TestCase
from json import dumps
from mock import patch
import bz2

from light.performance.results import PerformanceResult
from sample_data import RESULT1, RESULT2
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

    def readline(self):
        self._index += 1
        return self._data[self._index - 1]

    def __iter__(self):
        index = self._index
        self._index = len(self._data)
        return iter(self._data[index:])


class TestPerformanceResult(TestCase):
    """
    Tests for the PerformanceResult class.
    """
    def testEmptyFileValueError(self):
        """
        If an empty resultFile is given, a C{ValueError} must be raised when
        trying to read it.
        """
        mockOpener = mockOpen()
        with patch('__builtin__.open', mockOpener, create=True):
            error = "Result JSON file 'file.json' was empty."
            self.assertRaisesRegexp(
                ValueError, error, PerformanceResult, ['file.json'])

    def testNonJSONInput(self):
        """
        When given a file whose contents are not JSON, attempting to
        read the light matter results from it must raise a C{ValueError}.
        """
        mockOpener = mockOpen(read_data='not JSON\n')
        with patch('__builtin__.open', mockOpener, create=True):
            error = ("Content of file 'file.json' could not be converted to "
                     "JSON.")
            self.assertRaisesRegexp(
                ValueError, error, PerformanceResult, ['file.json'])

    def testCorrectJSONDictOneFile(self):
        """
        If one resultFile is given, the correct JSON dictionary must be made.
        """
        test = 'performance.perf_database.TestDatabase.testCreation'
        mockOpener = mockOpen(read_data=dumps(RESULT1) + '\n')
        with patch('__builtin__.open', mockOpener, create=True):
            performance = PerformanceResult(['file.json'])
            self.assertEqual([{1420819327.900306: {
                               'elapsed': 4.3869e-05,
                               'status': 'success',
                               }}], performance.result[test])

    def testCorrectJSONDictOneCompressedFile(self):
        """
        If one compressed resultFile is given, the correct JSON dictionary
        must be made.
        """
        result = BZ2([dumps(RESULT1) + '\n'])
        test = 'performance.perf_database.TestDatabase.testCreation'
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.return_value = result
            performance = PerformanceResult(['file.json.bz2'])
            self.assertEqual([{1420819327.900306: {
                               'elapsed': 4.3869e-05,
                               'status': 'success',
                               }}], performance.result[test])

    def testCorrectJSONDictTwoCompressedFiles(self):
        """
        If two compressed resultFiles are given, the correct JSON dictionary
        must be made.
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
        test = 'performance.perf_database.TestDatabase.testCreation'
        with patch.object(bz2, 'BZ2File') as mockMethod:
            mockMethod.side_effect = sideEffect.sideEffect
            performance = PerformanceResult(['file1.json.bz2',
                                             'file2.json.bz2'])
            self.assertEqual(2, len(performance.result[test]))

    def testShowAllTests(self):
        """
        ShowAllTests must return all test names correctly.
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
            performance = PerformanceResult(['file1.json.bz2',
                                             'file2.json.bz2'])
            allTests = performance.showAllTests()
            self.assertEqual(5, len(allTests))

    def testReturnResultNameNotPresent(self):
        """
        If returnResult is asked to return results from a non-existing test, a
        KeyError must be raised.
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
            performance = PerformanceResult(['file1.json.bz2',
                                             'file2.json.bz2'])
            error = 'weirdName not present in self.result'
            self.assertRaisesRegexp(
                KeyError, error, performance.returnResult, 'weirdName')

    def testReturnResultGoodName(self):
        """
        returnResult must return all results correctly, if asked to return
        results from a test present in the result.
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
            performance = PerformanceResult(['file1.json.bz2',
                                             'file2.json.bz2'])
            test = 'performance.perf_database.TestDatabase.testCreation'
            allResults = performance.returnResult(test)
            self.assertEqual(2, len(allResults))
            self.assertEqual('success',
                             allResults[0][1420819327.9003059864]['status'])
