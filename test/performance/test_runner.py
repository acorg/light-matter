from unittest import TestCase, TestSuite, TestResult

from light.performance.runner import PerformanceTestRunner, filterTestSuite


class TestRunner(TestCase):
    """
    Test the light.performance.PerformanceTestRunner class
    """

    def testVerbosityZero(self):
        """
        The verbosity attribute must be set to zero if args.hidePassFail is
        True.
        """
        class X(object):
            pass
        args = X()
        args.hidePassFail = True
        runner = PerformanceTestRunner(args)
        self.assertEqual(0, runner.verbosity)

    def testVerbosityOne(self):
        """
        The verbosity attribute must be set to one if args.hidePassFail is
        False.
        """
        class X(object):
            pass
        args = X()
        args.hidePassFail = False
        runner = PerformanceTestRunner(args)
        self.assertEqual(1, runner.verbosity)


class TestFilterTestSuite(TestCase):
    """
    Test the light.performance.filterTestSuite function.
    """
    def testFilterNothing(self):
        """
        Filtering an empty test suite must return an empty test suite.
        """
        self.assertEqual(0, filterTestSuite('', TestSuite()).countTestCases())

    def testFilterNoMatches(self):
        """
        Filtering a test suite on a name that doesn't match any tests should
        give a result with no tests.
        """

        class Dummy(TestCase):
            def runTest(self):
                pass

        orig = TestSuite()
        test = Dummy()
        orig.addTest(test)
        self.assertEqual(0, filterTestSuite('xxx', orig).countTestCases())

    def testFilterSomeMatches(self):
        """
        Filtering a test suite on a name that matches just one of two tests
        should give a result that contains just the matching test.
        """

        class Dummy1(TestCase):
            def runTest(self):
                pass

        class Dummy2(TestCase):
            def runTest(self):
                self.fail('Oops')

        orig = TestSuite()
        orig.addTests([Dummy1(), Dummy2()])
        new = filterTestSuite('Dummy1', orig)
        self.assertEqual(1, new.countTestCases())
        # Run the test suite to check that the Dummy2 test was not run.
        result = TestResult()
        new.run(result)
        self.assertTrue(result.wasSuccessful())

    def testFilterAllMatches(self):
        """
        Filtering a test suite on a name that matches both of two tests should
        give a result that contains both tests.
        """

        class Dummy1(TestCase):
            def runTest(self):
                pass

        class Dummy2(TestCase):
            def runTest(self):
                pass

        orig = TestSuite()
        orig.addTests([Dummy1(), Dummy2()])
        new = filterTestSuite('Dummy', orig)
        self.assertEqual(2, new.countTestCases())
