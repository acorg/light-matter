from unittest import TestCase
from io import StringIO
try:
    from ujson import loads
except ImportError:
    from json import loads

from dark.reads import AARead

from light.distance import scale
from light.features import Landmark, TrigPoint
from light.landmarks import AlphaHelix, BetaStrand
from light.trig import Peaks, Troughs
from light.checksum import Checksum
from light.parameters import Parameters
from light.backend import Backend
from light.reads import ScannedRead
from light.subject import SubjectStore


class TestBackend(TestCase):
    """
    Tests for the light.database.Backend class.
    """
    def testParametersAreStored(self):
        """
        The backend must call its super class so its parameters are stored.
        """
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend()
        be.configure(params)
        self.assertIs(params, be.params)

    def testInitialBackendIsEmpty(self):
        """
        The index must be empty if no reads have been added.
        """
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend()
        be.configure(params)
        self.assertEqual({}, be.d)

    def testHashWithFeatureOnLeft(self):
        """
        The database hash function must return the expected (negative offset)
        hash when the second feature is to the left of the first.
        """
        params = Parameters([], [])
        be = Backend()
        be.configure(params)
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 10)
        distanceMinus10 = str(scale(-10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual('A:B:' + distanceMinus10,
                         be.hash(landmark, trigPoint))

    def testHashWithFeatureOnRight(self):
        """
        The database hash function must return the expected (positive offset)
        hash when the second feature is to the right of the first.
        """
        params = Parameters([], [])
        be = Backend()
        be.configure(params)
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 30)
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual('A:B:' + distance10, be.hash(landmark, trigPoint))

    def testHashWithFeatureOnLeftAndNonDefaultDistanceBase(self):
        """
        The database hash function must return the expected hash when the
        database has a non-default distance base and the second feature is to
        the left of the first.
        """
        params = Parameters([], [], distanceBase=1.5)
        be = Backend()
        be.configure(params)
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 10)
        distanceMinus10 = str(scale(-10, 1.5))
        self.assertEqual('A:B:' + distanceMinus10,
                         be.hash(landmark, trigPoint))

    def testHashWithFeatureOnRightAndNonDefaultDistanceBase(self):
        """
        The database hash function must return the expected hash when the
        database has a non-default distance base and the second feature is to
        the right of the first.
        """
        params = Parameters([], [], distanceBase=1.5)
        be = Backend()
        be.configure(params)
        landmark = Landmark('name', 'A', 20, 0)
        trigPoint = TrigPoint('name', 'B', 30)
        distance10 = str(scale(10, 1.5))
        self.assertEqual('A:B:' + distance10, be.hash(landmark, trigPoint))

    def testHashWithSymbolDetail(self):
        """
        The database hash function must return the expected value when the
        landmark it is passed has a repeat count.
        """
        params = Parameters([], [])
        be = Backend()
        be.configure(params)
        landmark = Landmark('name', 'A', 20, 0, 5)
        trigPoint = TrigPoint('name', 'B', 30)
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual('A5:B:' + distance10, be.hash(landmark, trigPoint))

    def testScan(self):
        """
        The scan method must return a scanned subject.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        scannedSubject = be.scan(subject)
        self.assertIsInstance(scannedSubject, ScannedRead)

    def testGetScannedPairs(self):
        """
        The getSequencePairs method must return pairs of
        (landmark, trigPoints).
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix], [Peaks], distanceBase=1.0)
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        scannedSubject = be.scan(subject)
        pairs = list(be.getScannedPairs(scannedSubject))
        # First pair.
        landmark, trigPoint = pairs[0]
        self.assertEqual(Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL,
                                  0, 9, 2), landmark)
        self.assertEqual(TrigPoint(Peaks.NAME, Peaks.SYMBOL, 10), trigPoint)
        # Second pair.
        landmark, trigPoint = pairs[1]
        self.assertEqual(Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL,
                                  0, 9, 2), landmark)
        self.assertEqual(TrigPoint(Peaks.NAME, Peaks.SYMBOL, 13), trigPoint)
        self.assertEqual(2, len(pairs))

    def testCollectReadHashesWithNoHashes(self):
        """
        The getHashes method must return a dict keyed by (landmark, trigPoints)
        hash with values containing the read offsets. The result should be
        empty if there are no landmarks in the read.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        query = AARead('query', 'AAA')
        scannedQuery = be.scan(query)
        hashCount = be.getHashes(scannedQuery)
        self.assertEqual({}, hashCount)

    def testCollectReadHashesWithOneLandmark(self):
        """
        The getHashes method must return a dict keyed by (landmark, trigPoints)
        hash with values containing the read offsets. The result should be
        empty if there is only one landmark in the read.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        query = AARead('query', 'FRRRFRRRF')
        scannedQuery = be.scan(query)
        hashCount = be.getHashes(scannedQuery)
        self.assertEqual({}, hashCount)

    def testCollectReadHashes(self):
        """
        The getHashes method must return a dict keyed by (landmark, trigPoints)
        hash with values containing the read offsets.
        """
        params = Parameters([AlphaHelix], [Peaks], distanceBase=1.0)
        be = Backend()
        be.configure(params)
        query = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFASAASA')
        scannedQuery = be.scan(query)
        hashCount = be.getHashes(scannedQuery)
        helixAt0 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 0, 9, 2)
        helixAt15 = Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 15, 9, 2)
        peakAt10 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 10)
        peakAt13 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 13)
        peakAt25 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 25)
        peakAt28 = TrigPoint(Peaks.NAME, Peaks.SYMBOL, 28)
        self.assertEqual(
            {
                'A2:P:28': [[helixAt0, peakAt28]],
                'A2:P:25': [[helixAt0, peakAt25]],
                'A2:P:13': [[helixAt0, peakAt13], [helixAt15, peakAt28]],
                'A2:P:10': [[helixAt0, peakAt10], [helixAt15, peakAt25]],
                'A2:P:-5': [[helixAt15, peakAt10]],
                'A2:P:-2': [[helixAt15, peakAt13]],
                'A2:A2:15': [[helixAt0, helixAt15]],
            }, hashCount)

    def testAddSubjectReturnsCorrectResult(self):
        """
        If one subject is added, addSubject must return whether the subject
        already existed, the index ('0' in this case) of the added subject,
        and the backend name.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        subject = AARead('id', 'FRRRFRRRF')
        preExisting, subjectIndex, hashCount = be.addSubject(subject, '0')
        self.assertFalse(preExisting)
        self.assertEqual('0', subjectIndex)
        self.assertEqual(0, hashCount)

    def testAddSameSubjectIncreasesBackendSize(self):
        """
        If an identical subject is added multiple times, the backend size
        does not increase, because the backend subject store detect duplicates.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRF'), '0')
        self.assertEqual(1, be.subjectCount())
        be.addSubject(AARead('id', 'FRRRFRRRF'), '0')
        self.assertEqual(1, be.subjectCount())

    def testOneReadOneLandmark(self):
        """
        If one subject is added but it only has one landmark, nothing is added
        to the backend.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRF'), '0')
        self.assertEqual({}, be.d)

    def testOneReadTwoLandmarks(self):
        """
        If one subject is added and it has two landmarks, one key is added
        to the backend.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        be.addSubject(
            AARead('id', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '0')
        distance23 = str(scale(23, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:A3:' + distance23: {'0': [[0, 9, 23, 13]]},
            },
            be.d)

    def testTwoReadsTwoLandmarksLimitZeroPairsPerLandmark(self):
        """
        If two identical reads are added, both with two landmarks, no keys
        will be added to the dictionary if limitPerLandmark is zero.
        """
        params = Parameters([AlphaHelix], [], limitPerLandmark=0)
        be = Backend()
        be.configure(params)
        be.addSubject(
            AARead('id1', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '0')
        be.addSubject(
            AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '1')
        self.assertEqual({}, be.d)

    def testTwoReadsTwoLandmarksDifferentOffsets(self):
        """
        If two subjects are added, both with two landmarks separated by the
        same distance, only one key is added to the backend and both reads are
        listed in the dictionary values for the key.

        Note that A3:A2:-23 is not added to the backend since that would be
        redundant (it's the same two landmarks, with the same separation,
        just with the sign changed).
        """
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        be.addSubject(
            AARead('id1', 'AFRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '0')
        be.addSubject(
            AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '1')
        distance23 = str(scale(23, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:A3:' + distance23: {
                    '0': [[1, 9, 24, 13]],
                    '1': [[0, 9, 23, 13]],
                },
            },
            be.d)

    def testOneReadOneLandmarkOnePeak(self):
        """
        If one subject is added and it has one landmark and one peak, one pair
        is added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASA'), '0')
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [[0, 9, 10, 1]]},
            },
            be.d)

    def testOneReadOneLandmarkOnePeakDistanceBase(self):
        """
        If a non-default distanceBase of 2.0 is used, the right distance needs
        to be calculated. In this case, the offsets are 10 AA apart, and the
        distanceBase scaling will change that to a 3 (since int(log base 2 10)
        = 3), though we don't test the 3 value explicitly since that may change
        if we ever change the scale function. That's desirable, but we already
        have tests in test_distance.py that will break in that case.
        """
        distanceBase = 2.0
        params = Parameters([AlphaHelix], [Peaks], distanceBase=distanceBase)
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASA'), '0')
        distance10 = str(scale(10, distanceBase))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [[0, 9, 10, 1]]},
            },
            be.d)

    def testOneReadOneLandmarkOnePeakNoTrigFinders(self):
        """
        If one subject is added and it has one landmark and one peak, but no
        trig finders are in use, nothing is added to the backend.
        """
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASA'), '0')
        self.assertEqual({}, be.d)

    def testOneReadOneLandmarkTwoPeaks(self):
        """
        If one subject is added and it has one landmark and two peaks, two
        pairs are added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance13 = str(scale(13, Parameters.DEFAULT_DISTANCE_BASE))
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [[0, 9, 13, 1]]},
                'A2:P:' + distance10: {'0': [[0, 9, 10, 1]]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksLimitOnePairPerLandmark(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        limit of one pair per landmarks is imposed, only one key is added to
        the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], limitPerLandmark=1)
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [[0, 9, 10, 1]]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksSevereMaxDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        severe maximum distance is imposed, no keys are added to
        the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], maxDistance=1)
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        self.assertEqual({}, be.d)

    def testOneReadOneLandmarkTwoPeaksIntermediateMaxDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        maximum distance is imposed that makes one of the peaks too far
        away, only one key is added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], maxDistance=11)
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [[0, 9, 10, 1]]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksLargeMaxDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, and a
        maximum distance is imposed that is greater than the distance to the
        peaks, two keys are added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], maxDistance=15)
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance13 = str(scale(13, Parameters.DEFAULT_DISTANCE_BASE))
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [[0, 9, 13, 1]]},
                'A2:P:' + distance10: {'0': [[0, 9, 10, 1]]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksPermissiveMinDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        permissive minimum distance is imposed, all keys are added to
        the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], minDistance=1)
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance13 = str(scale(13, Parameters.DEFAULT_DISTANCE_BASE))
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [[0, 9, 13, 1]]},
                'A2:P:' + distance10: {'0': [[0, 9, 10, 1]]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksIntermediateMinDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but an
        intermediate minimum distance is imposed, only the key for the pair
        that exceeds the minimum distance is added to the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], minDistance=11)
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        distance13 = str(scale(13, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance13: {'0': [[0, 9, 13, 1]]},
            },
            be.d)

    def testOneReadOneLandmarkTwoPeaksSevereMinDistance(self):
        """
        If one subject is added and it has one landmark and two peaks, but a
        severe minimum distance is imposed, no keys are added to
        the backend.
        """
        params = Parameters([AlphaHelix], [Peaks], minDistance=100)
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        self.assertEqual({}, be.d)

    def testMultipleSubjectOffsets(self):
        """
        If one subject is added and it has one landmark and one peak separated
        by 10 bases and then, later in the subject, the same pair with the
        same separation, one key must be added to the backend and it
        should have two subject offsets.  Note that minDistance and
        maxDistance are used to discard the matches some longer and shorter
        distance pairs that only have one subject offset (i.e., that only
        appear in the subject once).
        """
        seq = 'FRRRFRRRFASA'
        params = Parameters([AlphaHelix], [Peaks], minDistance=5,
                            maxDistance=10)
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', seq + seq), '0')
        distance10 = str(scale(10, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:P:' + distance10: {'0': [[0, 9, 10, 1], [12, 9, 22, 1]]},
            },
            be.d)

    def testTwoReadsTwoLandmarksSameOffsets(self):
        """
        If two identical reads are added, both with two landmarks at the same
        offsets, only one key is added to the backend and both reads are
        listed in the dictionary values for the key.

        Note that A3:A2:-23 is not added to the backend since that would be
        redundant (it's the same two landmarks, with the same separation,
        just with the sign changed).
        """
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        be.addSubject(
            AARead('id1', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '0')
        be.addSubject(
            AARead('id2', 'FRRRFRRRFAAAAAAAAAAAAAAFRRRFRRRFRRRF'), '1')
        distance23 = str(scale(23, Parameters.DEFAULT_DISTANCE_BASE))
        self.assertEqual(
            {
                'A2:A3:' + distance23: {
                    '0': [[0, 9, 23, 13]],
                    '1': [[0, 9, 23, 13]],
                },
            },
            be.d)

    def testFindNoMatch(self):
        """
        A query against an empty backend must produce no results.
        """
        subject = AARead('subject', 'FRRRFRRRFASAASA')
        query = AARead('query', 'FRRR')
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(query)

        self.assertEqual({}, matches)
        self.assertEqual(0, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindOneMatchingHashInOneLocation(self):
        """
        One matching subject with one matching hash (that occurs in one
        location) must be found correctly.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=11)
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(query)

        self.assertEqual(
            {
                '0': [{
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 11),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 11),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                    'queryLandmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                }]
            },
            matches)
        self.assertEqual(1, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindOneMatchingHashButSubjectExcluded(self):
        """
        One matching subject with one matching hash (that occurs in one
        location) must not be returned if a subjectIndices argument that
        excludes it is passed to find.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=11)
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(query,
                                                        subjectIndices=set())

        self.assertEqual({}, matches)
        self.assertEqual(1, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindOneMatchingHashInTwoLocations(self):
        """
        One matching subject with one matching hash (that occurs in two
        locations) must be found correctly.
        """
        subject = AARead('subject', 'AFRRRFRRRFASAASAVVVVVVASAVVVASA')
        query = AARead('query', 'FRRRFRRRFASAASAFRRRFRRRFFRRRFRRRFFRRRFRRRF')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks])
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(query)

        self.assertEqual(
            {
                '0': [{
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 11),
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 10)
                }, {
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 14),
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 1, 9, 2),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 13)
                }]
            }, matches)

        self.assertEqual(14, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindWithNonMatchingHashes(self):
        """
        Non-matching hashes must be found correctly when storeFullAnalysis is
        passed to find() as True.
        """
        subject = AARead('subject', 'F')
        query = AARead('query', 'AFRRRFRRRFASAASAVV')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks])
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(query, True)

        self.assertEqual({}, matches)
        self.assertEqual(2, hashCount)
        self.assertEqual(
            {
                'A2:P:10': [
                    [Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 1, 9, 2),
                     TrigPoint(Peaks.NAME, Peaks.SYMBOL, 11)]
                ],
                'A2:P:13': [
                    [Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 1, 9, 2),
                     TrigPoint(Peaks.NAME, Peaks.SYMBOL, 14)]
                ],
            },
            nonMatchingHashes)

    def testFindWithIdenticalNonMatchingHashes(self):
        """
        Identical on-matching hashes must be found correctly when
        storeFullAnalysis is passed to find() as True.
        """
        subject = AARead('subject', 'F')
        query = AARead('query', 'AFRRRFRRRFASAAAAAAAAAAAFRRRFRRRFASA')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks], maxDistance=10)
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(query, True)

        self.assertEqual({}, matches)
        self.assertEqual(2, hashCount)
        self.assertEqual(
            {
                'A2:P:10': [
                    [Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 1, 9, 2),
                     TrigPoint(Peaks.NAME, Peaks.SYMBOL, 11)],
                    [Landmark(AlphaHelix.NAME, AlphaHelix.SYMBOL, 23, 9, 2),
                     TrigPoint(Peaks.NAME, Peaks.SYMBOL, 33)]
                ],
            },
            nonMatchingHashes)

    def testFindNoneMatchingTooSmallDistance(self):
        """
        No matches should be found if the max distance is too small.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks], maxDistance=1)
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(query)

        self.assertEqual({}, matches)
        self.assertEqual(0, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindNoneMatchingNoTrigPoint(self):
        """
        No matches should be found if there is only one landmark and there are
        no trig point finders.
        """
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [])
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(query)

        self.assertEqual({}, matches)
        self.assertEqual(0, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testFindTwoMatchingInSameSubject(self):
        """
        Two matching hashes in the subject must be found correctly.
        """
        sequence = 'FRRRFRRRFASAASA'
        subject = AARead('subject', sequence)
        query = AARead('query', sequence)
        params = Parameters([AlphaHelix], [Peaks])
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        matches, hashCount, nonMatchingHashes = be.find(query)

        self.assertEqual(
            {
                '0': [{
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 10),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 10),
                },
                    {
                    'queryLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'subjectLandmark': Landmark('AlphaHelix', 'A', 0, 9, 2),
                    'subjectTrigPoint': TrigPoint('Peaks', 'P', 13),
                    'queryTrigPoint': TrigPoint('Peaks', 'P', 13),
                }]
            }, matches)
        self.assertEqual(2, hashCount)
        self.assertEqual({}, nonMatchingHashes)

    def testInitialChecksum(self):
        """
        The backend checksum must be set to the value passed to its
        __init__ method.
        """
        params = Parameters([], [])
        be = Backend()
        be.configure(params, 'backend', 10)
        self.assertEqual(10, be.checksum())

    def testChecksumAfterSubjectAdded(self):
        """
        The backend checksum must be as expected when a subject has been
        added to the backend.
        """
        params = Parameters([], [])
        be = Backend()
        be.configure(params, 'backend', 10)
        sequence = 'AFRRRFRRRFASAASA'
        subject = AARead('id', sequence)
        be.addSubject(subject, '0')

        expected = Checksum(10).update([
            'id',
            sequence,
        ])
        self.assertEqual(expected.value, be.checksum())

    def testRestoreInvalidJSON(self):
        """
        If a backend restore is attempted from a file that does not contain
        valid JSON, a ValueError error must be raised.
        """
        error = '^Expected object or value$'
        self.assertRaisesRegex(ValueError, error, Backend.restore,
                               StringIO('xxx'))

    def testSaveContentIncludesExpectedKeysAndValues(self):
        """
        When a backend saves, its JSON content must include the expected
        keys and values.
        """
        params = Parameters([], [], limitPerLandmark=16, maxDistance=17,
                            minDistance=18, distanceBase=19.0)
        be = Backend()
        be.configure(params, 'backend', 33)
        fp = StringIO()
        be.save(fp)
        fp.seek(0)

        Parameters.restore(fp)
        SubjectStore.restore(fp)
        state = loads(fp.readline()[:-1])

        # Keys
        self.assertEqual(
            set(['checksum', 'd', 'name', '_totalCoveredResidues']),
            set(state.keys()))

        # Values
        self.assertEqual(be.checksum(), state['checksum'])
        self.assertEqual({}, state['d'])
        self.assertEqual('backend', state['name'])
        self.assertEqual(0, state['_totalCoveredResidues'])

    def testSaveRestoreWithNonDefaultParameters(self):
        """
        When asked to save and then restore a backend with non-default
        parameters, a backend with the correct parameters must result.
        """
        params = Parameters([], [], limitPerLandmark=16, maxDistance=17,
                            minDistance=18, distanceBase=19.0)
        be = Backend()
        be.configure(params)
        fp = StringIO()
        be.save(fp)
        fp.seek(0)
        result = be.restore(fp)
        self.assertIs(None, params.compare(result.params))

    def testSaveRestoreNonEmpty(self):
        """
        When asked to save and then restore a non-empty backend, the correct
        backend must result.
        """
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        be = Backend()
        be.configure(params)
        be.addSubject(AARead('id', 'FRRRFRRRFASAASA'), '0')
        fp = StringIO()
        be.save(fp)
        fp.seek(0)
        result = Backend.restore(fp)
        self.assertEqual(be.subjectCount(), result.subjectCount())
        self.assertEqual(be.d, result.d)
        self.assertEqual(be.params.landmarkClasses,
                         result.params.landmarkClasses)
        self.assertEqual(be.params.limitPerLandmark,
                         result.params.limitPerLandmark)
        self.assertEqual(be.params.maxDistance, result.params.maxDistance)
        self.assertEqual(be.params.minDistance, result.params.minDistance)
        self.assertEqual(be.params.trigPointClasses,
                         result.params.trigPointClasses)
        self.assertEqual(be.checksum(), result.checksum())

    def testChecksumAfterSaveRestore(self):
        """
        A backend that has a sequence added to it, which is then saved and
        restored, and then has a second sequence is added to it must have the
        same checksum as a backend that simply has the two sequences added to
        it without interruption.
        """
        seq1 = 'FRRRFRRRFASAASA'
        seq2 = 'MMMMMMMMMFRRRFR'
        params1 = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        be1 = Backend()
        be1.configure(params1, 'name1', 0)
        be1.addSubject(AARead('id1', seq1), '0')
        fp = StringIO()
        be1.save(fp)
        fp.seek(0)
        be1 = Backend.restore(fp)
        be1.addSubject(AARead('id2', seq2), '1')

        params2 = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs])
        be2 = Backend()
        be2.configure(params2, 'name2', 0)
        be2.addSubject(AARead('id1', seq1), '0')
        be2.addSubject(AARead('id2', seq2), '1')

        self.assertEqual(be1.checksum(), be2.checksum())

    def testPrint(self):
        """
        The print_ function should produce the expected output.
        """
        subject = AARead('subject-id', 'FRRRFRRRFASAASA')
        params = Parameters([AlphaHelix, BetaStrand], [Peaks, Troughs],
                            limitPerLandmark=16, maxDistance=10, minDistance=0,
                            distanceBase=1)
        be = Backend()
        be.configure(params)
        be.addSubject(subject, '0')
        expected = (
            'Name: backend\n'
            'Hash count: 3\n'
            'Checksum: 20379718\n'
            'Subjects (with offsets) by hash:\n'
            '  A2:P:10\n'
            '    0 [[0, 9, 10, 1]]\n'
            '  A2:T:4\n'
            '    0 [[0, 9, 4, 1]]\n'
            '  A2:T:8\n'
            '    0 [[0, 9, 8, 1]]\n'
            'Landmark symbol counts:\n'
            '  AlphaHelix (A2): 3\n'
            'Trig point symbol counts:\n'
            '  Peaks (P): 1\n'
            '  Troughs (T): 2')
        self.assertEqual(expected, be.print_())
