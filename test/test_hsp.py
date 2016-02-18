from unittest import TestCase

from light.hsp import normalizeBin
from light.features import Landmark, TrigPoint


class Template(object):
    def __init__(self, template):
        template = template.split('\n')
        # Allow the first template line to be empty.
        if len(template[0]) == 0:
            template = template[1:]

        # Analyze the template hit (subject).
        self.hit = template[0].rstrip()
        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen) = \
            self._analyze(self.hit)
        origin = spacesLen
        self.matchLen = matchLen
        self.subjectLength = len(self.hit) - spacesLen
        self.subjectLandmarkOffset = leadingDotsLen
        self.subjectTrigPointOffset = leadingDotsLen + matchLen

        # Analyze the template read (query).
        self.read = template[1].rstrip()
        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen) = \
            self._analyze(self.read)
        assert self.matchLen == matchLen
        self.matchLen = matchLen
        self.queryLength = len(self.read) - spacesLen
        self.queryLandmarkOffset = leadingDotsLen
        self.queryTrigPointOffset = leadingDotsLen + matchLen

        # Analyze the template read result.
        self.readResult = template[2].rstrip()
        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen) = \
            self._analyze(self.readResult)
        assert self.matchLen == matchLen
        self.queryResultStart = spacesLen - origin
        self.queryResultLen = len(self.readResult) - spacesLen
        assert self.queryResultLen == self.queryLength
        self.queryResultLandmarkOffset = leadingDotsLen
        self.queryResultTrigPointOffset = leadingDotsLen + matchLen

    def _leadingCharMatchLen(self, str, chars=' '):
        return len(str) - len(str.lstrip(chars))

    def _analyze(self, str):
        offset = spacesLen = self._leadingCharMatchLen(str)
        leadingDotsLen = self._leadingCharMatchLen(str[offset:], '.')
        offset += leadingDotsLen
        assert str[offset] == '+', 'Oops: %r should be a "+"' % str[offset]
        matchLen = self._leadingCharMatchLen(str[offset:], '+')
        offset += matchLen
        trailingDotsLen = self._leadingCharMatchLen(str[offset:], '.')
        return (spacesLen, leadingDotsLen, matchLen, trailingDotsLen)

    def bin_(self):
        """
        Make a bin. To make life easy, assume that the bin only has one pair.
        """
        queryLandmark = Landmark('AlphaHelix', 'A', self.queryLandmarkOffset,
                                 1)
        queryTrigPoint = TrigPoint('Peaks', 'P', self.queryTrigPointOffset)
        subjectLandmark = Landmark('AlphaHelix', 'A',
                                   self.subjectLandmarkOffset, 1)
        subjectTrigPoint = TrigPoint('Peaks', 'P', self.subjectTrigPointOffset)
        return [{
            'queryLandmark': queryLandmark,
            'queryTrigPoint': queryTrigPoint,
            'subjectLandmark': subjectLandmark,
            'subjectTrigPoint': subjectTrigPoint,
        }]


class TestTemplate(TestCase):
    """
    Tests for our helper Template class.
    """

    def testIt(self):
        template = Template('''
                                             ....+++...........
                                               ..+++...............
                                  ...............+++..
        ''')
        self.assertEqual(3, template.matchLen)

        self.assertEqual(18, template.subjectLength)
        self.assertEqual(4, template.subjectLandmarkOffset)

        self.assertEqual(20, template.queryLength)
        self.assertEqual(2, template.queryLandmarkOffset)

        self.assertEqual(20, template.queryResultLen)
        self.assertEqual(-11, template.queryResultStart)
        self.assertEqual(15, template.queryResultLandmarkOffset)

        bin_ = template.bin_()
        self.assertEqual(4, bin_[0]['subjectLandmark'].offset)
        self.assertEqual(2, bin_[0]['queryLandmark'].offset)
        self.assertEqual(3, (bin_[0]['queryTrigPoint'].offset -
                         bin_[0]['queryLandmark'].offset))


class TestBin(TestCase):

    def check(self, templateStr):
        template = Template(templateStr)
        normalized = normalizeBin(template.bin_(), template.queryLength)
        self.assertEqual({
            'subjectBinStart': template.subjectLandmarkOffset,
            'subjectBinEnd': (template.subjectLandmarkOffset +
                              template.matchLen),
            'queryBinStart': template.queryLandmarkOffset,
            'queryBinEnd': template.queryLandmarkOffset + template.matchLen,
            'queryStartInSubject': template.queryResultStart,
            'queryEndInSubject': (template.queryResultStart +
                                  template.queryLength),
        }, normalized)

    def testIdentical1(self):
        self.check('''
                                    +
                                    +
                                    +
                  ''')

    def testIdentical2(self):
        self.check('''
                                    ....++++...
                                    ....++++...
                                    ....++++...
                  ''')

    def testHitExtendsLeft1(self):
        self.check('''
                                    ....++
                                        ++
                                        ++
                  ''')

    def testHitExtendsLeft2(self):
        self.check('''
                                    ....+++...
                                        +++...
                                        +++...
                  ''')

    def testHitExtendsRight1(self):
        self.check('''
                                        ++......
                                        ++...
                                        ++...
                  ''')

    def testHitExtendsRight2(self):
        self.check('''
                                        ..++++......
                                        ..++++...
                                        ..++++...
                  ''')

    def testHitExtendsBoth(self):
        self.check('''
                                    ....+++...........
                                      ..+++....
                                      ..+++....
                  ''')

    def testReadExtendsLeft1(self):
        self.check('''
                                        ++
                                    ....++
                                    ....++
                  ''')

    def testReadExtendsLeft2(self):
        self.check('''
                                        +++...
                                    ....+++...
                                    ....+++...
                  ''')

    def testReadExtendsRight1(self):
        self.check('''
                                        ++...
                                        ++......
                                        ++......
                  ''')

    def testReadExtendsRight2(self):
        self.check('''
                                        ..++++...
                                        ..++++......
                                        ..++++......
                  ''')

    def testReadExtendsBoth(self):
        self.check('''
                                      ..+++....
                                    ....+++...........
                                    ....+++...........
                  ''')

    def testHitExtendsLeftReadExtendsRight(self):
        self.check('''
                                    ....+++...........
                                      ..+++...............
                                      ..+++...............
                  ''')

    def testHitExtendsRightReadExtendsLeft(self):
        self.check('''
                                      ..+++...............
                                    ....+++...........
                                    ....+++...........
                  ''')
