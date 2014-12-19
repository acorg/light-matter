from unittest import TestCase

from light.hsp import normalizeHSP


class Template(object):
    def __init__(self, template):
        template = template.split('\n')
        # Allow the first template line to be empty.
        if len(template[0]) == 0:
            template = template[1:]

        # Analyze the template hit.
        self.hit = template[0].rstrip()
        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen) = \
            self._analyze(self.hit)
        origin = spacesLen
        self.matchLen = matchLen
        self.subjectLength = len(self.hit) - spacesLen
        self.hitMatchStart = leadingDotsLen

        # Analyze the template read.
        self.read = template[1].rstrip()
        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen) = \
            self._analyze(self.read)
        assert self.matchLen == matchLen
        self.readLen = len(self.read) - spacesLen
        self.readMatchStart = leadingDotsLen

        # Analyze the template read result.
        self.readResult = template[2].rstrip()
        (spacesLen, leadingDotsLen, matchLen, trailingDotsLen) = \
            self._analyze(self.readResult)
        assert self.matchLen == matchLen
        self.readResultStart = spacesLen - origin
        self.readResultLen = len(self.readResult) - spacesLen
        assert self.readResultLen == self.readLen
        self.readResultMatchStart = leadingDotsLen

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

    def hsp(self):
        """
        Make an HSP.
        """
        return {
            'readOffset': self.readMatchStart,
            'subjectOffset': self.hitMatchStart,
            'landmarkLength': self.matchLen,
        }


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
        self.assertEqual(4, template.hitMatchStart)

        self.assertEqual(20, template.readLen)
        self.assertEqual(2, template.readMatchStart)

        self.assertEqual(20, template.readResultLen)
        self.assertEqual(-11, template.readResultStart)
        self.assertEqual(15, template.readResultMatchStart)

        hsp = template.hsp()
        self.assertEqual(4, hsp['subjectOffset'])
        self.assertEqual(2, hsp['readOffset'])
        self.assertEqual(3, hsp['landmarkLength'])


class TestHSP(TestCase):

    def check(self, templateStr):
        template = Template(templateStr)
        normalized = normalizeHSP(template.hsp(), template.readLen)
        self.assertEqual({
            'subjectStart': template.hitMatchStart,
            'subjectEnd': template.hitMatchStart + template.matchLen,
            'readStart': template.readMatchStart,
            'readEnd': template.readMatchStart + template.matchLen,
            'readStartInSubject': template.readResultStart,
            'readEndInSubject': template.readResultStart + template.readLen,
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
