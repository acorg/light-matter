import sys
from operator import attrgetter
try:
    from ujson import dumps, loads
except ImportError:
    from json import dumps, loads

from light.checksum import Checksum
from light.landmarks import findLandmark, DEFAULT_LANDMARK_CLASSES
from light.trig import findTrigPoint, DEFAULT_TRIG_CLASSES
from light.score import MinHashesScore
from light.significance import HashFraction


class Parameters(object):
    """
    Hold a collection of database parameter settings.

    @param landmarkClasses: A C{list} of landmark finder classes, or C{None}
        to use the default set of landmark finder classes.
    @param trigPointClasses: A C{list} of trig point finder classes, or C{None}
        to use the default set of trig point finder classes.
    @param limitPerLandmark: An C{int} limit on the number of pairs to
        yield per landmark.
    @param maxDistance: The C{int} maximum distance permitted between
        yielded pairs.
    @param minDistance: The C{int} minimum distance permitted between
        yielded pairs.
    @param distanceBase: The distance between a landmark and a trig point is
        scaled to be its logarithm using this C{float} base. This reduces
        sensitivity to relatively small differences in distance.
    """

    PARAMS = ('landmarkClasses', 'trigPointClasses', 'limitPerLandmark',
              'maxDistance', 'minDistance', 'distanceBase')

    # Database construction and look-up defaults. See explanations in
    # docstring above.
    DEFAULT_LIMIT_PER_LANDMARK = 10
    DEFAULT_MAX_DISTANCE = 200
    DEFAULT_MIN_DISTANCE = 1
    DEFAULT_DISTANCE_BASE = 1.1

    # The default fraction of all (landmark, trig point) pairs for a
    # scannedRead that need to fall into the same offset delta histogram
    # bucket for that bucket to be considered a significant match with a
    # database title.
    DEFAULT_SIGNIFICANCE_FRACTION = 0.25
    DEFAULT_SCORE_METHOD = MinHashesScore.__name__
    DEFAULT_SIGNIFICANCE_METHOD = HashFraction.__name__

    def __init__(self, landmarkClasses, trigPointClasses,
                 limitPerLandmark=None, maxDistance=None, minDistance=None,
                 distanceBase=None):

        self.distanceBase = (
            self.DEFAULT_DISTANCE_BASE if distanceBase is None
            else distanceBase)

        if self.distanceBase <= 0:
            raise ValueError('distanceBase must be > 0.')

        self.landmarkClasses = (
            DEFAULT_LANDMARK_CLASSES if landmarkClasses is None
            else landmarkClasses)

        self.landmarkFinders = []
        for landmarkClass in self.landmarkClasses:
            self.landmarkFinders.append(landmarkClass(self.distanceBase))

        self.trigPointClasses = (
            DEFAULT_TRIG_CLASSES if trigPointClasses is None
            else trigPointClasses)

        self.trigPointFinders = []
        for trigPointClass in self.trigPointClasses:
            self.trigPointFinders.append(trigPointClass(self.distanceBase))

        self.limitPerLandmark = (
            self.DEFAULT_LIMIT_PER_LANDMARK if limitPerLandmark is None
            else limitPerLandmark)

        self.maxDistance = (
            self.DEFAULT_MAX_DISTANCE if maxDistance is None
            else maxDistance)

        self.minDistance = (
            self.DEFAULT_MIN_DISTANCE if minDistance is None
            else minDistance)

    @property
    def checksum(self):
        """
        Calculate a checksum based on the database finders (their names
        and symbols) and other parameters.
        """
        # Add landmark and trig point finders in sorted order so databases
        # with the same finders will have identical parameter checksums
        # (all else being equal).
        checksum = Checksum()
        key = attrgetter('NAME')
        landmarkFinders = sorted(self.landmarkFinders, key=key)
        trigPointFinders = sorted(self.trigPointFinders, key=key)
        checksum.update(
            [f.NAME for f in landmarkFinders] +
            [f.SYMBOL for f in landmarkFinders] +
            [f.NAME for f in trigPointFinders] +
            [f.SYMBOL for f in trigPointFinders] +
            list(map(str, (self.limitPerLandmark, self.maxDistance,
                           self.minDistance, self.distanceBase))))
        return checksum.checksum

    def save(self, fp=sys.stdout):
        """
        Save the parameters in JSON format.

        @param fp: A file pointer.
        @return: The C{fp} we were passed (this is useful in testing).
        """
        print(dumps({
            'landmarkClassNames': [cls.NAME for cls in self.landmarkClasses],
            'trigPointClassNames': [cls.NAME for cls in self.trigPointClasses],
            'limitPerLandmark': self.limitPerLandmark,
            'maxDistance': self.maxDistance,
            'minDistance': self.minDistance,
            'distanceBase': self.distanceBase,
        }), file=fp)

        return fp

    @staticmethod
    def restore(fp=sys.stdin):
        """
        Load parameters from a file.

        @param fp: A file pointer.
        @return: An instance of L{Parameters}.
        @raises ValueError: If a now non-existent landmark or trig point name
            is found in the saved parameters file. Or if valid JSON cannot be
            loaded from C{fp}.
        """
        state = loads(fp.readline()[:-1])

        landmarkClasses = []
        for landmarkClassName in state['landmarkClassNames']:
            cls = findLandmark(landmarkClassName)
            if cls:
                landmarkClasses.append(cls)
            else:
                raise ValueError(
                    'Could not find landscape finder class %s. Has that '
                    'class been renamed or removed?' % landmarkClassName)

        trigPointClasses = []
        for trigPointClassName in state['trigPointClassNames']:
            cls = findTrigPoint(trigPointClassName)
            if cls:
                trigPointClasses.append(cls)
            else:
                raise ValueError(
                    'Could not find trig point finder class %s. Has that '
                    'class been renamed or removed?' % trigPointClassName)

        return Parameters(landmarkClasses, trigPointClasses,
                          limitPerLandmark=state['limitPerLandmark'],
                          maxDistance=state['maxDistance'],
                          minDistance=state['minDistance'],
                          distanceBase=state['distanceBase'])

    def print_(self, fp=sys.stdout):
        """
        Print parameter values.
        """
        print('Parameters:', file=fp)
        if self.landmarkFinders:
            print('  Landmark finders:', file=fp)
            print('    ' + '\n    '.join(sorted(
                finder.NAME for finder in self.landmarkFinders)),
                file=fp)
        else:
            print('  Landmark finders: none', file=fp)

        if self.trigPointFinders:
            print('  Trig point finders:', file=fp)
            print('    ' + '\n    '.join(sorted(
                finder.NAME for finder in self.trigPointFinders)),
                file=fp)
        else:
            print('  Trig point finders: none', file=fp)

        print('  Limit per landmark:', self.limitPerLandmark, file=fp)
        print('  Max distance:', self.maxDistance, file=fp)
        print('  Min distance:', self.minDistance, file=fp)
        print('  Distance base:', self.distanceBase, file=fp)

    def compare(self, other):
        """
        Compare our parameters against another Parameters instance.

        @param other: A C{Parameters} instance.
        @return: A C{str} summary of the parameter differences if any, else
            C{None}.
        """
        err = []
        for param in self.PARAMS:
            ours = getattr(self, param)
            theirs = getattr(other, param)
            if ours != theirs:
                err.append('\tParam %r values %r and %r differ.' %
                           (param, ours, theirs))

        return 'Summary of differences:\n%s' % '\n'.join(err) if err else None
