from __future__ import print_function

import sys
from operator import attrgetter
try:
    from ujson import dumps, loads
except ImportError:
    from json import dumps, loads

from light.checksum import Checksum
from light.landmarks import (
    findLandmark, DEFAULT_LANDMARK_CLASSES, ALL_LANDMARK_CLASSES)
from light.trig import findTrigPoint, DEFAULT_TRIG_CLASSES, ALL_TRIG_CLASSES
from light.bin_score import MinHashesScore, ALL_BIN_SCORE_CLASSES
from light.significance import HashFraction, ALL_SIGNIFICANCE_CLASSES
from light.overall_score import BestBinScore, ALL_OVERALL_SCORE_CLASSES
from light.string import MultilineString


def parseWeights(weights):
    """
    Parse feature weights specification.

    @param weights: A C{list} of C{str}s. Each item is of the form
        'landmark 2'.
    @return: A C{dict} whose keys are landmark or trigPoint names and whose
        values are weights.
    """
    allWeights = FindParameters.DEFAULT_WEIGHTS.copy()

    for featureWeight in weights:
        feature, weight = featureWeight.split()
        allWeights[feature] = float(weight)

    return allWeights


class FindParameters(object):
    """
    Hold a collection of parameter settings used by the database find command.

    @param significanceMethod: The C{str} name of the method used to calculate
        which histogram bins are considered significant.
    @param significanceFraction: The C{float} fraction of all (landmark, trig
        point) pairs for a scannedRead that need to fall into the same
        histogram bucket for that bucket to be considered a significant match
        with a database title.
    @param scoreMethod: The C{str} name of the method used to calculate the
        score of a histogram bin which is considered significant.
    @param featureMatchScore: The C{float} contribution (usually positive) to
        a score when a feature in a query and subject are part of a match.
    @param featureMismatchScore: The C{float} contribution (usually negative)
        to a score when a feature in a query and subject are part of a match.
    @param deltaScale: A C{float}. The delta between the queryOffset and the
        subjectOffset of a matching pair is scaled by dividing the delta by the
        deltaScale to reduce sensitivity.
    """
    # The methods to be used to calculate match scores and whether matches
    # are significant (i.e., worth reporting).
    DEFAULT_SCORE_METHOD = MinHashesScore.__name__
    DEFAULT_SIGNIFICANCE_METHOD = HashFraction.__name__
    DEFAULT_OVERALL_SCORE_METHOD = BestBinScore.__name__

    # The default fraction of all (landmark, trig point) pairs for a
    # scannedRead that need to fall into the same offset delta histogram
    # bucket for that bucket to be considered a significant match with a
    # database title.
    DEFAULT_SIGNIFICANCE_FRACTION = 0.25

    # Reward and penalty score for feature match / mismatch in
    # FeatureMatchingScore scoring (see light/score.py).
    DEFAULT_FEATURE_MATCH_SCORE = 1.0
    DEFAULT_FEATURE_MISMATCH_SCORE = -1.0

    DEFAULT_DELTA_SCALE = 1.0

    DEFAULT_WEIGHTS = {
        'AlphaHelix': 1.0,
        'AlphaHelix_3_10': 1.0,
        'AlphaHelix_pi': 1.0,
        'BetaStrand': 1.0,
        'BetaTurn': 1.0,
        'AminoAcidsLm': 1.0,
        'GOR4AlphaHelix': 1.0,
        'GOR4BetaStrand': 1.0,
        'GOR4Coil': 1.0,
        'Prosite': 1.0,
        'Peaks': 1.0,
        'Troughs': 1.0,
        'AminoAcids': 1.0,
        'IndividualPeaks': 1.0,
        'IndividualTroughs': 1.0,
    }

    def __init__(self, significanceMethod=None, significanceFraction=None,
                 scoreMethod=None, overallScoreMethod=None,
                 featureMatchScore=None, featureMismatchScore=None,
                 weights=None, deltaScale=None):
        self.significanceMethod = (
            self.DEFAULT_SIGNIFICANCE_METHOD if significanceMethod is None
            else significanceMethod)

        self.significanceFraction = (
            self.DEFAULT_SIGNIFICANCE_FRACTION if significanceFraction is None
            else significanceFraction)

        self.scoreMethod = (
            self.DEFAULT_SCORE_METHOD if scoreMethod is None else scoreMethod)

        self.featureMatchScore = (
            self.DEFAULT_FEATURE_MATCH_SCORE if featureMatchScore is None
            else featureMatchScore)

        self.featureMismatchScore = (
            self.DEFAULT_FEATURE_MISMATCH_SCORE if featureMismatchScore is None
            else featureMismatchScore)

        self.weights = self.DEFAULT_WEIGHTS if weights is None else weights

        self.overallScoreMethod = (
            self.DEFAULT_OVERALL_SCORE_METHOD if overallScoreMethod is
            None else overallScoreMethod)

        self.deltaScale = (
            self.DEFAULT_DELTA_SCALE if deltaScale is None else deltaScale)

        if self.deltaScale <= 0:
            raise ValueError('deltaScale must be > 0.')

    @staticmethod
    def addArgsToParser(parser):
        """
        Add arguments for doing a database find to an argparse parser.

        @param parser: An C{argparse.ArgumentParser} instance.
        """
        parser.add_argument(
            '--significanceMethod',
            default=FindParameters.DEFAULT_SIGNIFICANCE_METHOD,
            choices=[cls.__name__ for cls in ALL_SIGNIFICANCE_CLASSES],
            help=('The name of the method used to calculate which histogram '
                  'bins are considered significant.'))

        parser.add_argument(
            '--significanceFraction', type=float,
            default=FindParameters.DEFAULT_SIGNIFICANCE_FRACTION,
            help=('The (float) fraction of all (landmark, trig point) pairs '
                  'for a scannedRead that need to fall into the same '
                  'histogram bucket for that bucket to be considered a '
                  'significant match with a database title.'))

        parser.add_argument(
            '--scoreMethod', default=FindParameters.DEFAULT_SCORE_METHOD,
            choices=[cls.__name__ for cls in ALL_BIN_SCORE_CLASSES],
            help=('The name of the method used to calculate the score of a '
                  'histogram bin which is considered significant.'))

        parser.add_argument(
            '--featureMatchScore', type=float,
            default=FindParameters.DEFAULT_FEATURE_MATCH_SCORE,
            help=('The contribution (usually positive) to a score when a '
                  'feature in a query and subject are part of a match.'))

        parser.add_argument(
            '--featureMismatchScore', type=float,
            default=FindParameters.DEFAULT_FEATURE_MISMATCH_SCORE,
            help=('The contribution (usually negative) to a score when a '
                  'feature in a query and subject are part of a match.'))

        parser.add_argument(
            '--weights', action='append',
            help=('A string with the landmark name as the first element and '
                  'the weight as the second, separated by a space.'))

        parser.add_argument(
            '--overallScoreMethod',
            default=FindParameters.DEFAULT_OVERALL_SCORE_METHOD,
            choices=[cls.__name__ for cls in ALL_OVERALL_SCORE_CLASSES],
            help=('The name of the method used to calculate the overall '
                  'score of all histogram bins.'))

        parser.add_argument(
            '--deltaScale', type=float,
            default=FindParameters.DEFAULT_DELTA_SCALE,
            help=('The delta between the query offset and the subject offset '
                  'is scaled by dividing the delta by the deltaScale to '
                  'reduce sensitivity.'))

    @classmethod
    def fromArgs(cls, args):
        """
        Return an instance of FindParameters built from argument values.

        @param args: The result of calling C{parse_args()} on an
            C{argparse.ArgumentParser} instance.
        """
        return cls(significanceMethod=args.significanceMethod,
                   significanceFraction=args.significanceFraction,
                   scoreMethod=args.scoreMethod,
                   featureMatchScore=args.featureMatchScore,
                   featureMismatchScore=args.featureMismatchScore,
                   weights=parseWeights(args.weights or {}),
                   overallScoreMethod=args.overallScoreMethod,
                   deltaScale=args.deltaScale)

    def print_(self, margin=''):
        """
        Print find parameter values.

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the parameters.
        """
        result = MultilineString(margin=margin)
        result.append('Find parameters:')
        result.indent()
        result.extend([
            'Significance method: %s' % self.significanceMethod,
            'Significance fraction: %f' % self.significanceFraction,
            'Score method: %s' % self.scoreMethod,
            'Feature match score: %f' % self.featureMatchScore,
            'Feature mismatch score: %f' % self.featureMismatchScore,
            'OverallScoreMethod: %s' % self.overallScoreMethod,
            'Delta scale: %f' % self.deltaScale,
            'Weights: '
        ])
        result.indent()
        for key in sorted(self.weights.keys()):
            result.append('%s: %f' % (key, self.weights[key]))
        return str(result)


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
    @param featureLengthBase: The length of a landmark is scaled to be its
        logarithm using this C{float} base, for the purpose of matching
        landmarks via hashes. This reduces sensitivity to relatively small
        differences in lengths.
    """

    PARAMS = ('landmarkClasses', 'trigPointClasses', 'limitPerLandmark',
              'maxDistance', 'minDistance', 'distanceBase',
              'featureLengthBase')

    # Database defaults (see explanations in the above docstring).
    DEFAULT_LIMIT_PER_LANDMARK = 10
    DEFAULT_MAX_DISTANCE = 200
    DEFAULT_MIN_DISTANCE = 1
    DEFAULT_DISTANCE_BASE = 1.1
    DEFAULT_FEATURE_LENGTH_BASE = 1.35

    def __init__(self, landmarkClasses, trigPointClasses,
                 limitPerLandmark=None, maxDistance=None, minDistance=None,
                 distanceBase=None, featureLengthBase=None):

        self.distanceBase = (
            self.DEFAULT_DISTANCE_BASE if distanceBase is None
            else distanceBase)

        if self.distanceBase <= 0:
            raise ValueError('distanceBase must be > 0.')

        self.featureLengthBase = (
            self.DEFAULT_FEATURE_LENGTH_BASE if featureLengthBase is None
            else featureLengthBase)

        if self.featureLengthBase <= 0:
            raise ValueError('featureLengthBase must be > 0.')

        self.landmarkClasses = (
            DEFAULT_LANDMARK_CLASSES if landmarkClasses is None
            else landmarkClasses)

        self.landmarkFinders = [landmarkClass(self.featureLengthBase)
                                for landmarkClass in self.landmarkClasses]

        self.trigPointClasses = (
            DEFAULT_TRIG_CLASSES if trigPointClasses is None
            else trigPointClasses)

        self.trigPointFinders = [trigPointClass(self.featureLengthBase)
                                 for trigPointClass in self.trigPointClasses]

        self.limitPerLandmark = (
            self.DEFAULT_LIMIT_PER_LANDMARK if limitPerLandmark is None
            else limitPerLandmark)

        self.maxDistance = (
            self.DEFAULT_MAX_DISTANCE if maxDistance is None else maxDistance)

        self.minDistance = (
            self.DEFAULT_MIN_DISTANCE if minDistance is None else minDistance)

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
                           self.minDistance, self.distanceBase,
                           self.featureLengthBase))))
        return checksum.value

    @staticmethod
    def addArgsToParser(parser):
        """
        Add database creation arguments to an argparse parser.

        @param parser: An C{argparse.ArgumentParser} instance.
        """
        parser.add_argument(
            '--landmark', action='append', dest='landmarkFinderNames',
            choices=sorted(cl.NAME for cl in ALL_LANDMARK_CLASSES),
            help=('The name of a landmark finder to use. May be specified '
                  'multiple times.'))

        parser.add_argument(
            '--trig', action='append', dest='trigFinderNames',
            choices=sorted(cl.NAME for cl in ALL_TRIG_CLASSES),
            help=('The name of a trig point finder to use. May be '
                  'specified multiple times.'))

        parser.add_argument(
            '--defaultLandmarks', action='store_true', default=False,
            help=('If specified, use the default landmark finders: %s' %
                  sorted(cl.NAME for cl in
                         DEFAULT_LANDMARK_CLASSES)))

        parser.add_argument(
            '--defaultTrigPoints', action='store_true', default=False,
            help=('If specified, use the default trig point finders: %s' %
                  sorted(cl.NAME for cl in DEFAULT_TRIG_CLASSES)))

        parser.add_argument(
            '--limitPerLandmark', type=int,
            default=Parameters.DEFAULT_LIMIT_PER_LANDMARK,
            help=('A limit on the number of pairs to yield per landmark '
                  'per read.'))

        parser.add_argument(
            '--maxDistance', type=int,
            default=Parameters.DEFAULT_MAX_DISTANCE,
            help='The maximum distance permitted between yielded pairs.')

        parser.add_argument(
            '--minDistance', type=int,
            default=Parameters.DEFAULT_MIN_DISTANCE,
            help='The minimum distance permitted between yielded pairs.')

        parser.add_argument(
            '--distanceBase', type=float,
            default=Parameters.DEFAULT_DISTANCE_BASE,
            help=('The distance between a landmark and a trig point is '
                  'scaled to be its logarithm using this base. This '
                  'reduces sensitivity to relatively small differences in '
                  'distance.'))

        parser.add_argument(
            '--featureLengthBase', type=float,
            default=Parameters.DEFAULT_FEATURE_LENGTH_BASE,
            help=('The length of landmarks, for the purpose of matching '
                  'hashes, is scaled to be the logarithm of the landmark '
                  'length using this base. This reduces sensitivity to '
                  'relatively small differences in landmark length.'))

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
            'featureLengthBase': self.featureLengthBase,
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

        return Parameters(
            landmarkClasses, trigPointClasses,
            limitPerLandmark=state['limitPerLandmark'],
            maxDistance=state['maxDistance'],
            minDistance=state['minDistance'],
            distanceBase=state['distanceBase'],
            featureLengthBase=state['featureLengthBase'])

    def print_(self, margin=''):
        """
        Print parameter values.

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the parameters.
        """
        result = MultilineString(margin=margin)
        append = result.append
        extend = result.extend

        append('Parameters:')
        result.indent()

        if self.landmarkFinders:
            append('Landmark finders:')
            result.indent()
            extend(sorted(finder.NAME for finder in self.landmarkFinders))
            result.outdent()
        else:
            append('Landmark finders: none')

        if self.trigPointFinders:
            append('Trig point finders:')
            result.indent()
            extend(sorted(finder.NAME for finder in self.trigPointFinders))
            result.outdent()
        else:
            append('Trig point finders: none')

        extend([
            'Limit per landmark: %d' % self.limitPerLandmark,
            'Max distance: %d' % self.maxDistance,
            'Min distance: %d' % self.minDistance,
            'Distance base: %f' % self.distanceBase,
            'Feature length base: %f' % self.featureLengthBase,
        ])

        return str(result)

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
