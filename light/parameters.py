from __future__ import print_function

import sys
try:
    from ujson import dumps, loads
except ImportError:
    from json import dumps, loads

from six import string_types

from light.checksum import Checksum
from light.landmarks import (
    findLandmark, findLandmarks, DEFAULT_LANDMARK_CLASSES,
    ALL_LANDMARK_CLASSES_EVEN_BAD_ONES)
from light.trig import (
    findTrigPoint, findTrigPoints, DEFAULT_TRIG_CLASSES,
    ALL_TRIG_CLASSES_EVEN_BAD_ONES)
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
    @param binScoreMethod: The C{str} name of the method used to calculate the
        score of a histogram bin which is considered significant.
    @param featureMatchScore: The C{float} contribution (usually positive) to
        a score when a feature in a query and subject are part of a match.
    @param featureMismatchScore: The C{float} contribution (usually negative)
        to a score when a feature in a query and subject are part of a match.
    @param deltaScale: A C{float}. The delta between the queryOffset and the
        subjectOffset of a matching pair is scaled by dividing the delta by the
        deltaScale to reduce sensitivity.
    @param weights: C{dict} where the keys are feature names and the values are
        C{float} weights that should be assigned to that feature in the
        WeightedFeatureAAScore calculation.
    @param deltaScale: A C{float}. The delta between the query offset and the
        subject offset is scaled by dividing the delta by the deltaScale to
        reduce sensitivity.
    """
    # The methods to be used to calculate match scores and whether matches
    # are significant (i.e., worth reporting).
    DEFAULT_BIN_SCORE_METHOD = MinHashesScore.__name__
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
                 binScoreMethod=None, overallScoreMethod=None,
                 featureMatchScore=None, featureMismatchScore=None,
                 weights=None, deltaScale=None):
        self.significanceMethod = (
            self.DEFAULT_SIGNIFICANCE_METHOD if significanceMethod is None
            else significanceMethod)

        self.significanceFraction = (
            self.DEFAULT_SIGNIFICANCE_FRACTION if significanceFraction is None
            else significanceFraction)

        self.binScoreMethod = (
            self.DEFAULT_BIN_SCORE_METHOD if binScoreMethod is None else
            binScoreMethod)

        self.overallScoreMethod = (
            self.DEFAULT_OVERALL_SCORE_METHOD if overallScoreMethod is
            None else overallScoreMethod)

        self.featureMatchScore = (
            self.DEFAULT_FEATURE_MATCH_SCORE if featureMatchScore is None
            else featureMatchScore)

        self.featureMismatchScore = (
            self.DEFAULT_FEATURE_MISMATCH_SCORE if featureMismatchScore is None
            else featureMismatchScore)

        self.weights = self.DEFAULT_WEIGHTS if weights is None else weights

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
            '--binScoreMethod',
            default=FindParameters.DEFAULT_BIN_SCORE_METHOD,
            choices=[cls.__name__ for cls in ALL_BIN_SCORE_CLASSES],
            help=('The name of the method used to calculate the score of a '
                  'histogram bin which is considered significant.'))

        parser.add_argument(
            '--overallScoreMethod',
            default=FindParameters.DEFAULT_OVERALL_SCORE_METHOD,
            choices=[cls.__name__ for cls in ALL_OVERALL_SCORE_CLASSES],
            help=('The name of the method used to calculate the overall '
                  'score of all histogram bins.'))

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
                   binScoreMethod=args.binScoreMethod,
                   featureMatchScore=args.featureMatchScore,
                   featureMismatchScore=args.featureMismatchScore,
                   weights=parseWeights(args.weights or {}),
                   overallScoreMethod=args.overallScoreMethod,
                   deltaScale=args.deltaScale)

    def print_(self, margin='', result=None):
        """
        Print find parameter values.

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @param result: A C{MultilineString} instance, or C{None} if a new
            C{MultilineString} should be created.
        @return: If C{result} was C{None}, return a C{str} representation of
            the parameters, else C{None}.
        """
        if result is None:
            result = MultilineString(margin=margin)
            returnNone = False
        else:
            returnNone = True

        result.append('Find parameters:')
        result.indent()
        result.extend([
            'Significance method: %s' % self.significanceMethod,
            'Significance fraction: %f' % self.significanceFraction,
            'Bin Score Method: %s' % self.binScoreMethod,
            'Feature match score: %f' % self.featureMatchScore,
            'Feature mismatch score: %f' % self.featureMismatchScore,
            'Overall Score Method: %s' % self.overallScoreMethod,
            'Delta scale: %f' % self.deltaScale,
            'Weights: '
        ])
        result.indent()
        for key in sorted(self.weights):
            result.append('%s: %f' % (key, self.weights[key]))

        result.outdent()
        result.outdent()

        if not returnNone:
            return str(result)


class DatabaseParameters(object):
    """
    Hold a collection of database parameter settings, including parameters
    used by feature finders.

    @param landmarks: Either C{None} (to use the default landmark finders) or
        a mixed C{list} of landmark finder classes or C{str} landmark finder
        names. To specify no landmark finders, pass an empty list.
    @param trigPoints: Either C{None} (to use the default trig point finders)
        or a mixed C{list} of trig point finder classes or C{str} trig point
        finder names. To specify no trig point finders, pass an empty list.
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
    @param randomLandmarkDensity: The C{float} density of random length-1
        landmarks to generate.
    @param randomTrigPointDensity: The C{float} density of random trig points
        to generate.
    @raise ValueError: If an unknown landmark or trig point finder name is
        given in C{landmarks} or C{trigPoints}.
    """

    # Database defaults (see explanations in the above docstring).
    DEFAULT_LIMIT_PER_LANDMARK = 10
    DEFAULT_MAX_DISTANCE = 200
    DEFAULT_MIN_DISTANCE = 1
    DEFAULT_DISTANCE_BASE = 1.1
    DEFAULT_FEATURE_LENGTH_BASE = 1.35
    DEFAULT_RANDOM_LANDMARK_DENSITY = 0.1
    DEFAULT_RANDOM_TRIG_POINT_DENSITY = 0.1

    def __init__(self, landmarks=None, trigPoints=None,
                 limitPerLandmark=None, maxDistance=None, minDistance=None,
                 distanceBase=None, featureLengthBase=None,
                 randomLandmarkDensity=None, randomTrigPointDensity=None):

        # First set the simple scalar parameters.
        self.limitPerLandmark = (
            self.DEFAULT_LIMIT_PER_LANDMARK if limitPerLandmark is None
            else limitPerLandmark)

        self.maxDistance = (
            self.DEFAULT_MAX_DISTANCE if maxDistance is None else maxDistance)

        self.minDistance = (
            self.DEFAULT_MIN_DISTANCE if minDistance is None else minDistance)

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

        self.randomLandmarkDensity = (
            self.DEFAULT_RANDOM_LANDMARK_DENSITY if
            randomLandmarkDensity is None else randomLandmarkDensity)

        self.randomTrigPointDensity = (
            self.DEFAULT_RANDOM_TRIG_POINT_DENSITY if
            randomTrigPointDensity is None else randomTrigPointDensity)

        if landmarks is None:
            landmarkClasses = DEFAULT_LANDMARK_CLASSES
        else:
            landmarkClasses = set()
            for landmark in landmarks:
                if isinstance(landmark, string_types):
                    cls = findLandmark(landmark)
                    if cls:
                        landmarkClasses.add(cls)
                    else:
                        raise ValueError(
                            'Could not find landmark finder class %r.'
                            % landmark)
                else:
                    # Assume this is already a landmark class.
                    landmarkClasses.add(landmark)

        if trigPoints is None:
            trigPointClasses = DEFAULT_TRIG_CLASSES
        else:
            trigPointClasses = set()
            for trigPoint in trigPoints:
                if isinstance(trigPoint, string_types):
                    cls = findTrigPoint(trigPoint)
                    if cls:
                        trigPointClasses.add(cls)
                    else:
                        raise ValueError(
                            'Could not find trig point finder class %r.'
                            % trigPoint)
                else:
                    # Assume this is already a trig point class.
                    trigPointClasses.add(trigPoint)

        # The finders instantiated here are not used in this file. They are
        # used by the backend. We make sorted lists of them so we're guaranteed
        # to always process them in the same order (in printing, in checksums,
        # etc).
        self.landmarkFinders = sorted(cls(self) for cls in landmarkClasses)
        self.trigPointFinders = sorted(cls(self) for cls in trigPointClasses)

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
        checksum.update(
            [f.NAME for f in self.landmarkFinders] +
            [f.SYMBOL for f in self.landmarkFinders] +
            [f.NAME for f in self.trigPointFinders] +
            [f.SYMBOL for f in self.trigPointFinders] +
            list(map(str, (self.limitPerLandmark, self.maxDistance,
                           self.minDistance, self.distanceBase,
                           self.featureLengthBase,
                           self.randomLandmarkDensity,
                           self.randomTrigPointDensity))))
        return checksum.value

    @staticmethod
    def addArgsToParser(parser):
        """
        Add database creation arguments to an argparse parser.

        @param parser: An C{argparse.ArgumentParser} instance.
        """
        parser.add_argument(
            '--landmark', action='append', dest='landmarks',
            choices=sorted(c.NAME for c in ALL_LANDMARK_CLASSES_EVEN_BAD_ONES),
            help=('The name of a landmark finder to use. May be specified '
                  'multiple times. If no landmark finders are '
                  'given (and --noLandmarks is not specified), the default '
                  'set of landmark finders (%s) will be used.' %
                  ', '.join(sorted(c.NAME for c in DEFAULT_LANDMARK_CLASSES))))

        parser.add_argument(
            '--noLandmarks', action='store_true', default=False,
            help='If specified, no landmark finders will be used.')

        parser.add_argument(
            '--trig', action='append', dest='trigPoints',
            choices=sorted(c.NAME for c in ALL_TRIG_CLASSES_EVEN_BAD_ONES),
            help=('The name of a trig point finder to use. May be '
                  'specified multiple times. If no trig point finders are '
                  'given (and --noTrigPoints is not specified), the default '
                  'set of trig point finders (%s) will be used.' %
                  ', '.join(sorted(c.NAME for c in DEFAULT_TRIG_CLASSES))))

        parser.add_argument(
            '--noTrigPoints', action='store_true', default=False,
            help='If specified, no trig point finders will be used.')

        parser.add_argument(
            '--limitPerLandmark', type=int,
            default=DatabaseParameters.DEFAULT_LIMIT_PER_LANDMARK,
            help=('A limit on the number of pairs to yield per landmark '
                  'per read.'))

        parser.add_argument(
            '--maxDistance', type=int,
            default=DatabaseParameters.DEFAULT_MAX_DISTANCE,
            help='The maximum distance permitted between yielded pairs.')

        parser.add_argument(
            '--minDistance', type=int,
            default=DatabaseParameters.DEFAULT_MIN_DISTANCE,
            help='The minimum distance permitted between yielded pairs.')

        parser.add_argument(
            '--distanceBase', type=float,
            default=DatabaseParameters.DEFAULT_DISTANCE_BASE,
            help=('The distance between a landmark and a trig point is '
                  'scaled to be its logarithm using this base. This '
                  'reduces sensitivity to relatively small differences in '
                  'distance.'))

        parser.add_argument(
            '--featureLengthBase', type=float,
            default=DatabaseParameters.DEFAULT_FEATURE_LENGTH_BASE,
            help=('The length of landmarks, for the purpose of matching '
                  'hashes, is scaled to be the logarithm of the landmark '
                  'length using this base. This reduces sensitivity to '
                  'relatively small differences in landmark length.'))

        parser.add_argument(
            '--randomLandmarkDensity', type=float,
            default=DatabaseParameters.DEFAULT_RANDOM_LANDMARK_DENSITY,
            help=('The (float) density of random length-1 landmarks to '
                  'generate.'))

        parser.add_argument(
            '--randomTrigPointDensity', type=float,
            default=DatabaseParameters.DEFAULT_RANDOM_TRIG_POINT_DENSITY,
            help=('The (float) density of random trig points to generate.'))

    @classmethod
    def fromArgs(cls, args):
        if args.noLandmarks:
            if args.landmarks:
                raise RuntimeError(
                    'You cannot use both --noLandmarks and --landmark.')
            else:
                landmarkClasses = []
        elif args.landmarks:
            landmarkClasses = findLandmarks(args.landmarks)
        else:
            landmarkClasses = DEFAULT_LANDMARK_CLASSES

        if args.noTrigPoints:
            if args.trigPoints:
                raise RuntimeError(
                    'You cannot use both --noTrigPoints and --trig.')
            else:
                trigPointClasses = []
        elif args.trigPoints:
            trigPointClasses = findTrigPoints(args.trigPoints)
        else:
            trigPointClasses = DEFAULT_TRIG_CLASSES

        # Flag an error if no landmark or trig point finders were given. We
        # test this here instead of in __init__ because we want to allow
        # programmers (e.g., writing tests) to have no finders. But a user
        # giving command line arguments should probably be warned /
        # restricted.
        if len(landmarkClasses) + len(trigPointClasses) == 0:
            raise RuntimeError(
                'Cannot create database parameters as no landmark or trig '
                'point finders were specified. Use --landmark and/or --trig.')

        return cls(landmarks=landmarkClasses, trigPoints=trigPointClasses,
                   limitPerLandmark=args.limitPerLandmark,
                   maxDistance=args.maxDistance,
                   minDistance=args.minDistance,
                   distanceBase=args.distanceBase,
                   featureLengthBase=args.featureLengthBase,
                   randomLandmarkDensity=args.randomLandmarkDensity,
                   randomTrigPointDensity=args.randomTrigPointDensity)

    def save(self, fp=sys.stdout):
        """
        Save the parameters in JSON format.

        @param fp: A file pointer.
        @return: The C{fp} we were passed (this is useful in testing).
        """
        state = {
            'landmarks': [finder.NAME for finder in self.landmarkFinders],
            'trigPoints': [finder.NAME for finder in self.trigPointFinders],
            'limitPerLandmark': self.limitPerLandmark,
            'maxDistance': self.maxDistance,
            'minDistance': self.minDistance,
            'distanceBase': self.distanceBase,
            'featureLengthBase': self.featureLengthBase,
            'randomLandmarkDensity': self.randomLandmarkDensity,
            'randomTrigPointDensity': self.randomTrigPointDensity,
        }

        print(dumps(state), file=fp)
        return fp

    @staticmethod
    def restore(fp=sys.stdin):
        """
        Load parameters from a file.

        @param fp: A file pointer.
        @return: An instance of L{DatabaseParameters}.
        @raises ValueError: If a now non-existent landmark or trig point name
            is found in the saved parameters file. Or if valid JSON cannot be
            loaded from C{fp}.
        """
        state = loads(fp.readline()[:-1])
        return DatabaseParameters(**state)

    def print_(self, margin='', result=None):
        """
        Print parameter values.

        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @param result: A C{MultilineString} instance, or C{None} if a new
            C{MultilineString} should be created.
        @return: If C{result} was C{None}, return a C{str} representation of
            the parameters, else C{None}.
        """
        if result is None:
            result = MultilineString(margin=margin)
            returnNone = False
        else:
            returnNone = True

        append = result.append

        append('Parameters:')
        result.indent()

        if self.landmarkFinders:
            append('Landmark finders:')
            result.indent()
            for finder in self.landmarkFinders:
                append(finder.__class__.__name__)
            result.outdent()
        else:
            append('Landmark finders: none')

        if self.trigPointFinders:
            append('Trig point finders:')
            result.indent()
            for finder in self.trigPointFinders:
                append(finder.__class__.__name__)
            result.outdent()
        else:
            append('Trig point finders: none')

        result.extend([
            'Limit per landmark: %d' % self.limitPerLandmark,
            'Max distance: %d' % self.maxDistance,
            'Min distance: %d' % self.minDistance,
            'Distance base: %f' % self.distanceBase,
            'Feature length base: %f' % self.featureLengthBase,
            'Random landmark density: %f' % self.randomLandmarkDensity,
            'Random trig point density: %f' % self.randomTrigPointDensity,
        ])

        result.outdent()

        if not returnNone:
            return str(result)

    def compare(self, other):
        """
        Compare our parameters against another C{DatabaseParameters} instance.

        @param other: A C{DatabaseParameters} instance.
        @return: A C{str} summary of the parameter differences if any, else
            C{None}.
        """
        err = MultilineString()
        err.append('Summary of differences:')
        err.indent()

        ourLandmarks = self.landmarkFinderNames()
        otherLandmarks = other.landmarkFinderNames()
        if ourLandmarks != otherLandmarks:
            err.append("Param 'landmarks' values %r and %r differ." %
                       (ourLandmarks, otherLandmarks))

        ourTrigPoints = self.trigPointFinderNames()
        otherTrigPoints = other.trigPointFinderNames()
        if ourTrigPoints != otherTrigPoints:
            err.append("Param 'trigPoints' values %r and %r differ." %
                       (ourTrigPoints, otherTrigPoints))

        for param in (
                'limitPerLandmark', 'maxDistance', 'minDistance',
                'distanceBase', 'featureLengthBase', 'randomLandmarkDensity',
                'randomTrigPointDensity'):
            ours = getattr(self, param)
            others = getattr(other, param)
            if ours != others:
                err.append('Param %r values %r and %r differ.' %
                           (param, ours, others))

        # We have an error if the multi-line string result has more than
        # the initial 'Summary of differences' heading line.
        return str(err) if err.lineCount() > 1 else None

    def landmarkFinderNames(self):
        """
        Get a list of the names of our landmark finders.

        @return: A C{list} of C{str} landmark finder names. These will be in
            the sorted order of landmark finders set in our __init__.
        """
        return [finder.__class__.NAME for finder in self.landmarkFinders]

    def trigPointFinderNames(self):
        """
        Get a list of the names of our trig point finders.

        @return: A C{list} of C{str} trig point finder names. These will be in
            the sorted order of trig point finders set in our __init__.
        """
        return [finder.__class__.NAME for finder in self.trigPointFinders]
