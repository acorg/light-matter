import re

from dark.reads import AARead

from light.bin_score import (
    NoneScore, MinHashesScore, FeatureMatchingScore, FeatureAAScore,
    WeightedFeatureAAScore)
from light.database import Database
from light.features import Landmark, TrigPoint
from light.histogram import Histogram
from light.landmarks import findLandmark
from light.parameters import DatabaseParameters, FindParameters
from light.trig import findTrigPoint

_NONWHITE_REGEXP = re.compile('\S+')


class _QueryOrSubject(object):
    """
    Hold information about a query or subject found in a template.

    @param templateList: A C{list} of C{str}s with the lines from the template.
        Each line has been stripped of whitespace on its right.
    """
    def __init__(self, templateList):
        # Some of the intermediate values computed below are added to self
        # so that their values can be tested.
        self.landmarks = set()
        self.trigPoints = set()
        self.indentLength, self.relevantLines = self.extractRelevantLines(
            templateList)
        self.pipe1, self.pipe2 = self.findPipes()
        self.noPipes = [line.replace('|', '') for line in self.relevantLines]
        self.read = AARead(self.TYPE, self.noPipes[0][self.indentLength:])
        self.nonPairFeatures = list(self.extractNonPairFeatures())

        # Extract and check all unpaired features and remember the names of the
        # landmarks and trig points that were mentioned.
        for featureName, feature in self.nonPairFeatures:
            if isinstance(feature, Landmark):
                self.landmarks.add(featureName)
            else:
                self.trigPoints.add(featureName)

        # Extract and check all paired features and remember the names of the
        # landmarks and trig points that were mentioned.
        self.pairedFeatures = list(self.extractPairedFeatures())
        for (landmarkName, landmark,
                secondFeatureName, secondFeature) in self.pairedFeatures:
            self.landmarks.add(landmarkName)

            if isinstance(secondFeature, Landmark):
                self.landmarks.add(secondFeatureName)
            else:
                self.trigPoints.add(secondFeatureName)

    def extractRelevantLines(self, templateList):
        """
        Find lines for a subject or query in a list of template strings.

        @param templateList: A C{list} of C{str}s with the lines from the
            template.
        @raise ValueError: If no line starting the subject or query group is
            found.
        @return: A 2-tuple containing
            1. an C{int} giving the length of the prefix of the start line up
               to the start indicator ('query' or 'subject') and all its
               subsequent whitespace.
            2. A C{list} of lines from the subject or query section of the
               template line list.
        """
        startRe = re.compile('^\s*' + self.TYPE + '\s+', re.IGNORECASE)
        stopRe = re.compile('^\s*' + self.OTHER_TYPE + '\s+', re.IGNORECASE)
        foundStart = False
        result = []

        for line in templateList:
            if foundStart:
                if stopRe.match(line):
                    break
                else:
                    result.append(line)
            else:
                match = startRe.match(line)
                if match:
                    foundStart = True
                    indentLength = match.end() - match.start()
                    result.append(line)

        if not foundStart:
            raise ValueError('Could not find %r line in template' % self.TYPE)

        return indentLength, result

    def findPipes(self):
        """
        Check that all relevant lines for the subject or query have two pipe
        ('|') characters, and that these occur at the same places in each line.

        @raise ValueError: If any line does not contain two pipe characters
            or if the pipes are not all in the same places.
        @return: A 2-tuple containing the C{int} offsets of the two pipes.
        """
        for i, line in enumerate(self.relevantLines):
            if line.count('|') != 2:
                raise ValueError('Template line %r in %s does not contain '
                                 'two pipe characters' % (line, self.TYPE))
            if i == 0:
                pipe1 = line.index('|')
                pipe2 = line.index('|', pipe1 + 1)
            else:
                check1 = line.index('|')
                check2 = line.index('|', pipe1 + 1)
                if check1 != pipe1 or check2 != pipe2:
                    raise ValueError(
                        'Template line %r in %s has pipe characters at '
                        'offsets %d and %d, but these do not match the '
                        'offsets (%d and %d) of the pipe characters found on '
                        'the first %s template line' %
                        (line, self.TYPE, check1, check2, pipe1, pipe2,
                         self.TYPE))

        return pipe1, pipe2

    def extractNonPairFeatures(self):
        """
        Find all non-pair features in the template.

        @raise ValueError: If a line contains an unknown landmark or trig
            point name or if the sequence for a feature does not match the
            overall full sequence given in the first line of the template
            for the subject/query.
        @return: A generator that yields 2-tuples containing:
            1. A C{str} feature name
            2. A C{light.features.Landmark} or a C{light.features.TrigPoint}
                instance
            as they are found in the template lines.
        """
        for line in self.noPipes[1:]:
            if line.find(',') != -1:
                # There is a comma on the line, so this is a pair of
                # features from the matched region. Ignore it as it will be
                # processed in extractPairedFeatures.
                continue

            featureName = line[:self.indentLength].strip()
            featureStr = line[self.indentLength:]

            landmark = findLandmark(featureName)

            if landmark is None:
                trigPoint = findTrigPoint(featureName)

                if trigPoint is None:
                    raise ValueError('Unknown feature name %r found in %s '
                                     'template' % (featureName, self.TYPE))

            for match in _NONWHITE_REGEXP.finditer(featureStr):
                # print('Found %s match %d to %d in %r' % (
                #     featureName, match.start(), match.end(), line))
                offset = match.start()
                length = match.end() - offset
                sequence = featureStr[offset:offset + length]
                if sequence != self.read.sequence[offset:offset + length]:
                    raise ValueError(
                        '%s feature sequence %r found in %s template (offset '
                        '%d, length %d) does not match the full sequence for '
                        'the %s at those offsets' %
                        (featureName, sequence, self.TYPE, offset, length,
                         self.TYPE))
                if landmark:
                    yield (featureName,
                           Landmark(landmark.NAME, landmark.SYMBOL, offset,
                                    length))
                else:
                    yield (featureName,
                           TrigPoint(trigPoint.NAME, trigPoint.SYMBOL, offset))

    def extractPairedFeatures(self):
        """
        Find all paired features in the template.

        @raise ValueError: If a line contains an unknown landmark or trig
            point name or doesn't have exactly one comma (separating the two
            names).
        @return: A generator that yields 4-tuples containing
            1. A C{str} landmark feature name
            2. A C{light.features.Landmark} instance
            3. A C{str} trig point feature name (this may actually be a
                landmark, but we typically call it a trig point).
            4. A C{light.features.TrigPoint} instance (may actuall be a
                C{light.features.Landmark} instance)
            as they are found in the template lines.
        """
        for line in self.noPipes[1:]:
            if line.find(',') == -1:
                # This line contains non-paired features (it has just one
                # feature name - implied by the lack of a comma). Ignore it
                # as it will be processed in extractNonPairFeatures.
                continue

            if line.count(',') != 1:
                raise ValueError('%s template line %r contains multiple '
                                 'commas' % (self.TYPE, line))

            landmarkName, secondFeatureName = map(
                str.strip, line[:self.indentLength].split(','))
            featureStr = line[self.indentLength:]

            landmark = findLandmark(landmarkName)

            # The first feature of a matched pair must be a landmark.
            if landmark is None:
                raise ValueError('Unknown landmark name %r found in %s '
                                 'template line %r' %
                                 (landmarkName, self.TYPE, line))

            # The second feature of a matched pair can be a trig point or
            # another landmark.
            landmark2 = findLandmark(secondFeatureName)

            if landmark2 is None:
                trigPoint = findTrigPoint(secondFeatureName)

                if trigPoint is None:
                    raise ValueError('Unknown feature name %r found in %s '
                                     'template' % (secondFeatureName,
                                                   self.TYPE))

            for count, match in enumerate(
                    _NONWHITE_REGEXP.finditer(featureStr)):
                if count > 1:
                    raise ValueError(
                        'More than two features found in matched region pair '
                        'line in %s template. Line was %r' % (self.TYPE, line))
                first = count == 0
                # print('Found %s/%s match %d to %d in %r' % (
                #     landmarkName, secondFeatureName, match.start(),
                #     match.end(), line))
                offset = match.start()
                length = match.end() - offset
                sequence = featureStr[offset:offset + length]
                if sequence != self.read.sequence[offset:offset + length]:
                    if first:
                        raise ValueError(
                            '%s feature sequence %r found in %s template '
                            '(offset %d, length %d) does not match the full '
                            'sequence for the %s at those offsets' %
                            (landmarkName, sequence, self.TYPE, offset, length,
                             self.TYPE))
                    else:
                        raise ValueError(
                            '%s feature sequence %r found in %s template '
                            '(offset %d, length %d) does not match the full '
                            'sequence for the %s at those offsets' %
                            (secondFeatureName, sequence, self.TYPE, offset,
                             length, self.TYPE))

                if first:
                    landmarkFeature = Landmark(landmark.NAME, landmark.SYMBOL,
                                               offset, length)
                else:
                    if landmark2:
                        secondFeature = Landmark(landmark2.NAME,
                                                 landmark2.SYMBOL, offset,
                                                 length)
                    else:
                        secondFeature = TrigPoint(trigPoint.NAME,
                                                  trigPoint.SYMBOL, offset)

                    yield (landmarkName, landmarkFeature,
                           secondFeatureName, secondFeature)


class Query(_QueryOrSubject):
    TYPE, OTHER_TYPE = 'query', 'subject'


class Subject(_QueryOrSubject):
    TYPE, OTHER_TYPE = 'subject', 'query'


class Template(object):
    """
    Parse an ASCII art picture of a light matter match and provide access to
    it.

    @param template: A C{str} template picture of the match.
    @raise ValueError: If the query and subject do not have the same number of
        paired features.
    """
    def __init__(self, template):
        self.template = self.templateToList(template)
        self.query = Query(self.template)
        self.subject = Subject(self.template)

        if len(self.query.pairedFeatures) != len(self.subject.pairedFeatures):
            raise ValueError(
                'The query and subject do not have the same number of paired '
                'features (%d != %d)' % (len(self.query.pairedFeatures),
                                         len(self.subject.pairedFeatures)))

        # Union the landmark and trig point names from the query and subject.
        self.landmarks = self.query.landmarks | self.subject.landmarks
        self.trigPoints = self.query.trigPoints | self.subject.trigPoints

        self.histogram = Histogram(1)

        for queryPair, subjectPair in zip(self.query.pairedFeatures,
                                          self.subject.pairedFeatures):
            _, queryLandmark, _, queryTrigPoint = queryPair
            _, subjectLandmark, _, subjectTrigPoint = subjectPair

            self.histogram.add(0, {
                'queryLandmark': queryLandmark,
                'queryTrigPoint': queryTrigPoint,
                'subjectLandmark': subjectLandmark,
                'subjectTrigPoint': subjectTrigPoint,
            })

        self.histogram.finalize()

    @staticmethod
    def templateToList(template):
        """
        Convert a picture to a list of trimmed non-blank lines.

        @param template: A C{str} template picture of the match.
        @return: A C{list} of \n separated non-blank lines from C{template}.
        """
        result = []
        whitespace = re.compile('^\s*$')
        for line in template.split('\n'):
            if whitespace.match(line) is None:
                result.append(line.rstrip())
        return result

    def calculateScore(self, dbParams=None, findParams=None):
        """
        Using a given scoring method, calculate the score of the alignment
        between the query and subject in the template.

        @param findParams: An instance of C{light.parameters.FindParameters} or
            C{None} to use default find parameters.
        @raises ValueError: If C{dbParams} is passed and the landmarks and
            trig points it specifies do not include all the landmarks and trig
            points named in the template. Of if the C{binScoreMethod} in
            C{findParams} is unknown.
        @return: A 2-tuple, being the result of calling the C{calculateScore}
            method of the C{binScoreMethod} class. The tuple contains a
            C{float} score of the bin and a C{dict} with the analysis leading
            to the score (see light/bin_score.py).
        """
        findParams = findParams or FindParameters()
        if dbParams is None:
            dbParams = DatabaseParameters(landmarks=self.landmarks,
                                          trigPoints=self.trigPoints)
        else:
            missing = self.landmarks - set(dbParams.landmarkFinderNames())
            if missing:
                raise ValueError(
                    'The template mentions landmark finders (%s) that are '
                    'not present in the passed DatabaseParameters instance' %
                    ', '.join(sorted(missing)))

            missing = self.trigPoints - set(dbParams.trigPointFinderNames())
            if missing:
                raise ValueError(
                    'The template mentions trig point finders (%s) that are '
                    'not present in the passed DatabaseParameters instance' %
                    ', '.join(sorted(missing)))

        database = Database(dbParams=dbParams)
        _, subjectIndex, subjectHashCount = database.addSubject(
            self.subject.read)
        dbSubject = database.getSubjectByIndex(subjectIndex)

        binScoreMethod = findParams.binScoreMethod
        if binScoreMethod == 'NoneScore':
            scorer = NoneScore()
        elif binScoreMethod == 'MinHashesScore':
            be = database._connector._backend
            queryHashCount = 0
            scannedQuery = be.scan(self.query.read)
            for hashInfo in be.getHashes(scannedQuery).values():
                queryHashCount += len(hashInfo)
            scorer = MinHashesScore(self.histogram,
                                    min(queryHashCount, subjectHashCount))
        elif binScoreMethod == 'FeatureMatchingScore':
            scorer = FeatureMatchingScore(
                self.histogram, self.query.read, dbSubject, dbParams,
                findParams)
        elif binScoreMethod == 'FeatureAAScore':
            scorer = FeatureAAScore(
                self.histogram, self.query.read, dbSubject, dbParams)
        elif binScoreMethod == 'WeightedFeatureAAScore':
            scorer = WeightedFeatureAAScore(
                self.histogram, self.query.read, dbSubject, dbParams,
                findParams.weights)
        else:
            raise ValueError('Unknown bin score method %r' % binScoreMethod)

        return scorer.calculateScore(0)
