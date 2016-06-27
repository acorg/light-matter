from scipy import stats

from dark.fasta import FastaReads
from dark.fasta_ss import SSFastaReads
from dark.reads import AAReadWithX, SSAAReadWithX

from light.backend import Backend
from light.database import DatabaseSpecifier
from light.parameters import FindParameters
from light.performance.affinity import affinityMatrix

from light.performance.data.pdb_2hla_a import QUERIES as HLA_Q
from light.performance.data.pdb_2hla_a import SUBJECTS as HLA_S
from light.performance.data.pdb_4mtp_a import QUERIES as MTP_Q
from light.performance.data.pdb_4mtp_a import SUBJECTS as MTP_S
from light.performance.data.polymerase import QUERIES as POLY_Q
from light.performance.data.polymerase import SUBJECTS as POLY_S
from light.performance.data.ha import QUERIES as HA_Q
from light.performance.data.ha import SUBJECTS as HA_S


def evaluateMatch(structureString, start, end, structureType):
    """
    Test if a match is correct. There are four scenarios:
    1) The alpha helix matches part of a sequence that's not an alpha helix.
        --> false positive
    2) The alpha helix matches part of a sequence that's an alpha helix. The
       alpha helix in the sequence doesn't extend to the left or right.
        --> true positive.
    3) The alpha helix matches part of a sequence that's an alpha helix. The
       alpha helix in the sequence extends to the left.
        --> false positive.
    4) The alpha helix matches part of a sequence that's an alpha helix. The
       alpha helix in the sequence extends to the right.
        --> true positive.

    @param structureString: A C{str} of a structure sequence.
    @param start: An C{int} start of the match.
    @param end: An C{int} end of the match.
    @param structureType: A C{str} letter of the structure type that should
        be evaluated. H: Alpha helix, G: Alpha helix 3 10, I: Alpha helix pi,
        E: Extended strand.

    @return: C{True} if the match is a true positive and C{False} if the match
        is a false positive.
    """
    assert 0 <= start < end

    if start > 0 and structureString[start - 1] == structureType:
        return False

    return all(structureString[i] == structureType for i in range(start, end))


def evaluateMatchNoPrefix(structureString, start, end, structureType):
    """
    Test if a match is correct. There are two scenarios:
    1) The alpha helix matches part of a sequence that's not an alpha helix.
        --> false positive
    2) The alpha helix matches part of a sequence that's an alpha helix.
        --> true positive.

    @param structureString: A C{str} of a structure sequence.
    @param start: An C{int} start of the match.
    @param end: An C{int} end of the match.
    @param structureType: A C{str} letter of the structure type that should
        be evaluated. H: Alpha helix, G: Alpha helix 3 10, I: Alpha helix pi,
        E: Extended strand.

    @return: C{True} if the match is a true positive and C{False} if the match
        is a false positive.
    """
    assert 0 <= start < end

    return all(structureString[i] == structureType for i in range(start, end))


class PdbSubsetStatistics(object):

    def __init__(self, fileToEvaluate, fileToEvaluateTpFp, pdbFile,
                 structureFile, structureType):
        """
        A class which provides functions to assess the performance of a subset
        of substrings.

        @param fileToEvaluate: a C{str} filename with substrings that should be
            evaluated.
        @param fileToEvaluateTpFp: a C{str} filename with substrings that
            should be evaluated and their rates of true and false positives.
        @param pdbFile: a C{str} filename of the pdb ss.txt file.
        @param structureFile: a C{str} filename of the structure strings
            extracted from the pdbFile.
        @param structureType: a C{str} name of the structure type that should
            be considered. Must be one of 'AlphaHelix', 'AlphaHelix_3_10',
            'AlphaHelix_pi', 'ExtendedStrand'.
        """
        assert structureType in (
            'AlphaHelix', 'AlphaHelix_3_10', 'AlphaHelix_pi',
            'ExtendedStrand'), (
            'structureType %s must be one of "AlphaHelix", "AlphaHelix_3_10", '
            '"AlphaHelix_pi", "ExtendedStrand"' % structureType)

        self.fileToEvaluate = fileToEvaluate
        self.fileToEvaluateTpFp = fileToEvaluateTpFp
        self.pdbFile = pdbFile
        self.structureFile = structureFile
        self.structureType = structureType
        self.findParams = FindParameters(significanceFraction=0.01,
                                         binScoreMethod='FeatureAAScore')

    def getTotalTpr(self):
        """
        Return the total true positive rate of the subset of substrings that is
        being evaluated.
        """
        totalTp = 0
        totalFp = 0
        totalLength = 0

        with open(self.fileToEvaluateTpFp) as fp:
            for line in fp:
                substring, stats = line.split('')
                tp, fp, r = stats.split(' ')
                totalTp += int(tp[:-1])
                totalFp += int(fp[:-1])
                totalLength += len(substring.strip())
        return totalTp / (totalTp + totalFp)

    def getFractionOfPdbCovered(self):
        """
        Return the fraction of reads in PDB that are covered by at least one
        substring in the subset of substrings that is being evaluated.
        """
        hit = 0
        total = 0

        db = DatabaseSpecifier().getDatabaseFromKeywords(
            trigPoints=[],
            landmarks=['AC ' + self.structureType],
            acExtendedStrandFilename=self.fileToEvaluate)

        backend = Backend()
        backend.configure(db.dbParams)

        for read in SSFastaReads(self.pdbFile, readClass=SSAAReadWithX,
                                 checkAlphabet=0):
            if len(read) > 3:
                total += 1
                scannedRead = backend.scan(read)
                if len(scannedRead.landmarks) > 0:
                    hit += 1

        return hit / total

    def getFractionOfStructuresCovered(self):
        """
        Return the fraction of know structures covered by at least one
        substring in the subset that is being evaluated.
        """
        hit = 0
        total = 0

        db = DatabaseSpecifier().getDatabaseFromKeywords(
            trigPoints=[],
            landmarks=['AC ' + self.structureType],
            acExtendedStrandFilename=self.fileToEvaluate)

        backend = Backend()
        backend.configure(db.dbParams)

        for read in FastaReads(self.structureFile, readClass=AAReadWithX,
                               checkAlphabet=0):
            if len(read) > 2:
                total += 1
                scannedRead = backend.scan(read)
                if len(scannedRead.landmarks) > 0:
                    hit += 1

        return hit / total

    def getCorrelation(self):
        """
        Compute the correlation between light matter scores for the perfect
        PDB finders and for the finders in the subset that are being evaluated.
        """
        result = {}
        datasets = {
            '2HLA': {
                'queries': HLA_Q,
                'subjects': HLA_S,
            },
            '4MTP': {
                'queries': MTP_Q,
                'subjects': MTP_S,
            },
            'Polymerase': {
                'queries': POLY_Q,
                'subjects': POLY_S,
            },
            'HA': {
                'queries': HA_Q,
                'subjects': HA_S,
            },
        }

        for data in datasets:
            pdbScores = []
            evaluateScores = []

            pdbMatrix = affinityMatrix(datasets[data]['queries'],
                                       subjects=datasets[data]['subjects'],
                                       symmetric=False, computeDiagonal=True,
                                       returnDict=False,
                                       findParams=self.findParams,
                                       landmarks=['PDB ' + self.structureType],
                                       trigPoints=[])

            for i, query in enumerate(datasets[data]['queries']):
                for j, subject in enumerate(datasets[data]['subjects']):
                    pdbScores.append(pdbMatrix[i][j])

            evaluateMatrix = affinityMatrix(
                datasets[data]['queries'], subjects=datasets[data]['subjects'],
                symmetric=False, computeDiagonal=True, returnDict=False,
                findParams=self.findParams,
                landmarks=['AC ' + self.structureType],
                trigPoints=[])

            for i, query in enumerate(datasets[data]['queries']):
                for j, subject in enumerate(datasets[data]['subjects']):
                    evaluateScores.append(evaluateMatrix[i][j])

            slope, intercept, rValue, pValue, se = stats.linregress(
                pdbScores, evaluateScores)

            result[data] = rValue

        return result

    def runAll(self):
        """
        Return the result of all evaluations.
        """
        totalTpr = self.getTotalTpr()
        fractionOfPdbCovered = self.getFractionOfPdbCovered()
        fractionOfStructuresCovered = self.getFractionOfStructuresCovered()
        correlation = self.getCorrelation()

        return (totalTpr, fractionOfPdbCovered, fractionOfStructuresCovered,
                correlation)
