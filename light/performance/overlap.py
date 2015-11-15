from collections import defaultdict
from Bio import SeqIO

from light.backend import Backend
from light.database import DatabaseSpecifier

from dark.reads import Reads, AARead


class SSAARead(AARead):
    """
    Hold information to work with AAReads that have secondary structure
    information attached to them.

    @param id: A C{str} describing the read.
    @param sequence: A C{str} of sequence information.
    @param structure: A C{str} of structure information.
    """
    def __init__(self, id, sequence, structure):
        super().__init__(id, sequence)
        self.structure = structure


class CalculateOverlap(object):
    """
    A class which calculates the overlap between the features found by our
    finders and the secondary structures found by DSSP. The secondary
    structures found by DSSP were downloaded from
    http://www.rcsb.org/pdb/files/ss.txt on the 11/11/2015.

    @param pdbFile: The C{str} filename of the file from pdb containing the
        sequence and structural information.
    """
    def __init__(self, pdbFile):
        # The pdb file is in fasta format. For each structure it contains the
        # amino acid sequence on one line and the predicted secondary structure
        # sequence on the other.
        self.SSAAReads = Reads()
        records = [record for record in SeqIO.parse(pdbFile, 'fasta')]
        for i in range(0, len(records), 2):
            record = records[i]
            read = SSAARead(record.id, record.seq,
                            records[i + 1].seq)
            self.SSAAReads.add(read)

    def calculateOverlap(self, printAll=False):
        """
        Calculate the feature overlap for a set of features.

        @param printAll: A C{bool} saying whether the information about the
            overlap in each read should be printed out.
        """
        allAASequenceFeatures = defaultdict(list)
        allSSSequenceFeatures = defaultdict(list)
        allIntersects = defaultdict(list)

        for read in self.SSAAReads:
            aaSeqFeatures, ssSeqFeatures, intersects = self.getFeatures(
                read, printAll=printAll)
            for feature, indices in aaSeqFeatures:
                allAASequenceFeatures[feature].extend(feature[indices])
            for feature, indices in ssSeqFeatures:
                allSSSequenceFeatures[feature].extend(feature[indices])
            for featureIntersect, indices in intersects:
                allIntersects[feature].extend(feature[indices])

        self.calculateFraction(allAASequenceFeatures, allSSSequenceFeatures,
                               allIntersects, print_=printAll)

    @staticmethod
    def getFeatures(ssAARead, printAll=False, **kwargs):
        """
        Extract the features from the sequence and the structural information.

        @param ssAARead: An C{SSAARead} instance.
        @param printAll: A C{bool} saying whether the information about the
            overlap in the ssAARead should be printed out.
        @param kwargs: See
            C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
            additional keywords, all of which are optional.
        """
        if 'landmarkNames' not in kwargs:
            kwargs['landmarkNames'] = ['AlphaHelix', 'AlphaHelix_pi',
                                       'AlphaHelix_3_10', 'BetaStrand',
                                       'GOR4AlphaHelix', 'GOR4BetaStrand']
        if 'trigPointNames' not in kwargs:
            kwargs['trigPointNames'] = []

        db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        backend = Backend()
        backend.configure(db.params)

        aaSequenceFeatures = defaultdict(set)
        ssSequenceFeatures = defaultdict(set)
        intersects = defaultdict(set)

        scannedAaSequence = backend.scan(ssAARead)

        # Get all offsets for each landmark separately.
        for landmark in scannedAaSequence.landmarks:
            aaSequenceFeatures[landmark.name].update(landmark.coveredOffsets())

        # Get all offsets for each secondary structure feature separately.
        for offset, structure in enumerate(ssAARead.structure):
            ssSequenceFeatures[structure].add(offset)

        # Keep a record of the intersects for each read, so we can later do the
        # calculation on the whole dataset.
        intersects['AlphaHelix_H'] = (
            aaSequenceFeatures['AlphaHelix'] | ssSequenceFeatures['H'])

        intersects['AlphaHelix_3_10_G'] = (
            aaSequenceFeatures['AlphaHelix_3_10'] |
            ssSequenceFeatures['G'])

        intersects['AlphaHelix_pi_I'] = (
            aaSequenceFeatures['AlphaHelix_pi'] | ssSequenceFeatures['I'])

        intersects['GOR4AlphaHelix_HGI'] = (
            aaSequenceFeatures['GOR4AlphaHelix'] | ssSequenceFeatures['H'] |
            ssSequenceFeatures['G'] | ssSequenceFeatures['I'])

        intersects['Betastrand_E'] = (
            aaSequenceFeatures['Betastrand'] | ssSequenceFeatures['E'])

        intersects['GOR4BetaStrand_E'] = (
            aaSequenceFeatures['GOR4BetaStrand'] | ssSequenceFeatures['E'])

        if printAll:
            CalculateOverlap.calculateFraction(aaSequenceFeatures,
                                               ssSequenceFeatures, intersects,
                                               print_=True)

        return aaSequenceFeatures, ssSequenceFeatures, intersects

    @staticmethod
    def calculateFraction(aaSequenceFeatureDict, ssSequenceFeatureDict,
                          intersectDict, print_=False):
        """
        Calculate the fraction of covered residues in features.

        @param aaSequenceFeatureDict: A C{dict} which for each feature contains
            the offsets that are covered by it.
        @param ssSequenceFeatureDict: A C{dict} which for each secondary
            structure feature contains the offsets that are covered by it.
        @param intersectDict: A C{dict} which contains the intersects of
            secondary structure features and the features from our finders.
        @param print_: A C{bool} indicating whether the results of the
            calculation should be printed.
        """
        try:
            alphaHelix = (len(aaSequenceFeatureDict['AlphaHelix']) /
                          len(intersectDict['AlphaHelix_H']))
        except ZeroDivisionError:
            alphaHelix = 0.0
        try:
            alphaHelix_3_10 = (len(aaSequenceFeatureDict['AlphaHelix_3_10']) /
                               len(intersectDict['AlphaHelix_3_10_G']))
        except ZeroDivisionError:
            alphaHelix_3_10 = 0.0
        try:
            alphaHelix_pi = (len(aaSequenceFeatureDict['AlphaHelix_pi']) /
                             len(intersectDict['AlphaHelix_pi_I']))
        except ZeroDivisionError:
            alphaHelix_pi = 0.0
        try:
            gor4AlphaHelix = (len(aaSequenceFeatureDict['GOR4AlphaHelix']) /
                              len(intersectDict['GOR4AlphaHelix_HGI']))
        except ZeroDivisionError:
            gor4AlphaHelix = 0.0
        try:
            betaStrand = (len(aaSequenceFeatureDict['Betastrand']) /
                          len(intersectDict['Betastrand_E']))
        except ZeroDivisionError:
            betaStrand = 0.0
        try:
            gor4BetaStrand = (len(aaSequenceFeatureDict['GOR4BetaStrand']) /
                              len(intersectDict['GOR4BetaStrand_E']))
        except ZeroDivisionError:
            gor4BetaStrand = 0.0

        if print_:
            print('Overlap in percent:\n'
                  'AlphaHelix with H_AlphaHelix: %f\n'
                  'AlphaHelix_3_10 with G_AlphaHelix_3_10: %f\n'
                  'AlphaHelix_pi with I_AlphaHelix_pi: %f\n'
                  'GOR4AlphaHelix with H_AlphaHelix, G_AlphaHelix_3_10, '
                  'I_AlphaHelix_pi: %f\n'
                  'BetaStrand with E_BetaStrand: %f\n'
                  'GOR4BetaStrand E_BetaStrand: %f\n' %
                  (alphaHelix, alphaHelix_3_10, alphaHelix_pi, gor4AlphaHelix,
                   betaStrand, gor4BetaStrand))
            print('Number of covered residues per finder:\n'
                  'AlphaHelix: %i\nAlphaHelix_3_10: %i\nAlphaHelix_pi: %i\n'
                  'GOR4AlphaHelix: %i\nBetaStrand: %i\nGOR4BetaStrand: %i\n'
                  'H_AlphaHelix: %i\nG_AlphaHelix_3_10: %i\n'
                  'I_AlphaHelix_pi: %i\n'
                  'E_BetaStrand: %i' % (
                      len(aaSequenceFeatureDict['AlphaHelix']),
                      len(aaSequenceFeatureDict['AlphaHelix_3_10']),
                      len(aaSequenceFeatureDict['AlphaHelix_pi']),
                      len(aaSequenceFeatureDict['GOR4AlphaHelix']),
                      len(aaSequenceFeatureDict['BetaStrand']),
                      len(aaSequenceFeatureDict['GOR4BetaStrand']),
                      len(aaSequenceFeatureDict['H']),
                      len(aaSequenceFeatureDict['G']),
                      len(aaSequenceFeatureDict['I']),
                      len(aaSequenceFeatureDict['E'])))

        return (alphaHelix, alphaHelix_pi, alphaHelix_3_10, gor4AlphaHelix,
                betaStrand, gor4BetaStrand)
