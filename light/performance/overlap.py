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
            read = SSAARead(record.id, str(record.seq),
                            str(records[i + 1].seq))
            self.SSAAReads.add(read)

    def calculateOverlap(self, printAll=False):
        """
        Calculate the feature overlap for a set of features.

        @param printAll: A C{bool} saying whether the information about the
            overlap in each read should be printed out.
        """
        allSequenceFeatures = defaultdict(list)
        allCommons = defaultdict(list)
        allTotals = defaultdict(list)

        for read in self.SSAAReads:
            sequenceFeatures, commons, totals = self.getFeatures(
                read, printAll=printAll)
            for feature, indices in sequenceFeatures.items():
                allSequenceFeatures[feature].extend(indices)
            for name, indices in commons.items():
                allCommons[name].extend(indices)
            for name, indices in totals.items():
                allTotals[name].extend(indices)

        return allSequenceFeatures, allCommons, allTotals

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
        names = ['AlphaHelix', 'AlphaHelix_3_10', 'AlphaHelix_pi',
                 'AminoAcidsLm', 'BetaStrand', 'BetaTurn',
                 'GOR4AlphaHelix', 'GOR4BetaStrand', 'GOR4Coil', 'Prosite',
                 'AminoAcids', 'IndividualPeaks', 'IndividualTroughs', 'Peaks',
                 'Troughs', 'H', 'G', 'I', 'E']
        kwargs['landmarkNames'] = ['AlphaHelix', 'AlphaHelix_3_10',
                                   'AlphaHelix_pi', 'AminoAcidsLm',
                                   'BetaStrand', 'BetaTurn', 'GOR4AlphaHelix',
                                   'GOR4BetaStrand', 'GOR4Coil', 'Prosite']
        kwargs['trigPointNames'] = ['AminoAcids', 'IndividualPeaks',
                                    'IndividualTroughs', 'Peaks', 'Troughs']
        db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        backend = Backend()
        backend.configure(db.params)

        sequenceFeatures = defaultdict(set)
        totals = defaultdict(set)
        commons = defaultdict(set)

        scannedAaSequence = backend.scan(ssAARead)

        # Get all offsets for each landmark separately.
        for landmark in scannedAaSequence.landmarks:
            sequenceFeatures[landmark.name].update(landmark.coveredOffsets())

        # Get all offsets for each secondary structure feature separately.
        for offset, structure in enumerate(ssAARead.structure):
            sequenceFeatures[structure].add(offset)

        for i, name1 in enumerate(names):
            for j, name2 in enumerate(names):
                if i == j:
                    pass
                else:
                    common = sequenceFeatures[name1] & sequenceFeatures[name2]
                    total = sequenceFeatures[name1] | sequenceFeatures[name2]

                    name = name1 + name2
                    totals[name] = total
                    commons[name] = common

        return sequenceFeatures, commons, totals
