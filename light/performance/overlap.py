import six
from collections import defaultdict
from Bio import SeqIO

from light.backend import Backend
from light.database import DatabaseSpecifier
from light.landmarks import ALL_LANDMARK_CLASSES
from light.trig import ALL_TRIG_CLASSES

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
        if six.PY3:
            super().__init__(id, sequence)
        else:
            AARead.__init__(self, id, sequence)
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

    def calculateOverlap(self):
        """
        Calculate the feature overlap for a set of features.

        """
        allSequenceFeatures = defaultdict(list)
        allCommons = defaultdict(list)
        allTotals = defaultdict(list)

        for i, read in enumerate(self.SSAAReads):
            sequenceFeatures, commons, totals = self.getFeatures(
                read)
            for feature, indices in sequenceFeatures.items():
                allSequenceFeatures[feature].extend(indices)
            for name, indices in commons.items():
                allCommons[name].extend(indices)
            for name, indices in totals.items():
                allTotals[name].extend(indices)

        return allSequenceFeatures, allCommons, allTotals

    @staticmethod
    def getFeatures(ssAARead, **kwargs):
        """
        Extract the features from the sequence and the structural information.
        For each finder return the offsets covered by that finder and for
        each finder combination, return the offsets in common and the offsets
        in total.

        @param ssAARead: An C{SSAARead} instance.
        @param kwargs: See
            C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
            additional keywords, all of which are optional.
        """
        names = ['AlphaHelix', 'AlphaHelix_3_10', 'AlphaHelix_pi',
                 'AminoAcidsLm', 'BetaStrand', 'BetaTurn',
                 'GOR4AlphaHelix', 'GOR4BetaStrand', 'GOR4Coil', 'Prosite',
                 'AminoAcids', 'IndividualPeaks', 'IndividualTroughs', 'Peaks',
                 'Troughs', 'H', 'G', 'I', 'E']
        if 'landmarkNames' not in kwargs:
            kwargs['landmarkNames'] = [c.NAME for c in ALL_LANDMARK_CLASSES]
        if 'trigPointNames' not in kwargs:
            kwargs['trigPointNames'] = [c.NAME for c in ALL_TRIG_CLASSES]

        db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        backend = Backend()
        backend.configure(db.params)

        sequenceFeatures = defaultdict(set)
        totals = defaultdict(set)
        commons = defaultdict(set)

        scannedAaSequence = backend.scan(ssAARead)

        # Get all offsets for each landmark and trig point separately.
        for landmark in scannedAaSequence.landmarks:
            sequenceFeatures[landmark.name].update(landmark.coveredOffsets())
        for trigPoint in scannedAaSequence.trigPoints:
            sequenceFeatures[trigPoint.name].update(trigPoint.coveredOffsets())

        # Get all offsets for each secondary structure feature separately.
        for offset, structure in enumerate(ssAARead.structure):
            sequenceFeatures[structure].add(offset)

        # Get the overlap between all features
        seen = set()
        for i, name1 in enumerate(names):
            for j, name2 in enumerate(names):
                if name2 not in seen:
                    common = sequenceFeatures[name1] & sequenceFeatures[name2]
                    total = sequenceFeatures[name1] | sequenceFeatures[name2]

                    name = name1 + '-' + name2
                    totals[name] = total
                    commons[name] = common
                    seen.add(name1)

        return sequenceFeatures, commons, totals
