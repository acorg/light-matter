from collections import defaultdict

from light.backend import Backend
from light.database import DatabaseSpecifier
from light.landmarks import ALL_LANDMARK_CLASSES
from light.trig import ALL_TRIG_CLASSES

from dark.fasta_ss import SSFastaReads


class CalculateOverlap(object):
    """
    Calculate the overlap between the features found by our finders and the
    secondary structures found by DSSP. The secondary structures found by
    DSSP were downloaded from http://www.rcsb.org/pdb/files/ss.txt on
    11/11/2015 (also: see the NOTE below).

    @param pdbFile: The C{str} filename of the file from PDB containing the
        sequence and structural information. See the file C{dark/fasta_ss.py}
        in the dark matter repository at https://github.com/acorg/dark-matter
        for information about the file format.
    @raise ValueError: If C{pdbFile} has an odd number of FASTA records.
    """
    def __init__(self, pdbFile):
        self._ssAAReads = SSFastaReads(pdbFile)

    def calculateOverlap(self):
        """
        Calculate the feature overlap for a set of reads.

        @return: A triple, with sets of offsets...
        """
        allSequenceFeatures = defaultdict(set)
        allCommons = defaultdict(set)
        allTotals = defaultdict(set)

        for i, read in enumerate(self._ssAAReads):
            sequenceFeatures, commons, totals = self.getFeatures(read)
            for feature, indices in sequenceFeatures.items():
                allSequenceFeatures[feature].update(indices)
            for name, indices in commons.items():
                allCommons[name].update(indices)
            for name, indices in totals.items():
                allTotals[name].update(indices)

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
        if 'landmarks' not in kwargs:
            kwargs['landmarks'] = [c.NAME for c in ALL_LANDMARK_CLASSES]
        if 'trigPoints' not in kwargs:
            kwargs['trigPoints'] = ([c.NAME for c in ALL_TRIG_CLASSES if
                                     c.NAME != 'Volume'])

        db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        backend = Backend()
        backend.configure(db.dbParams)

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

        # Get the overlap between all features.
        seen = set()
        for name1 in names:
            for name2 in names:
                if name2 not in seen:
                    common = sequenceFeatures[name1] & sequenceFeatures[name2]
                    total = sequenceFeatures[name1] | sequenceFeatures[name2]

                    name = name1 + '-' + name2
                    totals[name] = total
                    commons[name] = common
                    seen.add(name1)

        return sequenceFeatures, commons, totals
