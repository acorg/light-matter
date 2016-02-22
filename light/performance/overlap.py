from collections import defaultdict

from light.backend import Backend
from light.database import DatabaseSpecifier
from light.landmarks import ALL_LANDMARK_CLASSES
from light.trig import ALL_TRIG_CLASSES


class CalculateOverlap(object):
    """
    Calculate the overlap between the features found by our finders and the
    secondary structures found by DSSP. The secondary structures found by
    DSSP were downloaded from http://www.rcsb.org/pdb/files/ss.txt on
    11/11/2015.

    @param kwargs: See
        C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
        additional keywords, all of which are optional.
    """
    def __init__(self, **kwargs):
        # TODO: The following (up to the point of setting up the database)
        #       will go away if we do
        #       https://github.com/acorg/light-matter/issues/442
        self._names = [
            'AlphaHelix', 'AlphaHelix_3_10', 'AlphaHelix_pi',
            'AminoAcidsLm', 'BetaStrand', 'BetaTurn', 'GOR4AlphaHelix',
            'GOR4BetaStrand', 'GOR4Coil', 'Prosite', 'AminoAcids',
            'IndividualPeaks', 'IndividualTroughs', 'Peaks', 'Troughs',
            'H', 'G', 'I', 'E']
        if 'landmarks' not in kwargs:
            kwargs['landmarks'] = [c.NAME for c in ALL_LANDMARK_CLASSES]
        if 'trigPoints' not in kwargs:
            kwargs['trigPoints'] = ([c.NAME for c in ALL_TRIG_CLASSES if
                                     c.NAME != 'Volume'])

        db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        self._backend = Backend()
        self._backend.configure(db.dbParams)

    def getFeatures(self, ssAARead):
        """
        Extract the features from the sequence. Return information about the
        offsets covered by each feature as well as the intersection and union
        of offsets for each pair of features.

        @param ssAARead: An C{SSAARead} instance.
        @return: A triple of C{defaultdict(set)}s. These contain:
            1) The sequence features, keyed by C{str} feature name, each with
               a C{set} of C{int}s as a value, giving the offsets in
               C{ssAARead} where the feature was found.
            2) The intersection of offsets for each pair of feature finders.
               This is keyed by C{str} "name1-name2" with the names of the
               finders. The values are C{set}s of C{int}s, as in (1).
            3) The union of offsets for each pair of feature finders. This is
               keyed by C{str} "name1-name2" with the names of the finders.
               The values are C{set}s of C{int}s, as in (1).
        """

        features = defaultdict(set)
        intersection = defaultdict(set)
        union = defaultdict(set)

        scannedSequence = self._backend.scan(ssAARead)

        # Get all offsets for each landmark and trig point separately.
        for feature in scannedSequence.landmarks + scannedSequence.trigPoints:
            features[feature.name].update(feature.coveredOffsets())

        # Get all offsets for each secondary structure feature separately.
        #
        # TODO: This will go away if we do
        #       https://github.com/acorg/light-matter/issues/442
        for offset, structure in enumerate(ssAARead.structure):
            features[structure].add(offset)

        # Get the offset intersection and union of all pairs of features.
        for i, name1 in enumerate(self._names):
            for name2 in self._names[i + 1:]:
                key = frozenset((name1, name2))
                intersection[key] = features[name1] & features[name2]
                union[key] = features[name1] | features[name2]

        return features, intersection, union
