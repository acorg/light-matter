from collections import defaultdict

from light.backend import Backend
from light.database import DatabaseSpecifier
from light.landmarks import ALL_LANDMARK_CLASSES, DEV_LANDMARK_CLASSES
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
        # Set default landmark and trig point finders.
        if 'landmarks' not in kwargs:
            kwargs['landmarks'] = ALL_LANDMARK_CLASSES + [
                c for c in DEV_LANDMARK_CLASSES if
                c.NAME.startswith('PDB ')]
        if 'trigPoints' not in kwargs:
            kwargs['trigPoints'] = [c for c in ALL_TRIG_CLASSES if
                                    c.NAME != 'Volume']

        db = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)
        self._backend = Backend()
        self._backend.configure(db.dbParams)

        self._names = (db.dbParams.landmarkFinderNames() +
                       db.dbParams.trigPointFinderNames())

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

        # Get the offset intersection and union of all pairs of features.
        for i, name1 in enumerate(self._names):
            for name2 in self._names[i + 1:]:
                key = frozenset((name1, name2))
                intersection[key] = features[name1] & features[name2]
                union[key] = features[name1] | features[name2]

        return features, intersection, union
