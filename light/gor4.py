from gor4 import GOR4

_GOR4 = None


def predictions(sequence):
    """
    Takes an amino acid sequence and returns the secondary structure
    predictions using the GOR4 algorithm.

    @param sequence: A C{str} amino acid sequence.

    @return: A C{str} sequence of predicted secondary structures.
    """
    global _GOR4
    if _GOR4 is None:
        _GOR4 = GOR4()
    return _GOR4.predict(sequence)['predictions']
