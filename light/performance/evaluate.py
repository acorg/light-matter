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
