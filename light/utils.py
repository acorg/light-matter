HYDROPHOBIC = ('F', 'Y', 'W', 'H', 'K', 'T', 'C', 'G', 'A', 'V', 'I', 'L', 'M')
HYDROPHILIC = ('P', 'S', 'N', 'D', 'Q', 'E', 'R')


def convertAAToAAProperties(sequence, properties):
    """
    Takes an amino acid sequence, converts it to a sequence of properties,
    e.g. hydrophobic / hydrophilic.

    @param sequence: an amino acid sequence.
    @param properties: a list of two properties
    """
    assert len(properties) == 2, ('Two properties must be specified')

    propertiesSequence = ''
    for aa in sequence:
        if aa in properties[0]:
            propertiesSequence += 'I'
        elif aa in properties[1]:
            propertiesSequence += 'O'

    return propertiesSequence
