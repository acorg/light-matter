# coding=utf-8

from itertools import chain

from light.features import Landmark
from light.finder import Finder


class PredictedStructure(Finder):
    """
    A class for computing landmarks based on predicted secondary structures.
    The predicted secondary structure sequences were downloaded from
    http://www.rcsb.org/pdb/files/ss.txt on the 11/11/2015. The predictions are
    made using the DSSP algorithm.
    This finder is intended to test the performance of our algorithm if perfect
    secondary structure information is available to be used as landmarks.
    """
    NAME = 'PredictedStructure'
    SYMBOL = 'ST'

    def getStructureOffsets(self, ssSequence):
        """
        Takes a structure sequence, and returns a dictionary that contains the
        offsets for each structure.
        Possible structures are:
        H = α-helix
        B = residue in isolated β-bridge (part of a strand)
        E = extended strand, participates in β ladder (set of one or more
                consecutive bridges of identical type, continuous stretches of
                E are beta strands)
        G = 3-helix (Alpha helix 3 10)
        I = 5 helix (Alpha helix pi)
        T = hydrogen bonded turn (pieces of helix too short to be helix or
                turn)
        S = bend (region of high curvature)

        @param ssSequence: A C{str} structure sequence.

        @return: A C{dict} which for each structure contains a list of lists,
            where each of the lists has the start and stop offset of
            consecutive sets of structures.
        """
        previous = None
        for offset, symbol in enumerate(chain(ssSequence, [object()])):
            if previous is None:  # Start of string.
                previous = symbol
                startOffset = 0
            elif symbol != previous:
                yield (previous, startOffset, offset)
                previous = symbol
                startOffset = offset

    def find(self, read, structureNames=None):
        """
        A function that returns landmarks based on known secondary structures.

        @param read: An instance of C{dark.reads.SSAARead}.
        @param structureNames: A C{set} of structure names that should be
            considered.
        """
        structureNames = structureNames or {'H', 'G', 'I', 'E'}

        for symbol, start, end in self.getStructureOffsets(read.structure):
            if symbol in structureNames:
                yield Landmark(symbol, self.SYMBOL, start, end - start)
