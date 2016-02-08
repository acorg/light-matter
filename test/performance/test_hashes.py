from unittest import TestCase

from dark.reads import Reads

from ..sample_data import (MONKEYPOX, MUMMYPOX)
from light.performance.hashes import HashesString


class TestHashesString(TestCase):
    """
    Tests for the light.performance.hashes.HashesString class.
    """

    def testCorrectNumberOfSequences(self):
        """
        The right number of strings must be made (one for each sequence).
        """
        reads = Reads()
        reads.add(MONKEYPOX)
        reads.add(MUMMYPOX)
        hashesString = HashesString(reads, 0.0, landmarks=['AlphaHelix'])
        result = len(hashesString.hashString)
        self.assertEqual(2, result)

    def testCorrectStringWithoutCutoff(self):
        """
        If no cutoff is given, a string with the correct number of '1' must be
        made (note that we can't check the actual string, because the order of
        '1' and '0' is different for each run.
        """
        reads = Reads()
        reads.add(MONKEYPOX)
        reads.add(MUMMYPOX)
        hashesString = HashesString(reads, 0.0, landmarks=['AlphaHelix'],
                                    trigPoints=['AminoAcids'])
        ones = hashesString.hashString['Monkeypox'].count('1')
        self.assertEqual(3, ones)

    def testCorrectStringWithCutoff(self):
        """
        If a cutoff is given, a string with the correct number of '1' must be
        made (note that we can't check the actual string, because the order of
        '1' and '0' is different for each run.
        """
        reads = Reads()
        reads.add(MONKEYPOX)
        reads.add(MUMMYPOX)
        hashesString = HashesString(reads, 0.5, landmarks=['AlphaHelix'],
                                    trigPoints=['AminoAcids'])
        zeroes = hashesString.hashString['Monkeypox'].count('0')
        self.assertEqual(0, zeroes)
