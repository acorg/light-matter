from dark.reads import AARead
from dark.fasta import FastaReads

from light.database import DatabaseSpecifier


class HashesString(object):

    def __init__(self, sequences, cutoff, **kwargs):
        """
        A class to work with hashes.
        For a set of given sequences, find all hashes and for each sequence
        make a string of 1 or 0 denoting whether a hash is present in that
        sequence or not. Only include hashes if they occur in more than at '
        least a specified percentage of all given sequences.

        @param sequences: A C{str} filename with a fasta file of sequences to
            be used or a C{dark.reads.Reads} object.
        @param cutoff: A C{float} between 0.0 and 1.0 of the fraction of
            sequences in which a hash has to be present to be included in the
            final string.
        @param kwargs: See
            C{database.DatabaseSpecifier.getDatabaseFromKeywords} for
            additional keywords, all of which are optional.
        """
        if isinstance(sequences, str):
            reads = FastaReads(sequences, readClass=AARead)
        else:
            reads = sequences

        database = DatabaseSpecifier().getDatabaseFromKeywords(**kwargs)

        # Make a dictionary where the keys are the sequence ids and the value
        # is an orderedDict of hashes as returned from getHashes().
        hashes = {}
        for read in reads:
            scannedRead = database.scan(read)
            readHashes = database.getHashes(scannedRead)
            hashes[read.id] = readHashes

        # Make a list of all unique hashes that occur.
        totalHashesList = []
        for read in hashes:
            totalHashesList += list(hashes[read].keys())
        totalHashes = set(totalHashesList)

        # Make a dictionary where the key is a hash and the value is a list of
        # the reads in which the hash occurs.
        byHashes = {}
        for hash_ in totalHashes:
            viruses = []
            for readId in hashes:
                try:
                    hashes[readId][hash_]
                except KeyError:
                    continue
                viruses.append(readId)
            byHashes[hash_] = viruses

        # Make a dictionary where the key is a readId and the value is a string
        # of 1 and 0 denoting which hashes occur in which read.
        co = cutoff * len(reads)
        self.hashString = {read.id: '' for read in reads}

        for hash_ in byHashes:
            if len(byHashes[hash_]) > co:
                for virus in self.hashString:
                    if virus in byHashes[hash_]:
                        self.hashString[virus] += '1'
                    else:
                        self.hashString[virus] += '0'

    def print_(self):
        """
        Print the string in self.hashString.
        """
        for virus in self.hashString:
            print('%s\t%s' % (virus, self.hashString[virus]))
