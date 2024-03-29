from os.path import dirname, join

from .bitScores import BIT_SCORES

from light.performance import data

from dark.fasta_ss import SSFastaReads
from dark.reads import SSAAReadWithX

DATASET = 'pdb_2hla_a_against_polymerase'

_DIR = join(dirname(data.__file__), DATASET)

QUERIES = list(SSFastaReads(join(_DIR, 'queries.fasta'),
                            readClass=SSAAReadWithX))

SUBJECTS = list(SSFastaReads(join(_DIR, 'subjects.fasta'),
                             readClass=SSAAReadWithX))

_ = BIT_SCORES  # Keep pyflakes quiet.
