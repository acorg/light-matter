#!/usr/bin/env python

"""
See the 'EPILOG' variable below, or (better) run with --help for help.
"""

from __future__ import print_function

import sys

from dark.fasta import FastaReads
from dark.reads import AAReadWithX


DEFAULT_PDB_FILE = '/usr/local/dark-matter/seqs/ss-03032016.txt'
DEFAULT_EXECUTABLE_NAME = ('/usr/local/dark-matter/light-matter/light/'
                           'performance/bin/evaluate-helices.py')
DEFAULT_EMAIL = 'tcj25@cam.ac.uk'
DEFAULT_SEQUENCES_PER_FILE = 100

EPILOG = """Given a FASTA file argument, write out the following:

  1) Files named 0.fasta, 1.fasta, 2.fasta, etc. each containing a maximum
     number of sequences (given by --seqs-per-evaluation).

  2) An HTCondor job spec file 'job.htcondor' that can be given to
     condor_submit to have it run the evaluate_helices.py on all the small
     FASTA files.

NOTE: A 'redo.sh' script that can be used to re-submit sub-FASTA files on which
the initial processing failed, a 'finalize.sh' script that removes empty
error files and a process.sh script which calls the specified executable with
the specified arguments are available in light/performance/bin/htcondor.

NOTE: files are overwritten. It is suggested you run this command in an
empty directory.
"""


def splitFASTA(params):
    """
    Read the FASTA file named params['fastaFile'] and print out its
    sequences into files named 0.fasta, 1.fasta, etc. with
    params['seqsPerJob'] sequences per file.
    """
    assert params['fastaFile'][-1] == 'a', ('You must specify a file in ',
                                            'fasta-format that ends in '
                                            '.fasta')

    fileCount = count = seqCount = 0
    outfp = None

    reads = FastaReads(params['fastaFile'], readClass=AAReadWithX,
                       checkAlphabet=0)
    for read in reads:
        seqCount += 1
        if count == params['seqsPerJob']:
            outfp.close()
            count = 0
        if count == 0:
            outfp = open('%d.fasta' % fileCount, 'w')
            fileCount += 1
        count += 1
        print(read.toString(format_='fasta'), end='', file=outfp)
    outfp.close()
    return fileCount, seqCount


def printJobSpec(params):
    """
    Write out a job spec file for HTCondor to process all the small
    FASTA input files via evaluate_helices.py.
    """
    with open('job.htcondor', 'w') as outfp:
        outfp.write("""\
universe                  = vanilla
executable                = /usr/local/dark-matter/light-matter/light/\
performance/bin/htcondor/process.sh
should_transfer_files     = YES
when_to_transfer_output   = ON_EXIT
notify_user               = %(email)s
max_transfer_input_mb     = -1
max_transfer_output_mb    = -1
log                       = job.log

# Job summary:
#   FASTA input from %(fastaFile)s
#   %(sequenceCount)d sequences split into %(nJobs)d jobs of \
%(seqsPerJob)d sequences each.

arguments                 = $(Process) %(executableName)s %(pdbFile)s \
                                       %(evaluateNoPrefix)s %(structureType)s
input                     = $(Process).fasta
output                    = $(Process).done
error                     = $(Process).error

queue %(nJobs)d
""" % params)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description=('Given a FASTA file with helices and an ss.txt file, '
                     'write an HTCondor job spec for evaluate_helices.py'),
        epilog=EPILOG)
    parser.add_argument(
        'fasta', metavar='FASTA-file',
        help='the FASTA file of sequences to evaluate.')
    parser.add_argument(
        '--seqs-per-evaluation', metavar='N',
        type=int, default=DEFAULT_SEQUENCES_PER_FILE, dest='seqsPerJob',
        help=('the number (>0) of sequences to pass to evaluate_helices.py on '
              'each run.'))
    parser.add_argument(
        '--pdb-file', default=DEFAULT_PDB_FILE, dest='pdbFile',
        help='the structure database to run against.')
    parser.add_argument(
        '--email', metavar='name@host',
        default=DEFAULT_EMAIL, dest='email',
        help='the email address to send the job completed message to.')
    parser.add_argument(
        '--executable-name',
        default=DEFAULT_EXECUTABLE_NAME, dest='executableName',
        help='the name of the executable to run.')
    parser.add_argument(
        '--evaluate-no-prefix', type=bool, default=True,
        dest='evaluateNoPrefix',
        help='Whether prefixes should not be evaluated.')
    parser.add_argument(
        '--structureType', default='H', choices={'H', 'G', 'I', 'E'},
        help=('The type of structure that should be evaluated against. '
              'H: Alpha helix, G: Alpha helix 3 10, I: Alpha helix pi, I: '
              'Extended strand.'))

    args = parser.parse_args()

    if args.seqsPerJob < 1:
        parser.print_help()
        sys.exit(1)

    params = {
        'pdbFile': args.pdbFile,
        'email': args.email,
        'executableName': args.executableName,
        'fastaFile': args.fasta,
        'seqsPerJob': args.seqsPerJob,
        'evaluateNoPrefix': args.evaluateNoPrefix,
        'structureType': args.structureType,
    }
    params['nJobs'], params['sequenceCount'] = splitFASTA(params)
    printJobSpec(params)

    print(('%(sequenceCount)d sequences split into %(nJobs)d jobs of '
           '%(seqsPerJob)d sequences each.' % params))
