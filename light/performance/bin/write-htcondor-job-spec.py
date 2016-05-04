#!/usr/bin/env python

"""
See the 'EPILOG' variable below, or (better) run with --help for help.
"""

from __future__ import print_function

import os
import sys

from dark.fasta import FastaReads
from dark.reads import AAReadWithX


DEFAULT_STRUCTURE_DB = '/usr/local/dark-matter/seqs/ss-03032016.txt'
DEFAULT_EXECUTABLE_NAME = ('/usr/local/dark-matter/light-matter/light/'
                           'performance/bin/evaluate_helices.py')
DEFAULT_EMAIL = 'tcj25@cam.ac.uk'
DEFAULT_SEQUENCES_PER_FILE = 100

EPILOG = """Given a FASTA file argument, write out the following:

  1) Files named 0.fasta, 1.fasta, 2.fasta, etc. each containing a maximum
     number of sequences (given by --seqs-per-evaluation).

  2) An HTCondor job spec file 'job.htcondor' that can be given to
     condor_submit to have it run the evaluate_helices.py on all the small
     FASTA files.

  3) A 'redo.sh' script that can be used to re-submit sub-FASTA files on which
     the initial processing failed.

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

    reads = FastaReads(params['fastaFile'], readClass=AAReadWithX)
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
executable                = process.sh
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

arguments                 = $(Process) %(db)s
input                     = $(Process).fasta
output                    = $(Process).done
error                     = $(Process).error

queue %(nJobs)d
""" % params)


def printRedoScript(params):
    """
    Write out a shell script that writes a job spec file for HTCondor to
    process a single FASTA input file via evaluate_helices.py, runs that job,
    and removes the one-time spec file it wrote.
    """
    with open('redo.sh', 'w') as outfp:
        outfp.write("""\
#!/bin/sh -e

case $# in
    0) echo "Usage: `basename $0` jobid1, jobid2, ..." >&2; exit 1;;
esac

tmp=redo.tmp.$$
trap "rm -f $tmp" 0 1 2 3 15

    cat >$tmp <<EOF
universe                  = vanilla
executable                = process.sh
should_transfer_files     = YES
when_to_transfer_output   = ON_EXIT
notify_user               = %(email)s
max_transfer_input_mb     = -1
max_transfer_output_mb    = -1
log                       = job.log
EOF

for jobid in "$@"
do
    cat >>$tmp <<EOF

arguments                 = $jobid %(db)s
input                     = $jobid.fasta
output                    = $jobid.done
error                     = $jobid.error

queue
EOF
    rm -f $jobid.error $jobid.done
done

rm -f job.log

condor_submit $tmp
""" % params)

    # Make the script executable so we can run it.
    os.chmod('redo.sh', 0o755)


def printProcessScript(params):
    """
    Write out a simple process script to call evaluate_helices.py.
    """
    with open('process.sh', 'w') as outfp:
        outfp.write("""\
#!/bin/sh

DM=/usr/local/dark-matter
export PYTHONPATH=$DM/light-matter/

jobid=$1
shift

errs=$jobid.error

%(executableName)s $jobid.fasta %(db)s $jobid.out 2> $errs

if [ -s $errs ]
then
    echo "Completed WITH ERRORS ($errs) on `hostname` at `date`." > $jobid.done
else
    rm $errs
    echo "Completed on `hostname` at `date`." > $jobid.done
fi
""" % params)

    # Make the script executable so we can run it.
    os.chmod('process.sh', 0o755)


def printFinalizeScript(params):
    """
    Write out a script that removes zero-length error files.

    Note that we need bash in order to set the nullglob shell option. That
    prevents an error if there are no *.fasta files.
    """
    with open('finalize.sh', 'w') as outfp:
        outfp.write("""\
#!/usr/bin/env bash

shopt -s nullglob

# redo will hold the numbers of jobs that need re-running, if any.
redo=

for i in *.fasta
do
    n=`echo $i | cut -f1 -d.`
    error=$n.error
    result=$n.out

    if [ -f $error ]
    then
        if [ -s $error ]
        then
            echo "WARNING: $error is non-empty." >&2
        else
            rm -f $error
        fi
    fi

done
""" % params)

    # Make the script executable so we can run it.
    os.chmod('finalize.sh', 0o755)


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
        '--structure-db-name', metavar='structure-database-name',
        type=str, default=DEFAULT_STRUCTURE_DB, dest='db',
        help='the structure database to run against.')
    parser.add_argument(
        '--email', metavar='name@host',
        type=str, default=DEFAULT_EMAIL, dest='email',
        help='the email address to send the job completed message to.')
    parser.add_argument(
        '--executable-name',
        type=str, default=DEFAULT_EXECUTABLE_NAME, dest='executableName',
        help='the name of the executable to run.')

    args = parser.parse_args()

    if args.seqsPerJob < 1:
        parser.print_help()
        sys.exit(1)

    params = {
        'db': args.db,
        'email': args.email,
        'executableName': args.executableName,
        'fastaFile': args.fasta,
        'seqsPerJob': args.seqsPerJob,
    }
    params['nJobs'], params['sequenceCount'] = splitFASTA(params)
    printJobSpec(params)
    printProcessScript(params)
    printRedoScript(params)
    printFinalizeScript(params)

    print(('%(sequenceCount)d sequences split into %(nJobs)d jobs of '
           '%(seqsPerJob)d sequences each.' % params))
