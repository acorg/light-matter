#!/bin/sh -e

# Convert the PDB FASTA polymerase file into plain FASTA and BLASTp it
# against itself.
#
# This script assumes you have BLAST installed and also the dark-matter
# Python package (which contains filter-fasta.py).

tmp=`mktemp -d /tmp/lm-XXXXXX`
trap "rm -fr $tmp" 1 2 3 15

base=$(dirname $(python -c 'import light; print(light.__file__)'))
pdbFASTA=$base/performance/data/polymerase.fasta

test -f $fasta || {
    echo "$(basename $0): Required input FASTA file '$fasta' does not exist." >&2
    exit 1
}

PATH="$PATH:$base/performance/bin"

cd $tmp

# Convert the PDB FASTA into plain FASTA for BLASTing.
FASTA=polymerase.fasta

set -x # Show commands as we run them.
set +e # Don't exit on error.

filter-fasta.py --readClass fasta-ss --saveAs fasta < $pdbFASTA > $FASTA 2>filter-fasta.errs
if [ $? -ne 0 ]
then
    echo "Filtering PDB FASTA failed." >&2
    cat filter-fasta.errs >&2
    exit 2
fi

# Make a BLAST db. Save the stdout and if makeblastdb exits non-zero cat
# the output to stderr. We cannot allow the makeblastdb output to reach
# stdout as our stdout (produced by convert-blast-10-to-json.py below) must
# be JSON.
makeblastdb -in $FASTA -dbtype prot -out polymerase > makeblastdb.out
if [ $? -ne 0 ]
then
    echo "Making BLAST database failed." >&2
    cat makeblastdb.out >&2
    exit 3
fi

set -e # Exit on error.

# Run BLAST and turn the output into JSON.
blastp -query $FASTA -db polymerase -outfmt '10 qseqid sseqid bitscore' |
  convert-blast-10-to-json.py --fasta $FASTA

set +x # Don't show commands.

rm -fr $tmp
