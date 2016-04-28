#!/bin/sh -e

# Convert PDB FASTA query and subject files into plain FASTA and BLASTp them.
# Then take the bit scores and use convert-blast-10-to-python.py to write out
# Python with a BIT_SCORES variable that's a dict of dicts (keyed by query id
# then subject id) holding the bit scores.
#
# This script assumes you have BLAST installed and also the dark-matter
# Python package (which contains filter-fasta.py).
#
#   Usage: create-bit-scores-python.sh dataset-name
#
# Where dataset-name is a sub-directory of light/performance/data

tmp=`mktemp -d /tmp/lm-perf-blastp-XXXXXX`
trap "rm -fr $tmp" 1 2 3 15

test -f queries.fasta || {
    echo "$(basename $0): Required input FASTA file 'queries.fasta' does not exist." >&2
    exit 1
}

test -f subjects.fasta || {
    echo "$(basename $0): Required input FASTA file 'subjects.fasta' does not exist." >&2
    exit 1
}

base=$(dirname $(python -c 'import light; print(light.__file__)'))
PATH="$PATH:$base/performance/bin"

# Convert the PDB FASTA into plain FASTA for BLASTing.

set -x # Show commands as we run them.
set +e # Don't exit on error.

filter-fasta.py --readClass fasta-ss --saveAs fasta < queries.fasta > \
    $tmp/queries.fasta 2>$tmp/filter-fasta.errs
if [ $? -ne 0 ]
then
    echo "Filtering query PDB FASTA failed." >&2
    cat $tmp/filter-fasta.errs >&2
    exit 2
fi

filter-fasta.py --readClass fasta-ss --saveAs fasta < subjects.fasta \
    > $tmp/subjects.fasta 2>$tmp/filter-fasta.errs
if [ $? -ne 0 ]
then
    echo "Filtering subject PDB FASTA failed." >&2
    cat $tmp/filter-fasta.errs >&2
    exit 2
fi

# Note that once we cd to $tmp, queries.fasta and subjects.fasta are the
# regular (i.e., non-PDB) FASTA files we just created.
cd $tmp

# Make a BLAST db. Save the stdout and if makeblastdb exits non-zero cat
# the output to stderr. We cannot allow the makeblastdb output to reach
# stdout as our stdout (produced by convert-blast-10-to-json.py below) must
# be JSON.
makeblastdb -in subjects.fasta -dbtype prot -out db-$$ > makeblastdb.out
if [ $? -ne 0 ]
then
    echo "Making BLAST database failed." >&2
    cat makeblastdb.out >&2
    exit 3
fi

set -e # Exit on error.

# Run BLAST and turn the output into Python.
blastp -query queries.fasta -db db-$$ -outfmt '10 qseqid sseqid bitscore' |
  convert-blast-10-to-python.py --queries queries.fasta --subjects subjects.fasta

set +x # Don't show commands.

rm -fr $tmp
