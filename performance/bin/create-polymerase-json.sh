#!/bin/sh

tmp=`mktemp -d /tmp/lm-XXXXXX`
trap "rm -fr $tmp" 1 2 3 15

base=$(dirname $(dirname $(python -c 'import light; print(light.__file__)')))
fasta=$base/performance/database/polymerase-db.fasta

test -f $fasta || {
    echo "$(basename $0): Required input FASTA file '$fasta' does not exist." >&2
    exit 1
}

PATH="$PATH:$base/performance/bin"

cd $tmp

# Make a BLAST db. Save the stdout and if makeblastdb exits non-zero cat
# the output to stderr.
makeblastdb -in $fasta -dbtype prot -out polymerase > $tmp/makeblastdb.out
if [ $? -ne 0 ]
then
    echo "Making BLAST database failed." >&2
    cat $tmp/makeblastdb.out >&2
    exit 2
fi

# Run BLAST and turn the output into JSON.
blastp -query $fasta -db polymerase -outfmt '10 qseqid sseqid bitscore' |
  convert-blast-10-to-json.py --fasta $fasta

rm -fr $tmp
