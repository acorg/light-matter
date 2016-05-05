#!/usr/bin/env bash

# A script that removes zero-length error files.

# Note that we need bash in order to set the nullglob shell option. That
# prevents an error if there are no *.fasta files.

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
