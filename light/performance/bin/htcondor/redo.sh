#!/bin/sh -e

# A shell script that writes a job spec file for HTCondor to process a single
# FASTA input file via evaluate_helices.py, runs that job, and removes the
# one-time spec file it wrote.

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
notify_user               = tcj25@cam.ac.uk
max_transfer_input_mb     = -1
max_transfer_output_mb    = -1
log                       = job.log
EOF

for jobid in "$@"
do
    cat >>$tmp <<EOF

arguments                 = $jobid
input                     = $jobid.fasta
output                    = $jobid.done
error                     = $jobid.error

queue
EOF
    rm -f $jobid.error $jobid.done
done

rm -f job.log

condor_submit $tmp
