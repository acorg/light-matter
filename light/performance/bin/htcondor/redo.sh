#!/bin/sh -e

# A shell script that writes a job spec file for HTCondor to process multiple
# FASTA input files via evaluate-helices.py, runs these jobs, and removes the
# one-time spec file it wrote.

if [ $# -lt 4 ]
then
    echo "Usage: `basename $0` executableName pdbFile evaluateNoPrefix jobid" >&2
else
    executableName=$1
    pdbFile=$2
    evaluateNoPrefix=$3
    shift 3
    jobids="$@"
fi

tmp=redo.tmp.$$
trap "rm -f $tmp" 0 1 2 3 15

    cat >$tmp <<EOF
universe                  = vanilla
executable                = /usr/local/dark-matter/light-matter/light/\
                            performance/bin/htcondor/process.sh
should_transfer_files     = YES
when_to_transfer_output   = ON_EXIT
notify_user               = tcj25@cam.ac.uk
max_transfer_input_mb     = -1
max_transfer_output_mb    = -1
log                       = job.log
EOF

for jobid in $jobids
do
    cat >>$tmp <<EOF

arguments                 = $jobid $1 $2 $3
input                     = $jobid.fasta
output                    = $jobid.done
error                     = $jobid.error

queue
EOF
    rm -f $jobid.error $jobid.done
done

rm -f job.log

condor_submit $tmp
