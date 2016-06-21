#!/bin/sh

# A script which calls the function specified in executableName with the given
# arguments.


case $# in
    5)
        jobid=$1
        executableName=$2
        pdbFile=$3
        evaluateNoPrefix=$4
        structureType=$5
        ;;

    *)
        echo "Usage: `basename $0` jobid executableName pdbFile evaluateNoPrefix structureType" >&2
        exit 1
        ;;
esac


DM=/usr/local/dark-matter
export PYTHONPATH=$DM/light-matter/:$DM/dark-matter

errs=$jobid.error

$DM/virtualenv/bin/python $executableName --pdbFile $pdbFile \
--evaluateNoPrefix $evaluateNoPrefix --structureType $structureType \
< $jobid.fasta > $jobid.out 2> $errs

if [ -s $errs ]
then
    echo "Completed WITH ERRORS ($errs) on `hostname` at `date`." > $jobid.done
else
    rm $errs
    echo "Completed on `hostname` at `date`." > $jobid.done
fi
