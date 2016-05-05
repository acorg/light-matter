This directory contais files for running jobs using HTCondor on albertine.

`write-htcondor-job-spec.py` writes a jobspec file for HTCondor and splits a fasta file into smaller files containing a specified number of sequences.

`process.sh` calls the function specified in executableName with the given arguments.

`redo.sh` writes a jobspec file for HTCondor to process multiple FASTA input files via evaluate-helices.py, runs these jobs, and removes the one-time spec file it wrote.

`finalize.sh` removes zero-length error files.
