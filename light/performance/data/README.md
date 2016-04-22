# Files in this directory

## `polymerase.fasta`

A set of 22 polymerase sequences and their structures (from PDB), as used
in Cerny et al., 2014, 'Evolution of Tertiary Structure of Viral RNA
Dependent Polymerases.'

## `polymerase.json`

Is a JSON file with an object that holds pairwise BLAST bit scores for the
sequences (in `polymerase.fasta`) taken from Table 2 in the paper mentioned
above.

It is created automatically by the top-level `Makefile` (using
`light/performance/bin/create-polymerase-json.sh`, which calls `blastp` and
converts its output to JSON via
`light/performance/bin/convert-blast-10-to-json.py`).

Do not edit this file manually as your changes will be overwritten.
Arguably, the file should not be checked into version control, but leaving
it out would then mean that light matter users would need to have `blastp`
installed.