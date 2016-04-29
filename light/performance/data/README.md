# Files in this directory and its subdirectories

This directory contains one directory for each dataset that we're using for testing.

Each subdirectory contains the following files:

## `queries.fasta`

A fasta file with the query sequences to be used in the tests.

`polymerase`: A set of 22 polymerase sequences and their structures (from PDB), as used in Cerny et al., 2014, 'Evolution of Tertiary Structure of Viral RNA Dependent Polymerases.'

`pdb_2hla_a`: A set of 141 sequences which have structural similarity to 2HLA:A.

`pdb_4mtp_a`: A set of 215 sequences which have structural similarity to 4MTP:A.

## `subjects.fasta`

A fasta file with the subject sequences to be used in the tests.

`polymerase`: A set of 22 polymerase sequences and their structures (from PDB), as used in Cerny et al., 2014, 'Evolution of Tertiary Structure of Viral RNA Dependent Polymerases.'

`pdb_2hla_a`: The 2HLA:A sequence.

`pdb_4mtp_a`: The 4MTP:A sequence.

## `zScores.py`

A file containing a python dictionary which lists the DALI Z scores for all pairwise structure comparisons needed by the tests.

## `Makefile`

In order to run the tests, you also need a bitScores.py file for each dataset. This can be obtained by running `$ make performance-data` at the top level of the tree. This will create a files named `bitScores.py` in each directory which contains all the necessary bit scores. Note that you need BLASTP as well as dark-matter installed to be able to run this.
