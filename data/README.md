# PDB secondary structures

## PDB known secondary structure file

The file `pdb-20160303-ss.txt.bz2` is produced from a PDB secondary
structure file downloaded on 2016-03-03.

The source file, `ss.txt.gz` can be found
[on this page](http://www.rcsb.org/pdb/static.do?p=download/http/index.html)
or via this direct link:
[http://www.rcsb.org/pdb/files/ss.txt.gz](http://www.rcsb.org/pdb/files/ss.txt.gz). Once
downloaded, clean up the file to be suitable for light matter processing as
follows (this assumes you're using `bash` and running in the light matter
`data` directory and that you have the `light/performance/bin` and
the [dark matter](https://github.com/acorg/dark-matter/) `bin` directory in
your shell's `PATH`)

    $ gzcat ss.txt.gz | clean-pdb-ss-fasta.py | bzip2 > pdb-`date '+%Y%m%d'`-ss.txt.bz2

See the `light/bin/performance/clean-pdb-ss-fasta.py` script for details of
what's done in the clean up.

The PDB `ss.txt.gz` file changes frequently, so we have a specific dated
copy that we work on to ensure consistency. It can be updated periodically,
but in that case the true and false positive counts computed earlier will
be slightly inaccurate.

## PDB alpha helix substrings

The file `pdb-20160303-alpha-helix-substrings.bz2` contains all substrings
of all alpha helices found in the processed `ss.txt.gz` file just
described. Each line has 3 fields:

  * the alpha helix substring
  * the number of true positives (i.e., number of times the string appears
    in alpha helices in PDB)
  * the number of false positives (i.e., number of times the string appears
    in PDB where the amino acids are not all in a helix)

This file was created on multiple cores on a cluster, where each core
processes a subset of the alpha helix substrings. Here's how that was set
up and run (this assumes you have both the `light/performance/bin` and the
[dark matter](https://github.com/acorg/dark-matter/) `bin` directories in
your shell's `PATH` in order for it to find the commands below):

    $ bzcat pdb-20160303-ss.txt.bz2 |
    extract-structures-from-pdb-ss.py --featureType 'PDB AlphaHelix' |
    fasta-subsequences.py --pdb --minLength 4 |
    filter-fasta.py --readClass fasta-ss --saveAs fasta --removeDuplicates >
        pdb-alpha-helix-substrings-no-duplicates.fasta

Then transfer the `pdb-alpha-helix-substrings-no-duplicates.fasta` file to
the cluster and there run:

    $ light/performance/bin/htcondor/write-htcondor-job-spec.py \
        --seqs-per-evaluation 10000 pdb-alpha-helix-substrings-no-duplicates.fasta
    13206964 sequences split into 1321 jobs of 10000 sequences each.

That command creates an HTCondor job file, which you submit for processing:

    $ condor_submit job.htcondor

Once the job has completed, make the final
`pdb-20160303-alpha-helix-substrings.bz2` file from the individual job
`*.out` files via:

    $ cat *.out |
    awk '{printf "%.8f %s\n", $2 / ($2 + $3), $0}' |
    sort -n -r |
    cut -f2- -d' ' |
    bzip2 > pdb-20160303-alpha-helix-substrings.bz2

The output file is sorted by fraction of true positives (out of all
positives), then by sequence and the true and false positive counts.

## Selecting PDB alpha helix substrings for the Aho Corasick alpha helix finder

Given the above file, extract some of the best performing alpha helix
substrings for use by the Aho Corasick alpha helix finder:

    $ bzcat pdb-20160303-alpha-helix-substrings.bz2 |
    alpha-helix-substrings-for-aho-corasick.py --minTruePositives 20 \
        --minTruePositiveFraction 0.75 --printSummary >
        aho-corasick-alpha-helix-substrings-20-0.75
    Kept 183653 of 13206964 (1.39%) candidate substrings seen on input.
    12051251 substrings did not meet the minimum true positive requirement (20).
    283137 substrings did not have a high enough true positive fraction (0.750000).
    688923 substrings were inferior to (at least) one of their own substrings.

If you want to see the true/false positive counts and the true positive /
all positives fraction in the selected subset, run
`alpha-helix-substrings-for-aho-corasick.py` with the `--printCounts`
option.  The Aho Corasick alpha helix substring finder expects just the
alpha helix substrings to be on its input, so a file with counts in it
cannot be used unless the counts are stripped out (via e.g., `awk '{print
$1}'`).

Note that the current repo substring file
`aho-corasick-alpha-helix-substrings-20-0.75` has not been evaluated so we
don't know how well it performs yet.

# Prosite

The current version of the prosite database is 20.119, downloaded from
ftp://ftp.expasy.org/databases/prosite/prosite.dat on Nov 11, 2015.

## Updating Prosite data

The prosite database can be updated by running

    $ make update-prosite

in the top-level light-matter repository.  This will download the latest
prosite database and convert it to our much-reduced JSON representation.

Note that if the prosite database is unchanged, the new database will have
the same version number as the old and the JSON output should be identical,
too.

If a new version is available, it will appear in this (`data`) directory
with a name like `prosite-20.119.json`. You'll then need to do the following.

### Update the db in git

    $ git rm data/prosite-AA.BBB.json  # Where AA.BBB is the old version.
    $ git add data/prosite-CC.DDD.json # Where CC.DDD is the new version.

### Update the version number that is loaded by the finder

Then edit `light/landmarks/prosite.py` and change the version number of
the database in the `_DB_FILE` variable to match the new version.

### Update the number of patterns in the test

In `test/landmarks/test_prosite.py` there is a `testDefaultDatabase` test
that knows the number of patterns in the prosite database. You will
probably need to change this if the database has changed.

### Run the tests

The test suite should now pass cleanly:

    $ make check

# aho-corasick-alpha-helix-prefixes-*

These files contain alpha helix prefixes to be used by the `AC_AlphaHelix`
finder. The number at the end of the filename corresponds to the cut-off
used to determine which helices should be in the file. E.g., all helices in
`aho-corasick-alpha-helix-prefixes-1` have a true positive to false
positive ratio of 1 or more.

# PDB structures by individual year and cumulatively

The file `pdb-structures-by-year.txt` contains a listing of which PDB id
entered PDB according to year. The first field of each line is a year and
subsequent fields are space-separated PDB sequence ids.

## PDB by year

This file can be processed (by
`../light/performance/bin/split-pdb-ss-by-category.py`) to produce
individual FASTA files giving the sequences that entered PDB by year.

E.g.,

```sh
$ bzcat pdb-20160303-ss.txt.bz2 |
  ../light/performance/bin/split-pdb-ss-by-category.py --prefix pdb- \
  --categories pdb-structures-by-year.txt
```

Which will produce files named `pdb-1976.fasta`, `pdb-1977.fasta`,
`pdb-1978.fasta`, etc.

Note that when you do the above, 445 sequences in `pdb-20160303-ss.txt.bz2`
will not be in any category. We have not investigated all of these, but in
one case here's what is happening. `5DTG` entered PDB and its sequence and
structure are still in the PDB `ss.txt` file.  But at some point PDB
realized it was the same as `5HCU` or they renamed it to `5HCU` for some
reason. As a result, there is only an entry in `pdb-structures-by-year.txt`
for `5HCU`. So when we split `pdb-20160303-ss.txt.bz2` by year we drop
`5HCU`.

## PDB by year, cumulatively

The `pdb-structures-by-year.txt` file can also be used to produce
cumulative FASTA containing all PDB sequences over time. This is done by
first making a cumulative category file:

```sh
$ ../light/performance/bin/convert-pdb-id-by-category-to-cumulative.py < \
    pdb-structures-by-year.txt > pdb-structures-by-year-cumulative.txt
```

and then using `../light/performance/bin/split-pdb-ss-by-category.py` (just
as above) to produce FASTA files with cumulative sets of sequences:

```sh
$ bzcat pdb-20160303-ss.txt.bz2 |
  ../light/performance/bin/split-pdb-ss-by-category.py --prefix pdb-cumulative- \
  --categories pdb-structures-by-year-cumulative.txt
```

which will produce files named `pdb-cumulative-1976-1976.fasta`,
`pdb-cumulative-1976-1977.fasta`, `pdb-cumulative-1976-1978.fasta`, etc.
