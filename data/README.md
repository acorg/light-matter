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

These files contain alpha helix prefixes to be used by the `AC_AlphaHelix` finder. The number at the end of the filename corresponds to the cut-off used to determine which helices should be in the file. E.g., all helices in aho-corasick-alpha-helix-prefixes-1 have a true positive to false positive ratio of 1 or more.
