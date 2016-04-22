This directory contains everything to do with performance testing.

Files and directories are as follows:

* `*.py`: Python used in performance evaluation and testing.
* `bin`: Scripts for running performance testing and creating data for testing.
* `data`: Input data files for performance tests.
* `database`: Historical input FASTA data for performance tests (to be removed once
    `bin/performance-test.py` has its tests either removed or added to `test`).
* `read`: Historical input FASTA data for performance tests (to be removed once
    `bin/performance-test.py` has its tests either removed or added to `test`).
* `test`: Performance tests, in files whose names start with `perf_`.
