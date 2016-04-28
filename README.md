**A collection of python tools for matching sequences based on predicted
secondary structures.**

For definitions of terms, look [here](doc/definitions.md).

For an explanation of the database and find parameters, look
[here](doc/parameters.md).

## Developers

You'll need Python 2 or 3, and the corresponding requirements in
`requirements-2.txt` or `requirements-3.txt`.

To run the performance tools, you'll need BLAST installed (so `blastp` is
available to compute bit scores for comparison).

### Creating needed data files

Run

```sh
$ make performance-data
```
