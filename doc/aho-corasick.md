It's possible to save RAM used by the
[Aho Corasick implementation](https://github.com/WojciechMula/pyahocorasick)
by compiling it so it uses `bytes` not `unicode` to store the characters at
each trie node.

This can be done by commenting out the `('AHOCORASICK_UNICODE', '')` line
in the packages's `setup.py` file (currently
[here](https://github.com/WojciechMula/pyahocorasick/blob/master/setup.py#L25)).

If you want to take that approach, `git clone` the repo, manually edit
`setup.py` as above, and then run `python setup.py install` to install the
Aho Corasick package manually.

Make sure you do this *after* you've run `pip install -r requirements.txt`
(to install all light matter requirements) so you don't accidentally
re-install a unicode version of the module over the bytes version you
manually installed.
