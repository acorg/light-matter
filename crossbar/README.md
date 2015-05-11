# Realm specification

The WAMP realm we're using is called `light-matter`. That value appears in
multiple places, all of which *must* agree:

1. In `auth.py`. Note that this is a python 2 file, so we don't assume we
   can import the realm from `light.wamp` (which is python 3).
1. In `.crossbar/config.json`.
1. In the `light-matter` repo in `light/wamp.py`
1. In the `light-matter` repo tests in `test/test_wamp.py`

So unless you have a good reason to change the name of the realm, you'll
probably want to leave it alone.
