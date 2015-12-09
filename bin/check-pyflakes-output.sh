#!/bin/sh

# Read the output from pyflakes (as run by make pyflakes or make lint in
# the top-level directory) on stdin and check that it looks ok. Show the
# output and exit non-zero if not (causing make to throw an error).
#
# This is necessary because we support Python 2 and 3 and the latter
# introduces some new keywords (async, await, yield from) that trip up the
# Python 2 version of pyflakes. So we check to see if the pyflakes in our
# PATH can handle the new keywords or not. If not, we check the pyflakes
# output (as received on stdin) to see it's exactly the error we
# expected. If so, we exit 0.
#
# In any other case, if we receive any input, cat it and exit non-zero.

out=/tmp/`basename $0`-out-$$
expected=/tmp/`basename $0`-expected-$$
trap "rm -f $out $expected" 1 2 3

cat > $out

if echo 'def x(): yield from [1, 2, 3]' | pyflakes 2>/dev/null
then
    # Python 3. We shouldn't have any output.
    if test -s $out
    then
        cat $out
        rm $out
        exit 1
    else
        # No output, all's well.
        rm $out
        exit 0
    fi
else
    # Python 2. If the output is exactly what we expect, ignore it and exit 0.
    # Else, cat the output and exit non-zero.
    cat >$expected <<EOF
./light/autobahn/backend.py:58:26: invalid syntax
                yield from self.registerAPIMethods()
                         ^
./light/autobahn/client.py:26:34: invalid syntax
            paramsStr = yield from self.call('parameters')
                                 ^
./light/autobahn/database.py:37:26: invalid syntax
                yield from self._connector.addBackend(sessionId)
                         ^
./light/autobahn/find-matches.py:18:31: invalid syntax
            result = yield from self.call('find', what)
                              ^
./light/autobahn/shutdown.py:25:22: invalid syntax
            yield from self.call('shutdown',
                     ^
./light/connector_wamp.py:120:49: invalid syntax
        name, checksum, subjectCount = yield from self._component.call(
                                                ^
EOF
    if cmp -s $out $expected
    then
        # Output as expected.
        rm -f $out $expected
        exit 0
    else
        # Output not as expected. Display it and exit non-zero. Note that
        # this includes the case where pyflakes for Python 2 has been
        # "fixed" to allow the Python 3 keywords. In that case there will
        # be no output (so it's not an error) but we'll still exit
        # non-zero.  Although we could test for that case, it's better to
        # issue a warning and fail because that will tell us we can remove
        # this ugly hack...
        if test -s $out
        then
            cat $out
        else
            echo "No output from pyflakes under Python 2. Is pyflakes fixed?" >&2
        fi
        rm -f $out $expected
        exit 1
    fi
fi
