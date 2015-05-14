#
# Provide a Python2 dynamic authentication service for WAMP.
#
# All our other WAMP components are written in Python3 and use asyncio. For
# now though, crossbar.io only runs under Twisted in Python2 (because
# Twisted is not yet quite Python3 ready). The authentication component
# must be internal to the WAMP router and so must also be in Python2. When
# Twisted works on Python3 (any day now), it will be possible to run
# crossbar.io in Python3 also and we can adjust the code below to use
# asyncio like the rest of our WAMP components.
#

from __future__ import print_function

import sys

try:
    from twisted.internet.defer import inlineCallbacks
except ImportError as e:
    print('Could not import Twisted (%s). Are you in the right virtualenv?' %
          e, file=sys.stderr)
    sys.exit(1)

from autobahn.twisted.wamp import ApplicationSession
from autobahn.wamp.exception import ApplicationError


class AppSession(ApplicationSession):

    # Note that the realm CANNOT be changed unless you make corresponding
    # changes elsewhere. See the explanation in README.md
    REALM = 'light-matter'

    # TODO: Get the secrets from the environment.
    USER_DB = {
        'frontend': {
            'role': 'frontend',
            'secret': 'secret2',
        },
        'backend': {
            'role': 'backend',
            'secret': 'secret2',
        },
        'shutdown': {
            'role': 'shutdown',
            'secret': 'secret2',
        },
        'find': {
            'role': 'find',
            'secret': 'secret2',
        },
    }

    @inlineCallbacks
    def onJoin(self, details):

        def authenticate(realm, authid, details):
            if realm != self.REALM:
                raise ApplicationError('unknown_realm',
                                       'Unknown realm %r' % realm)

            try:
                return self.USER_DB[authid]
            except KeyError:
                raise ApplicationError('unknown_user',
                                       'Unknown user %r' % authid)

        try:
            yield self.register(authenticate, 'authenticate')
        except Exception as e:
            print('Could not register custom WAMP-CRA authenticator: %s' % e)
        else:
            print('Custom WAMP-CRA authenticator registered.')
