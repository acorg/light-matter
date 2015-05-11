import asyncio
from autobahn.asyncio import wamp
from autobahn.wamp.auth import compute_wcs

from light.wamp import DEFAULT_REALM, DEFAULT_AUTH_METHOD


class Component(wamp.ApplicationSession):

    def onConnect(self):
        print('%r client connected to router' % self.NAME)
        self.join(DEFAULT_REALM, [DEFAULT_AUTH_METHOD], self.NAME)

    def onChallenge(self, challenge):
        if challenge.method == DEFAULT_AUTH_METHOD:
            signature = compute_wcs(u'secret2'.encode('utf8'),
                                    challenge.extra['challenge'])
            return signature.decode('ascii')
        raise ValueError('Cannot handle authmethod %r' % challenge.method)

    def onLeave(self, details):
        print('Leaving.')
        self.disconnect()

    def onDisconnect(self):
        print('Disconnecting.')
        asyncio.get_event_loop().stop()
