import logging
import asyncio
from autobahn.asyncio.wamp import ApplicationSession
from autobahn.wamp.auth import compute_wcs

from light.wamp import DEFAULT_REALM, DEFAULT_AUTH_METHOD


class Component(ApplicationSession):

    def onConnect(self):
        logging.info('%r client connected.', self.NAME)
        self.join(DEFAULT_REALM, [DEFAULT_AUTH_METHOD], self.NAME)

    def onChallenge(self, challenge):
        if challenge.method == DEFAULT_AUTH_METHOD:
            signature = compute_wcs(
                u'secret2'.encode('utf8'),
                challenge.extra['challenge'].encode('utf8'))
            return signature.decode('ascii')
        raise ValueError('Unknown authentication method %r' % challenge.method)

    def onLeave(self, details):
        super().onLeave(details)
        logging.info('%r client leaving.', self.NAME)

    def onDisconnect(self):
        logging.info('%r client disconnected.', self.NAME)
        asyncio.get_event_loop().stop()
