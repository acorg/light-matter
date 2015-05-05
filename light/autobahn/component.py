import asyncio
from autobahn.asyncio import wamp
from autobahn.wamp.auth import compute_wcs


class Component(wamp.ApplicationSession):

    def onConnect(self):
        print('%r client connected to router' % self.NAME)
        self.join(u'realm1', [u'wampcra'], self.NAME)

    def onChallenge(self, challenge):
        if challenge.method == u'wampcra':
            signature = compute_wcs(u'secret2'.encode('utf8'),
                                    challenge.extra['challenge'])
            return signature.decode('ascii')
        else:
            raise Exception("don't know how to handle authmethod {}".format(
                challenge.method))

    def onLeave(self, details):
        print('Leaving.')
        self.disconnect()

    def onDisconnect(self):
        print('Disconnecting.')
        asyncio.get_event_loop().stop()
