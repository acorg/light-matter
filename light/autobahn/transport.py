"""
NOTE: The TransportFactory class is currently unused. To be removed?
"""

from autobahn.asyncio.wamp import ApplicationSessionFactory
from autobahn.wamp import types
from autobahn.wamp.serializer import MsgPackSerializer, JsonSerializer

from light.wamp import DEFAULT_REALM, DEFAULT_URL


class TransportFactory(object):
    """
    Provide a transport factory.

    @param sessionComponent: A C{light.autobahn.component} instance.
    @param wampURL: The C{str} URL that wamp websocket connections should be
        made to.
    """
    def __init__(self, sessionComponent, wampURL=DEFAULT_URL):
        self.sessionFactory = ApplicationSessionFactory(
            types.ComponentConfig(realm=DEFAULT_REALM))
        self.sessionFactory.session = sessionComponent
        self.wampURL = wampURL

    def factory(self, transportType='webSocket', debug=False):
        """
        Create a transport factory.

        @param transportType: A C{str} transport type, either 'websocket',
            'rawsocket-json', or 'rawsocket-msgpack'.
        @param debug: A C{bool}. If C{True}, debugging information will be
            printed.
        @raises ValueError: If an unknown transportType is given.
        """

        if transportType == 'websocket':
            # Create a WAMP-over-WebSocket transport client factory.
            from autobahn.asyncio.websocket import WampWebSocketClientFactory
            transportFactory = WampWebSocketClientFactory(
                self.sessionFactory, url=self.wampURL, debug_wamp=debug)
            transportFactory.setProtocolOptions(failByDrop=False)

        elif transportType in ['rawsocket-json', 'rawsocket-msgpack']:
            from autobahn.asyncio.rawsocket import WampRawSocketClientFactory
            # Create a WAMP-over-RawSocket transport client factory.
            serializer = (MsgPackSerializer()
                          if transportType == 'rawsocket-msgpack'
                          else JsonSerializer())

            transportFactory = WampRawSocketClientFactory(
                self.sessionFactory, serializer, debug=debug)

        else:
            raise ValueError('Unknown transportType %r.' % transportType)

        return transportFactory
