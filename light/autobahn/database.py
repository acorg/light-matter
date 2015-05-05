import asyncio

from light.autobahn.component import Component


class DatabaseComponent(Component):

    NAME = 'database'

    @asyncio.coroutine
    def onJoin(self, details):
        self.sessionId = details.session

        print('database client joined', details)

        # Register our find command.
        def find(msg):
            print('find called with {}'.format(msg))
            return {}

        reg = yield from self.register(find, 'find-{}'.format(self.sessionId))
        print("Registered 'addSubject': {}".format(reg))

        # Register our quit command.
        def quit():
            print('quit called')
            asyncio.get_event_loop().call_soon(lambda: self.leave('goodbye!'))
            print('leave returned')

        reg = yield from self.register(quit, 'quit-{}'.format(self.sessionId))
        print("Registered 'quit': {}".format(reg))


if __name__ == '__main__':

    import argparse

    try:
        import asyncio
    except ImportError:
        # Trollius >= 0.3 was renamed
        import trollius as asyncio

    # parse command line arguments
    ##
    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--debug', action='store_true',
                        help='Enable debug output.')

    parser.add_argument('-c', '--component',
                        help="Start WAMP client with this application component, e.g. 'timeservice.TimeServiceFrontend'")

    parser.add_argument('-r', '--realm', default='realm1',
                        help='The WAMP realm to start the component in (if any).')

    parser.add_argument('--host', default='127.0.0.1',
                        help='IP or hostname to connect to.')

    parser.add_argument('--port', type=int, default=8080,
                        help='TCP port to connect to.')

    parser.add_argument('--transport', choices=['websocket', 'rawsocket-json', 'rawsocket-msgpack'], default='websocket',
                        help='WAMP transport type')

    parser.add_argument('--url', default='ws://127.0.0.1:8080/ws',
                        help='The WebSocket URL to connect to, e.g. ws://127.0.0.1:8080/ws.')

    args = parser.parse_args()

    # create a WAMP application session factory
    ##
    from autobahn.asyncio.wamp import ApplicationSessionFactory
    from autobahn.wamp import types
    session_factory = ApplicationSessionFactory(types.ComponentConfig(realm=args.realm))

    session_factory.session = DatabaseComponent

    if args.transport == 'websocket':

        # create a WAMP-over-WebSocket transport client factory
        ##
        from autobahn.asyncio.websocket import WampWebSocketClientFactory
        transport_factory = WampWebSocketClientFactory(session_factory, url=args.url, debug_wamp=args.debug)
        transport_factory.setProtocolOptions(failByDrop=False)

    elif args.transport in ['rawsocket-json', 'rawsocket-msgpack']:

        # create a WAMP-over-RawSocket transport client factory
        ##
        if args.transport == 'rawsocket-msgpack':
            from autobahn.wamp.serializer import MsgPackSerializer
            serializer = MsgPackSerializer()
        elif args.transport == 'rawsocket-json':
            from autobahn.wamp.serializer import JsonSerializer
            serializer = JsonSerializer()
        else:
            raise Exception('should not arrive here')

        from autobahn.asyncio.rawsocket import WampRawSocketClientFactory
        transport_factory = WampRawSocketClientFactory(session_factory, serializer, debug=args.debug)

    else:
        raise Exception('logic error')

    # start the client
    loop = asyncio.get_event_loop()
    coro = loop.create_connection(transport_factory, args.host, args.port)
    loop.run_until_complete(coro)

    try:
        loop.run_forever()
    except KeyboardInterrupt:
        pass
