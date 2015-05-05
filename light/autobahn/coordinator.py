import asyncio
import argparse

from light.autobahn.component import Component


class CoordinatorComponent(Component):

    NAME = 'coordinator'

    @asyncio.coroutine
    def onJoin(self, details):
        self.sessions = set()

        print('coordinator joined. details are {}.'.format(details))

        # Subscribe to events for clients joining the router.
        def clientJoin(details):
            # print('client joined', details)
            if details[b'authid'] == 'database':
                self.sessions.add(details[b'session'])
                print('Added database client, with session id',
                      details[b'session'])
            else:
                print('Client type %s joined.' % details[b'authid'])

        yield from self.subscribe(clientJoin, 'wamp.session.on_join')
        print("Subscribed to 'wamp.session.on_join'")

        # Subscribe to events for clients leaving the router.
        def clientLeave(session):
            print('wamp.session.on_leave event, session = {}'.format(session))
            self.sessions.discard(session)

        yield from self.subscribe(clientLeave, 'wamp.session.on_leave')
        print("Subscribed to 'wamp.session.on_leave'.")

        def addSubject(msg):
            # TODO: Pick a random database and give it the new subject.
            print('addSubject received msg {}'.format(msg))
            return True

        yield from self.register(addSubject, 'addSubject')
        print("Registered 'addSubject' method.")

        # Look for existing database sessions and add them to our set of
        # sessions. Such sessions will exist if they were started before us
        # or in case we have been restarted.
        sessions = yield from self.call('wamp.session.list')
        for sessionId in sessions:
            session = yield from self.call('wamp.session.get', sessionId)
            # from pprint import pprint
            # print('GOT SESSION')
            # pprint(session)
            if session[b'authid'] == 'database':
                self.sessions.add(sessionId)

        @asyncio.coroutine
        def quit():
            print('Quit received')
            result = yield from asyncio.gather(
                self.call('quit-%d' % session) for session in self.sessions)
            self.sessions = set()
            asyncio.get_event_loop().call_soon(lambda: self.leave('goodbye!'))
            return result

        yield from self.register(quit, 'quit')
        print("Registered 'quit' method.")

        @asyncio.coroutine
        def find(msg):
            print('Find received msg {}'.format(msg))
            calls = [self.call('find-%d' % session, msg)
                     for session in self.sessions]
            result = yield from asyncio.gather(*calls)
            return result

        yield from self.register(find, 'find')
        print("Registered 'find' method.")


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-d', '--debug', action='store_true',
                        help='Enable debug output.')

    parser.add_argument(
        '-r', '--realm', default='realm1',
        help='The WAMP realm to start the component in (if any).')

    parser.add_argument(
        '--host', default='127.0.0.1',
        help='IP or hostname to connect to.')

    parser.add_argument(
        '--port', type=int, default=8080,
        help='TCP port to connect to.')

    parser.add_argument(
        '--transport',
        choices=['websocket', 'rawsocket-json', 'rawsocket-msgpack'],
        default='websocket', help='WAMP transport type')

    parser.add_argument(
        '--url', default='ws://127.0.0.1:8080/ws',
        help='The WebSocket URL to connect to, e.g. ws://127.0.0.1:8080/ws.')

    args = parser.parse_args()

    # create a WAMP application session factory
    ##
    from autobahn.asyncio.wamp import ApplicationSessionFactory
    from autobahn.wamp import types
    session_factory = ApplicationSessionFactory(types.ComponentConfig(realm=args.realm))

    session_factory.session = CoordinatorComponent

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
