DEFAULT_HOST = '127.0.0.1'
DEFAULT_PORT = 8080
# TODO: Use wss for secure connections.
DEFAULT_URL = 'ws://%s:%d/ws' % (DEFAULT_HOST, DEFAULT_PORT)
DEFAULT_REALM = 'light-matter'
DEFAULT_TRANSPORT_TYPE = 'websocket'
DEFAULT_AUTH_METHOD = 'wampcra'


def addArgsToParser(parser):
    """
    Add standard WAMP arguments to an argparse parser.

    @param parser: An C{argparse.ArgumentParser} instance.
    """
    parser.add_argument('--debugWamp', action='store_true', default=False,
                        help='Enable debug output.')

    parser.add_argument('--realm', default=DEFAULT_REALM,
                        help='The WAMP realm to join.')

    parser.add_argument('--host', default=DEFAULT_HOST,
                        help='IP or hostname to connect to.')

    parser.add_argument('--port', type=int, default=DEFAULT_PORT,
                        help='TCP port to connect to.')

    parser.add_argument('--wampUrl', default=DEFAULT_URL,
                        help='The WAMP router URL to connect to.')

    parser.add_argument('--transport', default=DEFAULT_TRANSPORT_TYPE,
                        choices=['websocket', 'rawsocket-json',
                                 'rawsocket-msgpack'],
                        help='WAMP transport type')
