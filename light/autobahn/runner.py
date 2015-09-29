import txaio
txaio.use_asyncio()

# import signal

from autobahn.wamp.types import ComponentConfig
from autobahn.websocket.protocol import parseWsUrl
from autobahn.asyncio.websocket import WampWebSocketClientFactory

import asyncio


class ApplicationRunner:
    """
    This class is a convenience tool mainly for development and quick hosting
    of WAMP application components.

    It can host a WAMP application component in a WAMP-over-WebSocket client
    connecting to a WAMP router.
    """

    def __init__(self, url, realm, extra=None, serializers=None,
                 debug=False, debug_wamp=False, debug_app=False,
                 ssl=None):
        """
        :param url: The WebSocket URL of the WAMP router to connect to (e.g.
            `ws://somehost.com:8090/somepath`)
        :type url: unicode

        :param realm: The WAMP realm to join the application session to.
        :type realm: unicode

        :param extra: Optional extra configuration to forward to the
            application component.
        :type extra: dict

        :param serializers: A list of WAMP serializers to use (or None for
            default serializers). Serializers must implement
            :class:`autobahn.wamp.interfaces.ISerializer`.
        :type serializers: list

        :param debug: Turn on low-level debugging.
        :type debug: bool

        :param debug_wamp: Turn on WAMP-level debugging.
        :type debug_wamp: bool

        :param debug_app: Turn on app-level debugging.
        :type debug_app: bool

        :param ssl: An (optional) SSL context instance or a bool. See
           the documentation for the `loop.create_connection` asyncio
           method, to which this value is passed as the ``ssl=``
           kwarg.
        :type ssl: :class:`ssl.SSLContext` or bool
        """
        self.url = url
        self.realm = realm
        self.extra = extra or {}
        self.debug = debug
        self.debug_wamp = debug_wamp
        self.debug_app = debug_app
        self.make = None
        self.serializers = serializers
        self.ssl = ssl

    def run(self, make):
        """
        Run the application component.

        :param make: A factory that produces instances of
            :class:`autobahn.asyncio.wamp.ApplicationSession`
            when called with an instance of
            :class:`autobahn.wamp.types.ComponentConfig`.
        :type make: callable
        """
        # 1) factory for use ApplicationSession
        def create():
            cfg = ComponentConfig(self.realm, self.extra)
            try:
                session = make(cfg)
            except Exception as e:
                # the app component could not be created .. fatal
                print(e)
                asyncio.get_event_loop().stop()
            else:
                session.debug_app = self.debug_app
                return session

        isSecure, host, port, resource, path, params = parseWsUrl(self.url)

        if self.ssl is None:
            ssl = isSecure
        else:
            if self.ssl and not isSecure:
                raise RuntimeError(
                    'ssl argument value passed to %s conflicts with the "ws:" '
                    'prefix of the url argument. Did you mean to use "wss:"?' %
                    self.__class__.__name__)
            ssl = self.ssl

        # 2) create a WAMP-over-WebSocket transport client factory
        transport_factory = WampWebSocketClientFactory(
            create, url=self.url, serializers=self.serializers,
            debug=self.debug, debug_wamp=self.debug_wamp)

        # 3) start the client
        loop = asyncio.get_event_loop()
        # txaio.use_asyncio()
        # txaio.config.loop = loop
        coro = loop.create_connection(transport_factory, host, port, ssl=ssl)
        (transport, protocol) = loop.run_until_complete(coro)
        # loop.add_signal_handler(signal.SIGTERM, loop.stop)
