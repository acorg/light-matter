from io import StringIO
import logging
import asyncio

from light.autobahn.component import Component

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)


class DatabaseComponent(Component):

    NAME = 'database'

    @asyncio.coroutine
    def onJoin(self, details):
        print('joined.')
        self._database = self.config.extra['database']
        # We will talk to our database's connector directly to add/remove
        # backends. That way the database doesn't need to know anything
        # about WAMMP components connecting / disconnecting, etc.
        self._connector = self._database._connector
        self._connector.setComponent(self)

        @asyncio.coroutine
        def clientJoin(details):
            """
            A new client has connected to the WAMP router.

            @param details: The details of the client.
            """
            clientName = details['authid']
            sessionId = details['session']
            logging.info('WAMP %r client (with session id %s) joined.',
                         clientName, sessionId)
            if clientName == 'backend':
                yield from self._connector.addBackend(sessionId)

        # Subscribe to events for clients joining the router.
        self._joinSubscription = yield from self.subscribe(
            clientJoin, 'wamp.session.on_join')
        logging.info('Subscribed to wamp.session.on_join')

        def clientLeave(sessionId):
            """
            A client has disconnected to the WAMP router.

            @param session: The client session id.
            """
            self._connector.removeBackend(sessionId)

        # Subscribe to events for clients leaving the router.
        self._leaveSubscription = yield from self.subscribe(
            clientLeave, 'wamp.session.on_leave')
        logging.info('Subscribed to wamp.session.on_leave')

        # Look for existing backend sessions and tell the connector about
        # them.  Such sessions will exist if they were started before us or
        # if we exited and have been restarted.
        sessions = yield from self.call('wamp.session.list')
        for sessionId in sessions:
            session = yield from self.call('wamp.session.get', sessionId)
            if session['authid'] == 'backend':
                logging.info('Found existing backend session %s', sessionId)
                yield from self._connector.addBackend(sessionId)

        @asyncio.coroutine
        def shutdown(save, filePrefix):
            """
            Shut down the database, possibly saving state to a file.

            @param save: A C{bool}, if C{True} the database state will be
                saved to a file whose name starts with C{filePrefix}.
            @param filePrefix: A C{str} file name prefix to use if C{save} is
                C{True}.
            """
            logging.info('Shutdown received')
            # Unsubscribe from events for clients joining/leaving the
            # router. This shouldn't be needed, but due to a small bug in
            # the WAMP protocol implementation an exception is raised when
            # a Leave event (for ourselves) arrives after we call leave().
            # See https://github.com/tavendo/AutobahnPython/issues/459
            yield from self._leaveSubscription.unsubscribe()
            yield from self._joinSubscription.unsubscribe()

            self._database.shutdown(save, filePrefix)
            # Arrange to call self.leave in 0.1s time. Using self.leave
            # directly here or using call_soon causes the response not to
            # make it back to the WAMP client that called our shutdown
            # method. So we allow a little time for things to settle down.
            asyncio.get_event_loop().call_later(0.1, self.leave)

        yield from self.register(shutdown, 'shutdown')
        logging.info('Registered shutdown method.')

        def parameters():
            """
            Provide our database's parameters.

            @return: A C{str} representation of our parameters.
            """
            return self._database.params.save(StringIO()).getvalue()

        yield from self.register(parameters, 'parameters')
        logging.info('Registered parameters method.')

        # Most of our implementation comes directly from our database (which
        # has a WAMP connector that talks to WAMP backends).
        for method in ('addSubject', 'find', 'getIndexBySubject',
                       'getSubjectByIndex', 'getSubjects', 'subjectCount',
                       'hashCount', 'totalResidues', 'totalCoveredResidues',
                       'checksum', 'print_'):
            yield from self.register(getattr(self._database, method), method)
            logging.info('Registered %s method.', method)
