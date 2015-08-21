import logging

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)

import asyncio

from light.autobahn.component import Component


class DatabaseComponent(Component):

    NAME = 'database'

    @asyncio.coroutine
    def onJoin(self, details):
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
        yield from self.subscribe(clientJoin, 'wamp.session.on_join')
        logging.info('Subscribed to wamp.session.on_join')

        def clientLeave(sessionId):
            """
            A client has disconnected to the WAMP router.

            @param session: The client session id.
            """
            self._connector.removeBackend(sessionId)

        # Subscribe to events for clients leaving the router.
        yield from self.subscribe(clientLeave, 'wamp.session.on_leave')
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

        # @asyncio.coroutine
        def shutdown(noSave, filePrefix):
            logging.info('Shutdown received')
            self._database.shutdown(noSave, filePrefix)
            # self.leave()
            asyncio.get_event_loop().call_soon(lambda: self.leave('goodbye!'))

        yield from self.register(shutdown, 'shutdown')
        logging.info('Registered shutdown method.')

        # Most of our implementation comes directly from our database (which
        # has a WAMP connector that talks to WAMP backends).
        for method in ('addSubject', 'find', 'getIndexBySubject',
                       'getSubjectByIndex', 'getSubjects', 'subjectCount',
                       'hashCount', 'totalResidues', 'totalCoveredResidues',
                       'checksum'):
            yield from self.register(getattr(self._database, method), method)
            logging.info('Registered %s method.', method)
