import asyncio
import argparse

from autobahn.asyncio.wamp import ApplicationRunner

from light.autobahn.component import Component
from light.wamp import addArgsToParser


class CoordinatorComponent(Component):

    NAME = 'coordinator'

    @asyncio.coroutine
    def onJoin(self, details):
        self.sessions = set()

        # Subscribe to events for clients joining the router.
        def clientJoin(details):
            if details[b'authid'] == 'database':
                self.sessions.add(details[b'session'])
                print('Added database client, with session id',
                      details[b'session'])
            else:
                print('Non-database %r client joined.' % details[b'authid'])

        yield from self.subscribe(clientJoin, 'wamp.session.on_join')
        print("Subscribed to 'wamp.session.on_join'")

        # Subscribe to events for clients leaving the router.
        def clientLeave(session):
            print('wamp.session.on_leave event, session = %s' % session)
            self.sessions.discard(session)

        yield from self.subscribe(clientLeave, 'wamp.session.on_leave')
        print("Subscribed to 'wamp.session.on_leave'.")

        def addSubject(msg):
            # TODO: Pick a random database and give it the new subject.
            print('addSubject received msg %r' % msg)
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
        def shutdown():
            print('Shutdown received')
            if self.sessions:
                calls = [self.call('shutdown-%d' % session)
                         for session in self.sessions]
                result = yield from asyncio.gather(*calls)
            else:
                result = None
            print('shutdown result is', result)
            asyncio.get_event_loop().call_soon(lambda: self.leave('goodbye!'))
            self.sessions = set()
            # self.leave('goodbye!')
            return result

        yield from self.register(shutdown, 'shutdown')
        print("Registered 'shutdown' method.")

        @asyncio.coroutine
        def find(msg):
            print('Find received msg %r' % msg)
            calls = [self.call('find-%d' % session, msg)
                     for session in self.sessions]
            result = yield from asyncio.gather(*calls)
            return result

        yield from self.register(find, 'find')
        print("Registered 'find' method.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    addArgsToParser(parser)
    args = parser.parse_args()
    runner = ApplicationRunner(args.url, args.realm)
    runner.run(CoordinatorComponent)
