import asyncio
import argparse

from autobahn.asyncio.wamp import ApplicationRunner

from light.autobahn.component import Component
from light.wamp import addArgsToParser
# from light.database import Backend


class BackendComponent(Component):

    NAME = 'backend'

    @asyncio.coroutine
    def onJoin(self, details):
        self.sessionId = details.session

        print('database client joined', details)

        # Register our find command.
        def find(msg):
            print('find called with %r' % msg)
            return {}

        reg = yield from self.register(find, 'find-%s' % self.sessionId)
        print("Registered 'addSubject':", reg)

        # Register our shutdown command.
        def shutdown():
            print('shutdown called')
            self.leave('goodbye!')
            print('leave returned')

        reg = yield from self.register(shutdown,
                                       'shutdown-%s' % self.sessionId)
        print("Registered 'shutdown':", reg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    addArgsToParser(parser)
    args = parser.parse_args()
    runner = ApplicationRunner(args.url, args.realm)
    runner.run(BackendComponent)
