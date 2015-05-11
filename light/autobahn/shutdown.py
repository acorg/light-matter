import asyncio
import argparse

from autobahn.wamp.exception import ApplicationError
from autobahn.asyncio.wamp import ApplicationRunner

from light.autobahn.component import Component
from light.wamp import addArgsToParser


class ShutdownComponent(Component):

    NAME = 'shutdown'

    @asyncio.coroutine
    def onJoin(self, details):
        try:
            result = yield from self.call('shutdown')
        except ApplicationError as e:
            if e.error == 'wamp.error.no_such_procedure':
                print('Could not shutdown. It looks like no coordinator '
                      'worker is running on the router as the \'shutdown\' '
                      'procedure has not been registered.')
            else:
                print('Error', e.error)
        else:
            print('Shutdown result %r' % (result,))

        self.leave()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    addArgsToParser(parser)
    args = parser.parse_args()
    runner = ApplicationRunner(args.url, args.realm)
    runner.run(ShutdownComponent)
