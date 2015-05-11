import asyncio
import argparse

from autobahn.asyncio.wamp import ApplicationRunner

from light.autobahn.component import Component
from light.wamp import addArgsToParser


class FindComponent(Component):

    NAME = 'find'

    @asyncio.coroutine
    def onJoin(self, details):
        what = self.config.extra['what']
        try:
            result = yield from self.call('find', what)
        except Exception as e:
            print('Error:', e)
        else:
            print('Find %r = %r' % (what, result))

        self.leave()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--what', required=True,
                        help='The thing to find matches against.')
    addArgsToParser(parser)
    args = parser.parse_args()
    runner = ApplicationRunner(args.url, args.realm,
                               extra=dict(what=args.what))
    runner.run(FindComponent)
