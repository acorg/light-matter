import asyncio
import logging

from autobahn.wamp.exception import ApplicationError

from light.autobahn.component import Component


class ShutdownComponent(Component):
    """
    A WAMP component for shutting down a Database.
    """

    NAME = 'shutdown'

    @asyncio.coroutine
    def onJoin(self, details):
        """
        Shut down the WAMP database server.

        @param details: A C{dict} of WAMP connection details.
        """
        try:
            logging.info('Calling shutdown')
            yield from self.call('shutdown',
                                 self.config.extra['noSave'],
                                 self.config.extra['filePrefix'])
            logging.info('Called shutdown')
        except ApplicationError as e:
            if e.error == 'wamp.error.no_such_procedure':
                logging.warning(
                    'Could not shutdown. It looks like no Database '
                    'component is running on the WAMP router as the '
                    '"shutdown" procedure has not been registered.')
            else:
                logging.warning('Error trying to shutdown:', e.error)
        else:
            logging.info('Shutdown successful.')

        self.leave()
