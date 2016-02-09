from six import StringIO

import asyncio

from autobahn.wamp.exception import ApplicationError

from light.autobahn.component import Component
from light.parameters import DatabaseParameters
from light.exceptions import WampDbOffline


class ClientComponent(Component):

    NAME = 'client'

    @asyncio.coroutine
    def onJoin(self, details):
        """
        Configure the WAMP client.

        @param details: A C{dict} of WAMP connection details.
        """
        from light.database_wamp import WampDatabaseClient
        future = self.config.extra['future']
        try:
            paramsStr = yield from self.call('parameters')
        except ApplicationError as error:
            if error.error == ApplicationError.NO_SUCH_PROCEDURE:
                # Looks like the WAMP database is not online.
                future.set_exception(WampDbOffline())
            else:
                raise
        else:
            params = DatabaseParameters.restore(StringIO(paramsStr))
            future.set_result(WampDatabaseClient(params, self))
