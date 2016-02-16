import os
import logging

import asyncio
from six import StringIO

from light.autobahn.component import Component
from light.backend import Backend
from light.parameters import DatabaseParameters

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)


class BackendComponent(Component):

    NAME = 'backend'

    @asyncio.coroutine
    def onJoin(self, details):
        self._sessionId = details.session
        logging.info('Joined with WAMP server, session id %s', self._sessionId)
        self._ApiConfigured = False
        args = self.config.extra['args']

        if args.filePrefix:
            saveFile = args.filePrefix + Backend.SAVE_SUFFIX
            if os.path.exists(saveFile):
                self._backend = Backend.restore(saveFile)
        else:
            # Note that the backend will be configured by the WAMP connector
            # via a call to our configure method (below).
            self._backend = Backend()

        @asyncio.coroutine
        def configure(paramsStr, suggestedName, suggestedChecksum):
            """
            Configure our backend and register its API methods.

            @param paramsStr: The C{str} 'save' of a C{DatabaseParameters}
                instance.
            @param suggestedName: The C{str} suggested name for this backend.
                If the backend has already been configured (from a file
                restore) with a different name, the suggested name is ignored.
            @param suggestedChecksum: The C{int} suggested checksum for the
                backend, or C{None} if there is no initial value.
            @return: A 2-tuple consisting of the C{str} name of the backend and
                its checksum. These will either be the suggested values or
                those that were already in use (if the backend was already
                configured).
            """
            fp = StringIO(paramsStr)
            dbParams = DatabaseParameters.restore(fp)
            name, checksum, subjectCount = self._backend.configure(
                dbParams, suggestedName, suggestedChecksum)

            if not self._ApiConfigured:
                self._ApiConfigured = True
                yield from self.registerAPIMethods()

            return name, checksum, subjectCount

        # Register our configure command.
        yield from self.register(configure, 'configure-%s' % self._sessionId)
        logging.info('Registered configure method.')

        # Register our shutdown command.
        def shutdown(save, filePrefix):
            """
            Shut down the backend.

            @param save: If C{True}, save the backend state.
            @param filePrefix: When saving, use this C{str} as a file name
                prefix.
            """
            logging.info('Shutdown called.')
            self._backend.shutdown(save, filePrefix)
            self.leave('goodbye!')

        yield from self.register(shutdown, 'shutdown-%s' % self._sessionId)
        logging.info('Registered shutdown method.')

    @asyncio.coroutine
    def registerAPIMethods(self):
        """
        Register our API methods with the router.
        """
        yield from self.register(self.find, 'find-%s' % self._sessionId)
        logging.info('Registered find method.')

        # Most of our implementation comes directly from our backend.
        for method in ('addSubject', 'getIndexBySubject', 'getSubjectByIndex',
                       'getSubjects', 'subjectCount', 'hashCount',
                       'totalResidues', 'totalCoveredResidues', 'checksum'):
            yield from self.register(getattr(self._backend, method),
                                     '%s-%s' % (method, self._sessionId))
            logging.info('Registered %s method.', method)

    def find(self, read, storeFullAnalysis=False, subjectIndices=None):
        """
        Check which database sequences a read matches.

        @param read: A C{dark.read.AARead} instance.
        @param storeFullAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @param subjectIndices: A C{set} of subject indices, or C{None}. If a
            set is passed, only subject indices in the set will be returned
            in the results. If C{None}, all matching subject indices are
            returned.
        @return: The result of calling 'find' on our backend.
        """
        return self._backend.find(read, storeFullAnalysis, subjectIndices)
