import logging

from asyncio import gather
from collections import defaultdict
from io import StringIO
from random import uniform, choice

try:
    from ujson import dump, loads
except ImportError:
    from json import dump, loads

from Bio.File import as_handle

from light.backend import Backend
from light.checksum import Checksum
from light.parameters import Parameters
from light.result import Result
from light.exceptions import BackendException
from light.subject import Subject, SubjectStore
from light.string import MultilineString

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                    level=logging.INFO)


class SimpleConnector:
    """
    Provide a simple in-memory connection between a Database and a single
    Backend.

    @param params: A C{Parameters} instance.
    @param backend: A C{Backend} instance, or C{None} if a backend should be
        created.
    @param filePrefix: Either a C{str} file name prefix to use as a default
        when saving or C{None} if no default save file is needed.
    """

    SAVE_SUFFIX = '.lmco'

    def __init__(self, params, backend=None, filePrefix=None):
        self.params = params
        if backend:
            self._backend = backend
        else:
            self._backend = Backend(filePrefix=filePrefix)
            self._backend.configure(params)
        self._filePrefix = filePrefix

        # Most of our implementation comes directly from our backend.
        for method in ('addSubject', 'getIndexBySubject', 'getSubjectByIndex',
                       'getSubjects', 'subjectCount', 'hashCount',
                       'totalResidues', 'totalCoveredResidues', 'checksum'):
            setattr(self, method, getattr(self._backend, method))

    def find(self, read, significanceMethod=None, scoreMethod=None,
             significanceFraction=None, storeFullAnalysis=False):
        """
        Check which database sequences a read matches.

        @param read: A C{dark.read.AARead} instance.
        @param significanceMethod: The name of the method used to calculate
            which histogram bins are considered significant.
        @param scoreMethod: The C{str} name of the method used to calculate the
            score of a bin which is considered significant.
        @param significanceFraction: The C{float} fraction of all (landmark,
            trig point) pairs for a scannedRead that need to fall into the
            same histogram bucket for that bucket to be considered a
            significant match with a database title.
        @param storeFullAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @return: A C{light.result.Result} instance.
        """
        if significanceMethod is None:
            significanceMethod = self.params.DEFAULT_SIGNIFICANCE_METHOD
        if scoreMethod is None:
            scoreMethod = self.params.DEFAULT_SCORE_METHOD
        if significanceFraction is None:
            significanceFraction = self.params.DEFAULT_SIGNIFICANCE_FRACTION

        matches, hashCount, nonMatchingHashes = self._backend.find(
            read, significanceMethod, scoreMethod, significanceFraction,
            storeFullAnalysis)

        return Result(read, self, matches, hashCount,
                      significanceMethod, scoreMethod, significanceFraction,
                      nonMatchingHashes, storeFullAnalysis=storeFullAnalysis)

    def shutdown(self, save, filePrefix):
        """
        Shut down the connector.

        @param save: If C{True}, save the connector state.
        @param filePrefix: When saving, use this C{str} as a file name prefix.
        """
        self._backend.shutdown(save, filePrefix)

        if save:
            self.save(filePrefix)

    def save(self, fpOrFilePrefix=None):
        """
        Save state to a file.

        @param fpOrFilePrefix: A file pointer, or the C{str} prefix of a file
            name, or C{None}. If a C{str}, self.SAVE_SUFFIX is appended to get
            the full file name. If C{None}, self._filePrefix will be used as a
            file prefix unless it is also C{None}.
        @raises ValueError: If C{fpOrFilePrefix} and C{self._filePrefix} are
            both C{None}
        """
        if isinstance(fpOrFilePrefix, str):
            saveFile = fpOrFilePrefix + self.SAVE_SUFFIX
        elif fpOrFilePrefix is None:
            if self._filePrefix is None:
                raise ValueError('save must be given an argument (or the '
                                 'database must have been restored from a '
                                 'file).')
            else:
                saveFile = self._filePrefix + self.SAVE_SUFFIX
        else:
            saveFile = fpOrFilePrefix

        with as_handle(saveFile, 'w') as fp:
            self.params.save(fp)

        self._backend.save(fpOrFilePrefix)

    @classmethod
    def restore(cls, fpOrFilePrefix):
        """
        Restore state from a file.

        @param fpOrFilePrefix: A file pointer or the C{str} prefix of a file
            name. If a C{str}, self.SAVE_SUFFIX is appended to get the full
            file name.
        @return: An instance of L{SimpleConnector}.
        @raises ValueError: If valid JSON cannot be loaded from C{fp}.
        """
        if isinstance(fpOrFilePrefix, str):
            saveFile = fpOrFilePrefix + cls.SAVE_SUFFIX
            filePrefix = fpOrFilePrefix
        else:
            saveFile = fpOrFilePrefix
            filePrefix = None

        with as_handle(saveFile) as fp:
            params = Parameters.restore(fp)

        return cls(params, backend=Backend.restore(fpOrFilePrefix),
                   filePrefix=filePrefix)

    def print_(self, printHashes=False, margin=''):
        """
        Print information about this connector.

        @param printHashes: If C{True}, print all hashes and associated
            subjects from the backend.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the parameters.
        """
        result = MultilineString(margin=margin)

        if printHashes:
            result.append('Backends:')
            result.append(self._backend.print_(margin=margin + '  '),
                          verbatim=True)

        return str(result)


class WampServerConnector:
    """
    Provide a connection between a persistent WAMP database component and a set
    of remote WAMP backend instances.

    @param params: A C{Parameters} instance.
    @param _id: A C{str} id, or C{None} to have an id assigned randomly.
    @param subjectStore: A C{SubjectStore} instance, or C{None} if a
        C{SubjectStore} should be created.
    @param checksum: An C{Checksum} instance, or C{None} if a new checksum
        should be created.
    @param disconnectedBackends: A C{dict} mapping backend names to checksum
        values. Or C{None} if no backends have previously been connected to
        this database server.
    @param filePrefix: Either a C{str} file name prefix to use as a default
        when saving or C{None} if no default save file is needed.
    """
    SAVE_SUFFIX = '.lmwsco'

    def __init__(self, params, _id=None, subjectStore=None, checksum=None,
                 disconnectedBackends=None, filePrefix=None):
        self.params = params
        # Assign ourselves a random id to make backend name collisions
        # highly unlikely. I.e., if an unknown backend connects, the name
        # it supplies should be very unlikely to be accepted.
        self._id = _id or '%020d' % int(uniform(1, 1e20))
        self._subjectStore = subjectStore or SubjectStore()
        for method in ('getIndexBySubject', 'getSubjectByIndex',
                       'getSubjects'):
            setattr(self, method, getattr(self._subjectStore, method))
        # Our checksum will be based on our parameters. The checksum should
        # not be based on our id as that's a random value that changes from
        # run to run. If we used it in the checksum then if we restarted a
        # connector and it began to talk to existing backends, it would
        # have a different checksum and errors would result.
        self._checksum = checksum or Checksum(params.checksum)
        # In self._disconnectedBackends and self._backends, keys are
        # backend names and values are Checksum instances.
        self._disconnectedBackends = disconnectedBackends or {}
        self._backends = {}
        self._filePrefix = filePrefix
        # _sessionIdToName maps ephemeral WAMP session ids to (connected)
        # backend names.
        self._sessionIdToName = {}
        # The WAMP server component. This will be set by the server
        # database component using our setComponent method.
        self._component = None

        # TODO: Remove.
        logging.info('Connector started. Disconnected backends: %r' %
                     self._disconnectedBackends)

    async def addBackend(self, sessionId):
        """
        Process a new backend connection. The connection may be from a new
        backend name, or a previously connected one (e.g., a backend that was
        connected before an earlier database shutdown, or one that was
        stopped and restarted).

        @param sessionId: An C{int} WAMP session id.
        @raise BackendException: If 1) a known backend connects but with a
            non-matching checksum, or 2) a known backend connects but it is
            marked as already connected.
        """
        logging.info('Adding backend with WAMP session id %s', sessionId)

        # Suggest a name and a checksum for the backend.
        backendCount = len(self._backends) + len(self._disconnectedBackends)
        suggestedName = 'backend-%s-%d' % (self._id, backendCount)

        # Sanity check: Make sure we've not used this name before.
        assert suggestedName not in self._backends
        assert suggestedName not in self._disconnectedBackends

        # Derive an initial checksum for the backend from the database
        # parameters and the backend name.
        suggestedChecksum = Backend.initialChecksum(self.params, suggestedName)

        # Get a string representation of the params to send over the network.
        paramsStr = self.params.save(StringIO()).getvalue()

        # Attempt to configure the backend. It may already be configured,
        # which is fine. We must do this before adding the backend name to
        # self._backends so there is no possibility we will try to use the
        # backend (in another request) before it has been configured.

        name, checksum, subjectCount = await self._component.call(
            'configure-%d' % sessionId, paramsStr, suggestedName,
            suggestedChecksum.value)

        if name == suggestedName:
            # The backend is new - it has accepted the name we suggested.
            # Store the checksum and associate the session id with the
            # backend.
            assert(subjectCount == 0)
            logging.info('Added new backend with name %s', name)
            self._checksum.update([name])
            self._backends[name] = {
                'subjectCount': 0,
                'checksum': suggestedChecksum,
            }
        else:
            # The backend is a just-restored instance that already had a
            # name, or is an instance that was already online (and
            # therefore already named) when the connector came online.
            # Check that we do not have its name in our (connected)
            # backends and that we do have it in our disconnected backends.
            if name in self._backends:
                # This backend is already connected! Or at least it was, or
                # some backend with the same name is or was.
                raise BackendException(
                    'Already connected backend %r has re-connected!' % name)

            if name in self._disconnectedBackends:
                logging.info('Backend %s (checksum %s) reconnected.', name,
                             checksum)
                # Make sure the checksum and subject count sent by the
                # backend match what we have on record.
                assert(self._disconnectedBackends[name]['checksum'].value ==
                       checksum)
                assert(self._disconnectedBackends[name]['subjectCount'] ==
                       subjectCount)

                self._backends[name] = self._disconnectedBackends[name]
                del self._disconnectedBackends[name]

                logging.info('Added re-connected backend with name %s', name)
            else:
                logging.warning('Apparently pre-existing backend %s (checksum '
                                '%d) has connected but is not present in our '
                                'known disconnected backends.', name, checksum)

        self._sessionIdToName[sessionId] = name

    def removeBackend(self, sessionId):
        """
        A backend has disconnected.

        @param sessionId: The C{int} session id of a WAMP client that
            disconneted. Note that this could be any kind of client, it may
            not be a backend.
        """
        try:
            name = self._sessionIdToName.pop(sessionId)
        except KeyError:
            # Not a session we know anything about. Ignore.
            logging.debug('Ignoring disconnect from unknown session %s.',
                          sessionId)
        else:
            try:
                self._disconnectedBackends[name] = self._backends.pop(name)
            except KeyError:
                logging.warning(
                    'Session %s named %r disconnected, but its name is not '
                    'present in our connected backends!' % (sessionId, name))
            else:
                logging.info('Session %s with name %r disconnected.',
                             sessionId, name)

    def setComponent(self, component):
        """
        Set the WAMP component to use for talking to the WAMP router.

        @param component: A C{light.autobahn.DatabaseComponent} instance.
        """
        self._component = component

    def subjectCount(self):
        """
        How many subjects are being held across all our backends?

        @return: An C{int} total number of subjects.
        """
        print('in server connector subjectCount')
        return len(self._subjectStore)

    async def hashCount(self):
        """
        How many hashes are stored in all our backends?

        TODO: This total count will be too large due to identical hashes
              being stored on separate backends. What to do? We could have
              the backends return their hashes and de-dup them here.

        @return: An C{int} total number of hashes.
        """
        print('IN connector hashcount.....................')
        calls = [self._component.call('hashCount-%d' % sessionId)
                 for sessionId in self._sessionIdToName]
        results = await gather(*calls)
        # coro = gather(*calls)
        # loop = asyncio.get_event_loop()
        # results = loop.run_until_complete(coro)
        return sum(results)

    def totalResidues(self):
        """
        How many AA residues are stored in all our backends?

        @return: An C{int} total number of AA residues.
        """
        return sum(len(s) for s in self._subjectStore.getSubjects())

    async def totalCoveredResidues(self):
        """
        How many AA residues are covered by landmarks across all our backends?

        @return: An C{int} total number of covered residues.
        """
        calls = [self._component.call('coveredResidues-%d' % sessionId)
                 for sessionId in self._sessionIdToName]
        results = await gather(*calls)
        return sum(results)

    def checksum(self):
        """
        Our checksum value.

        @return: The C{int} value of our checksum.
        """
        return self._checksum.value

    async def addSubject(self, subject):
        """
        Ask a backend to add a subject.

        @param subject: A C{dark.read.AARead} instance. The subject sequence
            is passed as a read instance even though in many cases it will not
            be an actual read from a sequencing run.
        @raises BackendException: if no backends are available.
        @return: A tuple of 1) a C{bool} to indicate whether the subject was
            already in the database, and 2) the C{str} subject index.
        """
        if not self._backends:
            raise BackendException('Database has no backends.')

        # We could check to see if self._disconnectedBackends is non-empty
        # and issue a warning if so. But there's no reason we can't add a
        # subject to one of the connected backends.

        subject = Subject(subject.id, subject.sequence, 0, subject.quality)
        preExisting, subjectIndex = self._subjectStore.add(subject)

        if preExisting:
            return True, subjectIndex

        # Pick a backend to add the subject to. Do this by choosing a
        # session id at random.
        sessionId = choice([self._sessionIdToName.keys()])
        preExisting, BackendSubjectIndex, hashCount = (
            await self._component.call(
                'addSubject-%s' % sessionId, subject.id, subject.sequence,
                subject.quality))

        # Sanity check that the backend had not seen the subject before and
        # that it used the subjectIndex we gave it.
        assert(preExisting is False)
        assert(subjectIndex == BackendSubjectIndex)

        # Store the hash count returned by the backend so that our subject
        # store is accurate.
        subject.hashCount = hashCount

        # Update our own checksum and that of the specific backend to which
        # the subject was added.
        name = self._sessionIdToName[sessionId]
        update = (subject.id, subject.sequence)
        self._checksum.update(update)
        self._backends[name]['checksum'].update(update)

        self._backends[name]['subjectCount'] += 1

        return False, subjectIndex

    async def find(self, read, significanceMethod=None, scoreMethod=None,
                   significanceFraction=None, storeFullAnalysis=False):
        """
        Check which database sequences a read matches.

        @param read: A C{dark.read.AARead} instance.
        @param significanceMethod: The name of the method used to calculate
            which histogram bins are considered significant.
        @param scoreMethod: The C{str} name of the method used to calculate the
            score of a bin which is considered significant.
        @param significanceFraction: The C{float} fraction of all (landmark,
            trig point) pairs for a scannedRead that need to fall into the
            same histogram bucket for that bucket to be considered a
            significant match with a database title.
        @param storeFullAnalysis: A C{bool}. If C{True} the intermediate
            significance analysis computed in the Result will be stored.
        @return: A C{light.result.Result} instance.
        """
        if significanceMethod is None:
            significanceMethod = self.params.DEFAULT_SIGNIFICANCE_METHOD
        if scoreMethod is None:
            scoreMethod = self.params.DEFAULT_SCORE_METHOD
        if significanceFraction is None:
            significanceFraction = self.params.DEFAULT_SIGNIFICANCE_FRACTION

        allMatches = defaultdict(list)
        allNonMatchingHashes = {}
        hashCount = 0

        calls = [
            self._component.call(
                'find-%d' % sessionId, read, significanceMethod, scoreMethod,
                significanceFraction, storeFullAnalysis)
            for sessionId in self._sessionIdToName]
        results = await gather(*calls)

        for matches, hashCount, nonMatchingHashes in results:
            # TODO: This is probably wrong...
            for nonMatchingHash in nonMatchingHashes:
                if nonMatchingHash not in allNonMatchingHashes:
                    allNonMatchingHashes[nonMatchingHash] = nonMatchingHashes[
                        nonMatchingHash]
            for subjectIndex in matches:
                # Make sure we have not have seen this subject before. If
                # we have, it would mean that two backends are reporting
                # results for the same subject. We may allow that later,
                # but for now it should be an error.
                intSubjectIndex = int(subjectIndex)
                assert(intSubjectIndex not in allMatches)
                allMatches[intSubjectIndex] = matches[subjectIndex]

        return Result(read, self, allMatches, hashCount,
                      significanceMethod, scoreMethod, significanceFraction,
                      allNonMatchingHashes,
                      storeFullAnalysis=storeFullAnalysis)

    async def shutdown(self, save, filePrefix):
        """
        Shut down the connector.

        @param save: If C{True}, save the connector state.
        @param filePrefix: When saving, use this C{str} as a file name prefix.
        """
        if self._backends:
            calls = [self.call('shutdown-%d' % sessionId, save, filePrefix)
                     for sessionId in self._sessionIdToName]
            await gather(*calls)

        if not save:
            self.save(filePrefix)

    def save(self, fpOrFilePrefix=None):
        """
        Save state to a file.

        @param fpOrFilePrefix: A file pointer, or the C{str} prefix of a file
            name, or C{None}. If a C{str}, self.SAVE_SUFFIX is appended to get
            the full file name. If C{None}, self._filePrefix will be used as a
            file prefix unless it is also C{None}.
        @raises ValueError: If C{fpOrFilePrefix} and C{self._filePrefix} are
            both C{None}
        @return: The result returned by the 'save' method of our backend.
        """
        if isinstance(fpOrFilePrefix, str):
            saveFile = fpOrFilePrefix + self.SAVE_SUFFIX
        elif fpOrFilePrefix is None:
            if self._filePrefix is None:
                raise ValueError('save must be given an argument (or the '
                                 'database must have been restored from a '
                                 'file).')
            else:
                saveFile = self._filePrefix + self.SAVE_SUFFIX
        else:
            saveFile = fpOrFilePrefix

        # Save the checksums and subject counts from all known backends
        # (whether or not they have connected to us yet. I.e., we might be
        # shutting down before some backends that were previously connected
        # to us have re-connected).
        disconnectedBackends = {}
        for backends in self._backends, self._disconnectedBackends:
            for name, backendInfo in backends.items():
                disconnectedBackends[name] = {
                    'checksum': backendInfo['checksum'].value,
                    'subjectCount': backendInfo['subjectCount'],
                }

        state = {
            'checksum': self.checksum(),
            'disconnectedBackends': disconnectedBackends,
            'id': self._id,
        }

        with as_handle(saveFile, 'w') as fp:
            self.params.save(fp)
            dump(state, fp)
            fp.write('\n')

    @classmethod
    def restore(cls, fpOrFilePrefix):
        """
        Restore state from a file.

        @param fpOrFilePrefix: A file pointer or the C{str} prefix of a file
            name. If a C{str}, self.SAVE_SUFFIX is appended to get the full
            file name.
        @return: An instance of L{WampServerConnector}.
        @raises ValueError: If valid JSON cannot be loaded from C{fp}.
        """
        if isinstance(fpOrFilePrefix, str):
            saveFile = fpOrFilePrefix + cls.SAVE_SUFFIX
            filePrefix = fpOrFilePrefix
        else:
            saveFile = fpOrFilePrefix
            filePrefix = None

        with as_handle(saveFile) as fp:
            params = Parameters.restore(fp)
            state = loads(fp.readline()[:-1])

        disconnectedBackends = {}
        for name, backendInfo in state['disconnectedBackends'].items():
            disconnectedBackends[name] = {
                'checksum': Checksum(backendInfo['checksum']),
                'subjectCount': backendInfo['subjectCount'],
            }

        return cls(params, _id=state['id'],
                   checksum=Checksum(state['checksum']),
                   disconnectedBackends=disconnectedBackends,
                   filePrefix=filePrefix)

    async def print_(self, printHashes=False, margin=''):
        """
        Print information about this connector and its backends.

        @param printHashes: If C{True}, print all hashes and associated
            subjects from the backends.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the parameters.
        """
        finalResult = MultilineString(margin=margin)

        if printHashes:
            extend = finalResult.extend
            append = finalResult.append
            finalResult.append('Backends:')

            sessionIds = [self._sessionIdToName.keys()]
            names = [self._sessionIdToName[sessionId]
                     for sessionId in sessionIds]
            calls = [self.call('print_-%d' % sessionId, margin=margin + '    ')
                     for sessionId in sessionIds]
            results = await gather(*calls)

            for name, result in zip(names, results):
                finalResult.indent()
                extend([
                    'Backend name: %s' % name,
                    'Backend details:',
                ])
                append(result, verbatim=True)
                finalResult.outdent()

        return str(finalResult)
