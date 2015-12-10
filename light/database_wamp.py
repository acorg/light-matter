try:
    import asyncio
except ImportError:
    # Python 2
    import trollius as asyncio

from twisted.internet.error import ConnectionRefusedError

from light.autobahn.runner import ApplicationRunner
from light.exceptions import WampDbOffline
from light.autobahn.client import ClientComponent


def getWampClientDatabase(args):
    """
    Get a client to talk to a remote WAMP database.

    @param args: Command line arguments as returned by the C{argparse}
        C{parse_args} method.
    @return: An instance of C{light.database.WampDatabaseClient}, or
        C{None} if no WAMP database can be found.
    """
    loop = asyncio.get_event_loop()
    future = asyncio.Future()
    runner = ApplicationRunner(args.wampUrl, args.realm,
                               extra=dict(future=future))
    try:
        # This will raise if we fail to connect to the WAMP router.
        runner.run(ClientComponent)
    except ConnectionRefusedError:
        pass
    else:
        try:
            # This will raise if there is no 'parameters' method registered
            # with the router. That indicates that no WAMP database is
            # connected to the router. (The 'parameters' method is the first
            # thing called by a light.autobahn.client.ClientComponent
            # instance.)
            database = loop.run_until_complete(future)
        except WampDbOffline:
            pass
        else:
            return database


class WampDatabaseClient:
    """
    Provide an interface to a remote WAMP database.

    The API we expose must be the same as that provided by a Database
    instance. I.e., anyone using this class should be able to treat it
    just like a regular Database.

    @param params: A C{Parameters} instance.
    @param component: A C{ligh.autobahn.ClientComponent} instance for
        communicating with a remote WAMP-based database.
    """
    def __init__(self, params, component):
        self._component = component
        self.params = params

    def close(self):
        self._component.leave()

    @asyncio.coroutine
    def checksum(self):
        """
        Get the current checksum.

        @return: An C{int} checksum.
        """
        # TODO: This should use the 'await' keyword from python >= 3.3.
        # return (await self._component.call('checksum'))
        return self._component.call('checksum')

    def print_(self, printHashes=False, margin=''):
        """
        Print information about the database.

        @param printHashes: If C{True}, print all hashes and associated
            subjects.
        @param margin: A C{str} that should be inserted at the start of each
            line of output.
        @return: A C{str} representation of the database.
        """
        coro = self._component.call('print_', printHashes=printHashes,
                                    margin=margin)
        loop = asyncio.get_event_loop()
        return loop.run_until_complete(coro)

    def subjectCount(self):
        """
        How many subjects are stored in this backend?

        @return: An C{int} number of subjects.
        """
        coro = self._component.call('subjectCount')
        loop = asyncio.get_event_loop()
        subjectCount = loop.run_until_complete(coro)
        return subjectCount

    # @asyncio.coroutine
    def hashCount(self):
        """
        How many hashes are in the backend database?

        @return: An C{int} number of hashes.
        """
        coro = self._component.call('hashCount')
        loop = asyncio.get_event_loop()
        hashCount = loop.run_until_complete(coro)
        return hashCount

    @asyncio.coroutine
    def totalResidues(self):
        """
        How many AA residues are there in the subjects stored in this backend?

        @return: An C{int} number of residues.
        """
        # TODO: This should use the 'await' keyword from python >= 3.3.
        # return (await self._component.call('totalResidues'))
        return self._component.call('totalResidues')

    @asyncio.coroutine
    def totalCoveredResidues(self):
        """
        How many AA residues are covered by landmarks and trig points in the
        subjects stored in this backend?

        @return: The C{int} number of residues that are covered.
        """
        # TODO: This should use the 'await' keyword from python >= 3.3.
        # return (await self._component.call('totalCoveredResidues'))
        return self._component.call('totalCoveredResidues')

    @asyncio.coroutine
    def addSubject(self, subject, subjectIndex=None):
        """
        Examine a sequence for features and add its (landmark, trig point)
        pairs to the search dictionary.

        @param subject: A C{dark.read.AAReadWithX} instance. The subject
            sequence is passed as a read instance even though in many cases it
            will not be an actual read from a sequencing run.
        @param subjectIndex: A C{str} representing the index of the subject as
            known by the database front end. If C{None} the SubjectStore we
            call addSubject on will assign an index.
        @return: A tuple of 1) a C{bool} to indicate whether the subject was
            already in the database, 2) the C{str} subject index, and 3) the
            hash count for the subject.
        """
        pass
