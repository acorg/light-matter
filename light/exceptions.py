class BackendException(Exception):
    """
    Raised when there is a problem with a backend.
    """
    pass


class SubjectStoreException(Exception):
    """
    Raised when there is a problem with a SubjectStore.
    """
    pass


class WampDbOffline(Exception):
    """
    Raised when it appears a WAMP server database is offline.
    """
    pass
