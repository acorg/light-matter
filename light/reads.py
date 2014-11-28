class ScannedRead(object):
    """
    Hold information about a read that has been scanned for landmarks and trig
    points.

    @param read: A C{dark.reads.Read} instance.
    """

    def __init__(self, read):
        self.read = read
        self.landmarks = []
        self.trigPoints = []
