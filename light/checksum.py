from binascii import crc32


class Checksum(object):
    """
    Provide a simple CRC32 checksum with updating and equality checking.

    @param initialValue: The C{int} initial checksum value.
    """

    def __init__(self, initialValue=None):
        self._checksum = 0x0 if initialValue is None else initialValue

    def __eq__(self, other):
        return self._checksum == other._checksum

    def update(self, strings):
        """
        Update the checksum.

        @param strings: An iterable of strings to update the current checksum
            with.
        @return: self, to allow for chained calls.
        """
        text = b'\0'.join(map(lambda s: s.encode('UTF-8'), strings)) + b'\0'
        self._checksum = crc32(text, self._checksum) & 0xFFFFFFFF
        return self

    def _get(self):
        return self._checksum

    def _set(self, value):
        self._checksum = value

    checksum = property(_get, _set)
