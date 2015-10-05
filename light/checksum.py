from binascii import crc32


class Checksum(object):
    """
    Provide a simple CRC32 checksum with updating and equality checking.

    @param initialValue: The C{int} initial checksum value.
    """

    def __init__(self, initialValue=None):
        self._value = 0x0 if initialValue is None else initialValue

    def __eq__(self, other):
        return self._value == other._value

    def __repr__(self):
        return '<%s instance, value %d>' % (self.__class__.__name__,
                                            self._value)

    def update(self, strings):
        """
        Update the checksum value.

        @param strings: An iterable of strings to update the current checksum
            with.
        @return: self, to allow for chained calls.
        """
        text = b'\0'.join(map(lambda s: s.encode('UTF-8'), strings)) + b'\0'
        self._value = crc32(text, self._value) & 0xFFFFFFFF
        return self

    def _get(self):
        return self._value

    def _set(self, value):
        self._value = value

    value = property(_get, _set)
