import sys

_None = object()

if sys.version_info < (3, 4):
    def maxWithDefault(a, default=_None):
        try:
            return max(a)
        except ValueError:
            if default is _None:
                raise
            else:
                return default

    def minWithDefault(a, default=_None):
        try:
            return min(a)
        except ValueError:
            if default is _None:
                raise
            else:
                return default
else:
    maxWithDefault = max
    minWithDefault = min
