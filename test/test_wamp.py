from unittest import TestCase

from light.wamp import (
    DEFAULT_HOST, DEFAULT_PORT, DEFAULT_URL, DEFAULT_REALM,
    DEFAULT_TRANSPORT_TYPE, DEFAULT_AUTH_METHOD)


class TestWamp(TestCase):
    """
    Test the constant defaults in light.wamp
    """
    def testHost(self):
        """
        The host must be as expected.
        """
        self.assertEqual('127.0.0.1', DEFAULT_HOST)

    def testPort(self):
        """
        The port must be as expected.
        """
        self.assertEqual(8080, DEFAULT_PORT)

    def testURL(self):
        """
        The URL must be as expected.
        """
        self.assertEqual('ws://127.0.0.1:8080/ws', DEFAULT_URL)

    def testRealm(self):
        """
        The realm must be as expected.
        """
        self.assertEqual('light-matter', DEFAULT_REALM)

    def testTransportType(self):
        """
        The transport type must be as expected.
        """
        self.assertEqual('websocket', DEFAULT_TRANSPORT_TYPE)

    def testAuthMethod(self):
        """
        The authentication method must be as expected.
        """
        self.assertEqual('wampcra', DEFAULT_AUTH_METHOD)
