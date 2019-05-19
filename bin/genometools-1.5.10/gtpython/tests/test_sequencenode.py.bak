#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from gt import SequenceNode, Str


class SequenceNodeTestCase(unittest.TestCase):

    def setUp(self):
        self.fn = SequenceNode.create_new("testdesc", "AGATATAGA")

    def test_repr(self):
        self.assertEqual(str(self.fn),
                         'SequenceNode(start=0, end=0, seqid="testdesc")')

    def test_get_sequence(self):
        self.assertEqual(self.fn.get_sequence(), "AGATATAGA")

    def test_get_sequence_length(self):
        self.assertEqual(len(self.fn.get_sequence()),
                         self.fn.get_sequence_length())

    def test_get_description(self):
        self.assertEqual(self.fn.get_description(), 'testdesc')

if __name__ == "__main__":
    unittest.main()
