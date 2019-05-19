#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from gt import CommentNode, Str


class CommentNodeTestCase(unittest.TestCase):

    def setUp(self):
        self.fn = CommentNode.create_new("testcomment")
        self.fn2 = CommentNode.create_new(333)

    def test_repr(self):
        self.assertEqual(str(self.fn),
                         'CommentNode(start=0, end=0, seqid="None")')

    def test_get_comment(self):
        self.assertEqual(self.fn.get_comment(), 'testcomment')
        self.assertEqual(self.fn2.get_comment(), '333')

if __name__ == "__main__":
    unittest.main()
