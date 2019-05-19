#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from gt import FeatureNode, CommentNode, SequenceNode, RegionNode, \
    CustomVisitor, MetaNode, EOFNode, FeatureNodeIteratorDepthFirst, \
    GTError


class TestVisitor(CustomVisitor):

    def __init__(self):
        CustomVisitor.__init__(self)
        self.sn = None
        self.rn = None
        self.cn = None
        self.en = None
        self.mn = None

    def visit_feature_node(self, fn):
        new_child = FeatureNode.create_new(
            fn.get_seqid(), "bar", 100, 1000, "+")
        fn.add_child(new_child)

    def visit_region_node(self, rn):
        self.rn = rn

    def visit_comment_node(self, cn):
        self.cn = cn

    def visit_sequence_node(self, sn):
        self.sn = sn

    def visit_meta_node(self, mn):
        self.mn = mn

    def visit_eof_node(self, en):
        self.en = en


class ErrorTestVisitor(CustomVisitor):

    def __init__(self):
        CustomVisitor.__init__(self)

    def visit_feature_node(self, fn):
        raise GTError


class CustomVisitorTestCase(unittest.TestCase):

    def setUp(self):
        self.fn = FeatureNode.create_new("foo", "gene", 100, 10000, "+")
        self.cn = CommentNode.create_new("comment")
        self.rn = RegionNode.create_new("foo", 100, 2000)
        self.sn = SequenceNode.create_new("foo", "AGATATAGA")
        self.en = EOFNode.create_new()
        self.mn = MetaNode.create_new("foo", "bar")
        self.tv = TestVisitor()
        self.etv = ErrorTestVisitor()

    def test_accept_feature_node(self):
        self.fn.accept(self.tv)
        dfi = FeatureNodeIteratorDepthFirst(self.fn)
        f = dfi.next()
        nodes = []
        while f:
            nodes.append(f)
            f = dfi.next()
        self.assertEqual(nodes[1].get_type(), "bar")

        self.assertRaises(GTError, self.fn.accept, self.etv)
        self.cn.accept(self.etv)
        self.sn.accept(self.etv)
        self.rn.accept(self.etv)

    def test_accept_sequence_node(self):
        self.assertNotEqual(self.tv.sn, self.sn)
        self.sn.accept(self.tv)
        self.assertEqual(self.tv.sn, self.sn)

    def test_accept_region_node(self):
        self.assertNotEqual(self.tv.rn, self.rn)
        self.rn.accept(self.tv)
        self.assertEqual(self.tv.rn, self.rn)

    def test_accept_comment_node(self):
        self.assertNotEqual(self.tv.cn, self.cn)
        self.cn.accept(self.tv)
        self.assertEqual(self.tv.cn, self.cn)

    def test_accept_eof_node(self):
        self.assertNotEqual(self.tv.en, self.en)
        self.en.accept(self.tv)
        self.assertEqual(self.tv.en, self.en)

    def test_accept_meta_node(self):
        self.assertNotEqual(self.tv.mn, self.mn)
        self.mn.accept(self.tv)
        self.assertEqual(self.tv.mn, self.mn)

if __name__ == "__main__":
    unittest.main()
