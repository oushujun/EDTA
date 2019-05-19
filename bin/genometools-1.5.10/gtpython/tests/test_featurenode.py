#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from gt import FeatureNode, FeatureNodeIteratorDepthFirst, GenomeNode


class FeatureNodeTestCase(unittest.TestCase):

    def setUp(self):
        self.fn = FeatureNode.create_new("test", "type", 100, 500, "+")

    def test_repr(self):
        self.assertEqual(str(self.fn),
                         'FeatureNode(start=100, end=500, seqid="test")')

    def test_score(self):
        fn = self.fn
        self.assertTrue(not fn.score_is_defined())

        fn.set_score(2)
        self.assertTrue(fn.score_is_defined())
        self.assertEqual(2, fn.get_score())

        fn.unset_score()
        self.assertTrue(not fn.score_is_defined())

    def test_type(self):
        fn = self.fn
        self.assertTrue(not fn.has_type("foo"))
        self.assertTrue(fn.has_type("type"))

    def test_strand(self):
        fn = self.fn
        self.assertEqual(fn.get_strand(), "+")

    def test_seqid(self):
        fn = self.fn
        self.assertEqual(fn.seqid, "test")

    def test_start_end(self):
        fn = self.fn
        self.assertEqual(fn.start, 100)
        self.assertEqual(fn.end, 500)

    def test_attributes(self):
        fn = self.fn
        fn.add_attribute("test", "testval")
        fn.add_attribute("test2", "testval2")

        self.assertTrue("test" in fn.attribs)
        self.assertTrue("test2" in fn.attribs)

        nattrs = 0
        for (tag, val) in fn.each_attribute():
            self.assertEqual(val, fn.get_attribute(tag))
            nattrs += 1

        self.assertEqual(nattrs, 2)


class TestFeatureNodeChildren(unittest.TestCase):

    def setUp(self):
        self.fn = FeatureNode.create_new("test", "type", 100, 500, "+")
        self.fn2 = FeatureNode.create_new("test", "type2", 200, 300, "+")
        self.fn.add_child(self.fn2)

    def test_phase(self):
        fn = self.fn
        self.assertEqual(fn.get_phase(), 3)

        fn.set_phase(0)
        self.assertEqual(fn.get_phase(), 0)

    def test_fni(self):
        fn = self.fn
        fni = FeatureNodeIteratorDepthFirst(fn)
        num_features = 0
        tfn = fni.next()
        while tfn:
            tfn = fni.next()
            num_features += 1
        self.assertEqual(num_features, 2)

        fn3 = FeatureNode.create_new("test", "type3", 250, 300, "+")
        fn.add_child(fn3)
        fni = FeatureNodeIteratorDepthFirst(fn)

        num_features = 0
        tfn = fni.next()
        while tfn:
            num_features += 1
            tfn = fni.next()
        self.assertEqual(num_features, 3)

    def test_iterator(self):
        fn = self.fn
        fn3 = FeatureNode.create_new("test", "type3", 250, 300, "+")
        fn4 = FeatureNode.create_new("test", "type4", 250, 300, "+")
        fn.add_child(fn3)
        fn.add_child(fn4)
        # try object as iterator
        types = []
        for i, f in enumerate(fn):
            types.append(f.type)
        self.assertEqual(types, ["type", "type2", "type3", "type4"], types)
        self.assertTrue(i == 3, i)
        # try iterator method
        types = []
        for i, f in enumerate(fn.traverse_dfs()):
            types.append(f.type)
        self.assertEqual(types, ["type", "type2", "type3", "type4"], types)
        self.assertTrue(i == 3, i)
        # try callable object as iterator
        types = []
        for i, f in enumerate(fn(method="depth_first")):
            types.append(f.type)
        self.assertEqual(types, ["type", "type2", "type3", "type4"], types)
        self.assertTrue(i == 3, i)
        # direct
        types = []
        for i, f in enumerate(fn(method="direct")):
            types.append(f.type)
        self.assertEqual(types, ["type2", "type3", "type4"], types)
        self.assertTrue(i == 2, i)
        types = []
        for i, f in enumerate(fn.traverse_direct()):
            types.append(f.type)
        self.assertEqual(types, ["type2", "type3", "type4"], types)
        self.assertTrue(i == 2, i)
        fn.depth_first = False
        types = []
        for i, f in enumerate(fn):
            types.append(f.type)
        self.assertEqual(types, ["type2", "type3", "type4"], types)
        self.assertTrue(i == 2, i)


class TestFeatureNodeProperties(unittest.TestCase):

    def setUp(self):
        self.fn = FeatureNode.create_new("test", "type", 100, 500, "+")

    def test_strand(self):
        fn = self.fn
        self.assertEqual("+", fn.strand)
        fn.strand = "-"
        self.assertEqual("-", fn.strand)

    def test_score(self):
        fn = self.fn
        self.assertTrue(not fn.score_is_defined())
        fn.score = 2
        self.assertTrue(fn.score_is_defined())
        self.assertEqual(2, fn.get_score())
        self.assertEqual(2, fn.score)

        fn.set_score(4)

        self.assertEqual(2, fn.score)

    def test_range(self):
        fn = self.fn
        self.assertEqual((100, 500), fn.range)

    def test_conversion(self):
        fn = self.fn
        g = GenomeNode.create_from_ptr(fn.gn, True)
        self.assertEqual((100, 500), g.range)

        f2 = FeatureNode.create_from_ptr(g.gn, True)
        self.assertEqual((100, 500), f2.range)


if __name__ == "__main__":
    unittest.main()
