#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import gt
import os

op = os.path
datadir = op.abspath(op.join(op.dirname(__file__), "..", "..",
                             "testdata"))


class StreamTest(unittest.TestCase):

    def setUp(self):
        self.gff_file = op.join(datadir, "U89959_sas.gff3")
        self.ins = gt.GFF3InStream(self.gff_file)

    def test_pull(self):
        add_introns_stream = gt.AddIntronsStream(self.ins)
        fi = gt.FeatureIndexMemory()
        gt.FeatureStream(add_introns_stream, fi).pull()

        self.assertTrue('1877523' in fi.get_seqids())


class TestDuplicateStream(unittest.TestCase):

    def setUp(self):
        self.gff_file = op.join(datadir, "addintrons.gff3")
        self.ins = gt.GFF3InStream(self.gff_file)

    def test_dup(self):
        fi = gt.FeatureIndexMemory()
        gt.FeatureStream(gt.DuplicateFeatureStream(self.ins, "intron", "exon"),
                         fi).pull()

        f = fi.get_features_for_seqid('ctg123')
        dfi = gt.FeatureNodeIteratorDepthFirst(f[0])
        f = dfi.next()
        types = set([])
        while f:
            types.update([f.type])
            f = dfi.next()
        self.assertTrue('intron' in types, types)


class TestMergeStream(unittest.TestCase):

    def setUp(self):
        self.gff_file = op.join(datadir, "mergefeat.gff3")
        self.ins = gt.GFF3InStream(self.gff_file)

    def test_merge(self):
        fi = gt.FeatureIndexMemory()
        gt.FeatureStream(gt.MergeFeatureStream(self.ins), fi).pull()

        f = fi.get_features_for_seqid('seq1')
        self.assertEqual(len(f), 1)
        dfi = gt.FeatureNodeIteratorDirect(f[0])
        sub = dfi.next()
        self.assertEqual(None, dfi.next())

        self.assertEqual(sub.start, 1000)
        self.assertEqual(sub.end, 10000)


class TestInterFeat(unittest.TestCase):

    def setUp(self):
        self.gff_file = op.join(datadir, "addintrons.gff3")
        self.ins = gt.GFF3InStream(self.gff_file)

    def test_inter(self):
        fi = gt.FeatureIndexMemory()
        gt.FeatureStream(gt.InterFeatureStream(self.ins, "exon", "intron"),
                         fi).pull()

        f = fi.get_features_for_seqid('ctg123')
        dfi = gt.FeatureNodeIteratorDepthFirst(f[0])
        f = dfi.next()
        types = set([])
        while f:
            types.update([f.type])
            f = dfi.next()
        self.assertTrue('intron' in types, types)


class TestCustomExample(unittest.TestCase):

    def setUp(self):
        self.gff_file = op.join(datadir, "eden.gff3")
        self.ins = gt.GFF3InStream(self.gff_file)

    def test_inter(self):
        fi = gt.FeatureIndexMemory()
        gt.FeatureStream(gt.CustomStreamExample(self.ins),
                         fi).pull()

        f = fi.get_features_for_seqid('ctg123')
        dfi = gt.FeatureNodeIteratorDepthFirst(f[0])
        f = dfi.next()
        types = set([])
        while f:
            types.update([f.type])
            f = dfi.next()
        self.assertTrue('bar' in types, types)

if __name__ == "__main__":
    unittest.main()
