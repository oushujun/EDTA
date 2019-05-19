#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
import gt
import os

op = os.path
datadir = op.abspath(op.join(op.dirname(__file__), "..", "..",
                             "testdata"))

"""
###
1877523 gth gene    25221   26100   1   +   .   ID=gene12;Target=8690053 1 521 +
1877523 gth exon    25221   25310   1   +   .   Parent=gene12
1877523 gth five_prime_cis_splice_site  25311   25312   0.993   +   .   Parent=gene12
1877523 gth three_prime_cis_splice_site 25508   25509   0.439   +   .   Parent=gene12
1877523 gth exon    25510   25626   1   +   .   Parent=gene12
1877523 gth five_prime_cis_splice_site  25627   25628   0.998   +   .   Parent=gene12
1877523 gth three_prime_cis_splice_site 25711   25712   0.996   +   .   Parent=gene12
1877523 gth exon    25713   25841   1   +   .   Parent=gene12
1877523 gth five_prime_cis_splice_site  25842   25843   0.63    +   .   Parent=gene12
1877523 gth three_prime_cis_splice_site 25914   25915   0.926   +   .   Parent=gene12
1877523 gth exon    25916   26100   1   +   .   Parent=gene12
###
"""


class FeatureNodeIteratorTest(unittest.TestCase):

    def setUp(self):
        fi = gt.FeatureIndexMemory()
        self.gff_file = op.join(datadir, "U89959_sas.gff3")
        fi.add_gff3file(self.gff_file)
        self.fi = fi
        self.feature = fi.get_features_for_range(25000, 26000, '1877523')[0]

    def test_depth_first(self):
        dfi = gt.FeatureNodeIteratorDepthFirst(self.feature)

        found = dfi.next()
        self.assertEqual(found.get_attribute("ID"), "gene12")
        self.assertEqual(found.type, "gene")
        found = dfi.next()
        self.assertEqual(found.type, 'exon')
        found = dfi.next()
        self.assertEqual(found.type, 'five_prime_cis_splice_site')
        found = dfi.next()
        found = dfi.next()
        self.assertEqual(found.type, 'exon')

    def test_direct(self):
        di = gt.FeatureNodeIteratorDirect(self.feature)
        found = di.next()
        types = {}
        while found:
            types[found.type] = 1
            found = di.next()

        #self.assertTrue('gene' in types)

        self.assertTrue('exon' in types)
        self.assertTrue('five_prime_cis_splice_site' in types)
        self.assertTrue('three_prime_cis_splice_site' in types)

        self.assertEqual(found, None)


if __name__ == "__main__":
    unittest.main()
