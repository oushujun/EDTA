#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2010 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

from gt.core import *
import unittest
import tempfile
import os
import string


class AlphabetTest(unittest.TestCase):

    def setUp(self):
        self.dnaseqfile = tempfile.NamedTemporaryFile(mode="w", delete=False)
        self.dnaseqfile.write(">seq1\nagtccagctgtcagctagcgggcccgatgatatttt")
        self.aaseqfile = tempfile.NamedTemporaryFile(mode="w", delete=False)
        self.aaseqfile.write(">seq1\nMVHFTAEEKAAVTSLWSKMNVEEAGGEALG")
        self.dnaal1file = tempfile.NamedTemporaryFile(mode="w", suffix=".al1",
                                                      delete=False)
        self.dnaal1file.write("aA\ncC\ngG\ntTuU\nnsywrkvbdhmNSYWRKVBDHM\n")
        self.aaal1file = tempfile.NamedTemporaryFile(mode="w", suffix=".al1",
                                                     delete=False)
        self.aaal1file.write(
            "L\nV\nI\nF\nK\nR\nE\nD\nA\nG\nS\nT\nN\nQ\nY\nW\nP\nH\nM\nC\nXUBZO*- x\n")
        self.dnaseqfile.close()
        self.aaseqfile.close()
        self.dnaal1file.close()
        self.aaal1file.close()

    def tearDown(self):
        os.unlink(self.dnaseqfile.name)
        os.unlink(self.aaseqfile.name)
        os.unlink(self.dnaal1file.name)
        os.unlink(self.aaal1file.name)

    def dna_encodedecode(self, a_dna):
        # test DNA alphabet
        dnachars = {0: ('A', 'a'),   1: ('C', 'c'), 2: ('G', 'g'),
                    3: ('T', 't'),  254: ('N', 'n')}

        # encoding should yield the correct code independent of case
        for idx, c in dnachars.items():
            self.assertEqual(a_dna.encode(c[0]), idx)
            self.assertEqual(a_dna.decode(a_dna.encode(c[0])), c[1])
            self.assertEqual(a_dna.encode(c[1]), idx)
            self.assertEqual(a_dna.decode(a_dna.encode(c[1])), c[1])
        # invalid input should raise an exception
        self.assertRaises(GTError, a_dna.encode, 'Z')
        teststr = "agtccagctgtcagctagcgggcccgatgatatttt"
        # encode string 'by hand'
        encodedstr = [a_dna.encode(c) for c in teststr]
        # ensure that string is ok
        self.assertEqual(a_dna.decode_seq(encodedstr), teststr)

    def protein_encodedecode(self, a_protein):
        # test AA alphabet
        aas = "LVIFKREDAGSTNQYWPHMC"
        aachars = {}
        for idx, c in enumerate(aas):
            aachars[idx] = (c.upper(), c.lower())
        # print aachars
        wcs = ('X', 'x', 'U', 'u', 'B', 'b', 'Z', 'z', 'O', '*', '-')
        # encoding should yield the correct code independent of case
        for idx, c in aachars.items():
            self.assertEqual(a_protein.encode(c[0]), idx)
            self.assertEqual(a_protein.decode(a_protein.encode(c[0])), c[0])
        # invalid input should raise an exception
        self.assertRaises(GTError, a_protein.encode, '&')
        teststr = "MVHFTAEEKAAVTSLWSKMNVEEAGGEALG"
        # encode string 'by hand'
        encodedstr = [a_protein.encode(c) for c in teststr]
        # ensure that string is ok
        self.assertEqual(a_protein.decode_seq(encodedstr), teststr)

    def test_dna_newcreated(self):
        a_dna = Alphabet.create_dna()
        self.dna_encodedecode(a_dna)

    def test_dna_fromsequence(self):
        a_dna = Alphabet.create_from_sequence([self.dnaseqfile.name])
        self.dna_encodedecode(a_dna)

    def test_dna_fromsequence_fail(self):
        self.assertRaises(TypeError, Alphabet.create_from_sequence,
                          self.dnaseqfile.name)

    def test_dna_fromal1file(self):
        indexname = self.dnaal1file.name.rsplit('.', 1)[0]
        a_dna = Alphabet.create_from_file(indexname)
        self.dna_encodedecode(a_dna)

    def test_dna_fromal1file_fail(self):
        self.assertRaises(IOError, Alphabet.create_from_file, "nonexisting")

    def test_aa_newcreated(self):
        a_protein = Alphabet.create_protein()
        self.protein_encodedecode(a_protein)

    def test_aa_isdna(self):
        a_protein = Alphabet.create_protein()
        self.assertTrue(a_protein.is_protein())
        self.assertFalse(a_protein.is_dna())

    def test_isdna(self):
        a_protein = Alphabet.create_dna()
        self.assertFalse(a_protein.is_protein())
        self.assertTrue(a_protein.is_dna())

    def test_aa_fromsequence(self):
        a_protein = Alphabet.create_from_sequence([self.aaseqfile.name])
        self.protein_encodedecode(a_protein)

    def test_aa_fromsequence_fail(self):
        self.assertRaises(TypeError, Alphabet.create_from_sequence,
                          self.aaseqfile.name)

    def test_aa_fromsequence_fail_2(self):
        self.assertRaises(IOError, Alphabet.create_from_sequence,
                          ["nonexisting"])

    def test_aa_fromal1file(self):
        indexname = self.aaal1file.name.rsplit('.', 1)[0]
        a_protein = Alphabet.create_from_file(indexname)
        self.protein_encodedecode(a_protein)

    def test_aa_fromal1file_fail(self):
        self.assertRaises(IOError, Alphabet.create_from_file, "nonexisting")

if __name__ == "__main__":
    unittest.main()
