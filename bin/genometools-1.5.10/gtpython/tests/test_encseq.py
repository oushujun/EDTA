#!/usr/bin/python
# -*- coding: utf-8 -*-

from gt.core.encseq import *
from gt.core.error import GTError
from gt.core import readmode
import unittest
import tempfile
import sys
import os
import string


class EncodedsequenceTest(unittest.TestCase):

    def setUp(self):
        self.dnaseqfile = tempfile.NamedTemporaryFile(mode="w", delete=False)
        self.dnaseqfile.close()
        self.dnafname = self.dnaseqfile.name
        self.dnafile = open(self.dnafname, "w")
        self.dseq1 = "agtccagctgtcagctagcgggcccgatgatatttt"
        self.dseq2 = "gtgctgtac"
        self.dseq3 = "gtacagcac"
        self.dseq4 = "aaaatatcatcgggcccgctagctgacagctggact"
        self.dnafile.write(">seq1\n" + self.dseq1 + "\n")
        self.dnafile.write(">seq2\n" + self.dseq2 + "\n")

        self.aaseqfile = tempfile.NamedTemporaryFile(mode="w", delete=False)
        self.aaseqfile.close()
        self.aafname = self.aaseqfile.name
        self.aafile = open(self.aafname, "w")
        self.aaseq1 = "MVHFTAEEKAAVTSLWSKMNVEEAGGEALG"
        self.aaseq2 = "KMNAVE"
        self.aafile.write(">seq1\n" + self.aaseq1 + "\n")
        self.aafile.write(">seq2\n" + self.aaseq2 + "\n")
        self.dnafile.close()
        self.aafile.close()
        self.idxsuffixes = ['esq', 'des', 'ssp', 'sds']

    def tearDown(self):
        pass  # os.unlink(self.dnafname)
        # os.unlink(self.aafname)

    def create_es(self, indexname):
        ee = EncseqEncoder()
        return ee.encode([self.dnafname], indexname)

    def create_es_protein(self, indexname):
        ee = EncseqEncoder()
        return ee.encode([self.aafname], indexname)

    def create_mem(self):
        a = Alphabet.create_dna()
        eb = EncseqBuilder(a)
        eb.enable_description_support()
        eb.enable_multiseq_support()
        eb.add_string(self.dseq1, 'seq1')
        eb.add_string(self.dseq2, 'seq2')
        return eb.build()

    def create_mem_protein(self):
        a = Alphabet.create_protein()
        eb = EncseqBuilder(a)
        eb.enable_description_support()
        eb.enable_multiseq_support()
        eb.add_string(self.aaseq1, 'seq1')
        eb.add_string(self.aaseq2, 'seq2')
        return eb.build()

    def delete_idx(self, indexname):
        for suf in self.idxsuffixes:
            if os.path.isfile(indexname + "." + suf):
                os.unlink(indexname + "." + suf)

    def test_create_new(self):
        val = self.create_es("foo")
        for suf in self.idxsuffixes:
            self.assertTrue(os.path.isfile("foo." + suf))
        self.delete_idx("foo")

    def test_create_mapped(self):
        self.create_es("foo_mapped")
        el = EncseqLoader()
        es = el.load("foo_mapped")
        self.assertNotEqual(es, None)
        self.delete_idx("foo_mapped")

    def test_map_fail(self):
        el = EncseqLoader()
        self.assertRaises(IOError, el.load, "foo_fail")

    def test_dna(self):
        self.create_es("foo")
        el = EncseqLoader()
        es = el.load("foo")
        self.run_test_descriptions(es)
        self.run_test_get_encoded_char(
            es, self.dseq1, self.dseq2, self.dseq3, self.dseq4)
        self.run_test_num_seqs(es)
        self.run_test_num_files(es)
        self.run_test_seq_length(es)
        self.run_test_seq_startpos(es)
        self.run_test_seq_substr_encoded(
            es, self.dseq1, self.dseq2, self.dseq3, self.dseq4)
        self.run_test_seq_substr_plain(
            es, self.dseq1, self.dseq2, self.dseq3, self.dseq4)
        self.run_test_seq_substr_sequential(
            es, self.dseq1, self.dseq2, self.dseq3, self.dseq4)
        self.run_test_total_length(es)
        self.delete_idx("foo")
        es = self.create_mem()
        self.run_test_descriptions(es)
        self.run_test_get_encoded_char(
            es, self.dseq1, self.dseq2, self.dseq3, self.dseq4)
        self.run_test_num_seqs(es)
        self.run_test_num_files_mem(es)
        self.run_test_seq_length(es)
        self.run_test_seq_substr_encoded(
            es, self.dseq1, self.dseq2, self.dseq3, self.dseq4)
        self.run_test_seq_substr_plain(
            es, self.dseq1, self.dseq2, self.dseq3, self.dseq4)
        self.run_test_seq_substr_sequential(
            es, self.dseq1, self.dseq2, self.dseq3, self.dseq4)
        self.run_test_total_length(es)

    def test_protein(self):
        self.create_es_protein("foo")
        el = EncseqLoader()
        es = el.load("foo")
        self.run_test_descriptions(es)
        self.run_test_get_encoded_char(
            es, self.aaseq1, self.aaseq2, None, None)
        self.run_test_num_seqs(es)
        self.run_test_num_files(es)
        self.run_test_seq_startpos_protein(es)
        self.run_test_seq_length_protein(es)
        self.run_test_file_length_protein(es)
        self.run_test_seq_substr_encoded(
            es, self.aaseq1, self.aaseq2, None, None)
        self.run_test_seq_substr_plain(
            es, self.aaseq1, self.aaseq2, None, None)
        self.run_test_seq_substr_sequential(
            es, self.aaseq1, self.aaseq2, None, None)
        self.run_test_total_length_protein(es)
        self.delete_idx("foo")
        es = self.create_mem_protein()
        self.run_test_descriptions(es)
        self.run_test_get_encoded_char(
            es, self.aaseq1, self.aaseq2, None, None)
        self.run_test_num_seqs(es)
        self.run_test_num_files_mem(es)
        self.run_test_seq_length_protein(es)
        self.run_test_seq_substr_encoded(
            es, self.aaseq1, self.aaseq2, None, None)
        self.run_test_seq_substr_plain(
            es, self.aaseq1, self.aaseq2, None, None)
        self.run_test_seq_substr_sequential(
            es, self.aaseq1, self.aaseq2, None, None)
        self.run_test_total_length_protein(es)

    def run_test_num_seqs(self, es):
        self.assertEqual(es.num_of_sequences(), 2)
        if es.alphabet().is_dna():
            es.mirror()
            self.assertEqual(es.num_of_sequences(), 4)
            es.unmirror()

    def run_test_num_files(self, es):
        self.assertEqual(es.num_of_files(), 1)
        if es.alphabet().is_dna():
            es.mirror()
            self.assertEqual(es.num_of_files(), 1)
            es.unmirror()

    def run_test_num_files_mem(self, es):
        self.assertEqual(es.num_of_files(), 1)
        if es.alphabet().is_dna():
            es.mirror()
            self.assertEqual(es.num_of_files(), 1)
            es.unmirror()

    def run_test_descriptions(self, es):
        self.assertRaises(GTError, es.description, 2)
        self.assertEqual(es.description(0), "seq1")
        self.assertEqual(es.description(1), "seq2")
        if es.alphabet().is_dna():
            es.mirror()
            self.assertRaises(GTError, es.description, 5)
            self.assertEqual(es.description(2), "seq2")
            self.assertEqual(es.description(3), "seq1")
            es.unmirror()

    def run_test_total_length(self, es):
        self.assertEqual(es.total_length(), 46)
        if es.alphabet().is_dna():
            es.mirror()
            self.assertEqual(es.total_length(), 93)
            es.unmirror()

    def run_test_total_length_protein(self, es):
        self.assertEqual(es.total_length(), 37)

    def run_test_get_encoded_char(self, es, seq1, seq2, seq3, seq4):
        a = es.alphabet()
        for i, c in enumerate(seq1):
            encchar = es.get_encoded_char(i, readmode.FORWARD)
            self.assertEqual(a.decode(encchar), c)
        for i, c in enumerate(seq2[::-1]):
            encchar = es.get_encoded_char(i, readmode.REVERSE)
            self.assertEqual(a.decode(encchar), c)
        if es.alphabet().is_dna():
            es.mirror()
            for i, c in enumerate(seq3):
                encchar = es.get_encoded_char(47 + i, readmode.FORWARD)
                self.assertEqual(a.decode(encchar), c)
            for i, c in enumerate(seq4[::-1]):
                encchar = es.get_encoded_char(i, readmode.REVERSE)
                self.assertEqual(a.decode(encchar), c)
            es.unmirror()

    def run_test_seq_startpos(self, es):
        self.assertEqual(es.seqstartpos(0), 0)
        self.assertEqual(es.seqstartpos(1), 37)

    def run_test_seq_startpos_protein(self, es):
        self.assertEqual(es.seqstartpos(0), 0)
        self.assertEqual(es.seqstartpos(1), 31)

    def run_test_seq_length(self, es):
        self.assertEqual(es.seqlength(0), 36)
        self.assertEqual(es.seqlength(1), 9)

    def run_test_file_length(self, es):
        self.assertEqual(es.effective_filelength(0), 46)

    def run_test_seq_length_protein(self, es):
        self.assertEqual(es.seqlength(0), 30)
        self.assertEqual(es.seqlength(1), 6)

    def run_test_file_length_protein(self, es):
        self.assertEqual(es.effective_filelength(0), 37)

    def run_test_seq_substr_encoded(self, es, seq1, seq2, seq3, seq4):
        start = 3
        end = 13
        res = es.extract_encoded(start, end)
        a = es.alphabet()
        for i in range(start, end):
            self.assertEqual(a.decode(res[i - start]), seq1[i])
        start = 0
        end = 5
        ssp = es.seqstartpos(1)
        res = es.extract_encoded(ssp + start, ssp + end)
        for i in range(start, end):
            self.assertEqual(a.decode(res[i - start]), seq2[i])
        if es.alphabet().is_dna():
            es.mirror()
            start = 3
            end = 8
            ssp = es.seqstartpos(2)
            res = es.extract_encoded(ssp + start, ssp + end)
            for i in range(start, end):
                self.assertEqual(a.decode(res[i - start]), seq3[i])
            start = 0
            end = 5
            ssp = es.seqstartpos(3)
            res = es.extract_encoded(ssp + start, ssp + end)
            for i in range(start, end):
                self.assertEqual(a.decode(res[i - start]), seq4[i])
            es.unmirror()

    def run_test_seq_substr_plain(self, es, seq1, seq2, seq3, seq4):
        start = 3
        end = 13
        self.assertEqual(es.extract_decoded(start, end), seq1[start:end + 1])
        start = 0
        end = 5
        ssp = es.seqstartpos(1)
        self.assertEqual(es.extract_decoded(ssp + start, ssp + end),
                         seq2[start:end + 1])
        if es.alphabet().is_dna():
            es.mirror()
            start = 3
            end = 8
            ssp = es.seqstartpos(2)
            res = es.extract_encoded(ssp + start, ssp + end)
            for i in range(start, end):
                self.assertEqual(es.extract_decoded(ssp + start, ssp + end),
                                 seq3[start:end + 1])
            start = 0
            end = 5
            ssp = es.seqstartpos(3)
            res = es.extract_encoded(ssp + start, ssp + end)
            for i in range(start, end):
                self.assertEqual(es.extract_decoded(ssp + start, ssp + end),
                                 seq4[start:end + 1])
            es.unmirror()

    def run_test_seq_substr_sequential(self, es, seq1, seq2, seq3, seq4):
        start = 3
        end = 13
        er = es.create_reader_with_readmode(readmode.FORWARD, start)
        a = es.alphabet()
        for i in range(start, end):
            self.assertEqual(a.decode(er.next_encoded_char()), seq1[i])
        start = es.seqstartpos(1)
        end = start + 5
        er = es.create_reader_with_readmode(readmode.FORWARD, start)
        for i in range(start - start, end - start):
            self.assertEqual(a.decode(er.next_encoded_char()), seq2[i])
        if es.alphabet().is_dna():
            es.mirror()
            start = es.seqstartpos(2)
            end = start + 5
            er = es.create_reader_with_readmode(readmode.FORWARD, start)
            for i in range(start - start, end - start):
                self.assertEqual(a.decode(er.next_encoded_char()), seq3[i])
            start = es.seqstartpos(3)
            end = start + 5
            er = es.create_reader_with_readmode(readmode.FORWARD, start)
            for i in range(start - start, end - start):
                self.assertEqual(a.decode(er.next_encoded_char()), seq4[i])
            es.unmirror()

if __name__ == "__main__":
    unittest.main()
