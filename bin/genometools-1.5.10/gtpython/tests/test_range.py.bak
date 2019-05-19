#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest
from gt import Range, GTError


class RangeTestCase(unittest.TestCase):

    def test_newrng(self):
        rng = Range(10, 100)
        self.assertEqual(rng.start, 10)
        self.assertEqual(rng.end, 100)
        rng = Range(10, 10)
        rng = Range(0, 0)
        rng = Range()
        self.assertEqual(rng.start, 0)
        self.assertEqual(rng.end, 0)

    def test_invalidrng(self):
        self.assertRaises(GTError, Range, 100, 10)
        self.assertRaises(GTError, Range, -10, 100)
        self.assertRaises(GTError, Range, -10, -100)
        self.assertRaises(GTError, Range, 10, -100)

    def test_setgetrng(self):
        rng = Range(10, 100)

        # change start/end values
        rng.start = 20
        self.assertEqual(rng.start, 20)
        self.assertEqual(rng.end, 100)
        rng.end = 200
        self.assertEqual(rng.start, 20)
        self.assertEqual(rng.end, 200)

        # change start/end to invalid value
        def setend(rng, val):
            rng.end = val

        def setstart(rng, val):
            rng.start = val
        # 20 > 10
        self.assertRaises(GTError, setend, rng, 10)
        self.assertEqual(rng.start, 20)
        self.assertEqual(rng.end, 200)
        # 300 > 200
        self.assertRaises(GTError, setstart, rng, 300)
        self.assertEqual(rng.start, 20)
        self.assertEqual(rng.end, 200)

if __name__ == "__main__":
    unittest.main()
