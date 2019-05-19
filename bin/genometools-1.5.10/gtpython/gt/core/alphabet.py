#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Daniel Standage <daniel.standage@gmail.com>
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

from gt.dlload import gtlib
from gt.core.error import Error, gterror
from gt.core.gtstr import Str
from gt.core.str_array import StrArray
from ctypes import c_ulong, c_uint, c_int, c_char_p, c_void_p, c_char, \
    c_ubyte, c_int, c_bool, POINTER
import os.path


class Alphabet:

    def create_dna():
        return Alphabet(gtlib.gt_alphabet_new_dna(), True)

    create_dna = staticmethod(create_dna)

    def create_protein():
        return Alphabet(gtlib.gt_alphabet_new_protein(), True)

    create_protein = staticmethod(create_protein)

    def create_empty():
        return Alphabet(gtlib.gt_alphabet_new_empty(), True)

    create_empty = staticmethod(create_empty)

    def create_from_file(indexname):
        if not os.path.exists(indexname + ".al1"):
            raise IOError("file not found: %s" % (indexname + ".al1"))
        e = Error()
        a_ptr = gtlib.gt_alphabet_new_from_file(
            str(indexname).encode("UTF-8"), e._as_parameter_)
        a = Alphabet(a_ptr, True)
        if a == None:
            gterror(e)
        else:
            return a

    create_from_file = staticmethod(create_from_file)

    def create_from_sequence(seqfilenames):
        if not isinstance(seqfilenames, list):
            raise TypeError("argument must be a list of strings")
        e = Error()
        sa = StrArray()
        for f in seqfilenames:
            if not os.path.exists(f):
                raise IOError("file not found: %s" % f)
            sa.add(str(f))
        a_ptr = gtlib.gt_alphabet_new_from_sequence(sa._as_parameter_,
                                                    e._as_parameter_)
        a = Alphabet(a_ptr, True)
        if a == None:
            gterror(e)
        else:
            return a

    create_from_sequence = staticmethod(create_from_sequence)

    def __init__(self, alpha, own=False):
        if own:
            self.alpha = alpha
        else:
            self.alpha = gtlib.gt_alphabet_ref(alpha)
        self._as_parameter_ = self.alpha

    def __del__(self):
        try:
            gtlib.gt_alphabet_delete(self.slpha)
        except AttributeError:
            pass

    def from_param(cls, obj):
        if not isinstance(obj, Alphabet):
            raise TypeError("argument must be an Alphabet")
        return obj._as_parameter_

    from_param = classmethod(from_param)

    def size(self):
        return gtlib.gt_alphabet_size(self.alpha)

    def num_of_chars(self):
        return gtlib.gt_alphabet_num_of_chars(self.alpha)

    def encode(self, plainchar):
        if not self.valid_input(plainchar):
            gterror("No valid input character: %c" % plainchar)
        return gtlib.gt_alphabet_encode(self.alpha, plainchar.encode("UTF-8"))

    def decode(self, encodedchar):
        return gtlib.gt_alphabet_decode(self.alpha, encodedchar).decode("UTF-8")

    def decode_seq(self, encoded):
        encstr = (c_ubyte * len(encoded))(*encoded)
        str_p = gtlib.gt_alphabet_decode_seq_to_str(self.alpha, encstr,
                                                    len(encstr))
        s = Str(str_p)
        return str(s)

    def valid_input(self, plainchar):
        return gtlib.gt_alphabet_valid_input(self.alpha,
                                             plainchar.encode("UTF-8"))

    def is_dna(self):
        return gtlib.gt_alphabet_is_dna(self.alpha)

    def is_protein(self):
        return gtlib.gt_alphabet_is_protein(self.alpha)

    def register(cls, gtlib):
        gtlib.gt_alphabet_new_from_sequence.restype = c_void_p
        gtlib.gt_alphabet_new_from_sequence.argtypes = [c_void_p, c_void_p]
        gtlib.gt_alphabet_new_dna.restype = c_void_p
        gtlib.gt_alphabet_new_dna.argtypes = []
        gtlib.gt_alphabet_new_protein.restype = c_void_p
        gtlib.gt_alphabet_new_protein.argtypes = []
        gtlib.gt_alphabet_new_empty.restype = c_void_p
        gtlib.gt_alphabet_new_empty.argtypes = []
        gtlib.gt_alphabet_new_from_file.restype = c_void_p
        gtlib.gt_alphabet_new_from_file.argtypes = [c_char_p, c_void_p]
        gtlib.gt_alphabet_ref.restype = c_void_p
        gtlib.gt_alphabet_ref.argtypes = [c_void_p]
        gtlib.gt_alphabet_size.restype = c_uint
        gtlib.gt_alphabet_size.argtypes = [c_void_p]
        gtlib.gt_alphabet_num_of_chars.restype = c_uint
        gtlib.gt_alphabet_num_of_chars.argtypes = [c_void_p]
        gtlib.gt_alphabet_encode.restype = c_ubyte
        gtlib.gt_alphabet_encode.argtypes = [c_void_p, c_char]
        gtlib.gt_alphabet_decode.restype = c_char
        gtlib.gt_alphabet_decode.argtypes = [c_void_p, c_ubyte]
        gtlib.gt_alphabet_decode_seq_to_str.restype = c_void_p
        gtlib.gt_alphabet_decode_seq_to_str.argtypes = [c_void_p,
                                                        POINTER(c_ubyte), c_ulong]
        gtlib.gt_alphabet_valid_input.restype = c_int
        gtlib.gt_alphabet_valid_input.argtypes = [c_void_p, c_char]
        gtlib.gt_alphabet_is_dna.restype = c_bool
        gtlib.gt_alphabet_is_dna.argtypes = [c_void_p]
        gtlib.gt_alphabet_is_protein.restype = c_bool
        gtlib.gt_alphabet_is_protein.argtypes = [c_void_p]
        gtlib.gt_alphabet_delete.restype = None
        gtlib.gt_alphabet_delete.argtypes = [c_void_p]
        gtlib.gt_alphabet_ref.restype = c_void_p
        gtlib.gt_alphabet_ref.argtypes = [c_void_p]

    register = classmethod(register)
