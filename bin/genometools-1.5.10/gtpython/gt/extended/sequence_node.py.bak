#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2014 Daniel Standage <daniel.standage@gmail.com>
# Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
# Copyright (c) 2009 Center for Bioinformatics, University of Hamburg
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
from gt.core.gtstr import Str
from gt.extended.genome_node import GenomeNode


class SequenceNode(GenomeNode):

    def __init__(self):
        pass

    @classmethod
    def create_new(cls, description, sequence):
        seq_str = Str(sequence.encode('UTF-8'))
        description = description.encode('UTF-8')
        fn = gtlib.gt_sequence_node_new(description, seq_str._as_parameter_)
        n = cls.create_from_ptr(fn, True)
        return n

    def get_description(self):
        return gtlib.gt_sequence_node_get_description(self.gn).decode("UTF-8")

    def get_sequence(self):
        seq = gtlib.gt_sequence_node_get_sequence(self.gn)
        if seq:
            return seq.decode('UTF-8')
        else:
            return ""

    def get_sequence_length(self):
        return gtlib.gt_sequence_node_get_sequence_length(self.gn)

    def register(cls, gtlib):
        from ctypes import c_char_p,  c_void_p, c_ulong
        gtlib.gt_sequence_node_new.restype = c_void_p
        gtlib.gt_sequence_node_new.argtypes = [c_char_p, c_void_p]
        gtlib.gt_sequence_node_get_description.restype = c_char_p
        gtlib.gt_sequence_node_get_description.argtypes = [c_void_p]
        gtlib.gt_sequence_node_get_sequence.restype = c_char_p
        gtlib.gt_sequence_node_get_sequence.argtypes = [c_void_p]
        gtlib.gt_sequence_node_get_sequence_length.restype = c_ulong
        gtlib.gt_sequence_node_get_sequence_length.argtypes = [c_void_p]

    register = classmethod(register)
