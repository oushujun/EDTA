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
from gt.extended.genome_node import GenomeNode
from gt.core.gtstr import Str


class RegionNode(GenomeNode):

    def __init__(self):
        pass

    @classmethod
    def create_new(cls, seqid, start, end):
        if start > end:
            raise Error("start > end")
        seq_str = Str(str(seqid))
        fn = gtlib.gt_region_node_new(seq_str._as_parameter_, start, end)
        n = cls.create_from_ptr(fn, True)
        return n

    def register(cls, gtlib):
        from ctypes import c_ulong, c_void_p
        gtlib.gt_region_node_new.restype = c_void_p
        gtlib.gt_region_node_new.argtypes = [c_void_p, c_ulong, c_ulong]

    register = classmethod(register)
